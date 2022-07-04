###################################################################
## Cellular automata lattice model of human spread through Sahul ##
## SINGLE-SCENARIO AVERAGES (ITERATED)                           ##
## July 2022                                                      ##
## CJA Bradshaw                                                  ##
###################################################################

## remove everything
rm(list = ls())

## libraries
library(sp)
library(rgdal)
library(raster) # note: I have not yet moved over to terra to replace raster
library(oceanmap)
library(insol)
library(OceanView)
library(abind)
library(pracma)
library(binford)
library(rgl)
library(scatterplot3d) 
library(spatstat)
library(spatialEco)
library(SpatialPack)

## functions
# gradient function for immigration
# pI ~ a / (1 + (b * exp(-c * Rk)))
aI <- 0.95; bI <- 5000; cI <- 3
pI.func <- function(Rk) {
  pI <- aI / (1 + (bI * exp(-cI * Rk)))
  return(pI)}
pI.func(4)
xI <- seq(1,10,0.01)
yI <- pI.func(xI)

# gradient function for emigration
aE <- 1; bE <- -3.2
pE.func <- function(Rk) {
  pE <- aE * exp(bE * Rk)
  return(pE)}
pE.func(0.5)
xE <- seq(0.01,1,0.01)
yE <- pE.func(xE)
pE.out <- data.frame(xE,yE)

# stochastic beta sampler (single sample)
stoch.beta.func <- function(mu, var) {
  Sx <- rbeta(length(mu), (((1 - mu) / var - 1 / mu) * mu ^ 2), ((((1 - mu) / var - 1 / mu) * mu ^ 2)*(1 / mu - 1)))
  return(params=Sx)
}

# stochastic beta sampler (n samples)
stoch.n.beta.func <- function(n, mu, var) {
  Sx <- rbeta(n, (((1 - mu) / var - 1 / mu) * mu ^ 2), ((((1 - mu) / var - 1 / mu) * mu ^ 2)*(1 / mu - 1)))
  return(params=Sx)
}

# dynamics model function
Nproj.func <- function(Nt, rm, K) {
  Nt1 <- round(Nt * exp(rm*(1-(Nt/K))), 0)
  return(Nt1)
}

# rescale a range
rscale <- function (x, nx1, nx2, minx, maxx) {
  nx = nx1 + (nx2 - nx1) * (x - minx)/(maxx - minx)
  return(nx)
}

# matrix rotation
rot.mat <- function(x) t(apply(x, 2, rev))

# matrix poisson resampler
rpois.fun <- function(x,y,M) {
  rpois(1,M[x,y])
}
rpois.vec.fun <- Vectorize(rpois.fun,vectorize.args = c('x','y'))

## list coordinates to xyz
coordlist2xyz <- function (list) {
  rl <- length(list[[1]]); cl <- length(list[[2]])
  coords <- c(NA,NA)
  for (r in 1:rl) {
    for (c in 1:cl) {
      coords <- rbind(coords, c(list[[1]][r],list[[2]][c]))
    }
  }
  coords <- coords[-1,]
  return(coordxyz=coords)
}


####################################################
## set grids
####################################################

## relative density, carrying capacity & feedbacks
## NPP (LOVECLIM)
npp <- read.table("data/NppSahul(0-140ka).csv", header=T, sep=",") # 0.5 deg lat resolution
not.na <- which(is.na(npp[,3:dim(npp)[2]]) == F, arr.ind=T)
upper.row <- as.numeric(not.na[1,1])
lower.row <- as.numeric(not.na[dim(not.na)[1],1])
min.lat <- max(npp[not.na[,1], 1])  
max.lat <- min(npp[not.na[,1], 1])
min.lon <- min(npp[not.na[,1], 2])
max.lon <- max(npp[not.na[,1], 2])

sahul.sub <- rep(0, dim(npp)[1])
for (n in 1:dim(npp)[1]) {
  sahul.sub[n] <- ifelse(npp[n,1] <= min.lat & npp[n,1] >= max.lat & npp[n,2] >= min.lon & npp[n,2] <= max.lon, 1, 0)
}  
sah.keep <- which(sahul.sub == 1)
npp.sah <- npp[sah.keep,]

## import sea level & palaeo-lakes layer (1 = land; 2 = palaeo-lake)
sll <- read.table("data/SeaLevelSahul(0-140ka)&LakeC.csv", header=T, sep=",") # 0.5 deg lat resolution
sll.sah <- sll[sah.keep,]

## import ruggedness (average elevational slope) per cell
rug <- read.table("data/RuggednessSahul(0-140ka).csv", header=T, sep=",") # 0.5 deg lat resolution
sah.keep <- which(sahul.sub == 1)
rug.sah <- rug[sah.keep,]

## set base human Leslie matrix from Bradshaw et al. (2014; PNAS)
source("scripts/matrixOperators.r")

# Siler hazard h(x) (Gurven et al. 2007)
# average hunter-gatherer
a1 <- 0.422 # initial infant mortality rate (also known as αt)
b1 <- 1.131 # rate of mortality decline (also known as bt)
a2 <- 0.013 # age-independent mortality (exogenous mortality due to environment); also known as ct
a3 <- 1.47e-04 # initial adult mortality rate (also known as βt)
b3 <- 0.086 # rate of mortality increase
longev <- 80
x <- seq(0,longev,1) # age vector
h.x <- a1 * exp(-b1*x) + a2 + a3 * exp(b3 * x) # Siler's hazard model
l.x <- exp((-a1/b1) * (1 - exp(-b1*x))) * exp(-a2 * x) * exp(a3/b3 * (1 - exp(b3 * x))) # Siler's survival (proportion surviving) model

l.inf <- exp(-a1/b1) # survival at infinite time
T.m <- 1/b1 # time constant at which maturity is approached
h.m <- a2 # hazard for mature animals
l.m <- exp(-a2*x) # survival
h.s <- a3*exp(b3*x) # hazard for senescence
l.s <- exp((a3/b3)*(1 - exp(b3*x))) # survival for senescence
f.x <- a3*exp(b3*x)*exp((a3/b3)/(1-exp(b3*x))) # probability density function
(log(a3) - log(a1)) / a3
T.s <- (1/b3) # modal survival time

## survival
init.pop <- 10000
lx <- round(init.pop*l.x,0)
len.lx <- length(lx)
dx <- lx[1:(len.lx-1)]-lx[2:len.lx]
qx <- dx/lx[1:(length(lx)-1)]
Sx <- 1 - qx
sx <- lx[2:len.lx]/lx[1:(len.lx-1)]
mx <- 1 - sx
Lx <- (lx[1:(len.lx-1)] + lx[2:len.lx])/2
ex <- rev(cumsum(rev(Lx)))/lx[-len.lx]
ex.avg <- ex + x[-len.lx]

# set SD for Sx
Sx.sd <- 0.05 # can set to any value

# fertility (Walker et al. 2006)
primiparity.walker <- c(17.7,18.7,19.5,18.5,18.5,18.7,25.7,19,20.5,18.8,17.8,18.6,22.2,17,16.2,18.4)
prim.mean <- round(mean(primiparity.walker),0)
prim.lo <- round(quantile(primiparity.walker,probs=0.025),0)
prim.hi <- round(quantile(primiparity.walker,probs=0.975),0)
dat.world13 <- read.table("world2013lifetable.csv", header=T, sep=",")
fert.world13 <- dat.world13$m.f
fert.trunc <- fert.world13[1:(longev+1)]
pfert.trunc <- fert.trunc/sum(fert.trunc)
fert.bentley <- 4.69/2 # Bentley 1985 for !Kung
fert.vec <- fert.bentley * pfert.trunc

## construct matrix
stages <- len.lx
popmat <- matrix(0,nrow=stages,ncol=stages)
colnames(popmat) <- x
rownames(popmat) <- x

## populate matrix
popmat[1,] <- fert.vec
diag(popmat[2:stages,]) <- Sx
popmat[stages,stages] <- 0 # Sx[stages-1]
popmat.orig <- popmat ## save original matrix

## matrix properties
r.ann <- max.r(popmat) # rate of population change, 1-yr
#stable.stage.dist(popmat) ## stable stage distribution
R.val(popmat,stages) # reproductive value
gen.l <- G.val(popmat,stages) # mean generation length

## for r.max (set Sx=1)
Sx.1 <- Sx
Sx.1[] <- 1
popmat.max <- popmat.orig
diag(popmat.max[2:stages,]) <- Sx.1
popmat.max[stages,stages] <- 0 # Sx[stages-1]
max.lambda(popmat.max) ## 1-yr lambda
rm.ann <- max.r(popmat.max) # rate of population change, 1-yr

#stable.stage.dist(popmat) ## stable stage distribution
R.val(popmat,stages) # reproductive value
gen.l <- G.val(popmat,stages) # mean generation length

### population dynamics parameters
# dynamical model
# Nt+1 = Nt * exp(rm*(1-(Nt/K))) - (E - I)
lambda.ann <- exp(r.ann) # annual lambda
r.max.NEE <- 2 * log(lambda.ann^gen.l) # human rmax at generational scale (arbitrarily double)
lambda.max.ann <- exp(rm.ann)
rm.max.NEE <- log(lambda.max.ann^gen.l) # human rmax at generational scale (from Sx=1 Leslie matrix)

# Cole's allometric calculation (high)
alpha.Ab <- 15
a.Cole <- -0.16
a.lo.Cole <- -0.41
a.up.Cole <- 0.10
a.sd.Cole <- mean(c((a.Cole - a.lo.Cole)/1.96, (a.up.Cole - a.Cole)/1.96))
b.lo.Cole <- -1.2
b.up.Cole <- -0.79
b.Cole <- -0.99
b.sd.Cole <- mean(c((b.Cole - b.lo.Cole)/1.96, (b.up.Cole - b.Cole)/1.96))
r.max.Cole <- 10^(a.Cole + b.Cole*log10(alpha.Ab)) # from Hone et al. 2003-JApplEcol
r.max.lo.Cole <- 10^(a.lo.Cole + b.lo.Cole*log10(alpha.Ab))
r.max.up.Cole <- 10^(a.up.Cole + b.up.Cole*log10(alpha.Ab))

r.max.up.gen.Cole <- log((exp(r.max.up.Cole))^gen.l)
r.max.gen.Cole <- log((exp(r.max.Cole))^gen.l)
r.max.lo.gen.Cole <- log((exp(r.max.lo.Cole))^gen.l)
r.max.gen.Cole.sd <- mean(c((r.max.gen.Cole - r.max.lo.gen.Cole)/1.96, (r.max.up.gen.Cole - r.max.gen.Cole)/1.96))

## dispersal calculations
# natal dispersal distance function (from Sutherland et al. 2000 Conserv Biol)
# median
a.mid <- 1.45
b.mid <- 0.54
a.lo <- a.mid - 1.05
a.up <- a.mid + 1.05
b.lo <- b.mid - 0.01
b.up <- b.mid + 0.01
M <- 50 # average mass of hunter-gatherer adult
D.med.lo <- a.lo*M^b.lo
D.med.up <- a.up*M^b.up
D.vec <- 1:round((10*111/2), 0)
Pr.Dmed.lo <- exp(-D.vec/D.med.lo) 
Pr.Dmed.up <- exp(-D.vec/D.med.up)

# max
A.mid <- 3.31
B.mid <- 0.65
A.lo <- A.mid - 1.17
A.up <- A.mid + 1.17
B.lo <- B.mid - 0.05
B.up <- B.mid + 0.05
D.max.lo <- A.lo*M^B.lo
D.max.up <- A.up*M^B.up
Pr.Dmx.lo <- exp(-D.vec/D.max.lo) 
Pr.Dmx.up <- exp(-D.vec/D.max.up)

cellD <- D.vec/(111/2)
disp.max.out <- data.frame(cellD,Pr.Dmx.lo,Pr.Dmx.up)

## Hiscock rainfall-territory size relationsip (2008)
terr.rain <- read.table("data/territory.rainfall.Hiscock.csv", header=T, sep=",")
fit.terr.rain <- lm(log10(terr.rain$terrkm2) ~ log10(terr.rain$annrainmm))
fit.dispkm.rain <- lm(log10(terr.rain$terr.rkm) ~ log10(terr.rain$annrainmm))

## make rainfall relative to minimum
fit.dispkm.rain <- lm(log10(terr.rain$terr.rkm) ~ log10(terr.rain$annrainmm/min(terr.rain$annrainmm)))
yDmaxup <- (log10(D.max.up) - as.numeric(coef(fit.dispkm.rain)[1]))/as.numeric(coef(fit.dispkm.rain)[2])
yDminlo <- (log10(D.med.lo) - as.numeric(coef(fit.dispkm.rain)[1]))/as.numeric(coef(fit.dispkm.rain)[2])

hiscock.out <- data.frame(terr.rain$annrainmm/min(terr.rain$annrainmm), log10(terr.rain$terr.rkm))

## Binford's environmental and hunter-gatherer frames of reference (Binford 2001)
## to estimate effect of rugosity on dispersal
binforddat = NULL
binforddat$annual_move <- LRB$dismov
binforddat$annual_rain <- LRB$bio.12
binforddat$ID <- seq(1:length(binforddat$annual_rain))
binford.dat <- as.data.frame(binforddat)
binford.dat$altitude_max <- LRB$h25
binford.dat$altitude_min <- LRB$l25
binford.dat$altitude_dif <- binford.dat$altitude_max - binford.dat$altitude_min
binford.dat <- binford.dat[complete.cases(binford.dat),]
binford.dat <- binford.dat [binford.dat$altitude_dif>0,]  # remove the one case with negative altitude difference

# cube root
binford.dat$annual_move_tr <- binford.dat$annual_move^(1/3)
binford.dat$annual_rain_tr <- binford.dat$annual_rain^(1/3)
binford.dat$altitude_dif_tr <- binford.dat$altitude_dif^(1/3)

res <- lm(annual_move_tr ~ + annual_rain + altitude_dif_tr, data = binford.dat)

# 3D scatterplot with droplines
s3d <- scatterplot3d(binford.dat$annual_rain, binford.dat$altitude_dif_tr, binford.dat$annual_move_tr, pch=16, highlight.3d=F,
                     type="h", main="", angle=210, asp=4, grid=T, 
                     xlab="annual rainfall (mm)", ylab="cube root altitude difference (Δm)", zlab="cube root annual movement (km)")
s3d$plane3d(res, draw_polygon = T)
altdiff.vec <- seq(from=range(binford.dat$altitude_dif)[1], to=range(binford.dat$altitude_dif)[2], by=1)
altdiff.st <- altdiff.vec / max(altdiff.vec)
annmov.pred <- (coef(res)[1] + (altdiff.vec^(1/3) * coef(res)[3]^3) + (mean(binford.dat$annual_rain, na.rm=T) * coef(res)[2]))^3
annmov.st <- annmov.pred / max(annmov.pred)
altmov.dat <- data.frame(altdiff.st, annmov.st)

# exponential decay function fit
# y=a+b(x)^(1/3)
param.init <- c(1, -0.01)
fit.expd <- nls(annmov.st ~ a + b*(altdiff.st)^(1/3), 
                data = altmov.dat,
                algorithm = "port",
                start = c(a = param.init[1], b = param.init[2]),
                trace = TRUE,      
                nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
altdiff.st.vec <- seq(0,1,0.01)
pred.annmov.st <- as.numeric(coef(fit.expd)[1]) + as.numeric(coef(fit.expd)[2]) * (altdiff.st.vec)^(1/3)

## distance to water (based on modern-day water distribution; Damien O'grady)
dir.tmp <- getwd()
rastlist <- list.files(path = dir.tmp, pattern='.tif$', all.files=TRUE, full.names=FALSE)
allrasters <- lapply(rastlist, raster)

d2w <- read.table("Distance2Freshwater.csv", header=T, sep=",") # 0.5 deg lat resolution / units in dd
d2w.sah <- d2w[sah.keep,]


## data from Fagan & Holmes to estimate > mortality rate < MVP size
max.r.decl <- c(-3.24, -1.09, -1.96, -2.28, -0.69, -0.68, -1.88, -1.76, -2.05)
max.lam.decl <- exp(max.r.decl)
1-mean(max.lam.decl)

## set sub-sampled area characteristics
# boxed areas
CEN.crds <- rbind(coordlist2xyz(list(41:60, 16:48)),
                  coordlist2xyz(list(51:59, 49:66)))
TE.crds <- rbind(coordlist2xyz(list(25:37, 36:54)),
                  coordlist2xyz(list(28:36, 26:35)),
                  coordlist2xyz(list(32:45, 55:63)))
E.crds <- coordlist2xyz(list(44:67, 71:87))
NGH.crds <- rbind(coordlist2xyz(list(7:9, 53:59)),
                  coordlist2xyz(list(10, 59)),
                  coordlist2xyz(list(9:10, 60:64)),
                  coordlist2xyz(list(11, 64)),
                  coordlist2xyz(list(10:11, 65:66)),
                  coordlist2xyz(list(12, 66)),
                  coordlist2xyz(list(11:12, 67:68)),
                  coordlist2xyz(list(12:13, 69:70)),
                  coordlist2xyz(list(13:15, 71:72)),
                  coordlist2xyz(list(15:16, 73:74)),
                  coordlist2xyz(list(17, 74)),
                  coordlist2xyz(list(16:17, 75)))
CEN.cells <- dim(CEN.crds)[1]
TE.cells <- dim(TE.crds)[1]
E.cells <- dim(E.crds)[1]
NGH.cells <- dim(NGH.crds)[1]

# bioregions
PILGASMUR.crds <- rbind(coordlist2xyz(list(43:47, 11:22)), 
                        coordlist2xyz(list(50:55, 12:21)), 
                        coordlist2xyz(list(51:60, 22:25)),
                        coordlist2xyz(list(56:58, 15:21)),
                        coordlist2xyz(list(48:49, 11:21)), 
                        coordlist2xyz(list(42, 12:22)), 
                        coordlist2xyz(list(41, 16:20)), 
                        coordlist2xyz(list(59, 18:21)), 
                        coordlist2xyz(list(61, 23:27)), 
                        coordlist2xyz(list(50, 23:25)), 
                        coordlist2xyz(list(56, 13:14)), 
                        coordlist2xyz(list(60, 21:22)), 
                        coordlist2xyz(list(51:55, 11)), 
                        coordlist2xyz(list(51:54, 26)), 
                        coordlist2xyz(list(41, 14)), 
                        coordlist2xyz(list(40, 18)), 
                        coordlist2xyz(list(52, 27)), 
                        coordlist2xyz(list(60, 26))) # Pilbara + Gascoyne + Murchison
WAR.crds <- rbind(coordlist2xyz(list(69, 11:13)),
                    coordlist2xyz(list(70, 12:15)),
                    coordlist2xyz(list(68, 10)),
                    coordlist2xyz(list(68, 12))) # Warren
CER.crds <- rbind(coordlist2xyz(list(51:52, 34:40)), 
                  coordlist2xyz(list(52:53, 42:44)),
                  coordlist2xyz(list(48:49, 36:37)),
                  coordlist2xyz(list(54, 45:47)), 
                  coordlist2xyz(list(55, 39:40)), 
                  coordlist2xyz(list(56, 40:41)), 
                  coordlist2xyz(list(51:52, 41)), 
                  coordlist2xyz(list(47, 38)), 
                  coordlist2xyz(list(49, 35)), 
                  coordlist2xyz(list(52, 33)), 
                  coordlist2xyz(list(53, 38)), 
                  coordlist2xyz(list(53, 40)), 
                  coordlist2xyz(list(53, 45)), 
                  coordlist2xyz(list(55, 47))) # Central Ranges
OVP.crds <- rbind(coordlist2xyz(list(33:37, 37:39)), 
                  coordlist2xyz(list(34:36, 40:42)),
                  coordlist2xyz(list(38, 30:36)), 
                  coordlist2xyz(list(37, 31:32)), 
                  coordlist2xyz(list(37, 35:36)), 
                  coordlist2xyz(list(35:36, 36)), 
                  coordlist2xyz(list(31:32, 41)), 
                  coordlist2xyz(list(32:33, 42)), 
                  coordlist2xyz(list(30:35, 43))) # Ord Victoria Plain
ARP.crds <- rbind(coordlist2xyz(list(25:27, 46:47)),
                  coordlist2xyz(list(25, 48:49)), 
                  coordlist2xyz(list(24, 48)), 
                  coordlist2xyz(list(27, 45)), 
                  coordlist2xyz(list(28, 46))) # Arnhem Plateau
GUPEIU.crds <- rbind(coordlist2xyz(list(34:39, 62:69)),
                    coordlist2xyz(list(36:38, 59:61)), 
                    coordlist2xyz(list(34:37, 57:58)), 
                    coordlist2xyz(list(32:33, 63:67)), 
                    coordlist2xyz(list(38:39, 70:72)),
                    coordlist2xyz(list(32:36, 70)),
                    coordlist2xyz(list(36:37, 71)), 
                    coordlist2xyz(list(41:43, 62)), 
                    coordlist2xyz(list(39:41, 61)), 
                    coordlist2xyz(list(33:34, 56)), 
                    coordlist2xyz(list(31, 63:64)), 
                    coordlist2xyz(list(33, 68:69)), 
                    coordlist2xyz(list(41, 64:65)), 
                    coordlist2xyz(list(41, 68:69)), 
                    coordlist2xyz(list(40, 71:72)), 
                    coordlist2xyz(list(35, 59)), 
                    coordlist2xyz(list(39, 60)), 
                    coordlist2xyz(list(32, 69)), 
                    coordlist2xyz(list(34, 71)), 
                    coordlist2xyz(list(39, 73)), 
                    coordlist2xyz(list(41, 73))) # Gulf Plains + Einasleigh Uplands
BNSNS.crds <- rbind(coordlist2xyz(list(66:73, 74:79)),
                    coordlist2xyz(list(48:56, 78:81)), 
                    coordlist2xyz(list(50:53, 74:77)), 
                    coordlist2xyz(list(61:65, 78:81)),
                    coordlist2xyz(list(58:60, 80:82)), 
                    coordlist2xyz(list(52:57, 82:83)), 
                    coordlist2xyz(list(66:70, 80:81)), 
                    coordlist2xyz(list(75:76, 71:73)), 
                    coordlist2xyz(list(74:75, 74:75)), 
                    coordlist2xyz(list(49, 71:76)), 
                    coordlist2xyz(list(50, 72:73)), 
                    coordlist2xyz(list(47, 80:81)),
                    coordlist2xyz(list(57, 78:79)), 
                    coordlist2xyz(list(77, 72:73)), 
                    coordlist2xyz(list(55:56, 84)), 
                    coordlist2xyz(list(61:63, 82)), 
                    coordlist2xyz(list(65:68, 82)), 
                    coordlist2xyz(list(64:65, 77)),
                    coordlist2xyz(list(69:72, 73)), 
                    coordlist2xyz(list(46, 79)), 
                    coordlist2xyz(list(47, 78)), 
                    coordlist2xyz(list(48, 72)), 
                    coordlist2xyz(list(48, 75)), 
                    coordlist2xyz(list(48, 82)), 
                    coordlist2xyz(list(51, 73)), 
                    coordlist2xyz(list(51, 82)), 
                    coordlist2xyz(list(54, 77)), 
                    coordlist2xyz(list(57, 81)), 
                    coordlist2xyz(list(58, 78)), 
                    coordlist2xyz(list(58, 83)), 
                    coordlist2xyz(list(59, 75)), 
                    coordlist2xyz(list(59, 77)), 
                    coordlist2xyz(list(60, 74)), 
                    coordlist2xyz(list(61, 77)), 
                    coordlist2xyz(list(65, 75)), 
                    coordlist2xyz(list(66, 83)), 
                    coordlist2xyz(list(73, 72)), 
                    coordlist2xyz(list(74, 71)), 
                    coordlist2xyz(list(74, 76)), 
                    coordlist2xyz(list(74, 78)), 
                    coordlist2xyz(list(75, 70)), 
                    coordlist2xyz(list(77, 67))) # Brigalow Belt South + Nandewar + Sydney Basin + NSW South Western Slopes + South Eastern Highlands
MDDBHC.crds <- rbind(coordlist2xyz(list(65:73, 61:63)), 
                      coordlist2xyz(list(65:73, 65:66)), 
                      coordlist2xyz(list(62:63, 60:63)), 
                      coordlist2xyz(list(60:63, 64:65)), 
                      coordlist2xyz(list(64:67, 67:68)), 
                      coordlist2xyz(list(66:71, 60)), 
                      coordlist2xyz(list(67:68, 59)), 
                      coordlist2xyz(list(70:72, 59)), 
                      coordlist2xyz(list(69:73, 64)), 
                      coordlist2xyz(list(61:62, 66)), 
                      coordlist2xyz(list(64:67, 71)), 
                      coordlist2xyz(list(64, 62:63)), 
                      coordlist2xyz(list(74, 62:65)), 
                      coordlist2xyz(list(65, 69:70)), 
                      coordlist2xyz(list(71, 58)), 
                      coordlist2xyz(list(61, 60)), 
                      coordlist2xyz(list(66, 64)), 
                      coordlist2xyz(list(68, 67)), 
                      coordlist2xyz(list(71, 67)), 
                      coordlist2xyz(list(63, 69)), 
                      coordlist2xyz(list(62, 70)), 
                      coordlist2xyz(list(66, 72))) # Murray Darling Depression + Broken Hill Complex
TCH.crds <- rbind(coordlist2xyz(list(83:84, 71:72)),  
                  coordlist2xyz(list(84, 73:74)),  
                  coordlist2xyz(list(85, 71))) # Tasmanian Central Highlands

PILGASMUR.cells <- dim(PILGASMUR.crds)[1]                          
WAR.cells <- dim(WAR.crds)[1]
CER.cells <- dim(CER.crds)[1]
OVP.cells <- dim(OVP.crds)[1]
ARP.cells <- dim(ARP.crds)[1]
GUPEIU.cells <- dim(GUPEIU.crds)[1]
BNSNS.cells <- dim(BNSNS.crds)[1]
MDDBHC.cells <- dim(MDDBHC.crds)[1]
TCH.cells <- dim(TCH.crds)[1]

    ########################################################
    ########################################################
    ########################################################
    ## start diffusion model - iterate for average output ##
    ########################################################
    ########################################################
    ########################################################
    
    ### date of first colonisation?
    entry.date <- 75000
    
    ### how many generations to run?
    gen.run <- 2670 # 2680 # number of generations to run
    round(gen.run*gen.l, 0) # how many years is that?
    entry.date-round(gen.run*gen.l,0) # this means will run until year ...
    
    ### choose linear relationship between NPP and K, or 180 deg-rotated parabolic relationship
    #K.NPP <- "linear"
    K.NPP <- "rotated parabolic"
    #K.NPP <- "rotated quadratic yield density"
    
    ## K reduction scalar
    #modify.K <- "yes"
    modify.K <- "no"
    if (modify.K == "yes") {
      Kreduce <- 0.75 # if yes, by how much?
    } 
    
    ### for unknown SDs, choose % of mean (i.e., for M.cat.sd, pmov.sd, rm.max.NEE.sd)
    #SD.prop.xbar <- 0.05 # 5%
    SD.prop.xbar <- 0.10 # 10%
    
    # impose higher extinction probability below MVP size
    small.pop.ext <- "yes"
    #small.pop.ext <- "no"
    MVP.thresh <- 100 # increase mortality in cells with N < MVP.thresh
    ltMVP.red <- 0.2 # mean (beta resampled) additional mortality expected for cells meeting criteria
    
    ### choose low (2*NEE estimate) or high (generationally scaled Cole's estimate from alpha) rmax
    rmax.sel <- "low" # more defensible
    #rmax.sel <- "high" # seems unrealistically high
    
    ### if a lake is present, how much to reduce K in that grid (0 = 0 K; 1 = full K)?
    lake.red <- 0.01 # cannot be zero
    
    ### max long-distance dispersal modifier
    ldp.mult <- 1 # modify to deviate from theoretical expectations (0 - ∞)
    
    ## modifier of maximum ruggedness effect on movement
    rugmov.mod <- 1
    
    ## set resistance of landscape for distance-to-water function
    watmovresist <- 3 # from Saltré et al. 2009 Ecol Model
    
    ### apply catastrophe function (Reed et al. 2005)?
    catastrophe <- "yes"
    #catastrophe <- "no"
    pop.adjacency <- 0 # to account for 'population' area of influence, in incrementing neighbourhood (1 = immediate neighbours; 2 = 2 cells away in every direction, ...)
    cat.pr <- 0.14/((2*pop.adjacency+1)^2) # probability of catastrophe per generation (Reed et al. 2003) = 0.14, modified by adjacency from above
    M.cat <- 0.75 # mean mortality during a catastrophe event
    M.cat.sd <- SD.prop.xbar*M.cat # sd mortality during a catastrophe event
    
    ### are catastrophe's spatially clustered?
    spatial.cluster <- "yes"
    #spatial.cluster <- "no"
    
    # generate a random point pattern (Thomas cluster process)
    rpp.scale <- 0.015 # controls intensity of clustering (lower values = more clustering)
    kappa.mod <- 1 # modifies Thomas cluster process kappa upward or downward; ~ simulates changes to Pr(catastrophe)
    
    ### print each iteration's map (progression)?
    #print.map <- "yes"
    print.map <- "no"
    
    ### save each iteration's map as a jpg?
    save.map <- "no"
    #save.map <- "yes"
    
    # save output images?
    #save.outputs <- "yes"
    save.outputs <- "no"
    
    ### pick entry cell
    #start.col1 <- 45; start.row1 <- 1 # north of Bird's head, N Guinea
    #start.col1 <- 40; start.row1 <- 4 # west of Bird's head, N Guinea
    #start.col1 <- 38; start.row1 <- 21 # N Sahul shelf
    start.col1 <- 31; start.row1 <- 27 # S Sahul shelf
    #start.col1 <- 8; start.row1 <- 58 # SWA
    #start.col1 <- 87; start.row1 <- 59 # SNSW
    
    if (entry.date > 75000 & start.col1 == 31) {
      start.col1 <- 32}
    
    ### add secondary colonisation event?
    #second.col <- "no"
    second.col <- "yes"
    
    ### lag between first and second colonisation events
    lag.l <- 1000 # lag (between 1st & 2nd colonisation events) length?
    lag.n <- 2 # how many lag increments?
    start.time2 <- ifelse((lag.n * round(lag.l/gen.l)) == 0, 1, (lag.n * round(lag.l/gen.l))) # generations after first colonisation event (increments of ~ lag years)
    
    if (second.col == "yes") {
      #start.col2 <- 38; start.row2 <- 21 # N Sahul shelf
      #start.col2 <- 31; start.row2 <- 27 # S Sahul shelf
      start.col2 <- 40; start.row2 <- 4 # west of Bird's head, N Guinea
    }  
    
    if (entry.date > 75000 & start.col2 == 31) {
      start.col2 <- 32}
    
    # add this to jpg file names for scenario description
    name.add <- paste(entry.date/1000, "ka.", gen.run, "gen.", ifelse(K.NPP=="rotated quadratic yield density", "rqydK", ifelse(K.NPP=="linear", "linK", "parK")), ".rmax", ifelse(rmax.sel == "low","-lo","-hi"), ifelse(ldp.mult != 1, paste(".ldispmod", ldp.mult, sep=""), ""), ifelse(catastrophe=="yes",".CAT", ".NOCAT"), M.cat*100, ifelse(spatial.cluster=="yes",paste("cl",rpp.scale,sep=""),""), ".neigh-adj", pop.adjacency, ".", "intro", ifelse(second.col=="no",1,2), ifelse(start.row1==1, "nBH", ifelse(start.row1==4, "wBH", "SSh")), ifelse(second.col=="yes", ifelse(second.col=="yes" & start.row2==21,"-SSh", ifelse(start.row2==4, "-wBH", "-nBH")), ""), ifelse(second.col=="yes", paste("-lag", round(lag.l*lag.n/gen.l), "g", sep=""), ""), sep="") 
    
    ### minimum viable population size to consider a cell 'colonised' (for first-appearance map)
    min.pop.size <- 100
    
    ### start founding population in 1 cell
    N.found.mod <- 1
    N.found.lo <- 1300*N.found.mod; N.found.up <- 1500*N.found.mod
    
    ## estimate SD of rmax according to SD proportions et above
    rm.max.NEE.sd <- SD.prop.xbar * rm.max.NEE
    
    ### dispersal parameters
    pmov.mean <- 1/3 # mu for beta function indicating proportion of people immigrating/emigrating
    pmov.sd <- SD.prop.xbar*pmov.mean # sd for beta function
    # stoch.beta.func(pmov.mean, pmov.sd)
    # long.disp <- 5 # Poisson lambda distance (in cell units) for random long-range dispersal // no longer used
    # long.disp.mu <- 0.04 # beta mu for long-term dispersal probability & proportion dispersing // no longer used //
    # long.disp.sd <- 0.01 # beta sd for above // no longer used //
    # stoch.beta.func(long.disp.mu, long.disp.sd)
    NK.emig.min <- 0.3; NK.emig.max <- 0.7 # N/K when pmov.mean emigrates
    
    # update ruggedness movement function with rugmov.mod
    rugmovmod.func <- function(x) {
      rugmovmod <- as.numeric(coef(fit.expd)[1]) + rugmov.mod*as.numeric(coef(fit.expd)[2]) * (x)^(1/3)
      return(rugmovmod=rugmovmod)
    }
    
    # npp @ entry date ka
    sub.entry <- which(colnames(npp.sah) == paste("X",entry.date,sep=""))
    npp.sah.entry <- npp.sah[,c(1,2,sub.entry)]

    # plot raster
    coordinates(npp.sah.entry) = ~ Lon + Lat
    proj4string(npp.sah.entry)=CRS("+proj=longlat +datum=WGS84") # set it to lat-long
    gridded(npp.sah.entry) = TRUE
    npp.entry = raster(npp.sah.entry)
    lim.exts <- 5

    # transform to array
    lz <- dim(npp.sah)[2] - 2
    npp.array <- array(data=NA, dim=c(dim(raster2matrix(npp.entry)),lz))
    for (k in 3:(lz+2)) {
      npp.sah.k <- npp.sah[,c(1,2,k)] 
      coordinates(npp.sah.k) = ~ Lon + Lat
      proj4string(npp.sah.k)=CRS("+proj=longlat +datum=WGS84") # set it to lat-long
      gridded(npp.sah.k) = TRUE
      npp.k = raster(npp.sah.k)
      npp.array[,,k-2] <- raster2matrix(npp.k)
    }

    ## calculate all Ks as relative to current
    npp.sah.rel <- npp.array
    for (z in 1:lz) {
      npp.sah.rel[,,z] <- npp.array[,,z] / npp.array[,,3]
    }
    npp.sah.rel[,,3] <- npp.array[,,1]
    
    # npp to K
    hum.dens.med <- 6.022271e-02
    hum.dens.lq <- 3.213640e-02
    hum.dens.uq <- 1.439484e-01
    hum.dens.max <- 1.152206e+00
    hum.dens.min <- 1.751882e-02
    cell.area <- (111.12/2)*(111.12/2) # km2
    
    # create vector of K reduction scalars for projection interval
    if (modify.K == "yes") {
      Kreduce.vec <- stoch.n.beta.func(lz, Kreduce, 0.05*Kreduce)
    }
    if (modify.K == "no") {
      Kreduce <- rep(1,lz)
    } 
    
    # modify underlying K magnitude by modifying NPP across the board
    K.array <- npp.sah.rel
    for (z in 1:lz) {
      K.array[,,z] <- rscale(npp.array[,,z], round(hum.dens.min*cell.area, 0), round(hum.dens.max*cell.area, 0), min(npp.array[,,z], na.rm=T), max(npp.array[,,z], na.rm=T))
    }

    # 180-deg rotated parabola
    # y = a(x - h)^2 + k
    # h = median NPP; k = max K; a = negative for 180 flipped
    k.Kmax <- max(K.array, na.rm=T)/2
    h.NPPmed <- mean(npp.array, na.rm=T)
    h.NPPmed <- mean(range(npp.array, na.rm=T))
    Kmin <- min(K.array, na.rm=T)
    NPP.seq <- seq(min(npp.array, na.rm=T), max(npp.array, na.rm=T), 0.01)
    K.parab.pred <- (-3 * (NPP.seq - h.NPPmed)^2) + k.Kmax
    K.parab.pred.rescale <- rscale(K.parab.pred, round(hum.dens.min*cell.area, 0), 0.5*round(hum.dens.max*cell.area, 0), min(K.parab.pred), max(K.parab.pred))

    # slow exponential increase combined with peak
    # reciprocal quadratic yield density
    # y = x/(a + b*x + c*x^2)
    # y = K, x = NPP
    a.rqyd <- 200; b.rqyd <- 0.60; c.rqyd <- 0.2
    K.rqyd.pred <- a.rqyd * exp(-(NPP.seq-b.rqyd)^2/(2*c.rqyd^2))
    K.rqyd.pred.rescale <- rscale(K.rqyd.pred, round(hum.dens.min*cell.area, 0), 0.5*round(hum.dens.max*cell.area, 0), min(K.rqyd.pred), max(K.rqyd.pred))

    K.lin.x <- c(min(npp.array, na.rm=T), max(npp.array, na.rm=T))
    K.lin.y <- c(min(K.array, na.rm=T), max(K.array, na.rm=T))
    fit.K.lin <- lm(K.lin.y ~ K.lin.x)
    K.lin.pred <- as.numeric(coef(fit.K.lin)[1]) + as.numeric(coef(fit.K.lin)[2])*NPP.seq

    # rotated parabolic
    K.array.parab <- (-3 * (npp.array - h.NPPmed)^2) + k.Kmax
    K.array.parab.rescale <- K.array.parab
    for (z in 1:lz) {
      K.array.parab.rescale[,,z] <- rscale(K.array.parab[,,z], round(hum.dens.min*cell.area, 0), round(hum.dens.max*cell.area, 0), min(K.array.parab[,,z], na.rm=T), max(K.array.parab[,,z], na.rm=T))
    }

    # reciprocal quadratic yield density
    K.array.rqyd <- a.rqyd * exp(-(npp.array - b.rqyd)^2 / (2*c.rqyd^2))
    K.array.rqyd.rescale <- K.array.rqyd
    for (z in 1:lz) {
      K.array.rqyd.rescale[,,z] <- rscale(K.array.rqyd[,,z], round(hum.dens.min*cell.area, 0), round(hum.dens.max*cell.area, 0), min(K.array.rqyd[,,z], na.rm=T), max(K.array.rqyd[,,z], na.rm=T))
    }

    # rescale so that parabolic total K = linear total K
    hist.K.array <- hist(K.array, br=12)
    hist.K.array.dat <- data.frame(hist.K.array$mids, hist.K.array$density)
    
    # rescale K.array.parab.rescale to same sum as K.array (distribution of Ks = same total)
    K.array.parab.rescale2 <- K.array.parab.rescale / (sum(K.array.parab.rescale, na.rm=T)/sum(K.array, na.rm=T))
    sum(K.array.parab.rescale2, na.rm=T)
    
    K.parab.pred.rescale2 <- K.parab.pred.rescale / (sum(K.array.parab.rescale, na.rm=T)/sum(K.array, na.rm=T))
    hist.K.parab.pred.rescale2 <- hist(K.parab.pred.rescale2,br=12)
    hist.K.parab.pred.rescale2.dat <- data.frame(hist.K.parab.pred.rescale2$mids, hist.K.parab.pred.rescale2$density)

    # rescale K.array.parab.rescale to same sum as K.array (distribution of Ks = same total)
    K.array.rqyd.rescale2 <- K.array.rqyd.rescale / (sum(K.array.rqyd.rescale, na.rm=T)/sum(K.array, na.rm=T))
    sum(K.array.rqyd.rescale2, na.rm=T)
    
    K.rqyd.pred.rescale2 <- K.rqyd.pred.rescale / (sum(K.array.rqyd.rescale, na.rm=T)/sum(K.array, na.rm=T))
    hist.K.rqyd.pred.rescale2 <- hist(K.rqyd.pred.rescale2,br=12)
    hist.K.rqyd.pred.rescale2.dat <- data.frame(hist.K.rqyd.pred.rescale2$mids, hist.K.rqyd.pred.rescale2$density)
    
    NPP.K.out <- data.frame(NPP.seq,K.lin.pred,K.parab.pred.rescale2,K.rqyd.pred.rescale2)
    colnames(NPP.K.out) <- c("NPP","K.lin","K.para","K.rqyd")

    # rotate matrix -90 & renumber from oldest to youngest
    if (K.NPP == "linear") {
      K.array.use <- K.array
    }
    if (K.NPP == "rotated parabolic") {
      K.array.use <- K.array.parab.rescale
    }
    if (K.NPP == "rotated quadratic yield density") {
      K.array.use <- K.array.rqyd.rescale
    }
    
    K.rot.array <- array(data=NA, c(dim(K.array.use)[2], dim(K.array.use)[1], lz))
    for (z in 1:lz) {
      K.rot.array[,,z] <- apply(t(K.array.use[,,142-z]),2,rev)
    }

    if (modify.K == "yes") {
      for (z in 1:lz) {
        K.rot.array[,,z] <- K.rot.array[,,z] * Kreduce.vec[z]
      }
    }
    
    ## block out Indonesia & never-connected islands (make NA)
    K.rot.array[1:20, 1:39, ] <- NA
    K.rot.array[6:20, 40:44, ] <- NA
    K.rot.array[9:17, 45:46, ] <- NA
    K.rot.array[21:22, 24:33, ] <- NA
    K.rot.array[27:28, 24:26, ] <- NA
    K.rot.array[34:35, 77:83, ] <- NA
    K.rot.array[22:23, 85:87, ] <- NA
    K.rot.array[18, 84:85, ] <- NA
    K.rot.array[9:12, 77:86, ] <- NA
    K.rot.array[2:6, 70:87, ] <- NA
    K.rot.array[20, 71, ] <- NA
    K.rot.array[2, 52, ] <- NA
    
    # block passage to Tasmania (70 to 67K; 60 to 46K; 43-42K cannot cross)
    K.rot.array[80, 68:77, c(71:101)] <- NA
    K.rot.array[79, 68:77, c(71:101)] <- NA
    
    # interpolate between 1000-year NPP intervals per human generation
    # interpolate Ks at gen.l intervals between 1000-yr slices
    mill.gen.div <- round(1000/gen.l, 0)
    subtentry <- dim(K.rot.array)[3] - (entry.date/1000)
    Kentry.array <- K.rot.array[,,subtentry:(dim(K.rot.array)[3])]

    K.start <- entry.date
    K.end <- 0
    yr.proj.vec <- seq(K.start, K.end, -round((1000/mill.gen.div),0))
    Kentry.interp.array <- array(data=0, dim=c(dim(Kentry.array[,,1])[1],dim(Kentry.array[,,1])[2],length(yr.proj.vec)))
    mill.vec <- seq(K.start,K.end,-1000)
    
    for (i in 1:dim(Kentry.array)[1]) {
      for (j in 1:dim(Kentry.array)[2]) {
        if (length(which(is.na(Kentry.array[i,j,]))==T) <  (dim(Kentry.array)[3]-1))
          Kentry.interp.array[i,j,] <- approx(mill.vec, Kentry.array[i,j,], xout=yr.proj.vec, method="linear")$y
        else {
          Kentry.interp.array[i,j,] <- rep(NA, length(yr.proj.vec))
        }
      }
    }

    # transform sea level & palaeo-lakes file (sll) to an array
    lz <- dim(sll.sah)[2] - 2
    sll.array <- array(data=NA, dim=c(dim(raster2matrix(npp.entry)),lz))
    for (k in 3:(lz+2)) {
      sll.sah.k <- sll.sah[,c(1,2,k)] 
      coordinates(sll.sah.k) = ~ Lon + Lat
      proj4string(sll.sah.k)=CRS("+proj=longlat +datum=WGS84") # set it to lat-long
      gridded(sll.sah.k) = TRUE
      sll.k = raster(sll.sah.k)
      sll.array[,,k-2] <- raster2matrix(sll.k)
    }

    # rotate matrix -90 & renumber from oldest to youngest
    sll.rot.array <- array(data=NA, c(dim(sll.array)[2], dim(sll.array)[1], lz))
    for (z in 1:lz) {
      sll.rot.array[,,z] <- apply(t(sll.array[,,142-z]),2,rev)
    }

    ## copy values between 1000-yr intervals
    sllentry.array <- sll.rot.array[,,subtentry:(dim(sll.rot.array)[3])]
    sllentry.copy.array <- array(data=0, dim=c(dim(sllentry.array[,,1])[1],dim(sllentry.array[,,1])[2],length(yr.proj.vec)))
    
    for (i in 1:dim(sllentry.array)[1]) {
      for (j in 1:dim(sllentry.array)[2]) {
        if (length(which(is.na(sllentry.array[i,j,]))==T) < (dim(sllentry.array)[3]-1))
          sllentry.copy.array[i,j,] <- approx(mill.vec, sllentry.array[i,j,], xout=yr.proj.vec, method="linear")$y
        else {
          sllentry.copy.array[i,j,] <- rep(NA, length(yr.proj.vec))
        }
      }
    }

    # transform distance to water file (d2w) to an array
    lz <- dim(d2w.sah)[2] - 2
    d2w.array <- array(data=NA, dim=c(dim(raster2matrix(npp.entry)),lz))
    for (k in 3:(lz+2)) {
      d2w.sah.k <- d2w.sah[,c(1,2,k)] 
      coordinates(d2w.sah.k) = ~ Lon + Lat
      proj4string(d2w.sah.k)=CRS("+proj=longlat +datum=WGS84") # set it to lat-long
      gridded(d2w.sah.k) = TRUE
      d2w.k = raster(d2w.sah.k)
      d2w.array[,,k-2] <- raster2matrix(d2w.k)
    }
    # just use matrix
    d2w.mat <- t(as.matrix(d2w.array[,,1]))
    
    # transform ruggedness file to array
    lz <- dim(rug.sah)[2] - 2
    rug.array <- array(data=NA, dim=c(dim(raster2matrix(npp.entry)),lz))
    for (k in 3:(lz+2)) {
      rug.sah.k <- rug.sah[,c(1,2,k)] 
      coordinates(rug.sah.k) = ~ Lon + Lat
      proj4string(rug.sah.k)=CRS("+proj=longlat +datum=WGS84") # set it to lat-long
      gridded(rug.sah.k) = TRUE
      rug.k = raster(rug.sah.k)
      rug.array[,,k-2] <- raster2matrix(rug.k)
    }
   
    # rescale rugosity from 0-1
    rug.array.rescale <- rug.array
    for (z in 1:lz) {
      rug.array.rescale[,,z] <- rscale(rug.array[,,z], 0, 1, min(rug.array[,,z], na.rm=T), max(rug.array[,,z], na.rm=T))
    }

    # rotate matrix -90 & renumber from oldest to youngest
    rug.rot.array <- array(data=NA, c(dim(rug.array.rescale)[2], dim(rug.array.rescale)[1], lz))
    for (z in 1:lz) {
      rug.rot.array[,,z] <- apply(t(rug.array.rescale[,,142-z]),2,rev)
    }

    ## interpolate between 1000-yr intervals
    rugentry.array <- rug.rot.array[,,subtentry:(dim(rug.rot.array)[3])]
    rugentry.interp.array <- array(data=0, dim=c(dim(rugentry.array[,,1])[1],dim(rugentry.array[,,1])[2],length(yr.proj.vec)))
    
    for (i in 1:dim(rugentry.array)[1]) {
      for (j in 1:dim(rugentry.array)[2]) {
        if (length(which(is.na(rugentry.array[i,j,]))==T) <  (dim(rugentry.array)[3]-1))
          rugentry.interp.array[i,j,] <- approx(mill.vec, rugentry.array[i,j,], xout=yr.proj.vec, method="linear")$y
        else {
          rugentry.interp.array[i,j,] <- rep(NA, length(yr.proj.vec))
        }
      }
    }

    # spatial clustering of catastraphe events controlling parameters
    kappa.mult <- 0.9
    cellslo <- 1
    cellshi <- 3772
    kappa.mult.up <- 1.2*kappa.mod
    kappa.mult.lo <- 0.3*kappa.mod
    kappa.rep <- seq(kappa.mult.up,kappa.mult.lo,-(kappa.mult.up-kappa.mult.lo)/cellshi)
    cells.rep <- seq(cellslo,cellshi, (cellshi-cellslo)/cellshi)
    kappa.fit <- lm(kappa.rep ~ cells.rep)
    kappaP.func <- function(cells.occ) {
      kappa.pred <- (as.numeric(coef(kappa.fit)[1])) + as.numeric(coef(kappa.fit)[2])*cells.occ
      return(kappa.pred)
    }
    rpp.mu.mult <- 0.6 # this, with the fluctuating kappa multiplier, keeps overall mean proportion of cells experiencing catastrophic mortality ~ 0.14
    
    proc.sim.start <- proc.time()
    reps <- 100 # number of times to redo spread model to find average values
    first.arrive.array <- array(data=NA, dim=c(dim(rugentry.interp.array[,,1]), reps))
    N.mat <- NG.N.mat <- AUS.N.mat <- CEN.N.mat <- TE.N.mat <- E.N.mat <- NGH.N.mat <-
      PILGASMUR.N.mat <- WAR.N.mat <- CER.N.mat <- OVP.N.mat <- ARP.N.mat <- GUPEIU.N.mat <- BNSNS.N.mat <- MDDBHC.N.mat <- TCH.N.mat <- 
      matrix(data=NA,nrow=reps,ncol=(gen.run+1))
    CEN.D.mat <- TE.D.mat <- E.D.mat <- NGH.D.mat <- PILGASMUR.D.mat <- WAR.D.mat <- CER.D.mat <- OVP.D.mat <- ARP.D.mat <- 
      GUPEIU.D.mat <- BNSNS.D.mat <- MDDBHC.D.mat <- TCH.D.mat <- 
      matrix(data=NA,nrow=reps,ncol=(gen.run+1))
    
    ALL.mode.imm.dir.arr <- MIS3.mode.imm.dir.arr <- MIS2.mode.imm.dir.arr <- LGM.mode.imm.dir.arr <- LGMc.mode.imm.dir.arr <- 
      LGMu.mode.imm.dir.arr <- HOL.mode.imm.dir.arr <- array(data=NA, dim=c(dim(Kentry.interp.array[,,1]), reps))
    
    L26.mode.imm.dir.arr <- L25.mode.imm.dir.arr <- L24.mode.imm.dir.arr <- L23.mode.imm.dir.arr <- L22.mode.imm.dir.arr <-
      L21.mode.imm.dir.arr <- L20.mode.imm.dir.arr <- L19.mode.imm.dir.arr <- L18.mode.imm.dir.arr <- L17.mode.imm.dir.arr <-
      L16.mode.imm.dir.arr <- L15.mode.imm.dir.arr <- L14.mode.imm.dir.arr <- L13.mode.imm.dir.arr <- L12.mode.imm.dir.arr <-
      L11.mode.imm.dir.arr <- L10.mode.imm.dir.arr <- L09.mode.imm.dir.arr <- L08.mode.imm.dir.arr <- L07.mode.imm.dir.arr <-
      array(data=NA, dim=c(dim(Kentry.interp.array[,,1]), reps))
    
    L2624.mode.imm.dir.arr <- L2422.mode.imm.dir.arr <- L2220.mode.imm.dir.arr <- L2018.mode.imm.dir.arr <- L1816.mode.imm.dir.arr <-
      L1614.mode.imm.dir.arr <- L1412.mode.imm.dir.arr <- L1210.mode.imm.dir.arr <- L1008.mode.imm.dir.arr <- L0806.mode.imm.dir.arr <-
      array(data=NA, dim=c(dim(Kentry.interp.array[,,1]), reps))
    
    maxt <- dim(Kentry.interp.array[,,1:(gen.run+1)])[3] # 75 ka
    MIS3.rnge <- c(964, 1373); MIS3.t.en <- maxt - MIS3.rnge[1]; MIS3.t.st <- maxt - MIS3.rnge[2]
    MIS2.rnge <- c(1648, 2271); MIS2.t.en <- maxt - MIS2.rnge[1]; MIS2.t.st <- maxt - MIS2.rnge[2]
    LGM.rnge <- c(1863, 2007); LGM.t.en <- maxt - LGM.rnge[1]; LGM.t.st <- maxt - LGM.rnge[2]
    LGMc.rnge <- c(1762, 1974); LGMc.t.en <- maxt - LGMc.rnge[1]; LGMc.t.st <- maxt - LGMc.rnge[2] # 25800-19900 
    LGMu.rnge <- c(1562, 2131); LGMu.t.en <- maxt - LGMu.rnge[1]; LGMu.t.st <- maxt - LGMu.rnge[2] # 15500-31400
    HOL.rnge <- c(2007, maxt); HOL.t.en <- maxt - HOL.rnge[1]; HOL.t.st <- maxt - HOL.rnge[2]
    
    L26.rnge <- c(1755,1791); L26.t.en <- L26.rnge[2]; L26.t.st <- L26.rnge[1]
    L25.rnge <- c(1791, 1827); L25.t.en <- L25.rnge[2]; L25.t.st <- L25.rnge[1]
    L24.rnge <- c(1827, 1862); L24.t.en <- L24.rnge[2]; L24.t.st <- L24.rnge[1]
    L23.rnge <- c(1862, 1898); L23.t.en <- L23.rnge[2]; L23.t.st <- L23.rnge[1]
    L22.rnge <- c(1898, 1934); L22.t.en <- L22.rnge[2]; L22.t.st <- L22.rnge[1]
    L21.rnge <- c(1934, 1970); L21.t.en <- L21.rnge[2]; L21.t.st <- L21.rnge[1]
    L20.rnge <- c(1970, 2006); L20.t.en <- L20.rnge[2]; L20.t.st <- L20.rnge[1]
    L19.rnge <- c(2006, 2041); L19.t.en <- L19.rnge[2]; L19.t.st <- L19.rnge[1]
    L18.rnge <- c(2041, 2077); L18.t.en <- L18.rnge[2]; L18.t.st <- L18.rnge[1]
    L17.rnge <- c(2077, 2113); L17.t.en <- L17.rnge[2]; L17.t.st <- L17.rnge[1]
    L16.rnge <- c(2113, 2149); L16.t.en <- L16.rnge[2]; L16.t.st <- L16.rnge[1]
    L15.rnge <- c(2149, 2185); L15.t.en <- L15.rnge[2]; L15.t.st <- L15.rnge[1]
    L14.rnge <- c(2185, 2220); L14.t.en <- L14.rnge[2]; L14.t.st <- L14.rnge[1]
    L13.rnge <- c(2220, 2256); L13.t.en <- L13.rnge[2]; L13.t.st <- L13.rnge[1]
    L12.rnge <- c(2256, 2292); L12.t.en <- L12.rnge[2]; L12.t.st <- L12.rnge[1]
    L11.rnge <- c(2292, 2328); L11.t.en <- L11.rnge[2]; L11.t.st <- L11.rnge[1]
    L10.rnge <- c(2328, 2364); L10.t.en <- L10.rnge[2]; L10.t.st <- L10.rnge[1]
    L09.rnge <- c(2364, 2400); L09.t.en <- L09.rnge[2]; L09.t.st <- L09.rnge[1]
    L08.rnge <- c(2400, 2435); L08.t.en <- L08.rnge[2]; L08.t.st <- L08.rnge[1]
    L07.rnge <- c(2435, 2471); L07.t.en <- L07.rnge[2]; L07.t.st <- L07.rnge[1]
    
    ##########################
    ## start of m reps loop ##
    ##########################
    
    for (m in 1:reps) {
      
      ## parallel set-up
      # cores <- 12
      # getDoParWorkers()
      # getDoParName()
      # cl <- makeCluster(cores, type = "SOCK")
      # registerDoSNOW(cl)
      # getDoParWorkers()
      # getDoParName()
      # getDoParVersion()
      
      # set up NA array
      NA.array <- Kentry.interp.array * 0
      NA.array <- NA.array[,,1:(gen.run+1)]
      
      # set direction codes
      dir.vec <- c("NW", "N", "NE", "W", "E", "SW", "S", "SE")
      
      # generations to run
      gen.run * gen.l # numbers years to run
      
      # set up N array
      array.N <- Kentry.interp.array * 0
      array.N <- array.N[,,1:(gen.run+1)]
      
      # set up direction array (dominant direction of influx)
      dir.array <- array.N
      array.N[start.row1, start.col1, 1] <- round(runif(1, N.found.lo, N.found.up), 0)
      
      # log10 relative K array
      Kentry.interp.lrel.array <- log10(Kentry.interp.array / min(Kentry.interp.array, na.rm=T)) # log10 relative K
      
      ## assume same slope between relative NPP and radius movement
      disp.npp.slope <- as.numeric(coefficients((fit.dispkm.rain))[2])
      disp.npp.int <- as.numeric(coefficients((fit.dispkm.rain))[1])
      disp.npp.max.int <- 1.6646371
      
      i.rows <- dim(array.N[,,1])[1]
      j.cols <- dim(array.N[,,1])[2]
      z.layers <- dim(array.N)[3]
      
      # storage vectors
      N.vec <- NG.N.vec <- AUS.N.vec <- CEN.N.vec <- TE.N.vec <- E.N.vec <- NGH.N.vec <- 
        PILGASMUR.N.vec <- WAR.N.vec <- CER.N.vec <- OVP.N.vec <- ARP.N.vec <- GUPEIU.N.vec <- BNSNS.N.vec <- MDDBHC.N.vec <- TCH.N.vec <- 
        poparea.vec <- pc.complete <- cat.pr.est <- rep(0,z.layers)
      CEN.D.vec <- TE.D.vec <- E.D.vec <- NGH.D.vec <- 
        PILGASMUR.D.vec <- WAR.D.vec <- CER.D.vec <- OVP.D.vec <- ARP.D.vec <- GUPEIU.D.vec <- BNSNS.D.vec <- MDDBHC.D.vec <- TCH.D.vec <- rep(0,z.layers)
      N.vec[1] <- AUS.N.vec <- array.N[start.row1,start.col1,1]
      if (second.col == "yes") {
        N.vec[start.time2] <- array.N[start.row1,start.col1,start.time2]
        NG.N.vec[start.time2] <- array.N[start.row1,start.col1,start.time2]
      }
      poparea.vec[1] <- cell.area/1000
      
      #############################
      ## project
      for (t in 1:(z.layers-1)) {
        
        # Poisson-resampled K matrix
        Kentry.interp.poiss <- outer(1:nrow(round(Kentry.interp.array[,,t],0)), 1:ncol(round(Kentry.interp.array[,,t],0)), rpois.vec.fun, round(Kentry.interp.array[,,t],0))
        
        # reduce Ks if lake present
        Kentry.interp.pois <- Kentry.interp.poiss
        for (i in 1:i.rows) { # i rows
          for (j in 1:j.cols) { # j columns
            Kentry.interp.pois[i,j] <- ifelse(sllentry.copy.array[i,j,t] > 1, lake.red * Kentry.interp.poiss[i,j], Kentry.interp.poiss[i,j])
          }
        }
        
        # step through array in t for immigration & emigration
        for (i in 1:i.rows) { # i rows
          for (j in 1:j.cols) { # j columns
            
            if (second.col == "yes") {
              if (t == start.time2) {
                array.N[start.row2, start.col2, start.time2] <- round(runif(1, N.found.lo, N.found.up), 0)
              }
            }
            
            # set cell-cell gradients relative to focal cell
            NW.RK <- ifelse(i > 1 & j > 1, (Kentry.interp.pois[i,j] / ifelse(length(Kentry.interp.pois[i-1,j-1]) > 0, Kentry.interp.pois[i-1,j-1], NA)), NA)  # NW
            N.RK <- ifelse(i > 1, (Kentry.interp.pois[i,j] / ifelse(length(Kentry.interp.pois[i-1,j]) > 0, Kentry.interp.pois[i-1,j], NA)), NA) # N
            NE.RK <- ifelse(i > 1 & j < j.cols, (Kentry.interp.pois[i,j] / ifelse(length(Kentry.interp.pois[i-1,j+1]) > 0, Kentry.interp.pois[i-1,j+1], NA)), NA)  # NE
            W.RK <- ifelse(j > 1, (Kentry.interp.pois[i,j] / ifelse(length(Kentry.interp.pois[i,j-1]) > 0, Kentry.interp.pois[i,j-1], NA)), NA)    # W
            E.RK <- ifelse(j < j.cols, (Kentry.interp.pois[i,j] / ifelse(length(Kentry.interp.pois[i,j+1]) > 0, Kentry.interp.pois[i,j+1], NA)), NA)      # E
            SW.RK <- ifelse(i < i.rows & j > 1, (Kentry.interp.pois[i,j] / ifelse(length(Kentry.interp.pois[i+1,j-1]) > 0, Kentry.interp.pois[i+1,j-1], NA)), NA)  # SW
            S.RK <- ifelse(i < i.rows, (Kentry.interp.pois[i,j] / ifelse(length(Kentry.interp.pois[i+1,j]) > 0, Kentry.interp.pois[i+1,j], NA)), NA)      # S
            SE.RK <- ifelse(i < i.rows & j < j.cols, (Kentry.interp.pois[i,j] / ifelse(length(Kentry.interp.pois[i+1,j+1]) > 0, Kentry.interp.pois[i+1,j+1], NA)), NA)  # SE
            
            ## current population distances from K in each cell
            focal.dK1 <- array.N[i,j,t] / Kentry.interp.pois[i,j]
            NW.dK1 <- ifelse(i > 1 & j > 1, array.N[i-1,j-1,t] / Kentry.interp.pois[i-1,j-1], NA)
            N.dK1 <- ifelse(i > 1, array.N[i-1,j,t] / Kentry.interp.pois[i-1,j], NA)
            NE.dK1 <- ifelse(i > 1 & j < j.cols, array.N[i-1,j+1,t] / Kentry.interp.pois[i-1,j+1], NA)
            W.dK1 <- ifelse(j > 1, array.N[i,j-1,t] / Kentry.interp.pois[i,j-1], NA)
            E.dK1 <- ifelse(j < j.cols, array.N[i,j+1,t] / Kentry.interp.pois[i,j+1], NA)
            SW.dK1 <- ifelse(i < i.rows & j > 1, array.N[i+1,j-1,t] / Kentry.interp.pois[i+1,j-1], NA)
            S.dK1 <- ifelse(i < i.rows, array.N[i+1,j,t] / Kentry.interp.pois[i+1,j], NA)
            SE.dK1 <- ifelse(i < i.rows & j < j.cols, array.N[i+1,j+1,t] / Kentry.interp.pois[i+1,j+1], NA)
            
            focal.dK <- ifelse(focal.dK1 == 0 | is.na(focal.dK1) == T, 0, focal.dK1)
            NW.dK <- ifelse(length(NW.dK1) == 0 | is.na(NW.dK1) == T, 0, NW.dK1)
            N.dK <- ifelse(length(N.dK1) == 0 | is.na(N.dK1) == T, 0, N.dK1)
            NE.dK <- ifelse(length(NE.dK1) == 0 | is.na(NE.dK1) == T, 0, NE.dK1)
            W.dK <- ifelse(length(W.dK1) == 0 | is.na(W.dK1) == T, 0, W.dK1)
            E.dK <- ifelse(length(E.dK1) == 0 | is.na(E.dK1) == T, 0, E.dK1)
            SW.dK <- ifelse(length(SW.dK1) == 0 | is.na(SW.dK1) == T, 0, SW.dK1)
            S.dK <- ifelse(length(S.dK1) == 0 | is.na(S.dK1) == T, 0, S.dK1)
            SE.dK <- ifelse(length(SE.dK1) == 0 | is.na(SE.dK1) == T, 0, SE.dK1)
            
            # direction indices initialised as NA
            NWdir <- Ndir <- NEdir <- Wdir <- Edir <- SWdir <- Sdir <- SEdir <- NA
            
            # immigration into focal cell
            if (is.na(NW.RK) == F & NW.RK > 1 & NW.dK >= runif(1, min=NK.emig.min, max=NK.emig.max)) {
              NW.I <- (pI.func(NW.RK) * stoch.beta.func(pmov.mean, pmov.sd) * array.N[i-1,j-1,t] * (rugmovmod.func(rugentry.interp.array[i-1,j-1,t])))
              array.N[i,j,t] <- sum(c(array.N[i,j,t], NW.I), na.rm=T)
              array.N[i-1,j-1,t] <- sum(c(array.N[i-1,j-1,t], -NW.I), na.rm=T)
              NWdir <- ifelse(is.na(NW.RK) == F & NW.RK > 1, NW.I, NA)}
            if (is.na(N.RK) == F & N.RK > 1 & N.dK >= runif(1, min=NK.emig.min, max=NK.emig.max)) {
              N.I <- (pI.func(N.RK) * stoch.beta.func(pmov.mean, pmov.sd) * array.N[i-1,j,t] * (rugmovmod.func(rugentry.interp.array[i-1,j,t])))
              array.N[i,j,t] <- sum(c(array.N[i,j,t], N.I), na.rm=T)
              array.N[i-1,j,t] <- sum(c(array.N[i-1,j,t], -N.I), na.rm=T)
              Ndir <- ifelse(is.na(N.RK) == F & N.RK > 1, N.I, NA)}
            if (is.na(NE.RK) == F & NE.RK > 1 & NE.dK >= runif(1, min=NK.emig.min, max=NK.emig.max)) {
              NE.I <- (pI.func(NE.RK) * stoch.beta.func(pmov.mean, pmov.sd) * array.N[i-1,j+1,t] * (rugmovmod.func(rugentry.interp.array[i-1,j+1,t])))
              array.N[i,j,t] <- sum(c(array.N[i,j,t], NE.I), na.rm=T)
              array.N[i-1,j+1,t] <- sum(c(array.N[i-1,j+1,t], -NE.I), na.rm=T)
              NEdir <- ifelse(is.na(NE.RK) == F & NE.RK > 1, NE.I, NA)}
            if (is.na(W.RK) == F & W.RK > 1 & W.dK >= runif(1, min=NK.emig.min, max=NK.emig.max)) {
              W.I <- (pI.func(W.RK) * stoch.beta.func(pmov.mean, pmov.sd) * array.N[i,j-1,t] * (rugmovmod.func(rugentry.interp.array[i,j-1,t])))
              array.N[i,j,t] <- sum(c(array.N[i,j,t], W.I), na.rm=T)
              array.N[i,j-1,t] <- sum(c(array.N[i,j-1,t], -W.I), na.rm=T)
              Wdir <- ifelse(is.na(W.RK) == F & W.RK > 1, W.I, NA)}
            if (is.na(E.RK) == F & E.RK > 1 & E.dK >= runif(1, min=NK.emig.min, max=NK.emig.max)) {
              E.I <- (pI.func(E.RK) * stoch.beta.func(pmov.mean, pmov.sd) * array.N[i,j+1,t] * (rugmovmod.func(rugentry.interp.array[i,j+1,t])))
              array.N[i,j,t] <- sum(c(array.N[i,j,t], E.I), na.rm=T)
              array.N[i,j+1,t] <- sum(c(array.N[i,j+1,t], -E.I), na.rm=T)
              Edir <- ifelse(is.na(E.RK) == F & E.RK > 1, E.I, NA)}
            if (is.na(SW.RK) == F & SW.RK > 1 & SW.dK >= runif(1, min=NK.emig.min, max=NK.emig.max)) {
              SW.I <- (pI.func(SW.RK) * stoch.beta.func(pmov.mean, pmov.sd) * array.N[i+1,j-1,t] * (rugmovmod.func(rugentry.interp.array[i+1,j-1,t])))
              array.N[i,j,t] <- sum(c(array.N[i,j,t], SW.I), na.rm=T)
              array.N[i+1,j-1,t] <- sum(c(array.N[i+1,j-1,t], -SW.I), na.rm=T)
              SWdir <- ifelse(is.na(SW.RK) == F & SW.RK > 1, SW.I, NA)}
            if (is.na(S.RK) == F & S.RK > 1 & S.dK >= runif(1, min=NK.emig.min, max=NK.emig.max)) {
              S.I <- (pI.func(S.RK) * stoch.beta.func(pmov.mean, pmov.sd) * array.N[i+1,j,t] * (rugmovmod.func(rugentry.interp.array[i+1,j,t])))
              array.N[i,j,t] <- sum(c(array.N[i,j,t], S.I), na.rm=T)
              array.N[i+1,j,t] <- sum(c(array.N[i+1,j,t], -S.I), na.rm=T)
              Sdir <- ifelse(is.na(S.RK) == F & S.RK > 1, S.I, NA)}
            if (is.na(SE.RK) == F & SE.RK > 1 & SE.dK >= runif(1, min=NK.emig.min, max=NK.emig.max)) {
              SE.I <- (pI.func(SE.RK) * stoch.beta.func(pmov.mean, pmov.sd) * array.N[i+1,j+1,t] * (rugmovmod.func(rugentry.interp.array[i+1,j+1,t])))
              array.N[i,j,t] <- sum(c(array.N[i,j,t], SE.I), na.rm=T)
              array.N[i+1,j+1,t] <- sum(c(array.N[i+1,j+1,t], -SE.I), na.rm=T)
              SEdir <- ifelse(is.na(SE.RK) == F & SE.RK > 1, SE.I, NA)}
            
            # direction of dominant influx per time layer
            I.vec <- c(NWdir, Ndir, NEdir, Wdir, Edir, SWdir, Sdir, SEdir)
            dir.array[i,j,t] <- ifelse((length(which(I.vec > 1))) > 0, (dir.vec[max(which(I.vec > 1))]), NA)
            
            # emigration out of focal cell
            if ((ifelse(focal.dK >= runif(1, min=NK.emig.min, max=NK.emig.max), 1, 0)) == 1) {
              if (is.na(NW.RK) == F & NW.RK <= 1) {
                NW.E <- (pE.func(NW.RK) * stoch.beta.func(pmov.mean, pmov.sd) * array.N[i,j,t] * (rugmovmod.func(rugentry.interp.array[i,j,t])))
                array.N[i,j,t] <- sum(c(array.N[i,j,t], ifelse(NW.E < 0, 0, -NW.E)), na.rm=T)
                array.N[i-1,j-1,t] <- sum(c(array.N[i-1,j-1,t], ifelse(NW.E < 0, 0, NW.E)), na.rm=T)}
              if (is.na(N.RK) == F & N.RK <= 1) {
                N.E <- (pE.func(N.RK) * stoch.beta.func(pmov.mean, pmov.sd) * array.N[i,j,t] * (rugmovmod.func(rugentry.interp.array[i,j,t]))) - NW.E
                array.N[i,j,t] <- sum(c(array.N[i,j,t], ifelse(N.E < 0, 0, -N.E)), na.rm=T)
                array.N[i-1,j,t] <- sum(c(array.N[i-1,j,t], ifelse(N.E < 0, 0, N.E)), na.rm=T)}
              if (is.na(NE.RK) == F & NE.RK <= 1) {
                NE.E <- (pE.func(NE.RK) * stoch.beta.func(pmov.mean, pmov.sd) * array.N[i,j,t] * (rugmovmod.func(rugentry.interp.array[i,j,t]))) - NW.E - N.E
                array.N[i,j,t] <- sum(c(array.N[i,j,t], ifelse(NE.E < 0, 0, -NE.E)), na.rm=T)
                array.N[i-1,j+1,t] <- sum(c(array.N[i-1,j+1,t], ifelse(NE.E < 0, 0, NE.E)), na.rm=T)}
              if (is.na(W.RK) == F & W.RK <= 1) {
                W.E <- (pE.func(W.RK) * stoch.beta.func(pmov.mean, pmov.sd) * array.N[i,j,t] * (rugmovmod.func(rugentry.interp.array[i,j,t]))) - NW.E - N.E - NE.E
                array.N[i,j,t] <- sum(c(array.N[i,j,t], ifelse(W.E < 0, 0, -W.E)), na.rm=T)
                array.N[i,j-1,t] <- sum(c(array.N[i,j-1,t], ifelse(W.E < 0, 0, W.E)), na.rm=T)}
              if (is.na(E.RK) == F & E.RK <= 1) {
                E.E <- (pE.func(E.RK) * stoch.beta.func(pmov.mean, pmov.sd) * array.N[i,j,t] * (rugmovmod.func(rugentry.interp.array[i,j,t]))) - NW.E - N.E - NE.E - W.E
                array.N[i,j,t] <- sum(c(array.N[i,j,t], ifelse(E.E < 0, 0, -E.E)), na.rm=T)
                array.N[i,j+1,t] <- sum(c(array.N[i,j+1,t], ifelse(E.E < 0, 0, E.E)), na.rm=T)}
              if (is.na(SW.RK) == F & SW.RK <= 1) {
                SW.E <- (pE.func(SW.RK) * stoch.beta.func(pmov.mean, pmov.sd) * array.N[i,j,t] * (rugmovmod.func(rugentry.interp.array[i,j,t]))) - NW.E - N.E - NE.E - W.E - E.E
                array.N[i,j,t] <- sum(c(array.N[i,j,t], ifelse(SW.E < 0, 0, -SW.E)), na.rm=T)
                array.N[i+1,j-1,t] <- sum(c(array.N[i+1,j-1,t], ifelse(SW.E < 0, 0, SW.E)), na.rm=T)}
              if (is.na(S.RK) == F & S.RK <= 1) {
                S.E <- (pE.func(S.RK) * stoch.beta.func(pmov.mean, pmov.sd) * array.N[i,j,t] * (rugmovmod.func(rugentry.interp.array[i,j,t]))) - NW.E - N.E - NE.E - W.E - E.E - SW.E
                array.N[i,j,t] <- sum(c(array.N[i,j,t], ifelse(S.E < 0, 0, -S.E)), na.rm=T)
                array.N[i+1,j,t] <- sum(c(array.N[i+1,j,t], ifelse(S.E < 0, 0, S.E)), na.rm=T)}
              if (is.na(SE.RK) == F & SE.RK <= 1) {
                SE.E <- (pE.func(SE.RK) * stoch.beta.func(pmov.mean, pmov.sd) * array.N[i,j,t] * (rugmovmod.func(rugentry.interp.array[i,j,t]))) - NW.E - N.E - NE.E - W.E - E.E - SW.E - S.E
                array.N[i,j,t] <- sum(c(array.N[i,j,t], ifelse(SE.E < 0, 0, -SE.E)), na.rm=T)
                array.N[i+1,j+1,t] <- sum(c(array.N[i+1,j+1,t], ifelse(SE.E < 0, 0, SE.E)), na.rm=T)}
            }
            
            # reset emigration values to zero for next run
            NW.E <- N.E <- NE.E <- W.E <- E.E <- SW.E <- S.E <- SE.E <- 0
            
            # long-range dispersal
            disp.prob <- (exp(-D.vec/(ifelse(Kentry.interp.lrel.array[i, j, t] >= log10(6.706985), D.max.lo, 10^(disp.npp.max.int + as.numeric(disp.npp.slope) * Kentry.interp.lrel.array[i, j, t])))))
            if (is.na(disp.prob[1])==F) {
              disp.prob.ran <- disp.prob} 
            if (is.na(disp.prob[1])==T) {
              disp.prob.ran <- ((Pr.Dmx.up+Pr.Dmx.lo)/2)}
            
            cellD.move <- ldp.mult * ((sample(cellD, size=1, replace=F, prob=disp.prob.ran)) * (2*focal.dK)) # multiply probability upwards if closer to focal cell K
            
            # condition on distance 2 water, where
            # P(reach) = 1 – ((D/max(D))^J); D = distance to water; higher J = more difficult reach gridcell; max(D) = max distance to water
            if (is.na(d2w.mat[i,j]) == F & cellD.move < d2w.mat[i,j]) {
              P.reach <- 1 - ((cellD.move / d2w.mat[i,j])^(watmovresist))
              reach.trial <- rbinom(1, 1, prob=P.reach)
            }
            reach.succ <- ifelse(is.na(d2w.mat[i,j]) == F, reach.trial, 1)
            
            if ((round(cellD.move)) > 0 & reach.succ == 1) {
              dx <- sample(c(-1,1), 1) * rpois(1,(round(cellD.move)))
              dy <- sample(c(-1,1), 1) * rpois(1,(round(cellD.move)))
              if ((i+dy) > 0 & (i+dy) <= i.rows & (j+dx) > 0 & (j+dx) <= j.cols) {
                if (length(which(is.na(array.N[(i+sign(dy)):(i+dy), (j+sign(dx)):(j+dx), t]) == T)) == 0) { # if there is an NA cell anywhere in block between [i,j] & [i+dy,j+dx], don't let dispersal happen
                  N.disp <- round(array.N[i, j, t] * stoch.beta.func(pmov.mean/10, pmov.sd/10) * (rugmovmod.func(rugentry.interp.array[i,j,t])), 0) # number dispersing to new cell
                  array.N[(i+dy), (j+dx), t] <- array.N[(i+dy), (j+dx), t] + N.disp
                  array.N[i, j, t] <- ifelse((array.N[i, j, t] - N.disp) < 0, 0, (array.N[i, j, t] - N.disp))
                } # end if
              } # end if
            } # end if
            
          } # end i loop
        } # end j loop
        
        # remove negative values
        array.N[,,t] <- ifelse(array.N[,,t] < 0, 0, round(array.N[,,t], 0))
        
        # apply dynamical model after movements from previous step
        if (rmax.sel == "low") {
          array.N[,,t+1] <- Nproj.func(Nt=array.N[,,t], rm=rnorm(1, rm.max.NEE, rm.max.NEE.sd), K=Kentry.interp.pois)
        }
        if (rmax.sel == "high") {
          array.N[,,t+1] <- Nproj.func(Nt=array.N[,,t], rm=log((exp(10^(rnorm(1, a.Cole, a.sd.Cole) + rnorm(1, b.Cole, b.sd.Cole) * log10(alpha.Ab))))^gen.l), K=Kentry.interp.pois)
        }
        
        # apply catastrophe; resampling cat.pr (binomial) & M.cat mortality (beta)
        if (catastrophe == "yes") {
          mat.noNA.g0 <- which(is.na(array.N[,,t+1]) != T & array.N[,,t+1] > 0, arr.ind=F)
          mat.noNA.g0.sel <- rbinom(length(mat.noNA.g0), 1, cat.pr) # spatially random
          
          if (spatial.cluster == "yes" & length(mat.noNA.g0) > 1) {
            mat.noNA.g0.arr <- which(is.na(array.N[,,t+1]) != T & array.N[,,t+1] > 0, arr.ind=T)
            rows.ng0 <- sort(unique(mat.noNA.g0.arr[,1]))
            cols.ng0 <- sort(unique(mat.noNA.g0.arr[,2]))
            lrowcol.mn <- round(mean(length(rows.ng0),length(cols.ng0)), 0)
            mat.N.it <- array.N[,,t+1]
            Y <- rThomas(round(kappaP.func(length(which(array.N[,,t+1] > 0, arr.ind=F)))*lrowcol.mn, 0), rpp.scale, round(rpp.mu.mult*lrowcol.mn,0), drop=T, nsim = 1)
            Yrows.calc <- round(Y$y * length(rows.ng0), 0)
            Yrows <- Yrows.calc + min(rows.ng0) - 1
            Yrows <- ifelse((Yrows.calc + min(rows.ng0) - 1) == 0, 1, (Yrows.calc + min(rows.ng0) - 1)) 
            Ycols.calc <- round(Y$x * length(cols.ng0), 0)
            Ycols <- Yrows.calc + min(cols.ng0) - 1
            Ycols <- ifelse((Ycols.calc + min(cols.ng0) - 1) == 0, 1, (Ycols.calc + min(cols.ng0) - 1)) 
            Yrowscols <- unique(cbind(Yrows,Ycols))
            mat.N.it[Yrowscols] <- round((stoch.n.beta.func(length(mat.N.it[Yrowscols]), M.cat, M.cat.sd)) * (mat.N.it[Yrowscols]), 0)
            mat.noNA.g0.spat <- which(mat.N.it %in% na.omit(mat.N.it[Yrowscols])[na.omit(mat.N.it[Yrowscols]) > 0])
            cat.pr.est[t] <- length(mat.noNA.g0.spat)/length(mat.noNA.g0)
            cat.pr.mx <- rnorm(1, cat.pr, cat.pr*SD.prop.xbar)
            if (cat.pr.est[t] > cat.pr.mx & length(mat.noNA.g0.spat) > 0 & length(mat.noNA.g0) > 1 & spatial.cluster=="yes") {
              mat.noNA.g0.spat <- sample(mat.noNA.g0.spat, round(cat.pr.mx * length(mat.noNA.g0)), replace=F)
              cat.pr.est[t] <- length(mat.noNA.g0.spat)/dim(mat.noNA.g0.arr)[1]}
            if (length(mat.noNA.g0.spat) > 0 & length(mat.noNA.g0) > 1 & spatial.cluster=="yes") {
              array.N[,,t+1][mat.noNA.g0.spat] <- round((stoch.n.beta.func(length(mat.noNA.g0.spat), M.cat, M.cat.sd)) * ((array.N[,,t+1])[mat.noNA.g0.spat]), 0)}
            if (length(mat.noNA.g0.spat) == 0 & length(mat.noNA.g0) > 1 & spatial.cluster=="yes") {
              mat.noNA.g0.sel <- rbinom(length(mat.noNA.g0), 1, cat.pr) # spatially random
              array.N[,,t+1][mat.noNA.g0[which(mat.noNA.g0.sel > 0)]] <- round((stoch.n.beta.func(length(mat.noNA.g0[which(mat.noNA.g0.sel > 0)]), M.cat, M.cat.sd)) * ((array.N[,,t+1])[mat.noNA.g0[which(mat.noNA.g0.sel > 0)]]), 0)
              cat.pr.est[t] <- length(mat.noNA.g0[which(mat.noNA.g0.sel > 0)])/dim(mat.noNA.g0.arr)[1]}
          }
          if (length(mat.noNA.g0[-mat.noNA.g0.sel]) > 0 & length(mat.noNA.g0) > 1 & spatial.cluster=="no") {
            array.N[,,t+1][mat.noNA.g0[which(mat.noNA.g0.sel > 0)]] <- round((stoch.n.beta.func(length(mat.noNA.g0[which(mat.noNA.g0.sel > 0)]), M.cat, M.cat.sd)) * ((array.N[,,t+1])[mat.noNA.g0[which(mat.noNA.g0.sel > 0)]]), 0)}
        }
        
        # apply NA layer post-growth
        array.N[,,t+1] <- ifelse(is.na(NA.array[,,t+1]) == T, NA, array.N[,,t+1])
        
        # apply additional mortality for cells with N < MVP.thresh
        if (small.pop.ext == "yes") {
          ltMVPthresh.sub <- which(array.N[,,t+1] < MVP.thresh) # which cells in array.N[,,t] < MVP.thresh
          ltMVPred.vec <- stoch.n.beta.func(length(ltMVPthresh.sub), (1-ltMVP.red), 0.05*ltMVP.red) # assume 5% SD
          array.N[,,t+1][ltMVPthresh.sub] <- round(array.N[,,t+1][ltMVPthresh.sub] * ltMVPred.vec, 0)
        }
        
        # if extinct, restart colonisation at origin
        if ((length(which(is.na(array.N[,,t+1]) == T | array.N[,,t+1] < 1) == T)) == i.rows*j.cols) {
          array.N[start.row1,start.col1,t+1] <- array.N[start.row1,start.col1,1]
        }
        
        N.vec[t+1] <- sum(array.N[,,t+1], na.rm=T)
        NG.N.vec[t+1] <- sum((array.N[1:19, 40:69, t+1]), na.rm=T) + sum((array.N[1:24, 70:87, t+1]), na.rm=T)
        AUS.N.vec[t+1] <- sum(array.N[,, t+1], na.rm = T) - NG.N.vec[t+1]
        
        CEN.N.cell <- rep(NA,dim(CEN.crds)[1])
        for (i in 1:dim(CEN.crds)[1]) {
          CEN.N.cell[i] <- array.N[CEN.crds[i,1], CEN.crds[i,2], t+1]
        }
        TE.N.cell <- rep(NA,dim(TE.crds)[1])
        for (i in 1:dim(TE.crds)[1]) {
          TE.N.cell[i] <- array.N[TE.crds[i,1], TE.crds[i,2], t+1]
        }
        E.N.cell <- rep(NA,dim(E.crds)[1])
        for (i in 1:dim(E.crds)[1]) {
          E.N.cell[i] <- array.N[E.crds[i,1], E.crds[i,2], t+1]
        }
        NGH.N.cell <- rep(NA,dim(NGH.crds)[1])
        for (i in 1:dim(NGH.crds)[1]) {
          NGH.N.cell[i] <- array.N[NGH.crds[i,1], NGH.crds[i,2], t+1]
        }
        
        CEN.N.vec[t+1] <- sum(CEN.N.cell, na.rm=T)
        TE.N.vec[t+1] <- sum(TE.N.cell, na.rm=T)
        E.N.vec[t+1] <- sum(E.N.cell, na.rm=T)
        NGH.N.vec[t+1] <- sum(NGH.N.cell, na.rm=T)
        
        CEN.D.vec[t+1] <- CEN.N.vec[t+1]/(length(which(is.na(CEN.N.cell)==F)) * cell.area)
        TE.D.vec[t+1] <- TE.N.vec[t+1]/(length(which(is.na(TE.N.cell)==F)) * cell.area)
        E.D.vec[t+1] <- E.N.vec[t+1]/(length(which(is.na(E.N.cell)==F)) * cell.area)
        NGH.D.vec[t+1] <- NGH.N.vec[t+1]/(length(which(is.na(NGH.N.cell)==F)) * cell.area)
        
        # bioregions
        PILGASMUR.N.cell <- rep(NA,dim(PILGASMUR.crds)[1])
        for (i in 1:dim(PILGASMUR.crds)[1]) {
          PILGASMUR.N.cell[i] <- array.N[PILGASMUR.crds[i,1], PILGASMUR.crds[i,2], t+1]
        }
        WAR.N.cell <- rep(NA,dim(WAR.crds)[1])
        for (i in 1:dim(WAR.crds)[1]) {
          WAR.N.cell[i] <- array.N[WAR.crds[i,1], WAR.crds[i,2], t+1]
        }
        CER.N.cell <- rep(NA,dim(CER.crds)[1])
        for (i in 1:dim(CER.crds)[1]) {
          CER.N.cell[i] <- array.N[CER.crds[i,1], CER.crds[i,2], t+1]
        }
        OVP.N.cell <- rep(NA,dim(OVP.crds)[1])
        for (i in 1:dim(OVP.crds)[1]) {
          OVP.N.cell[i] <- array.N[OVP.crds[i,1], OVP.crds[i,2], t+1]
        }
        ARP.N.cell <- rep(NA,dim(ARP.crds)[1])
        for (i in 1:dim(ARP.crds)[1]) {
          ARP.N.cell[i] <- array.N[ARP.crds[i,1], ARP.crds[i,2], t+1]
        }
        GUPEIU.N.cell <- rep(NA,dim(GUPEIU.crds)[1])
        for (i in 1:dim(GUPEIU.crds)[1]) {
          GUPEIU.N.cell[i] <- array.N[GUPEIU.crds[i,1], GUPEIU.crds[i,2], t+1]
        }
        BNSNS.N.cell <- rep(NA,dim(BNSNS.crds)[1])
        for (i in 1:dim(BNSNS.crds)[1]) {
          BNSNS.N.cell[i] <- array.N[BNSNS.crds[i,1], BNSNS.crds[i,2], t+1]
        }
        MDDBHC.N.cell <- rep(NA,dim(MDDBHC.crds)[1])
        for (i in 1:dim(MDDBHC.crds)[1]) {
          MDDBHC.N.cell[i] <- array.N[MDDBHC.crds[i,1], MDDBHC.crds[i,2], t+1]
        }
        TCH.N.cell <- rep(NA,dim(TCH.crds)[1])
        for (i in 1:dim(TCH.crds)[1]) {
          TCH.N.cell[i] <- array.N[TCH.crds[i,1], TCH.crds[i,2], t+1]
        }
        
        PILGASMUR.N.vec[t+1] <- sum(PILGASMUR.N.cell, na.rm=T)
        WAR.N.vec[t+1] <- sum(WAR.N.cell, na.rm=T)
        CER.N.vec[t+1] <- sum(CER.N.cell, na.rm=T)
        OVP.N.vec[t+1] <- sum(OVP.N.cell, na.rm=T)
        ARP.N.vec[t+1] <- sum(ARP.N.cell, na.rm=T)
        GUPEIU.N.vec[t+1] <- sum(GUPEIU.N.cell, na.rm=T)
        BNSNS.N.vec[t+1] <- sum(BNSNS.N.cell, na.rm=T)
        MDDBHC.N.vec[t+1] <- sum(MDDBHC.N.cell, na.rm=T)
        TCH.N.vec[t+1] <- sum(TCH.N.cell, na.rm=T)
        
        PILGASMUR.D.vec[t+1] <- PILGASMUR.N.vec[t+1]/(length(which(is.na(PILGASMUR.N.cell)==F)) * cell.area)
        WAR.D.vec[t+1] <- WAR.N.vec[t+1]/(length(which(is.na(WAR.N.cell)==F)) * cell.area)
        CER.D.vec[t+1] <- CER.N.vec[t+1]/(length(which(is.na(CER.N.cell)==F)) * cell.area)
        OVP.D.vec[t+1] <- OVP.N.vec[t+1]/(length(which(is.na(OVP.N.cell)==F)) * cell.area)
        ARP.D.vec[t+1] <- ARP.N.vec[t+1]/(length(which(is.na(ARP.N.cell)==F)) * cell.area)
        GUPEIU.D.vec[t+1] <- GUPEIU.N.vec[t+1]/(length(which(is.na(GUPEIU.N.cell)==F)) * cell.area)
        BNSNS.D.vec[t+1] <- BNSNS.N.vec[t+1]/(length(which(is.na(BNSNS.N.cell)==F)) * cell.area)
        MDDBHC.D.vec[t+1] <- MDDBHC.N.vec[t+1]/(length(which(is.na(MDDBHC.N.cell)==F)) * cell.area)
        TCH.D.vec[t+1] <- TCH.N.vec[t+1]/(length(which(is.na(TCH.N.cell)==F)) * cell.area)
        
        poparea.vec[t+1] <- length(which(array.N[,,t+1] > 0)) * cell.area/1000
        pc.complete[t+1] <- (round(length(which(array.N[,,t+1] > 0, arr.ind=F)) / length(which(is.na(array.N[,,t+1]) != T, arr.ind=F)) * 99, 2))
        
      } # end t loop
      
      N.mat[m,] <- N.vec
      NG.N.mat[m,] <- NG.N.vec
      AUS.N.mat[m,] <- AUS.N.vec
      CEN.N.mat[m,] <- CEN.N.vec
      TE.N.mat[m,] <- TE.N.vec
      E.N.mat[m,] <- E.N.vec
      NGH.N.mat[m,] <- NGH.N.vec
      PILGASMUR.N.mat[m,] <- PILGASMUR.N.vec
      WAR.N.mat[m,] <- WAR.N.vec
      CER.N.mat[m,] <- CER.N.vec
      OVP.N.mat[m,] <- OVP.N.vec
      ARP.N.mat[m,] <- ARP.N.vec
      GUPEIU.N.mat[m,] <- GUPEIU.N.vec
      BNSNS.N.mat[m,] <- BNSNS.N.vec
      MDDBHC.N.mat[m,] <- MDDBHC.N.vec
      TCH.N.mat[m,] <- TCH.N.vec
      
      CEN.D.mat[m,] <- CEN.D.vec
      TE.D.mat[m,] <- TE.D.vec
      E.D.mat[m,] <- E.D.vec
      NGH.D.mat[m,] <- NGH.D.vec
      PILGASMUR.D.mat[m,] <- PILGASMUR.D.vec
      WAR.D.mat[m,] <- WAR.D.vec
      CER.D.mat[m,] <- CER.D.vec
      OVP.D.mat[m,] <- OVP.D.vec
      ARP.D.mat[m,] <- ARP.D.vec
      GUPEIU.D.mat[m,] <- GUPEIU.D.vec
      BNSNS.D.mat[m,] <- BNSNS.D.vec
      MDDBHC.D.mat[m,] <- MDDBHC.D.vec
      TCH.D.mat[m,] <- TCH.D.vec
      
      ## array-based summaries
      ## CV of N by cell
      cv.N.arr <- MIS3.cv.N <- MIS2.cv.N <- LGM.cv.N <- LGMc.cv.N <- LGMu.cv.N <- HOL.cv.N <- array(data=NA, dim=c(dim(array.N)[1:2], reps))
      array.N.nzero <- ifelse(array.N < 1, NA, array.N)
      cv.N.arr[,,m] <- apply(array.N.nzero, c(1,2), cv, na.rm=T) # CV of N over all time by cell
      #image(rot.mat(cv.N.arr[,,m]), col=rev(grey(1:100/100)))
      
      ## CV by epoch
      MIS3.array.N <- array.N[,, MIS3.t.st:MIS3.t.en]
      MIS3.array.N.nzero <- ifelse(MIS3.array.N < 1, NA, MIS3.array.N)
      MIS3.cv.N[,,m] <- apply(MIS3.array.N.nzero, c(1,2), cv, na.rm=T) # CV of N over all time by cell
      
      MIS2.array.N <- array.N[,, MIS2.t.st:MIS2.t.en] # 75
      MIS2.array.N.nzero <- ifelse(MIS2.array.N < 1, NA, MIS2.array.N)
      MIS2.cv.N[,,m] <- apply(MIS2.array.N.nzero, c(1,2), cv, na.rm=T) # CV of N over all time by cell
      
      LGM.array.N <- array.N[,, LGM.t.st:LGM.t.en] # 75
      LGM.array.N.nzero <- ifelse(LGM.array.N < 1, NA, LGM.array.N)
      LGM.cv.N[,,m] <- apply(LGM.array.N.nzero, c(1,2), cv, na.rm=T) # CV of N over all time by cell
      
      LGMc.array.N <- array.N[,, LGMc.t.st:LGMc.t.en] # 75
      LGMc.array.N.nzero <- ifelse(LGMc.array.N < 1, NA, LGMc.array.N)
      LGMc.cv.N[,,m] <- apply(LGMc.array.N.nzero, c(1,2), cv, na.rm=T) # CV of N over all time by cell
      
      LGMu.array.N <- array.N[,, LGMu.t.st:LGMu.t.en] # 75
      LGMu.array.N.nzero <- ifelse(LGMu.array.N < 1, NA, LGMu.array.N)
      LGMu.cv.N[,,m] <- apply(LGMu.array.N.nzero, c(1,2), cv, na.rm=T) # CV of N over all time by cell
      
      HOL.array.N <- array.N[,, HOL.t.st:HOL.t.en] # 75
      HOL.array.N.nzero <- ifelse(HOL.array.N < 1, NA, HOL.array.N)
      HOL.cv.N[,,m] <- apply(HOL.array.N.nzero, c(1,2), cv, na.rm=T) # CV of N over all time by cell
      
      # calculate modal direction array
      MIS3.dir.array <- dir.array[,,MIS3.t.st:MIS3.t.en]
      MIS2.dir.array <- dir.array[,,MIS2.t.st:MIS2.t.en]
      LGM.dir.array <- dir.array[,,LGM.t.st:LGM.t.en]
      LGMc.dir.array <- dir.array[,,LGMc.t.st:LGMc.t.en]
      LGMu.dir.array <- dir.array[,,LGMu.t.st:LGMu.t.en]
      HOL.dir.array <- dir.array[,,HOL.t.st:HOL.t.en]
      L26.dir.array <- dir.array[,,L26.t.st:L26.t.en]
      L25.dir.array <- dir.array[,,L25.t.st:L25.t.en]
      L24.dir.array <- dir.array[,,L24.t.st:L24.t.en]
      L23.dir.array <- dir.array[,,L23.t.st:L23.t.en]
      L22.dir.array <- dir.array[,,L22.t.st:L22.t.en]
      L21.dir.array <- dir.array[,,L21.t.st:L21.t.en]
      L20.dir.array <- dir.array[,,L20.t.st:L20.t.en]
      L19.dir.array <- dir.array[,,L19.t.st:L19.t.en]
      L18.dir.array <- dir.array[,,L18.t.st:L18.t.en]
      L17.dir.array <- dir.array[,,L17.t.st:L17.t.en]
      L16.dir.array <- dir.array[,,L16.t.st:L16.t.en]
      L15.dir.array <- dir.array[,,L15.t.st:L15.t.en]
      L14.dir.array <- dir.array[,,L14.t.st:L14.t.en]
      L13.dir.array <- dir.array[,,L13.t.st:L13.t.en]
      L12.dir.array <- dir.array[,,L12.t.st:L12.t.en]
      L11.dir.array <- dir.array[,,L11.t.st:L11.t.en]
      L10.dir.array <- dir.array[,,L10.t.st:L10.t.en]
      L09.dir.array <- dir.array[,,L09.t.st:L09.t.en]
      L08.dir.array <- dir.array[,,L08.t.st:L08.t.en]
      L07.dir.array <- dir.array[,,L07.t.st:L07.t.en]
      
      L2624.dir.array <- dir.array[,,L26.t.st:L25.t.en]
      L2422.dir.array <- dir.array[,,L24.t.st:L23.t.en]
      L2220.dir.array <- dir.array[,,L22.t.st:L21.t.en]
      L2018.dir.array <- dir.array[,,L20.t.st:L19.t.en]
      L1816.dir.array <- dir.array[,,L18.t.st:L17.t.en]
      L1614.dir.array <- dir.array[,,L16.t.st:L15.t.en]
      L1412.dir.array <- dir.array[,,L14.t.st:L13.t.en]
      L1210.dir.array <- dir.array[,,L12.t.st:L11.t.en]
      L1008.dir.array <- dir.array[,,L10.t.st:L09.t.en]
      L0806.dir.array <- dir.array[,,L08.t.st:L07.t.en]
      
      ALL.mode.imm.dir <- MIS3.mode.imm.dir <- MIS2.mode.imm.dir <- LGM.mode.imm.dir <- LGMc.mode.imm.dir <- LGMu.mode.imm.dir <- HOL.mode.imm.dir <-
        L26.mode.imm.dir <- L25.mode.imm.dir <- L24.mode.imm.dir <- L23.mode.imm.dir <- L22.mode.imm.dir <- L21.mode.imm.dir <- L20.mode.imm.dir <-
        L19.mode.imm.dir <- L18.mode.imm.dir <- L17.mode.imm.dir <- L16.mode.imm.dir <- L15.mode.imm.dir <- L14.mode.imm.dir <- L13.mode.imm.dir <-
        L12.mode.imm.dir <- L11.mode.imm.dir <- L10.mode.imm.dir <- L09.mode.imm.dir <- L08.mode.imm.dir <- L07.mode.imm.dir <-
        L2624.mode.imm.dir <- L2422.mode.imm.dir <- L2220.mode.imm.dir <- L2018.mode.imm.dir <- L1816.mode.imm.dir <- L1614.mode.imm.dir <- 
        L1412.mode.imm.dir <- L1210.mode.imm.dir <- L1008.mode.imm.dir <- L0806.mode.imm.dir <-
        matrix(data = NA, nrow = i.rows, ncol = j.cols)
      for (i in 1:i.rows) { # i rows
        for (j in 1:j.cols) { # j columns
          ALL.mode.imm.dir[i,j] <- ifelse((length(which(is.na(dir.array[i,j,]) == F)) > 0), (attr(sort(table(dir.array[i, j, which(is.na(dir.array[i,j,]) == F)]), decreasing=T), 'names')[1]), NA)
          MIS3.mode.imm.dir[i,j] <- ifelse((length(which(is.na(MIS3.dir.array[i,j,]) == F)) > 0), (attr(sort(table(MIS3.dir.array[i, j, which(is.na(MIS3.dir.array[i,j,]) == F)]), decreasing=T), 'names')[1]), NA)
          MIS2.mode.imm.dir[i,j] <- ifelse((length(which(is.na(MIS2.dir.array[i,j,]) == F)) > 0), (attr(sort(table(MIS2.dir.array[i, j, which(is.na(MIS2.dir.array[i,j,]) == F)]), decreasing=T), 'names')[1]), NA)
          LGM.mode.imm.dir[i,j] <- ifelse((length(which(is.na(LGM.dir.array[i,j,]) == F)) > 0), (attr(sort(table(LGM.dir.array[i, j, which(is.na(LGM.dir.array[i,j,]) == F)]), decreasing=T), 'names')[1]), NA)
          LGMc.mode.imm.dir[i,j] <- ifelse((length(which(is.na(LGMc.dir.array[i,j,]) == F)) > 0), (attr(sort(table(LGMc.dir.array[i, j, which(is.na(LGMc.dir.array[i,j,]) == F)]), decreasing=T), 'names')[1]), NA)
          LGMu.mode.imm.dir[i,j] <- ifelse((length(which(is.na(LGMu.dir.array[i,j,]) == F)) > 0), (attr(sort(table(LGMu.dir.array[i, j, which(is.na(LGMu.dir.array[i,j,]) == F)]), decreasing=T), 'names')[1]), NA)
          HOL.mode.imm.dir[i,j] <- ifelse((length(which(is.na(HOL.dir.array[i,j,]) == F)) > 0), (attr(sort(table(HOL.dir.array[i, j, which(is.na(HOL.dir.array[i,j,]) == F)]), decreasing=T), 'names')[1]), NA)
          L26.mode.imm.dir[i,j] <- ifelse((length(which(is.na(L26.dir.array[i,j,]) == F)) > 0), (attr(sort(table(L26.dir.array[i, j, which(is.na(L26.dir.array[i,j,]) == F)]), decreasing=T), 'names')[1]), NA)
          L25.mode.imm.dir[i,j] <- ifelse((length(which(is.na(L25.dir.array[i,j,]) == F)) > 0), (attr(sort(table(L25.dir.array[i, j, which(is.na(L25.dir.array[i,j,]) == F)]), decreasing=T), 'names')[1]), NA)
          L24.mode.imm.dir[i,j] <- ifelse((length(which(is.na(L24.dir.array[i,j,]) == F)) > 0), (attr(sort(table(L24.dir.array[i, j, which(is.na(L24.dir.array[i,j,]) == F)]), decreasing=T), 'names')[1]), NA)
          L23.mode.imm.dir[i,j] <- ifelse((length(which(is.na(L23.dir.array[i,j,]) == F)) > 0), (attr(sort(table(L23.dir.array[i, j, which(is.na(L23.dir.array[i,j,]) == F)]), decreasing=T), 'names')[1]), NA)
          L22.mode.imm.dir[i,j] <- ifelse((length(which(is.na(L22.dir.array[i,j,]) == F)) > 0), (attr(sort(table(L22.dir.array[i, j, which(is.na(L22.dir.array[i,j,]) == F)]), decreasing=T), 'names')[1]), NA)
          L21.mode.imm.dir[i,j] <- ifelse((length(which(is.na(L21.dir.array[i,j,]) == F)) > 0), (attr(sort(table(L21.dir.array[i, j, which(is.na(L21.dir.array[i,j,]) == F)]), decreasing=T), 'names')[1]), NA)
          L20.mode.imm.dir[i,j] <- ifelse((length(which(is.na(L20.dir.array[i,j,]) == F)) > 0), (attr(sort(table(L20.dir.array[i, j, which(is.na(L20.dir.array[i,j,]) == F)]), decreasing=T), 'names')[1]), NA)
          L19.mode.imm.dir[i,j] <- ifelse((length(which(is.na(L19.dir.array[i,j,]) == F)) > 0), (attr(sort(table(L19.dir.array[i, j, which(is.na(L19.dir.array[i,j,]) == F)]), decreasing=T), 'names')[1]), NA)
          L18.mode.imm.dir[i,j] <- ifelse((length(which(is.na(L18.dir.array[i,j,]) == F)) > 0), (attr(sort(table(L18.dir.array[i, j, which(is.na(L18.dir.array[i,j,]) == F)]), decreasing=T), 'names')[1]), NA)
          L17.mode.imm.dir[i,j] <- ifelse((length(which(is.na(L17.dir.array[i,j,]) == F)) > 0), (attr(sort(table(L17.dir.array[i, j, which(is.na(L17.dir.array[i,j,]) == F)]), decreasing=T), 'names')[1]), NA)
          L16.mode.imm.dir[i,j] <- ifelse((length(which(is.na(L16.dir.array[i,j,]) == F)) > 0), (attr(sort(table(L16.dir.array[i, j, which(is.na(L16.dir.array[i,j,]) == F)]), decreasing=T), 'names')[1]), NA)
          L15.mode.imm.dir[i,j] <- ifelse((length(which(is.na(L15.dir.array[i,j,]) == F)) > 0), (attr(sort(table(L15.dir.array[i, j, which(is.na(L15.dir.array[i,j,]) == F)]), decreasing=T), 'names')[1]), NA)
          L14.mode.imm.dir[i,j] <- ifelse((length(which(is.na(L14.dir.array[i,j,]) == F)) > 0), (attr(sort(table(L14.dir.array[i, j, which(is.na(L14.dir.array[i,j,]) == F)]), decreasing=T), 'names')[1]), NA)
          L13.mode.imm.dir[i,j] <- ifelse((length(which(is.na(L13.dir.array[i,j,]) == F)) > 0), (attr(sort(table(L13.dir.array[i, j, which(is.na(L13.dir.array[i,j,]) == F)]), decreasing=T), 'names')[1]), NA)
          L12.mode.imm.dir[i,j] <- ifelse((length(which(is.na(L12.dir.array[i,j,]) == F)) > 0), (attr(sort(table(L12.dir.array[i, j, which(is.na(L12.dir.array[i,j,]) == F)]), decreasing=T), 'names')[1]), NA)
          L11.mode.imm.dir[i,j] <- ifelse((length(which(is.na(L11.dir.array[i,j,]) == F)) > 0), (attr(sort(table(L11.dir.array[i, j, which(is.na(L11.dir.array[i,j,]) == F)]), decreasing=T), 'names')[1]), NA)
          L10.mode.imm.dir[i,j] <- ifelse((length(which(is.na(L10.dir.array[i,j,]) == F)) > 0), (attr(sort(table(L10.dir.array[i, j, which(is.na(L10.dir.array[i,j,]) == F)]), decreasing=T), 'names')[1]), NA)
          L09.mode.imm.dir[i,j] <- ifelse((length(which(is.na(L09.dir.array[i,j,]) == F)) > 0), (attr(sort(table(L09.dir.array[i, j, which(is.na(L09.dir.array[i,j,]) == F)]), decreasing=T), 'names')[1]), NA)
          L08.mode.imm.dir[i,j] <- ifelse((length(which(is.na(L08.dir.array[i,j,]) == F)) > 0), (attr(sort(table(L08.dir.array[i, j, which(is.na(L08.dir.array[i,j,]) == F)]), decreasing=T), 'names')[1]), NA)
          L07.mode.imm.dir[i,j] <- ifelse((length(which(is.na(L07.dir.array[i,j,]) == F)) > 0), (attr(sort(table(L07.dir.array[i, j, which(is.na(L07.dir.array[i,j,]) == F)]), decreasing=T), 'names')[1]), NA)
          
          L2624.mode.imm.dir[i,j] <- ifelse((length(which(is.na(L2624.dir.array[i,j,]) == F)) > 0), (attr(sort(table(L2624.dir.array[i, j, which(is.na(L2624.dir.array[i,j,]) == F)]), decreasing=T), 'names')[1]), NA)
          L2422.mode.imm.dir[i,j] <- ifelse((length(which(is.na(L2422.dir.array[i,j,]) == F)) > 0), (attr(sort(table(L2422.dir.array[i, j, which(is.na(L2422.dir.array[i,j,]) == F)]), decreasing=T), 'names')[1]), NA)
          L2220.mode.imm.dir[i,j] <- ifelse((length(which(is.na(L2220.dir.array[i,j,]) == F)) > 0), (attr(sort(table(L2220.dir.array[i, j, which(is.na(L2220.dir.array[i,j,]) == F)]), decreasing=T), 'names')[1]), NA)
          L2018.mode.imm.dir[i,j] <- ifelse((length(which(is.na(L2018.dir.array[i,j,]) == F)) > 0), (attr(sort(table(L2018.dir.array[i, j, which(is.na(L2018.dir.array[i,j,]) == F)]), decreasing=T), 'names')[1]), NA)
          L1816.mode.imm.dir[i,j] <- ifelse((length(which(is.na(L1816.dir.array[i,j,]) == F)) > 0), (attr(sort(table(L1816.dir.array[i, j, which(is.na(L1816.dir.array[i,j,]) == F)]), decreasing=T), 'names')[1]), NA)
          L1614.mode.imm.dir[i,j] <- ifelse((length(which(is.na(L1614.dir.array[i,j,]) == F)) > 0), (attr(sort(table(L1614.dir.array[i, j, which(is.na(L1614.dir.array[i,j,]) == F)]), decreasing=T), 'names')[1]), NA)
          L1412.mode.imm.dir[i,j] <- ifelse((length(which(is.na(L1412.dir.array[i,j,]) == F)) > 0), (attr(sort(table(L1412.dir.array[i, j, which(is.na(L1412.dir.array[i,j,]) == F)]), decreasing=T), 'names')[1]), NA)
          L1210.mode.imm.dir[i,j] <- ifelse((length(which(is.na(L1210.dir.array[i,j,]) == F)) > 0), (attr(sort(table(L1210.dir.array[i, j, which(is.na(L1210.dir.array[i,j,]) == F)]), decreasing=T), 'names')[1]), NA)
          L1008.mode.imm.dir[i,j] <- ifelse((length(which(is.na(L1008.dir.array[i,j,]) == F)) > 0), (attr(sort(table(L1008.dir.array[i, j, which(is.na(L1008.dir.array[i,j,]) == F)]), decreasing=T), 'names')[1]), NA)
          L0806.mode.imm.dir[i,j] <- ifelse((length(which(is.na(L0806.dir.array[i,j,]) == F)) > 0), (attr(sort(table(L0806.dir.array[i, j, which(is.na(L0806.dir.array[i,j,]) == F)]), decreasing=T), 'names')[1]), NA)
        }
      }
      
      ALL.mode.imm.dir.arr[,,m] <- ALL.mode.imm.dir
      MIS3.mode.imm.dir.arr[,,m] <- MIS3.mode.imm.dir
      MIS2.mode.imm.dir.arr[,,m] <- MIS2.mode.imm.dir
      LGM.mode.imm.dir.arr[,,m] <- LGM.mode.imm.dir
      LGMc.mode.imm.dir.arr[,,m] <- LGMc.mode.imm.dir
      LGMu.mode.imm.dir.arr[,,m] <- LGMu.mode.imm.dir
      HOL.mode.imm.dir.arr[,,m] <- HOL.mode.imm.dir
      L26.mode.imm.dir.arr[,,m] <- L26.mode.imm.dir
      L25.mode.imm.dir.arr[,,m] <- L25.mode.imm.dir
      L24.mode.imm.dir.arr[,,m] <- L24.mode.imm.dir
      L23.mode.imm.dir.arr[,,m] <- L23.mode.imm.dir
      L22.mode.imm.dir.arr[,,m] <- L22.mode.imm.dir
      L21.mode.imm.dir.arr[,,m] <- L21.mode.imm.dir
      L20.mode.imm.dir.arr[,,m] <- L20.mode.imm.dir
      L19.mode.imm.dir.arr[,,m] <- L19.mode.imm.dir
      L18.mode.imm.dir.arr[,,m] <- L18.mode.imm.dir
      L17.mode.imm.dir.arr[,,m] <- L17.mode.imm.dir
      L16.mode.imm.dir.arr[,,m] <- L16.mode.imm.dir
      L15.mode.imm.dir.arr[,,m] <- L15.mode.imm.dir
      L14.mode.imm.dir.arr[,,m] <- L14.mode.imm.dir
      L13.mode.imm.dir.arr[,,m] <- L13.mode.imm.dir
      L12.mode.imm.dir.arr[,,m] <- L12.mode.imm.dir
      L11.mode.imm.dir.arr[,,m] <- L11.mode.imm.dir
      L10.mode.imm.dir.arr[,,m] <- L10.mode.imm.dir
      L09.mode.imm.dir.arr[,,m] <- L09.mode.imm.dir
      L08.mode.imm.dir.arr[,,m] <- L08.mode.imm.dir
      L07.mode.imm.dir.arr[,,m] <- L07.mode.imm.dir
      
      L2624.mode.imm.dir.arr[,,m] <- L2624.mode.imm.dir
      L2422.mode.imm.dir.arr[,,m] <- L2422.mode.imm.dir
      L2220.mode.imm.dir.arr[,,m] <- L2220.mode.imm.dir
      L2018.mode.imm.dir.arr[,,m] <- L2018.mode.imm.dir
      L1816.mode.imm.dir.arr[,,m] <- L1816.mode.imm.dir
      L1614.mode.imm.dir.arr[,,m] <- L1614.mode.imm.dir
      L1412.mode.imm.dir.arr[,,m] <- L1412.mode.imm.dir
      L1210.mode.imm.dir.arr[,,m] <- L1210.mode.imm.dir
      L1008.mode.imm.dir.arr[,,m] <- L1008.mode.imm.dir
      L0806.mode.imm.dir.arr[,,m] <- L0806.mode.imm.dir
      
      
      print("*********************************************************************")
      print(paste("run = ", m, "; ", "time elapsed = ", round(as.numeric((proc.time() - proc.sim.start)[3])/60/60, 1), " hours", sep=""))
      print("*********************************************************************")
      
      
    } # end m reps loop
    
    cv.N.mat.mn <- apply(cv.N.arr, c(1,2), mean, na.rm=T)
    #image(rot.mat(cv.N.mat.mn), col=rev(grey(1:100/100)))
    
    MIS3.cv.N.mat.mn <-  apply(MIS3.cv.N, c(1,2), mean, na.rm=T)
    #image(rot.mat(MIS3.cv.N.mat.mn), col=rev(grey(1:100/100)))
    
    MIS2.cv.N.mat.mn <-  apply(MIS2.cv.N, c(1,2), mean, na.rm=T)
    #image(rot.mat(MIS2.cv.N.mat.mn), col=rev(grey(1:100/100)))
    
    LGM.cv.N.mat.mn <-  apply(LGM.cv.N, c(1,2), mean, na.rm=T)
    #image(rot.mat(LGM.cv.N.mat.mn), col=rev(grey(1:100/100)))
    
    LGMc.cv.N.mat.mn <-  apply(LGMc.cv.N, c(1,2), mean, na.rm=T)
    #image(rot.mat(LGMc.cv.N.mat.mn), col=rev(grey(1:100/100)))
    
    LGMu.cv.N.mat.mn <-  apply(LGMu.cv.N, c(1,2), mean, na.rm=T)
    #image(rot.mat(LGMu.cv.N.mat.mn), col=rev(grey(1:100/100)))
    
    HOL.cv.N.mat.mn <-  apply(HOL.cv.N, c(1,2), mean, na.rm=T)
    #image(rot.mat(HOL.cv.N.mat.mn), col=rev(grey(1:100/100)))
    
    ## create & export rasters
    lat.vec <- -seq(0.5,43,0.5)
    lon.vec <- seq(110.5,153.5,0.5)
    llat <- length(lat.vec)
    llon <- length(lon.vec)
    
    Nxyz <- MIS3xyz <- MIS2xyz <- LGMxyz <- LGMcxyz <- LGMuxyz <- HOLxyz <- matrix(data=NA,nrow=1,ncol=3)
    for (i in 1:llat) {
      for (j in 1:llon) {
        Nxyz <- rbind(Nxyz, c(lon.vec[j],lat.vec[i],cv.N.mat.mn[i,j]))
        MIS3xyz <- rbind(MIS3xyz, c(lon.vec[j],lat.vec[i],MIS3.cv.N.mat.mn[i,j]))
        MIS2xyz <- rbind(MIS2xyz, c(lon.vec[j],lat.vec[i],MIS2.cv.N.mat.mn[i,j]))
        LGMxyz <- rbind(LGMxyz, c(lon.vec[j],lat.vec[i],LGM.cv.N.mat.mn[i,j]))
        LGMcxyz <- rbind(LGMcxyz, c(lon.vec[j],lat.vec[i],LGMc.cv.N.mat.mn[i,j]))
        LGMuxyz <- rbind(LGMuxyz, c(lon.vec[j],lat.vec[i],LGMu.cv.N.mat.mn[i,j]))
        HOLxyz <- rbind(HOLxyz, c(lon.vec[j],lat.vec[i],HOL.cv.N.mat.mn[i,j]))
        
      }
    }
    
    Nxyz <- Nxyz[-1,]
    N.xyz <- as.data.frame(Nxyz)
    colnames(N.xyz) <- c("x","y","CV")
    head(N.xyz)
    N.rast <- rasterFromXYZ(N.xyz, crs=CRS("+proj=longlat +datum=WGS84"))
    plot(N.rast)
    writeRaster(N.rast, filename="NCV100ka2sn75.grd", format="raster", overwrite=T)
    
    MIS3xyz <- MIS3xyz[-1,]
    MIS3.xyz <- as.data.frame(MIS3xyz)
    colnames(MIS3.xyz) <- c("x","y","CV")
    head(MIS3.xyz)
    MIS3.rast <- rasterFromXYZ(MIS3.xyz, crs=CRS("+proj=longlat +datum=WGS84"))
    plot(MIS3.rast)
    writeRaster(MIS3.rast, filename="MIS3CV100ka2sn75.grd", format="raster", overwrite=T)
    
    MIS2xyz <- MIS2xyz[-1,]
    MIS2.xyz <- as.data.frame(MIS2xyz)
    colnames(MIS2.xyz) <- c("x","y","CV")
    head(MIS2.xyz)
    MIS2.rast <- rasterFromXYZ(MIS2.xyz, crs=CRS("+proj=longlat +datum=WGS84"))
    plot(MIS2.rast)
    writeRaster(MIS2.rast, filename="MIS2CV100ka2sn75.grd", format="raster", overwrite=T)
    
    LGMxyz <- LGMxyz[-1,]
    LGM.xyz <- as.data.frame(LGMxyz)
    colnames(LGM.xyz) <- c("x","y","CV")
    head(LGM.xyz)
    LGM.rast <- rasterFromXYZ(LGM.xyz, crs=CRS("+proj=longlat +datum=WGS84"))
    plot(LGM.rast)
    writeRaster(LGM.rast, filename="LGMCV100ka2sn75.grd", format="raster", overwrite=T)
    
    LGMcxyz <- LGMcxyz[-1,]
    LGMc.xyz <- as.data.frame(LGMcxyz)
    colnames(LGMc.xyz) <- c("x","y","CV")
    head(LGMc.xyz)
    LGMc.rast <- rasterFromXYZ(LGMc.xyz, crs=CRS("+proj=longlat +datum=WGS84"))
    plot(LGMc.rast)
    writeRaster(LGMc.rast, filename="LGMcCV100ka2sn75.grd", format="raster", overwrite=T)
    
    LGMuxyz <- LGMuxyz[-1,]
    LGMu.xyz <- as.data.frame(LGMuxyz)
    colnames(LGMu.xyz) <- c("x","y","CV")
    head(LGMu.xyz)
    LGMu.rast <- rasterFromXYZ(LGMu.xyz, crs=CRS("+proj=longlat +datum=WGS84"))
    plot(LGMu.rast)
    writeRaster(LGMu.rast, filename="LGMuCV100ka2sn75.grd", format="raster", overwrite=T)
    
    HOLxyz <- HOLxyz[-1,]
    HOL.xyz <- as.data.frame(HOLxyz)
    colnames(HOL.xyz) <- c("x","y","CV")
    head(HOL.xyz)
    HOL.rast <- rasterFromXYZ(HOL.xyz, crs=CRS("+proj=longlat +datum=WGS84"))
    plot(HOL.rast)
    writeRaster(HOL.rast, filename="HOLCV100ka2sn75.grd", format="raster", overwrite=T)
    
    ## output mean + CI for N trajectory
    # Sahul (all)
    N.mean.vec <- apply(N.mat, 2, median, na.rm=T)
    N.lo.vec <- apply(N.mat, 2, quantile, probs=0.025, na.rm=T)
    N.up.vec <- apply(N.mat, 2, quantile, probs=0.975, na.rm=T)
    
    # AUS
    N.AUS.mean.vec <- apply(AUS.N.mat, 2, median, na.rm=T)
    N.AUS.lo.vec <- apply(AUS.N.mat, 2, quantile, probs=0.025, na.rm=T)
    N.AUS.up.vec <- apply(AUS.N.mat, 2, quantile, probs=0.975, na.rm=T)
    
    # NG
    N.NG.mean.vec <- apply(NG.N.mat, 2, median, na.rm=T)
    N.NG.lo.vec <- apply(NG.N.mat, 2, quantile, probs=0.025, na.rm=T)
    N.NG.up.vec <- apply(NG.N.mat, 2, quantile, probs=0.975, na.rm=T)
    
    # CENTRE
    N.CEN.mean.vec <- apply(CEN.N.mat, 2, median, na.rm=T)
    N.CEN.lo.vec <- apply(CEN.N.mat, 2, quantile, probs=0.025, na.rm=T)
    N.CEN.up.vec <- apply(CEN.N.mat, 2, quantile, probs=0.975, na.rm=T)
    
    # TOP END
    N.TE.mean.vec <- apply(TE.N.mat, 2, median, na.rm=T)
    N.TE.lo.vec <- apply(TE.N.mat, 2, quantile, probs=0.025, na.rm=T)
    N.TE.up.vec <- apply(TE.N.mat, 2, quantile, probs=0.975, na.rm=T)
    
    # EASTERN
    N.E.mean.vec <- apply(E.N.mat, 2, median, na.rm=T)
    N.E.lo.vec <- apply(E.N.mat, 2, quantile, probs=0.025, na.rm=T)
    N.E.up.vec <- apply(E.N.mat, 2, quantile, probs=0.975, na.rm=T)
    
    # NEW GUINEA HIGHLANDS
    N.NGH.mean.vec <- apply(NGH.N.mat, 2, median, na.rm=T)
    N.NGH.lo.vec <- apply(NGH.N.mat, 2, quantile, probs=0.025, na.rm=T)
    N.NGH.up.vec <- apply(NGH.N.mat, 2, quantile, probs=0.975, na.rm=T)
    
    ## bioregions
    # Pilbara + Gascoyne + Murchison
    N.PILGASMUR.mean.vec <- apply(PILGASMUR.N.mat, 2, median, na.rm=T)
    N.PILGASMUR.lo.vec <- apply(PILGASMUR.N.mat, 2, quantile, probs=0.025, na.rm=T)
    N.PILGASMUR.up.vec <- apply(PILGASMUR.N.mat, 2, quantile, probs=0.975, na.rm=T)
    
    # Warren
    N.WAR.mean.vec <- apply(WAR.N.mat, 2, median, na.rm=T)
    N.WAR.lo.vec <- apply(WAR.N.mat, 2, quantile, probs=0.025, na.rm=T)
    N.WAR.up.vec <- apply(WAR.N.mat, 2, quantile, probs=0.975, na.rm=T)
    
    # Central Ranges
    N.CER.mean.vec <- apply(CER.N.mat, 2, median, na.rm=T)
    N.CER.lo.vec <- apply(CER.N.mat, 2, quantile, probs=0.025, na.rm=T)
    N.CER.up.vec <- apply(CER.N.mat, 2, quantile, probs=0.975, na.rm=T)
    
    # Ord Victoria Plain
    N.OVP.mean.vec <- apply(OVP.N.mat, 2, median, na.rm=T)
    N.OVP.lo.vec <- apply(OVP.N.mat, 2, quantile, probs=0.025, na.rm=T)
    N.OVP.up.vec <- apply(OVP.N.mat, 2, quantile, probs=0.975, na.rm=T)
    
    # Arnhem Plateau
    N.ARP.mean.vec <- apply(ARP.N.mat, 2, median, na.rm=T)
    N.ARP.lo.vec <- apply(ARP.N.mat, 2, quantile, probs=0.025, na.rm=T)
    N.ARP.up.vec <- apply(ARP.N.mat, 2, quantile, probs=0.975, na.rm=T)
    
    # Gulf Plains + Einasleigh Uplands
    N.GUPEIU.mean.vec <- apply(GUPEIU.N.mat, 2, median, na.rm=T)
    N.GUPEIU.lo.vec <- apply(GUPEIU.N.mat, 2, quantile, probs=0.025, na.rm=T)
    N.GUPEIU.up.vec <- apply(GUPEIU.N.mat, 2, quantile, probs=0.975, na.rm=T)
    
    # Brigalow Belt South + Nandewar + Sydney Basin + NSW South Western Slopes + South Eastern Highlands
    N.BNSNS.mean.vec <- apply(BNSNS.N.mat, 2, median, na.rm=T)
    N.BNSNS.lo.vec <- apply(BNSNS.N.mat, 2, quantile, probs=0.025, na.rm=T)
    N.BNSNS.up.vec <- apply(BNSNS.N.mat, 2, quantile, probs=0.975, na.rm=T)
    
    # Murray Darling Depression + Broken Hill Complex
    N.MDDBHC.mean.vec <- apply(MDDBHC.N.mat, 2, median, na.rm=T)
    N.MDDBHC.lo.vec <- apply(MDDBHC.N.mat, 2, quantile, probs=0.025, na.rm=T)
    N.MDDBHC.up.vec <- apply(MDDBHC.N.mat, 2, quantile, probs=0.975, na.rm=T)
    
    # Tasmanian Central Highlands
    N.TCH.mean.vec <- apply(TCH.N.mat, 2, median, na.rm=T)
    N.TCH.lo.vec <- apply(TCH.N.mat, 2, quantile, probs=0.025, na.rm=T)
    N.TCH.up.vec <- apply(TCH.N.mat, 2, quantile, probs=0.975, na.rm=T)
    
    ## output mean + CI for D trajectory
    # CENTRE
    D.CEN.mean.vec <- apply(CEN.D.mat, 2, median, na.rm=T)
    D.CEN.lo.vec <- apply(CEN.D.mat, 2, quantile, probs=0.025, na.rm=T)
    D.CEN.up.vec <- apply(CEN.D.mat, 2, quantile, probs=0.975, na.rm=T)
    
    # TOP END
    D.TE.mean.vec <- apply(TE.D.mat, 2, median, na.rm=T)
    D.TE.lo.vec <- apply(TE.D.mat, 2, quantile, probs=0.025, na.rm=T)
    D.TE.up.vec <- apply(TE.D.mat, 2, quantile, probs=0.975, na.rm=T)
    
    # EASTERN
    D.E.mean.vec <- apply(E.D.mat, 2, median, na.rm=T)
    D.E.lo.vec <- apply(E.D.mat, 2, quantile, probs=0.025, na.rm=T)
    D.E.up.vec <- apply(E.D.mat, 2, quantile, probs=0.975, na.rm=T)
    
    # NEW GUINEA HIGHLANDS
    D.NGH.mean.vec <- apply(NGH.D.mat, 2, median, na.rm=T)
    D.NGH.lo.vec <- apply(NGH.D.mat, 2, quantile, probs=0.025, na.rm=T)
    D.NGH.up.vec <- apply(NGH.D.mat, 2, quantile, probs=0.975, na.rm=T)
    
    ## bioregions
    # Pilbara + Gascoyne + Murchison
    D.PILGASMUR.mean.vec <- apply(PILGASMUR.D.mat, 2, median, na.rm=T)
    D.PILGASMUR.lo.vec <- apply(PILGASMUR.D.mat, 2, quantile, probs=0.025, na.rm=T)
    D.PILGASMUR.up.vec <- apply(PILGASMUR.D.mat, 2, quantile, probs=0.975, na.rm=T)
    
    # Warren
    D.WAR.mean.vec <- apply(WAR.D.mat, 2, median, na.rm=T)
    D.WAR.lo.vec <- apply(WAR.D.mat, 2, quantile, probs=0.025, na.rm=T)
    D.WAR.up.vec <- apply(WAR.D.mat, 2, quantile, probs=0.975, na.rm=T)
    
    # Central Ranges
    D.CER.mean.vec <- apply(CER.D.mat, 2, median, na.rm=T)
    D.CER.lo.vec <- apply(CER.D.mat, 2, quantile, probs=0.025, na.rm=T)
    D.CER.up.vec <- apply(CER.D.mat, 2, quantile, probs=0.975, na.rm=T)
    
    # Ord Victoria Plain
    D.OVP.mean.vec <- apply(OVP.D.mat, 2, median, na.rm=T)
    D.OVP.lo.vec <- apply(OVP.D.mat, 2, quantile, probs=0.025, na.rm=T)
    D.OVP.up.vec <- apply(OVP.D.mat, 2, quantile, probs=0.975, na.rm=T)
    
    # Arnhem Plateau
    D.ARP.mean.vec <- apply(ARP.D.mat, 2, median, na.rm=T)
    D.ARP.lo.vec <- apply(ARP.D.mat, 2, quantile, probs=0.025, na.rm=T)
    D.ARP.up.vec <- apply(ARP.D.mat, 2, quantile, probs=0.975, na.rm=T)
    
    # Gulf Plains + Einasleigh Uplands
    D.GUPEIU.mean.vec <- apply(GUPEIU.D.mat, 2, median, na.rm=T)
    D.GUPEIU.lo.vec <- apply(GUPEIU.D.mat, 2, quantile, probs=0.025, na.rm=T)
    D.GUPEIU.up.vec <- apply(GUPEIU.D.mat, 2, quantile, probs=0.975, na.rm=T)
    
    # Brigalow Belt South + Nandewar + Sydney Basin + NSW South Western Slopes + South Eastern Highlands
    D.BNSNS.mean.vec <- apply(BNSNS.D.mat, 2, median, na.rm=T)
    D.BNSNS.lo.vec <- apply(BNSNS.D.mat, 2, quantile, probs=0.025, na.rm=T)
    D.BNSNS.up.vec <- apply(BNSNS.D.mat, 2, quantile, probs=0.975, na.rm=T)
    
    # Murray Darling Depression + Broken Hill Complex
    D.MDDBHC.mean.vec <- apply(MDDBHC.D.mat, 2, median, na.rm=T)
    D.MDDBHC.lo.vec <- apply(MDDBHC.D.mat, 2, quantile, probs=0.025, na.rm=T)
    D.MDDBHC.up.vec <- apply(MDDBHC.D.mat, 2, quantile, probs=0.975, na.rm=T)
    
    # Tasmanian Central Highlands
    D.TCH.mean.vec <- apply(TCH.D.mat, 2, median, na.rm=T)
    D.TCH.lo.vec <- apply(TCH.D.mat, 2, quantile, probs=0.025, na.rm=T)
    D.TCH.up.vec <- apply(TCH.D.mat, 2, quantile, probs=0.975, na.rm=T)
    
    forgen.proj <- 1:(gen.run+1)
    ybp <- entry.date - round(forgen.proj*gen.l, 0)
    
    dat.N.out <- data.frame(ybp,N.mean.vec,N.up.vec,N.lo.vec,N.AUS.mean.vec,N.AUS.up.vec,N.AUS.lo.vec,
                            N.NG.mean.vec,N.NG.up.vec,N.NG.lo.vec, N.CEN.mean.vec, N.CEN.up.vec, N.CEN.lo.vec,
                            N.TE.mean.vec, N.TE.up.vec, N.TE.lo.vec, N.E.mean.vec, N.E.up.vec, N.E.lo.vec,
                            N.NGH.mean.vec, N.NGH.up.vec, N.NGH.lo.vec, N.PILGASMUR.mean.vec, N.PILGASMUR.up.vec, N.PILGASMUR.lo.vec,
                            N.WAR.mean.vec, N.WAR.up.vec, N.WAR.lo.vec, N.CER.mean.vec, N.CER.up.vec, N.CER.lo.vec,
                            N.OVP.mean.vec, N.OVP.up.vec, N.OVP.lo.vec, N.ARP.mean.vec, N.ARP.up.vec, N.ARP.lo.vec,
                            N.GUPEIU.mean.vec, N.GUPEIU.up.vec, N.GUPEIU.lo.vec, N.BNSNS.mean.vec, N.BNSNS.up.vec, N.BNSNS.lo.vec,
                            N.MDDBHC.mean.vec, N.MDDBHC.up.vec, N.MDDBHC.lo.vec, N.TCH.mean.vec, N.TCH.up.vec, N.TCH.lo.vec)
    
    colnames(dat.N.out) <- c("ybp","Nm","Nup","Nlo","NmAUS","NupAUS","NloAUS","NmNG","NupNG","NloNG","NmCEN","NupCEN","NloCEN",
                             "NmTE","NupTE","NloTE","NmE","NupE","NloE","NmNGH","NupNGH","NloNGH","NmPILGASMUR","NupPILGASMUR","NloPILGASMUR",
                             "NmWAR","NupWAR","NloWAR","NmCER","NupCER","NloCER","NmOVP","NupOVP","NloOVP","NmARP","NupARP","NloARP",
                             "NmGUPEIU","NupGUPEIU","NloGUPEIU","NmBNSNS","NupBNSNS","NloBNSNS","NmMDDBHC","NupMDDBHC","NloMDDBHC",
                             "NmTCH","NupTCH","NloTCH")
    
    dat.D.out <- data.frame(ybp,D.CEN.mean.vec, D.CEN.up.vec, D.CEN.lo.vec,
                            D.TE.mean.vec, D.TE.up.vec, D.TE.lo.vec, D.E.mean.vec, D.E.up.vec, D.E.lo.vec,
                            D.NGH.mean.vec, D.NGH.up.vec, D.NGH.lo.vec,D.PILGASMUR.mean.vec, D.PILGASMUR.up.vec, D.PILGASMUR.lo.vec,
                            D.WAR.mean.vec, D.WAR.up.vec, D.WAR.lo.vec,D.CER.mean.vec, D.CER.up.vec, D.CER.lo.vec,
                            D.OVP.mean.vec, D.OVP.up.vec, D.OVP.lo.vec,D.ARP.mean.vec, D.ARP.up.vec, D.ARP.lo.vec,
                            D.GUPEIU.mean.vec, D.GUPEIU.up.vec, D.GUPEIU.lo.vec, D.BNSNS.mean.vec, D.BNSNS.up.vec, D.BNSNS.lo.vec,
                            D.MDDBHC.mean.vec, D.MDDBHC.up.vec, D.MDDBHC.lo.vec, D.TCH.mean.vec, D.TCH.up.vec, D.TCH.lo.vec)
    
    colnames(dat.D.out) <- c("ybp","DmCEN","DupCEN","DloCEN",
                             "DmTE","DupTE","DloTE","DmE","DupE","DloE","DmNGH","DupNGH","DloNGH","DmPILGASMUR","DupPILGASMUR","DloPILGASMUR",
                             "DmWAR","DupWAR","DloWAR","DmCER","DupCER","DloCER","DmOVP","DupOVP","DloOVP","DmARP","DupARP","DloARP",
                             "DmGUPEIU","DupGUPEIU","DloGUPEIU","DmBNSNS","DupBNSNS","DloBNSNS","DmMDDBHC","DupMDDBHC","DloMDDBHC",
                             "DmTCH","DupTCH","DloTCH")
    
    datN.out.name <- paste("ScenarioMean100ka2-LT-Nproj-Nout-",name.add,".csv",sep="")
    datD.out.name <- paste("ScenarioMean100ka2-LT-Dproj-Nout-",name.add,".csv",sep="")
    write.table(dat.N.out,file=datN.out.name,sep=",", row.names = F, col.names = T)
    write.table(dat.D.out,file=datD.out.name,sep=",", row.names = F, col.names = T)
    
    save.image(paste("ScenarioMeanVecDir100ka2-LT-Nproj.",reps,"Reps.",name.add,".RData",sep=""))
    
    # mode of directions over m reps
    detach("package:pracma")  
    library(DescTools)
    ALL.dir.mode <- apply(ALL.mode.imm.dir.arr, c(1,2), Mode, na.rm=T)
    MIS3.dir.mode <- apply(MIS3.mode.imm.dir.arr, c(1,2), Mode, na.rm=T)
    MIS2.dir.mode <- apply(MIS2.mode.imm.dir.arr, c(1,2), Mode, na.rm=T)
    LGM.dir.mode <- apply(LGM.mode.imm.dir.arr, c(1,2), Mode, na.rm=T)
    LGMc.dir.mode <- apply(LGMc.mode.imm.dir.arr, c(1,2), Mode, na.rm=T)
    LGMu.dir.mode <- apply(LGMu.mode.imm.dir.arr, c(1,2), Mode, na.rm=T)
    HOL.dir.mode <- apply(HOL.mode.imm.dir.arr, c(1,2), Mode, na.rm=T)
    L26.dir.mode <- apply(L26.mode.imm.dir.arr, c(1,2), Mode, na.rm=T)
    L25.dir.mode <- apply(L25.mode.imm.dir.arr, c(1,2), Mode, na.rm=T)
    L24.dir.mode <- apply(L24.mode.imm.dir.arr, c(1,2), Mode, na.rm=T)
    L23.dir.mode <- apply(L23.mode.imm.dir.arr, c(1,2), Mode, na.rm=T)
    L22.dir.mode <- apply(L22.mode.imm.dir.arr, c(1,2), Mode, na.rm=T)
    L21.dir.mode <- apply(L21.mode.imm.dir.arr, c(1,2), Mode, na.rm=T)
    L20.dir.mode <- apply(L20.mode.imm.dir.arr, c(1,2), Mode, na.rm=T)
    L19.dir.mode <- apply(L19.mode.imm.dir.arr, c(1,2), Mode, na.rm=T)
    L18.dir.mode <- apply(L18.mode.imm.dir.arr, c(1,2), Mode, na.rm=T)
    L17.dir.mode <- apply(L17.mode.imm.dir.arr, c(1,2), Mode, na.rm=T)
    L16.dir.mode <- apply(L16.mode.imm.dir.arr, c(1,2), Mode, na.rm=T)
    L15.dir.mode <- apply(L15.mode.imm.dir.arr, c(1,2), Mode, na.rm=T)
    L14.dir.mode <- apply(L14.mode.imm.dir.arr, c(1,2), Mode, na.rm=T)
    L13.dir.mode <- apply(L13.mode.imm.dir.arr, c(1,2), Mode, na.rm=T)
    L12.dir.mode <- apply(L12.mode.imm.dir.arr, c(1,2), Mode, na.rm=T)
    L11.dir.mode <- apply(L11.mode.imm.dir.arr, c(1,2), Mode, na.rm=T)
    L10.dir.mode <- apply(L10.mode.imm.dir.arr, c(1,2), Mode, na.rm=T)
    L09.dir.mode <- apply(L09.mode.imm.dir.arr, c(1,2), Mode, na.rm=T)
    L08.dir.mode <- apply(L08.mode.imm.dir.arr, c(1,2), Mode, na.rm=T)
    L07.dir.mode <- apply(L07.mode.imm.dir.arr, c(1,2), Mode, na.rm=T)
    
    L2624.dir.mode <- apply(L2624.mode.imm.dir.arr, c(1,2), Mode, na.rm=T)
    L2422.dir.mode <- apply(L2422.mode.imm.dir.arr, c(1,2), Mode, na.rm=T)
    L2220.dir.mode <- apply(L2220.mode.imm.dir.arr, c(1,2), Mode, na.rm=T)
    L2018.dir.mode <- apply(L2018.mode.imm.dir.arr, c(1,2), Mode, na.rm=T)
    L1816.dir.mode <- apply(L1816.mode.imm.dir.arr, c(1,2), Mode, na.rm=T)
    L1614.dir.mode <- apply(L1614.mode.imm.dir.arr, c(1,2), Mode, na.rm=T)
    L1412.dir.mode <- apply(L1412.mode.imm.dir.arr, c(1,2), Mode, na.rm=T)
    L1210.dir.mode <- apply(L1210.mode.imm.dir.arr, c(1,2), Mode, na.rm=T)
    L1008.dir.mode <- apply(L1008.mode.imm.dir.arr, c(1,2), Mode, na.rm=T)
    L0806.dir.mode <- apply(L0806.mode.imm.dir.arr, c(1,2), Mode, na.rm=T)
    
    ALL.dir.xyz <- MIS3.dir.xyz <- MIS2.dir.xyz <- LGM.dir.xyz <- LGMc.dir.xyz <- LGMu.dir.xyz <- HOL.dir.xyz <- matrix(data=NA,nrow=1,ncol=3)
    for (i in 1:llat) {
      for (j in 1:llon) {
        ALL.dir.xyz <- rbind(ALL.dir.xyz, c(lon.vec[j],lat.vec[i],ALL.dir.mode[i,j]))
        MIS3.dir.xyz <- rbind(MIS3.dir.xyz, c(lon.vec[j],lat.vec[i],MIS3.dir.mode[i,j]))
        MIS2.dir.xyz <- rbind(MIS2.dir.xyz, c(lon.vec[j],lat.vec[i],MIS2.dir.mode[i,j]))
        LGM.dir.xyz <- rbind(LGM.dir.xyz, c(lon.vec[j],lat.vec[i],LGM.dir.mode[i,j]))
        LGMc.dir.xyz <- rbind(LGMc.dir.xyz, c(lon.vec[j],lat.vec[i],LGMc.dir.mode[i,j]))
        LGMu.dir.xyz <- rbind(LGMu.dir.xyz, c(lon.vec[j],lat.vec[i],LGMu.dir.mode[i,j]))
        HOL.dir.xyz <- rbind(HOL.dir.xyz, c(lon.vec[j],lat.vec[i],HOL.dir.mode[i,j]))
        
      }
    }
    
    L26.dir.xyz <- L25.dir.xyz <- L24.dir.xyz <- L23.dir.xyz <- L22.dir.xyz <- L21.dir.xyz <- L20.dir.xyz <- L19.dir.xyz <- L18.dir.xyz <- L17.dir.xyz <- 
      L16.dir.xyz <- L15.dir.xyz <- L14.dir.xyz <- L13.dir.xyz <- L12.dir.xyz <- L11.dir.xyz <- L10.dir.xyz <-
      L09.dir.xyz <- L08.dir.xyz <- L07.dir.xyz <- matrix(data=NA,nrow=1,ncol=3)
    for (i in 1:llat) {
      for (j in 1:llon) {
        L26.dir.xyz <- rbind(L26.dir.xyz, c(lon.vec[j],lat.vec[i],L26.dir.mode[i,j]))
        L25.dir.xyz <- rbind(L25.dir.xyz, c(lon.vec[j],lat.vec[i],L25.dir.mode[i,j]))
        L24.dir.xyz <- rbind(L24.dir.xyz, c(lon.vec[j],lat.vec[i],L24.dir.mode[i,j]))
        L23.dir.xyz <- rbind(L23.dir.xyz, c(lon.vec[j],lat.vec[i],L23.dir.mode[i,j]))
        L22.dir.xyz <- rbind(L22.dir.xyz, c(lon.vec[j],lat.vec[i],L22.dir.mode[i,j]))
        L21.dir.xyz <- rbind(L21.dir.xyz, c(lon.vec[j],lat.vec[i],L21.dir.mode[i,j]))
        L20.dir.xyz <- rbind(L20.dir.xyz, c(lon.vec[j],lat.vec[i],L20.dir.mode[i,j]))
        L19.dir.xyz <- rbind(L19.dir.xyz, c(lon.vec[j],lat.vec[i],L19.dir.mode[i,j]))
        L18.dir.xyz <- rbind(L18.dir.xyz, c(lon.vec[j],lat.vec[i],L18.dir.mode[i,j]))
        L17.dir.xyz <- rbind(L17.dir.xyz, c(lon.vec[j],lat.vec[i],L17.dir.mode[i,j]))
        L16.dir.xyz <- rbind(L16.dir.xyz, c(lon.vec[j],lat.vec[i],L16.dir.mode[i,j]))
        L15.dir.xyz <- rbind(L15.dir.xyz, c(lon.vec[j],lat.vec[i],L15.dir.mode[i,j]))
        L14.dir.xyz <- rbind(L14.dir.xyz, c(lon.vec[j],lat.vec[i],L14.dir.mode[i,j]))
        L13.dir.xyz <- rbind(L13.dir.xyz, c(lon.vec[j],lat.vec[i],L13.dir.mode[i,j]))
        L12.dir.xyz <- rbind(L12.dir.xyz, c(lon.vec[j],lat.vec[i],L12.dir.mode[i,j]))
        L11.dir.xyz <- rbind(L11.dir.xyz, c(lon.vec[j],lat.vec[i],L11.dir.mode[i,j]))
        L10.dir.xyz <- rbind(L10.dir.xyz, c(lon.vec[j],lat.vec[i],L10.dir.mode[i,j]))
        L09.dir.xyz <- rbind(L09.dir.xyz, c(lon.vec[j],lat.vec[i],L09.dir.mode[i,j]))
        L08.dir.xyz <- rbind(L08.dir.xyz, c(lon.vec[j],lat.vec[i],L08.dir.mode[i,j]))
        L07.dir.xyz <- rbind(L07.dir.xyz, c(lon.vec[j],lat.vec[i],L07.dir.mode[i,j]))
      }
    }
    
    L2624.dir.xyz <- L2422.dir.xyz <- L2220.dir.xyz <- L2018.dir.xyz <- L1816.dir.xyz <- L1614.dir.xyz <- L1412.dir.xyz <- L1210.dir.xyz <- L1008.dir.xyz <- L0806.dir.xyz <- 
      matrix(data=NA,nrow=1,ncol=3)
    for (i in 1:llat) {
      for (j in 1:llon) {
        L2624.dir.xyz <- rbind(L2624.dir.xyz, c(lon.vec[j],lat.vec[i],L2624.dir.mode[i,j]))
        L2422.dir.xyz <- rbind(L2422.dir.xyz, c(lon.vec[j],lat.vec[i],L2422.dir.mode[i,j]))
        L2220.dir.xyz <- rbind(L2220.dir.xyz, c(lon.vec[j],lat.vec[i],L2220.dir.mode[i,j]))
        L2018.dir.xyz <- rbind(L2018.dir.xyz, c(lon.vec[j],lat.vec[i],L2018.dir.mode[i,j]))
        L1816.dir.xyz <- rbind(L1816.dir.xyz, c(lon.vec[j],lat.vec[i],L1816.dir.mode[i,j]))
        L1614.dir.xyz <- rbind(L1614.dir.xyz, c(lon.vec[j],lat.vec[i],L1614.dir.mode[i,j]))
        L1412.dir.xyz <- rbind(L1412.dir.xyz, c(lon.vec[j],lat.vec[i],L1412.dir.mode[i,j]))
        L1210.dir.xyz <- rbind(L1210.dir.xyz, c(lon.vec[j],lat.vec[i],L1210.dir.mode[i,j]))
        L1008.dir.xyz <- rbind(L1008.dir.xyz, c(lon.vec[j],lat.vec[i],L1008.dir.mode[i,j]))
        L0806.dir.xyz <- rbind(L0806.dir.xyz, c(lon.vec[j],lat.vec[i],L0806.dir.mode[i,j]))
      }
    }
    
    ALL.dir.xyz <- ALL.dir.xyz[-1,]
    ALLdir.xyz <- as.data.frame(ALL.dir.xyz)
    colnames(ALLdir.xyz) <- c("x","y","dir")
    ALLdir.xyz$dirdeg <- ifelse(ALLdir.xyz$dir == "E", 0, NA)
    ALLdir.xyz$dirdeg <- ifelse(ALLdir.xyz$dir == "SE", 45, ALLdir.xyz$dirdeg)
    ALLdir.xyz$dirdeg <- ifelse(ALLdir.xyz$dir == "S", 90, ALLdir.xyz$dirdeg)
    ALLdir.xyz$dirdeg <- ifelse(ALLdir.xyz$dir == "SW", 135, ALLdir.xyz$dirdeg)
    ALLdir.xyz$dirdeg <- ifelse(ALLdir.xyz$dir == "W", 180, ALLdir.xyz$dirdeg)
    ALLdir.xyz$dirdeg <- ifelse(ALLdir.xyz$dir == "NW", 235, ALLdir.xyz$dirdeg)
    ALLdir.xyz$dirdeg <- ifelse(ALLdir.xyz$dir == "N", 270, ALLdir.xyz$dirdeg)
    ALLdir.xyz$dirdeg <- ifelse(ALLdir.xyz$dir == "NE", 315, ALLdir.xyz$dirdeg)
    
    ALLdir.dat <- data.frame(x=unlist(ALLdir.xyz$x), y=unlist(ALLdir.xyz$y), dirdeg=unlist(ALLdir.xyz$dirdeg))
    ALLdir.rast <- rasterFromXYZ(ALLdir.dat, crs=CRS("+proj=longlat +datum=WGS84"))
    writeRaster(ALLdir.rast, filename="ALLdir.grd", format="raster", overwrite=T)
    
    MIS3.dir.xyz <- MIS3.dir.xyz[-1,]
    MIS3dir.xyz <- as.data.frame(MIS3.dir.xyz)
    colnames(MIS3dir.xyz) <- c("x","y","dir")
    MIS3dir.xyz$dirdeg <- ifelse(MIS3dir.xyz$dir == "E", 0, NA)
    MIS3dir.xyz$dirdeg <- ifelse(MIS3dir.xyz$dir == "SE", 45, MIS3dir.xyz$dirdeg)
    MIS3dir.xyz$dirdeg <- ifelse(MIS3dir.xyz$dir == "S", 90, MIS3dir.xyz$dirdeg)
    MIS3dir.xyz$dirdeg <- ifelse(MIS3dir.xyz$dir == "SW", 135, MIS3dir.xyz$dirdeg)
    MIS3dir.xyz$dirdeg <- ifelse(MIS3dir.xyz$dir == "W", 180, MIS3dir.xyz$dirdeg)
    MIS3dir.xyz$dirdeg <- ifelse(MIS3dir.xyz$dir == "NW", 235, MIS3dir.xyz$dirdeg)
    MIS3dir.xyz$dirdeg <- ifelse(MIS3dir.xyz$dir == "N", 270, MIS3dir.xyz$dirdeg)
    MIS3dir.xyz$dirdeg <- ifelse(MIS3dir.xyz$dir == "NE", 315, MIS3dir.xyz$dirdeg)
    
    MIS3dir.dat <- data.frame(x=unlist(MIS3dir.xyz$x), y=unlist(MIS3dir.xyz$y), dirdeg=unlist(MIS3dir.xyz$dirdeg))
    MIS3dir.rast <- rasterFromXYZ(MIS3dir.dat, crs=CRS("+proj=longlat +datum=WGS84"))
    writeRaster(MIS3dir.rast, filename="MIS3dir.grd", format="raster", overwrite=T)
    
    MIS2.dir.xyz <- MIS2.dir.xyz[-1,]
    MIS2dir.xyz <- as.data.frame(MIS2.dir.xyz)
    colnames(MIS2dir.xyz) <- c("x","y","dir")
    MIS2dir.xyz$dirdeg <- ifelse(MIS2dir.xyz$dir == "E", 0, NA)
    MIS2dir.xyz$dirdeg <- ifelse(MIS2dir.xyz$dir == "SE", 45, MIS2dir.xyz$dirdeg)
    MIS2dir.xyz$dirdeg <- ifelse(MIS2dir.xyz$dir == "S", 90, MIS2dir.xyz$dirdeg)
    MIS2dir.xyz$dirdeg <- ifelse(MIS2dir.xyz$dir == "SW", 135, MIS2dir.xyz$dirdeg)
    MIS2dir.xyz$dirdeg <- ifelse(MIS2dir.xyz$dir == "W", 180, MIS2dir.xyz$dirdeg)
    MIS2dir.xyz$dirdeg <- ifelse(MIS2dir.xyz$dir == "NW", 235, MIS2dir.xyz$dirdeg)
    MIS2dir.xyz$dirdeg <- ifelse(MIS2dir.xyz$dir == "N", 270, MIS2dir.xyz$dirdeg)
    MIS2dir.xyz$dirdeg <- ifelse(MIS2dir.xyz$dir == "NE", 315, MIS2dir.xyz$dirdeg)
    
    MIS2dir.dat <- data.frame(x=unlist(MIS2dir.xyz$x), y=unlist(MIS2dir.xyz$y), dirdeg=unlist(MIS2dir.xyz$dirdeg))
    MIS2dir.rast <- rasterFromXYZ(MIS2dir.dat, crs=CRS("+proj=longlat +datum=WGS84"))
    writeRaster(MIS2dir.rast, filename="MIS2dir.grd", format="raster", overwrite=T)
    
    LGM.dir.xyz <- LGM.dir.xyz[-1,]
    LGMdir.xyz <- as.data.frame(LGM.dir.xyz)
    colnames(LGMdir.xyz) <- c("x","y","dir")
    LGMdir.xyz$dirdeg <- ifelse(LGMdir.xyz$dir == "E", 0, NA)
    LGMdir.xyz$dirdeg <- ifelse(LGMdir.xyz$dir == "SE", 45, LGMdir.xyz$dirdeg)
    LGMdir.xyz$dirdeg <- ifelse(LGMdir.xyz$dir == "S", 90, LGMdir.xyz$dirdeg)
    LGMdir.xyz$dirdeg <- ifelse(LGMdir.xyz$dir == "SW", 135, LGMdir.xyz$dirdeg)
    LGMdir.xyz$dirdeg <- ifelse(LGMdir.xyz$dir == "W", 180, LGMdir.xyz$dirdeg)
    LGMdir.xyz$dirdeg <- ifelse(LGMdir.xyz$dir == "NW", 235, LGMdir.xyz$dirdeg)
    LGMdir.xyz$dirdeg <- ifelse(LGMdir.xyz$dir == "N", 270, LGMdir.xyz$dirdeg)
    LGMdir.xyz$dirdeg <- ifelse(LGMdir.xyz$dir == "NE", 315, LGMdir.xyz$dirdeg)
    
    LGMdir.dat <- data.frame(x=unlist(LGMdir.xyz$x), y=unlist(LGMdir.xyz$y), dirdeg=unlist(LGMdir.xyz$dirdeg))
    LGMdir.rast <- rasterFromXYZ(LGMdir.dat, crs=CRS("+proj=longlat +datum=WGS84"))
    writeRaster(LGMdir.rast, filename="LGMdir.grd", format="raster", overwrite=T)
    
    LGMc.dir.xyz <- LGMc.dir.xyz[-1,]
    LGMcdir.xyz <- as.data.frame(LGMc.dir.xyz)
    colnames(LGMcdir.xyz) <- c("x","y","dir")
    LGMcdir.xyz$dirdeg <- ifelse(LGMcdir.xyz$dir == "E", 0, NA)
    LGMcdir.xyz$dirdeg <- ifelse(LGMcdir.xyz$dir == "SE", 45, LGMcdir.xyz$dirdeg)
    LGMcdir.xyz$dirdeg <- ifelse(LGMcdir.xyz$dir == "S", 90, LGMcdir.xyz$dirdeg)
    LGMcdir.xyz$dirdeg <- ifelse(LGMcdir.xyz$dir == "SW", 135, LGMcdir.xyz$dirdeg)
    LGMcdir.xyz$dirdeg <- ifelse(LGMcdir.xyz$dir == "W", 180, LGMcdir.xyz$dirdeg)
    LGMcdir.xyz$dirdeg <- ifelse(LGMcdir.xyz$dir == "NW", 235, LGMcdir.xyz$dirdeg)
    LGMcdir.xyz$dirdeg <- ifelse(LGMcdir.xyz$dir == "N", 270, LGMcdir.xyz$dirdeg)
    LGMcdir.xyz$dirdeg <- ifelse(LGMcdir.xyz$dir == "NE", 315, LGMcdir.xyz$dirdeg)
    
    LGMcdir.dat <- data.frame(x=unlist(LGMcdir.xyz$x), y=unlist(LGMcdir.xyz$y), dirdeg=unlist(LGMcdir.xyz$dirdeg))
    LGMcdir.rast <- rasterFromXYZ(LGMcdir.dat, crs=CRS("+proj=longlat +datum=WGS84"))
    writeRaster(LGMcdir.rast, filename="LGMcdir.grd", format="raster", overwrite=T)
    
    LGMu.dir.xyz <- LGMu.dir.xyz[-1,]
    LGMudir.xyz <- as.data.frame(LGMu.dir.xyz)
    colnames(LGMudir.xyz) <- c("x","y","dir")
    LGMudir.xyz$dirdeg <- ifelse(LGMudir.xyz$dir == "E", 0, NA)
    LGMudir.xyz$dirdeg <- ifelse(LGMudir.xyz$dir == "SE", 45, LGMudir.xyz$dirdeg)
    LGMudir.xyz$dirdeg <- ifelse(LGMudir.xyz$dir == "S", 90, LGMudir.xyz$dirdeg)
    LGMudir.xyz$dirdeg <- ifelse(LGMudir.xyz$dir == "SW", 135, LGMudir.xyz$dirdeg)
    LGMudir.xyz$dirdeg <- ifelse(LGMudir.xyz$dir == "W", 180, LGMudir.xyz$dirdeg)
    LGMudir.xyz$dirdeg <- ifelse(LGMudir.xyz$dir == "NW", 235, LGMudir.xyz$dirdeg)
    LGMudir.xyz$dirdeg <- ifelse(LGMudir.xyz$dir == "N", 270, LGMudir.xyz$dirdeg)
    LGMudir.xyz$dirdeg <- ifelse(LGMudir.xyz$dir == "NE", 315, LGMudir.xyz$dirdeg)
    
    LGMudir.dat <- data.frame(x=unlist(LGMudir.xyz$x), y=unlist(LGMudir.xyz$y), dirdeg=unlist(LGMudir.xyz$dirdeg))
    LGMudir.rast <- rasterFromXYZ(LGMudir.dat, crs=CRS("+proj=longlat +datum=WGS84"))
    writeRaster(LGMudir.rast, filename="LGMudir.grd", format="raster", overwrite=T)
    
    HOL.dir.xyz <- HOL.dir.xyz[-1,]
    HOLdir.xyz <- as.data.frame(HOL.dir.xyz)
    colnames(HOLdir.xyz) <- c("x","y","dir")
    HOLdir.xyz$dirdeg <- ifelse(HOLdir.xyz$dir == "E", 0, NA)
    HOLdir.xyz$dirdeg <- ifelse(HOLdir.xyz$dir == "SE", 45, HOLdir.xyz$dirdeg)
    HOLdir.xyz$dirdeg <- ifelse(HOLdir.xyz$dir == "S", 90, HOLdir.xyz$dirdeg)
    HOLdir.xyz$dirdeg <- ifelse(HOLdir.xyz$dir == "SW", 135, HOLdir.xyz$dirdeg)
    HOLdir.xyz$dirdeg <- ifelse(HOLdir.xyz$dir == "W", 180, HOLdir.xyz$dirdeg)
    HOLdir.xyz$dirdeg <- ifelse(HOLdir.xyz$dir == "NW", 235, HOLdir.xyz$dirdeg)
    HOLdir.xyz$dirdeg <- ifelse(HOLdir.xyz$dir == "N", 270, HOLdir.xyz$dirdeg)
    HOLdir.xyz$dirdeg <- ifelse(HOLdir.xyz$dir == "NE", 315, HOLdir.xyz$dirdeg)
    
    HOLdir.dat <- data.frame(x=unlist(HOLdir.xyz$x), y=unlist(HOLdir.xyz$y), dirdeg=unlist(HOLdir.xyz$dirdeg))
    HOLdir.rast <- rasterFromXYZ(HOLdir.dat, crs=CRS("+proj=longlat +datum=WGS84"))
    writeRaster(HOLdir.rast, filename="HOLdir.grd", format="raster", overwrite=T)
    
    L26.dir.xyz <- L26.dir.xyz[-1,]
    L26dir.xyz <- as.data.frame(L26.dir.xyz)
    colnames(L26dir.xyz) <- c("x","y","dir")
    L26dir.xyz$dirdeg <- ifelse(L26dir.xyz$dir == "E", 0, NA)
    L26dir.xyz$dirdeg <- ifelse(L26dir.xyz$dir == "SE", 45, L26dir.xyz$dirdeg)
    L26dir.xyz$dirdeg <- ifelse(L26dir.xyz$dir == "S", 90, L26dir.xyz$dirdeg)
    L26dir.xyz$dirdeg <- ifelse(L26dir.xyz$dir == "SW", 135, L26dir.xyz$dirdeg)
    L26dir.xyz$dirdeg <- ifelse(L26dir.xyz$dir == "W", 180, L26dir.xyz$dirdeg)
    L26dir.xyz$dirdeg <- ifelse(L26dir.xyz$dir == "NW", 235, L26dir.xyz$dirdeg)
    L26dir.xyz$dirdeg <- ifelse(L26dir.xyz$dir == "N", 270, L26dir.xyz$dirdeg)
    L26dir.xyz$dirdeg <- ifelse(L26dir.xyz$dir == "NE", 315, L26dir.xyz$dirdeg)
    
    L26dir.dat <- data.frame(x=unlist(L26dir.xyz$x), y=unlist(L26dir.xyz$y), dirdeg=unlist(L26dir.xyz$dirdeg))
    L26dir.rast <- rasterFromXYZ(L26dir.dat, crs=CRS("+proj=longlat +datum=WGS84"))
    writeRaster(L26dir.rast, filename="L26dir.grd", format="raster", overwrite=T)
    
    L25.dir.xyz <- L25.dir.xyz[-1,]
    L25dir.xyz <- as.data.frame(L25.dir.xyz)
    colnames(L25dir.xyz) <- c("x","y","dir")
    L25dir.xyz$dirdeg <- ifelse(L25dir.xyz$dir == "E", 0, NA)
    L25dir.xyz$dirdeg <- ifelse(L25dir.xyz$dir == "SE", 45, L25dir.xyz$dirdeg)
    L25dir.xyz$dirdeg <- ifelse(L25dir.xyz$dir == "S", 90, L25dir.xyz$dirdeg)
    L25dir.xyz$dirdeg <- ifelse(L25dir.xyz$dir == "SW", 135, L25dir.xyz$dirdeg)
    L25dir.xyz$dirdeg <- ifelse(L25dir.xyz$dir == "W", 180, L25dir.xyz$dirdeg)
    L25dir.xyz$dirdeg <- ifelse(L25dir.xyz$dir == "NW", 235, L25dir.xyz$dirdeg)
    L25dir.xyz$dirdeg <- ifelse(L25dir.xyz$dir == "N", 270, L25dir.xyz$dirdeg)
    L25dir.xyz$dirdeg <- ifelse(L25dir.xyz$dir == "NE", 315, L25dir.xyz$dirdeg)
    
    L25dir.dat <- data.frame(x=unlist(L25dir.xyz$x), y=unlist(L25dir.xyz$y), dirdeg=unlist(L25dir.xyz$dirdeg))
    L25dir.rast <- rasterFromXYZ(L25dir.dat, crs=CRS("+proj=longlat +datum=WGS84"))
    writeRaster(L25dir.rast, filename="L25dir.grd", format="raster", overwrite=T)
    
    L24.dir.xyz <- L24.dir.xyz[-1,]
    L24dir.xyz <- as.data.frame(L24.dir.xyz)
    colnames(L24dir.xyz) <- c("x","y","dir")
    L24dir.xyz$dirdeg <- ifelse(L24dir.xyz$dir == "E", 0, NA)
    L24dir.xyz$dirdeg <- ifelse(L24dir.xyz$dir == "SE", 45, L24dir.xyz$dirdeg)
    L24dir.xyz$dirdeg <- ifelse(L24dir.xyz$dir == "S", 90, L24dir.xyz$dirdeg)
    L24dir.xyz$dirdeg <- ifelse(L24dir.xyz$dir == "SW", 135, L24dir.xyz$dirdeg)
    L24dir.xyz$dirdeg <- ifelse(L24dir.xyz$dir == "W", 180, L24dir.xyz$dirdeg)
    L24dir.xyz$dirdeg <- ifelse(L24dir.xyz$dir == "NW", 235, L24dir.xyz$dirdeg)
    L24dir.xyz$dirdeg <- ifelse(L24dir.xyz$dir == "N", 270, L24dir.xyz$dirdeg)
    L24dir.xyz$dirdeg <- ifelse(L24dir.xyz$dir == "NE", 315, L24dir.xyz$dirdeg)
    
    L24dir.dat <- data.frame(x=unlist(L24dir.xyz$x), y=unlist(L24dir.xyz$y), dirdeg=unlist(L24dir.xyz$dirdeg))
    L24dir.rast <- rasterFromXYZ(L24dir.dat, crs=CRS("+proj=longlat +datum=WGS84"))
    writeRaster(L24dir.rast, filename="L24dir.grd", format="raster", overwrite=T)
    
    L23.dir.xyz <- L23.dir.xyz[-1,]
    L23dir.xyz <- as.data.frame(L23.dir.xyz)
    colnames(L23dir.xyz) <- c("x","y","dir")
    L23dir.xyz$dirdeg <- ifelse(L23dir.xyz$dir == "E", 0, NA)
    L23dir.xyz$dirdeg <- ifelse(L23dir.xyz$dir == "SE", 45, L23dir.xyz$dirdeg)
    L23dir.xyz$dirdeg <- ifelse(L23dir.xyz$dir == "S", 90, L23dir.xyz$dirdeg)
    L23dir.xyz$dirdeg <- ifelse(L23dir.xyz$dir == "SW", 135, L23dir.xyz$dirdeg)
    L23dir.xyz$dirdeg <- ifelse(L23dir.xyz$dir == "W", 180, L23dir.xyz$dirdeg)
    L23dir.xyz$dirdeg <- ifelse(L23dir.xyz$dir == "NW", 235, L23dir.xyz$dirdeg)
    L23dir.xyz$dirdeg <- ifelse(L23dir.xyz$dir == "N", 270, L23dir.xyz$dirdeg)
    L23dir.xyz$dirdeg <- ifelse(L23dir.xyz$dir == "NE", 315, L23dir.xyz$dirdeg)
    
    L23dir.dat <- data.frame(x=unlist(L23dir.xyz$x), y=unlist(L23dir.xyz$y), dirdeg=unlist(L23dir.xyz$dirdeg))
    L23dir.rast <- rasterFromXYZ(L23dir.dat, crs=CRS("+proj=longlat +datum=WGS84"))
    writeRaster(L23dir.rast, filename="L23dir.grd", format="raster", overwrite=T)
    
    L22.dir.xyz <- L22.dir.xyz[-1,]
    L22dir.xyz <- as.data.frame(L22.dir.xyz)
    colnames(L22dir.xyz) <- c("x","y","dir")
    L22dir.xyz$dirdeg <- ifelse(L22dir.xyz$dir == "E", 0, NA)
    L22dir.xyz$dirdeg <- ifelse(L22dir.xyz$dir == "SE", 45, L22dir.xyz$dirdeg)
    L22dir.xyz$dirdeg <- ifelse(L22dir.xyz$dir == "S", 90, L22dir.xyz$dirdeg)
    L22dir.xyz$dirdeg <- ifelse(L22dir.xyz$dir == "SW", 135, L22dir.xyz$dirdeg)
    L22dir.xyz$dirdeg <- ifelse(L22dir.xyz$dir == "W", 180, L22dir.xyz$dirdeg)
    L22dir.xyz$dirdeg <- ifelse(L22dir.xyz$dir == "NW", 235, L22dir.xyz$dirdeg)
    L22dir.xyz$dirdeg <- ifelse(L22dir.xyz$dir == "N", 270, L22dir.xyz$dirdeg)
    L22dir.xyz$dirdeg <- ifelse(L22dir.xyz$dir == "NE", 315, L22dir.xyz$dirdeg)
    
    L22dir.dat <- data.frame(x=unlist(L22dir.xyz$x), y=unlist(L22dir.xyz$y), dirdeg=unlist(L22dir.xyz$dirdeg))
    L22dir.rast <- rasterFromXYZ(L22dir.dat, crs=CRS("+proj=longlat +datum=WGS84"))
    writeRaster(L22dir.rast, filename="L22dir.grd", format="raster", overwrite=T)
    
    L21.dir.xyz <- L21.dir.xyz[-1,]
    L21dir.xyz <- as.data.frame(L21.dir.xyz)
    colnames(L21dir.xyz) <- c("x","y","dir")
    L21dir.xyz$dirdeg <- ifelse(L21dir.xyz$dir == "E", 0, NA)
    L21dir.xyz$dirdeg <- ifelse(L21dir.xyz$dir == "SE", 45, L21dir.xyz$dirdeg)
    L21dir.xyz$dirdeg <- ifelse(L21dir.xyz$dir == "S", 90, L21dir.xyz$dirdeg)
    L21dir.xyz$dirdeg <- ifelse(L21dir.xyz$dir == "SW", 135, L21dir.xyz$dirdeg)
    L21dir.xyz$dirdeg <- ifelse(L21dir.xyz$dir == "W", 180, L21dir.xyz$dirdeg)
    L21dir.xyz$dirdeg <- ifelse(L21dir.xyz$dir == "NW", 235, L21dir.xyz$dirdeg)
    L21dir.xyz$dirdeg <- ifelse(L21dir.xyz$dir == "N", 270, L21dir.xyz$dirdeg)
    L21dir.xyz$dirdeg <- ifelse(L21dir.xyz$dir == "NE", 315, L21dir.xyz$dirdeg)
    
    L21dir.dat <- data.frame(x=unlist(L21dir.xyz$x), y=unlist(L21dir.xyz$y), dirdeg=unlist(L21dir.xyz$dirdeg))
    L21dir.rast <- rasterFromXYZ(L21dir.dat, crs=CRS("+proj=longlat +datum=WGS84"))
    writeRaster(L21dir.rast, filename="L21dir.grd", format="raster", overwrite=T)
    
    L20.dir.xyz <- L20.dir.xyz[-1,]
    L20dir.xyz <- as.data.frame(L20.dir.xyz)
    colnames(L20dir.xyz) <- c("x","y","dir")
    L20dir.xyz$dirdeg <- ifelse(L20dir.xyz$dir == "E", 0, NA)
    L20dir.xyz$dirdeg <- ifelse(L20dir.xyz$dir == "SE", 45, L20dir.xyz$dirdeg)
    L20dir.xyz$dirdeg <- ifelse(L20dir.xyz$dir == "S", 90, L20dir.xyz$dirdeg)
    L20dir.xyz$dirdeg <- ifelse(L20dir.xyz$dir == "SW", 135, L20dir.xyz$dirdeg)
    L20dir.xyz$dirdeg <- ifelse(L20dir.xyz$dir == "W", 180, L20dir.xyz$dirdeg)
    L20dir.xyz$dirdeg <- ifelse(L20dir.xyz$dir == "NW", 235, L20dir.xyz$dirdeg)
    L20dir.xyz$dirdeg <- ifelse(L20dir.xyz$dir == "N", 270, L20dir.xyz$dirdeg)
    L20dir.xyz$dirdeg <- ifelse(L20dir.xyz$dir == "NE", 315, L20dir.xyz$dirdeg)
    
    L20dir.dat <- data.frame(x=unlist(L20dir.xyz$x), y=unlist(L20dir.xyz$y), dirdeg=unlist(L20dir.xyz$dirdeg))
    L20dir.rast <- rasterFromXYZ(L20dir.dat, crs=CRS("+proj=longlat +datum=WGS84"))
    writeRaster(L20dir.rast, filename="L20dir.grd", format="raster", overwrite=T)
    
    L19.dir.xyz <- L19.dir.xyz[-1,]
    L19dir.xyz <- as.data.frame(L19.dir.xyz)
    colnames(L19dir.xyz) <- c("x","y","dir")
    L19dir.xyz$dirdeg <- ifelse(L19dir.xyz$dir == "E", 0, NA)
    L19dir.xyz$dirdeg <- ifelse(L19dir.xyz$dir == "SE", 45, L19dir.xyz$dirdeg)
    L19dir.xyz$dirdeg <- ifelse(L19dir.xyz$dir == "S", 90, L19dir.xyz$dirdeg)
    L19dir.xyz$dirdeg <- ifelse(L19dir.xyz$dir == "SW", 135, L19dir.xyz$dirdeg)
    L19dir.xyz$dirdeg <- ifelse(L19dir.xyz$dir == "W", 180, L19dir.xyz$dirdeg)
    L19dir.xyz$dirdeg <- ifelse(L19dir.xyz$dir == "NW", 235, L19dir.xyz$dirdeg)
    L19dir.xyz$dirdeg <- ifelse(L19dir.xyz$dir == "N", 270, L19dir.xyz$dirdeg)
    L19dir.xyz$dirdeg <- ifelse(L19dir.xyz$dir == "NE", 315, L19dir.xyz$dirdeg)
    
    L19dir.dat <- data.frame(x=unlist(L19dir.xyz$x), y=unlist(L19dir.xyz$y), dirdeg=unlist(L19dir.xyz$dirdeg))
    L19dir.rast <- rasterFromXYZ(L19dir.dat, crs=CRS("+proj=longlat +datum=WGS84"))
    writeRaster(L19dir.rast, filename="L19dir.grd", format="raster", overwrite=T)
    
    L18.dir.xyz <- L18.dir.xyz[-1,]
    L18dir.xyz <- as.data.frame(L18.dir.xyz)
    colnames(L18dir.xyz) <- c("x","y","dir")
    L18dir.xyz$dirdeg <- ifelse(L18dir.xyz$dir == "E", 0, NA)
    L18dir.xyz$dirdeg <- ifelse(L18dir.xyz$dir == "SE", 45, L18dir.xyz$dirdeg)
    L18dir.xyz$dirdeg <- ifelse(L18dir.xyz$dir == "S", 90, L18dir.xyz$dirdeg)
    L18dir.xyz$dirdeg <- ifelse(L18dir.xyz$dir == "SW", 135, L18dir.xyz$dirdeg)
    L18dir.xyz$dirdeg <- ifelse(L18dir.xyz$dir == "W", 180, L18dir.xyz$dirdeg)
    L18dir.xyz$dirdeg <- ifelse(L18dir.xyz$dir == "NW", 235, L18dir.xyz$dirdeg)
    L18dir.xyz$dirdeg <- ifelse(L18dir.xyz$dir == "N", 270, L18dir.xyz$dirdeg)
    L18dir.xyz$dirdeg <- ifelse(L18dir.xyz$dir == "NE", 315, L18dir.xyz$dirdeg)
    
    L18dir.dat <- data.frame(x=unlist(L18dir.xyz$x), y=unlist(L18dir.xyz$y), dirdeg=unlist(L18dir.xyz$dirdeg))
    L18dir.rast <- rasterFromXYZ(L18dir.dat, crs=CRS("+proj=longlat +datum=WGS84"))
    writeRaster(L18dir.rast, filename="L18dir.grd", format="raster", overwrite=T)
    
    L17.dir.xyz <- L17.dir.xyz[-1,]
    L17dir.xyz <- as.data.frame(L17.dir.xyz)
    colnames(L17dir.xyz) <- c("x","y","dir")
    L17dir.xyz$dirdeg <- ifelse(L17dir.xyz$dir == "E", 0, NA)
    L17dir.xyz$dirdeg <- ifelse(L17dir.xyz$dir == "SE", 45, L17dir.xyz$dirdeg)
    L17dir.xyz$dirdeg <- ifelse(L17dir.xyz$dir == "S", 90, L17dir.xyz$dirdeg)
    L17dir.xyz$dirdeg <- ifelse(L17dir.xyz$dir == "SW", 135, L17dir.xyz$dirdeg)
    L17dir.xyz$dirdeg <- ifelse(L17dir.xyz$dir == "W", 180, L17dir.xyz$dirdeg)
    L17dir.xyz$dirdeg <- ifelse(L17dir.xyz$dir == "NW", 235, L17dir.xyz$dirdeg)
    L17dir.xyz$dirdeg <- ifelse(L17dir.xyz$dir == "N", 270, L17dir.xyz$dirdeg)
    L17dir.xyz$dirdeg <- ifelse(L17dir.xyz$dir == "NE", 315, L17dir.xyz$dirdeg)
    
    L17dir.dat <- data.frame(x=unlist(L17dir.xyz$x), y=unlist(L17dir.xyz$y), dirdeg=unlist(L17dir.xyz$dirdeg))
    L17dir.rast <- rasterFromXYZ(L17dir.dat, crs=CRS("+proj=longlat +datum=WGS84"))
    writeRaster(L17dir.rast, filename="L17dir.grd", format="raster", overwrite=T)
    
    L16.dir.xyz <- L16.dir.xyz[-1,]
    L16dir.xyz <- as.data.frame(L16.dir.xyz)
    colnames(L16dir.xyz) <- c("x","y","dir")
    L16dir.xyz$dirdeg <- ifelse(L16dir.xyz$dir == "E", 0, NA)
    L16dir.xyz$dirdeg <- ifelse(L16dir.xyz$dir == "SE", 45, L16dir.xyz$dirdeg)
    L16dir.xyz$dirdeg <- ifelse(L16dir.xyz$dir == "S", 90, L16dir.xyz$dirdeg)
    L16dir.xyz$dirdeg <- ifelse(L16dir.xyz$dir == "SW", 135, L16dir.xyz$dirdeg)
    L16dir.xyz$dirdeg <- ifelse(L16dir.xyz$dir == "W", 180, L16dir.xyz$dirdeg)
    L16dir.xyz$dirdeg <- ifelse(L16dir.xyz$dir == "NW", 235, L16dir.xyz$dirdeg)
    L16dir.xyz$dirdeg <- ifelse(L16dir.xyz$dir == "N", 270, L16dir.xyz$dirdeg)
    L16dir.xyz$dirdeg <- ifelse(L16dir.xyz$dir == "NE", 315, L16dir.xyz$dirdeg)
    
    L16dir.dat <- data.frame(x=unlist(L16dir.xyz$x), y=unlist(L16dir.xyz$y), dirdeg=unlist(L16dir.xyz$dirdeg))
    L16dir.rast <- rasterFromXYZ(L16dir.dat, crs=CRS("+proj=longlat +datum=WGS84"))
    writeRaster(L16dir.rast, filename="L16dir.grd", format="raster", overwrite=T)
    
    L15.dir.xyz <- L15.dir.xyz[-1,]
    L15dir.xyz <- as.data.frame(L15.dir.xyz)
    colnames(L15dir.xyz) <- c("x","y","dir")
    L15dir.xyz$dirdeg <- ifelse(L15dir.xyz$dir == "E", 0, NA)
    L15dir.xyz$dirdeg <- ifelse(L15dir.xyz$dir == "SE", 45, L15dir.xyz$dirdeg)
    L15dir.xyz$dirdeg <- ifelse(L15dir.xyz$dir == "S", 90, L15dir.xyz$dirdeg)
    L15dir.xyz$dirdeg <- ifelse(L15dir.xyz$dir == "SW", 135, L15dir.xyz$dirdeg)
    L15dir.xyz$dirdeg <- ifelse(L15dir.xyz$dir == "W", 180, L15dir.xyz$dirdeg)
    L15dir.xyz$dirdeg <- ifelse(L15dir.xyz$dir == "NW", 235, L15dir.xyz$dirdeg)
    L15dir.xyz$dirdeg <- ifelse(L15dir.xyz$dir == "N", 270, L15dir.xyz$dirdeg)
    L15dir.xyz$dirdeg <- ifelse(L15dir.xyz$dir == "NE", 315, L15dir.xyz$dirdeg)
    
    L15dir.dat <- data.frame(x=unlist(L15dir.xyz$x), y=unlist(L15dir.xyz$y), dirdeg=unlist(L15dir.xyz$dirdeg))
    L15dir.rast <- rasterFromXYZ(L15dir.dat, crs=CRS("+proj=longlat +datum=WGS84"))
    writeRaster(L15dir.rast, filename="L15dir.grd", format="raster", overwrite=T)
    
    L14.dir.xyz <- L14.dir.xyz[-1,]
    L14dir.xyz <- as.data.frame(L14.dir.xyz)
    colnames(L14dir.xyz) <- c("x","y","dir")
    L14dir.xyz$dirdeg <- ifelse(L14dir.xyz$dir == "E", 0, NA)
    L14dir.xyz$dirdeg <- ifelse(L14dir.xyz$dir == "SE", 45, L14dir.xyz$dirdeg)
    L14dir.xyz$dirdeg <- ifelse(L14dir.xyz$dir == "S", 90, L14dir.xyz$dirdeg)
    L14dir.xyz$dirdeg <- ifelse(L14dir.xyz$dir == "SW", 135, L14dir.xyz$dirdeg)
    L14dir.xyz$dirdeg <- ifelse(L14dir.xyz$dir == "W", 180, L14dir.xyz$dirdeg)
    L14dir.xyz$dirdeg <- ifelse(L14dir.xyz$dir == "NW", 235, L14dir.xyz$dirdeg)
    L14dir.xyz$dirdeg <- ifelse(L14dir.xyz$dir == "N", 270, L14dir.xyz$dirdeg)
    L14dir.xyz$dirdeg <- ifelse(L14dir.xyz$dir == "NE", 315, L14dir.xyz$dirdeg)
    
    L14dir.dat <- data.frame(x=unlist(L14dir.xyz$x), y=unlist(L14dir.xyz$y), dirdeg=unlist(L14dir.xyz$dirdeg))
    L14dir.rast <- rasterFromXYZ(L14dir.dat, crs=CRS("+proj=longlat +datum=WGS84"))
    writeRaster(L14dir.rast, filename="L14dir.grd", format="raster", overwrite=T)
    
    L13.dir.xyz <- L13.dir.xyz[-1,]
    L13dir.xyz <- as.data.frame(L13.dir.xyz)
    colnames(L13dir.xyz) <- c("x","y","dir")
    L13dir.xyz$dirdeg <- ifelse(L13dir.xyz$dir == "E", 0, NA)
    L13dir.xyz$dirdeg <- ifelse(L13dir.xyz$dir == "SE", 45, L13dir.xyz$dirdeg)
    L13dir.xyz$dirdeg <- ifelse(L13dir.xyz$dir == "S", 90, L13dir.xyz$dirdeg)
    L13dir.xyz$dirdeg <- ifelse(L13dir.xyz$dir == "SW", 135, L13dir.xyz$dirdeg)
    L13dir.xyz$dirdeg <- ifelse(L13dir.xyz$dir == "W", 180, L13dir.xyz$dirdeg)
    L13dir.xyz$dirdeg <- ifelse(L13dir.xyz$dir == "NW", 235, L13dir.xyz$dirdeg)
    L13dir.xyz$dirdeg <- ifelse(L13dir.xyz$dir == "N", 270, L13dir.xyz$dirdeg)
    L13dir.xyz$dirdeg <- ifelse(L13dir.xyz$dir == "NE", 315, L13dir.xyz$dirdeg)
    
    L13dir.dat <- data.frame(x=unlist(L13dir.xyz$x), y=unlist(L13dir.xyz$y), dirdeg=unlist(L13dir.xyz$dirdeg))
    L13dir.rast <- rasterFromXYZ(L13dir.dat, crs=CRS("+proj=longlat +datum=WGS84"))
    writeRaster(L13dir.rast, filename="L13dir.grd", format="raster", overwrite=T)
    
    L12.dir.xyz <- L12.dir.xyz[-1,]
    L12dir.xyz <- as.data.frame(L12.dir.xyz)
    colnames(L12dir.xyz) <- c("x","y","dir")
    L12dir.xyz$dirdeg <- ifelse(L12dir.xyz$dir == "E", 0, NA)
    L12dir.xyz$dirdeg <- ifelse(L12dir.xyz$dir == "SE", 45, L12dir.xyz$dirdeg)
    L12dir.xyz$dirdeg <- ifelse(L12dir.xyz$dir == "S", 90, L12dir.xyz$dirdeg)
    L12dir.xyz$dirdeg <- ifelse(L12dir.xyz$dir == "SW", 135, L12dir.xyz$dirdeg)
    L12dir.xyz$dirdeg <- ifelse(L12dir.xyz$dir == "W", 180, L12dir.xyz$dirdeg)
    L12dir.xyz$dirdeg <- ifelse(L12dir.xyz$dir == "NW", 235, L12dir.xyz$dirdeg)
    L12dir.xyz$dirdeg <- ifelse(L12dir.xyz$dir == "N", 270, L12dir.xyz$dirdeg)
    L12dir.xyz$dirdeg <- ifelse(L12dir.xyz$dir == "NE", 315, L12dir.xyz$dirdeg)
    
    L12dir.dat <- data.frame(x=unlist(L12dir.xyz$x), y=unlist(L12dir.xyz$y), dirdeg=unlist(L12dir.xyz$dirdeg))
    L12dir.rast <- rasterFromXYZ(L12dir.dat, crs=CRS("+proj=longlat +datum=WGS84"))
    writeRaster(L12dir.rast, filename="L12dir.grd", format="raster", overwrite=T)
    
    L11.dir.xyz <- L11.dir.xyz[-1,]
    L11dir.xyz <- as.data.frame(L11.dir.xyz)
    colnames(L11dir.xyz) <- c("x","y","dir")
    L11dir.xyz$dirdeg <- ifelse(L11dir.xyz$dir == "E", 0, NA)
    L11dir.xyz$dirdeg <- ifelse(L11dir.xyz$dir == "SE", 45, L11dir.xyz$dirdeg)
    L11dir.xyz$dirdeg <- ifelse(L11dir.xyz$dir == "S", 90, L11dir.xyz$dirdeg)
    L11dir.xyz$dirdeg <- ifelse(L11dir.xyz$dir == "SW", 135, L11dir.xyz$dirdeg)
    L11dir.xyz$dirdeg <- ifelse(L11dir.xyz$dir == "W", 180, L11dir.xyz$dirdeg)
    L11dir.xyz$dirdeg <- ifelse(L11dir.xyz$dir == "NW", 235, L11dir.xyz$dirdeg)
    L11dir.xyz$dirdeg <- ifelse(L11dir.xyz$dir == "N", 270, L11dir.xyz$dirdeg)
    L11dir.xyz$dirdeg <- ifelse(L11dir.xyz$dir == "NE", 315, L11dir.xyz$dirdeg)
    
    L11dir.dat <- data.frame(x=unlist(L11dir.xyz$x), y=unlist(L11dir.xyz$y), dirdeg=unlist(L11dir.xyz$dirdeg))
    L11dir.rast <- rasterFromXYZ(L11dir.dat, crs=CRS("+proj=longlat +datum=WGS84"))
    writeRaster(L11dir.rast, filename="L11dir.grd", format="raster", overwrite=T)
    
    L10.dir.xyz <- L10.dir.xyz[-1,]
    L10dir.xyz <- as.data.frame(L10.dir.xyz)
    colnames(L10dir.xyz) <- c("x","y","dir")
    L10dir.xyz$dirdeg <- ifelse(L10dir.xyz$dir == "E", 0, NA)
    L10dir.xyz$dirdeg <- ifelse(L10dir.xyz$dir == "SE", 45, L10dir.xyz$dirdeg)
    L10dir.xyz$dirdeg <- ifelse(L10dir.xyz$dir == "S", 90, L10dir.xyz$dirdeg)
    L10dir.xyz$dirdeg <- ifelse(L10dir.xyz$dir == "SW", 135, L10dir.xyz$dirdeg)
    L10dir.xyz$dirdeg <- ifelse(L10dir.xyz$dir == "W", 180, L10dir.xyz$dirdeg)
    L10dir.xyz$dirdeg <- ifelse(L10dir.xyz$dir == "NW", 235, L10dir.xyz$dirdeg)
    L10dir.xyz$dirdeg <- ifelse(L10dir.xyz$dir == "N", 270, L10dir.xyz$dirdeg)
    L10dir.xyz$dirdeg <- ifelse(L10dir.xyz$dir == "NE", 315, L10dir.xyz$dirdeg)
    
    L10dir.dat <- data.frame(x=unlist(L10dir.xyz$x), y=unlist(L10dir.xyz$y), dirdeg=unlist(L10dir.xyz$dirdeg))
    L10dir.rast <- rasterFromXYZ(L10dir.dat, crs=CRS("+proj=longlat +datum=WGS84"))
    writeRaster(L10dir.rast, filename="L10dir.grd", format="raster", overwrite=T)
    
    L09.dir.xyz <- L09.dir.xyz[-1,]
    L09dir.xyz <- as.data.frame(L09.dir.xyz)
    colnames(L09dir.xyz) <- c("x","y","dir")
    L09dir.xyz$dirdeg <- ifelse(L09dir.xyz$dir == "E", 0, NA)
    L09dir.xyz$dirdeg <- ifelse(L09dir.xyz$dir == "SE", 45, L09dir.xyz$dirdeg)
    L09dir.xyz$dirdeg <- ifelse(L09dir.xyz$dir == "S", 90, L09dir.xyz$dirdeg)
    L09dir.xyz$dirdeg <- ifelse(L09dir.xyz$dir == "SW", 135, L09dir.xyz$dirdeg)
    L09dir.xyz$dirdeg <- ifelse(L09dir.xyz$dir == "W", 180, L09dir.xyz$dirdeg)
    L09dir.xyz$dirdeg <- ifelse(L09dir.xyz$dir == "NW", 235, L09dir.xyz$dirdeg)
    L09dir.xyz$dirdeg <- ifelse(L09dir.xyz$dir == "N", 270, L09dir.xyz$dirdeg)
    L09dir.xyz$dirdeg <- ifelse(L09dir.xyz$dir == "NE", 315, L09dir.xyz$dirdeg)
    
    L09dir.dat <- data.frame(x=unlist(L09dir.xyz$x), y=unlist(L09dir.xyz$y), dirdeg=unlist(L09dir.xyz$dirdeg))
    L09dir.rast <- rasterFromXYZ(L09dir.dat, crs=CRS("+proj=longlat +datum=WGS84"))
    writeRaster(L09dir.rast, filename="L09dir.grd", format="raster", overwrite=T)
    
    L08.dir.xyz <- L08.dir.xyz[-1,]
    L08dir.xyz <- as.data.frame(L08.dir.xyz)
    colnames(L08dir.xyz) <- c("x","y","dir")
    L08dir.xyz$dirdeg <- ifelse(L08dir.xyz$dir == "E", 0, NA)
    L08dir.xyz$dirdeg <- ifelse(L08dir.xyz$dir == "SE", 45, L08dir.xyz$dirdeg)
    L08dir.xyz$dirdeg <- ifelse(L08dir.xyz$dir == "S", 90, L08dir.xyz$dirdeg)
    L08dir.xyz$dirdeg <- ifelse(L08dir.xyz$dir == "SW", 135, L08dir.xyz$dirdeg)
    L08dir.xyz$dirdeg <- ifelse(L08dir.xyz$dir == "W", 180, L08dir.xyz$dirdeg)
    L08dir.xyz$dirdeg <- ifelse(L08dir.xyz$dir == "NW", 235, L08dir.xyz$dirdeg)
    L08dir.xyz$dirdeg <- ifelse(L08dir.xyz$dir == "N", 270, L08dir.xyz$dirdeg)
    L08dir.xyz$dirdeg <- ifelse(L08dir.xyz$dir == "NE", 315, L08dir.xyz$dirdeg)
    
    L08dir.dat <- data.frame(x=unlist(L08dir.xyz$x), y=unlist(L08dir.xyz$y), dirdeg=unlist(L08dir.xyz$dirdeg))
    L08dir.rast <- rasterFromXYZ(L08dir.dat, crs=CRS("+proj=longlat +datum=WGS84"))
    writeRaster(L08dir.rast, filename="L08dir.grd", format="raster", overwrite=T)
    
    
    L2624.dir.xyz <- L2624.dir.xyz[-1,]
    L2624dir.xyz <- as.data.frame(L2624.dir.xyz)
    colnames(L2624dir.xyz) <- c("x","y","dir")
    L2624dir.xyz$dirdeg <- ifelse(L2624dir.xyz$dir == "E", 0, NA)
    L2624dir.xyz$dirdeg <- ifelse(L2624dir.xyz$dir == "SE", 45, L2624dir.xyz$dirdeg)
    L2624dir.xyz$dirdeg <- ifelse(L2624dir.xyz$dir == "S", 90, L2624dir.xyz$dirdeg)
    L2624dir.xyz$dirdeg <- ifelse(L2624dir.xyz$dir == "SW", 135, L2624dir.xyz$dirdeg)
    L2624dir.xyz$dirdeg <- ifelse(L2624dir.xyz$dir == "W", 180, L2624dir.xyz$dirdeg)
    L2624dir.xyz$dirdeg <- ifelse(L2624dir.xyz$dir == "NW", 235, L2624dir.xyz$dirdeg)
    L2624dir.xyz$dirdeg <- ifelse(L2624dir.xyz$dir == "N", 270, L2624dir.xyz$dirdeg)
    L2624dir.xyz$dirdeg <- ifelse(L2624dir.xyz$dir == "NE", 315, L2624dir.xyz$dirdeg)
    
    L2624dir.dat <- data.frame(x=unlist(L2624dir.xyz$x), y=unlist(L2624dir.xyz$y), dirdeg=unlist(L2624dir.xyz$dirdeg))
    L2624dir.rast <- rasterFromXYZ(L2624dir.dat, crs=CRS("+proj=longlat +datum=WGS84"))
    writeRaster(L2624dir.rast, filename="L2624dir.grd", format="raster", overwrite=T)
    
    L2422.dir.xyz <- L2422.dir.xyz[-1,]
    L2422dir.xyz <- as.data.frame(L2422.dir.xyz)
    colnames(L2422dir.xyz) <- c("x","y","dir")
    L2422dir.xyz$dirdeg <- ifelse(L2422dir.xyz$dir == "E", 0, NA)
    L2422dir.xyz$dirdeg <- ifelse(L2422dir.xyz$dir == "SE", 45, L2422dir.xyz$dirdeg)
    L2422dir.xyz$dirdeg <- ifelse(L2422dir.xyz$dir == "S", 90, L2422dir.xyz$dirdeg)
    L2422dir.xyz$dirdeg <- ifelse(L2422dir.xyz$dir == "SW", 135, L2422dir.xyz$dirdeg)
    L2422dir.xyz$dirdeg <- ifelse(L2422dir.xyz$dir == "W", 180, L2422dir.xyz$dirdeg)
    L2422dir.xyz$dirdeg <- ifelse(L2422dir.xyz$dir == "NW", 235, L2422dir.xyz$dirdeg)
    L2422dir.xyz$dirdeg <- ifelse(L2422dir.xyz$dir == "N", 270, L2422dir.xyz$dirdeg)
    L2422dir.xyz$dirdeg <- ifelse(L2422dir.xyz$dir == "NE", 315, L2422dir.xyz$dirdeg)
    
    L2422dir.dat <- data.frame(x=unlist(L2422dir.xyz$x), y=unlist(L2422dir.xyz$y), dirdeg=unlist(L2422dir.xyz$dirdeg))
    L2422dir.rast <- rasterFromXYZ(L2422dir.dat, crs=CRS("+proj=longlat +datum=WGS84"))
    writeRaster(L2422dir.rast, filename="L2422dir.grd", format="raster", overwrite=T)
    
    L2220.dir.xyz <- L2220.dir.xyz[-1,]
    L2220dir.xyz <- as.data.frame(L2220.dir.xyz)
    colnames(L2220dir.xyz) <- c("x","y","dir")
    L2220dir.xyz$dirdeg <- ifelse(L2220dir.xyz$dir == "E", 0, NA)
    L2220dir.xyz$dirdeg <- ifelse(L2220dir.xyz$dir == "SE", 45, L2220dir.xyz$dirdeg)
    L2220dir.xyz$dirdeg <- ifelse(L2220dir.xyz$dir == "S", 90, L2220dir.xyz$dirdeg)
    L2220dir.xyz$dirdeg <- ifelse(L2220dir.xyz$dir == "SW", 135, L2220dir.xyz$dirdeg)
    L2220dir.xyz$dirdeg <- ifelse(L2220dir.xyz$dir == "W", 180, L2220dir.xyz$dirdeg)
    L2220dir.xyz$dirdeg <- ifelse(L2220dir.xyz$dir == "NW", 235, L2220dir.xyz$dirdeg)
    L2220dir.xyz$dirdeg <- ifelse(L2220dir.xyz$dir == "N", 270, L2220dir.xyz$dirdeg)
    L2220dir.xyz$dirdeg <- ifelse(L2220dir.xyz$dir == "NE", 315, L2220dir.xyz$dirdeg)
    
    L2220dir.dat <- data.frame(x=unlist(L2220dir.xyz$x), y=unlist(L2220dir.xyz$y), dirdeg=unlist(L2220dir.xyz$dirdeg))
    L2220dir.rast <- rasterFromXYZ(L2220dir.dat, crs=CRS("+proj=longlat +datum=WGS84"))
    writeRaster(L2220dir.rast, filename="L2220dir.grd", format="raster", overwrite=T)
    
    L2018.dir.xyz <- L2018.dir.xyz[-1,]
    L2018dir.xyz <- as.data.frame(L2018.dir.xyz)
    colnames(L2018dir.xyz) <- c("x","y","dir")
    L2018dir.xyz$dirdeg <- ifelse(L2018dir.xyz$dir == "E", 0, NA)
    L2018dir.xyz$dirdeg <- ifelse(L2018dir.xyz$dir == "SE", 45, L2018dir.xyz$dirdeg)
    L2018dir.xyz$dirdeg <- ifelse(L2018dir.xyz$dir == "S", 90, L2018dir.xyz$dirdeg)
    L2018dir.xyz$dirdeg <- ifelse(L2018dir.xyz$dir == "SW", 135, L2018dir.xyz$dirdeg)
    L2018dir.xyz$dirdeg <- ifelse(L2018dir.xyz$dir == "W", 180, L2018dir.xyz$dirdeg)
    L2018dir.xyz$dirdeg <- ifelse(L2018dir.xyz$dir == "NW", 235, L2018dir.xyz$dirdeg)
    L2018dir.xyz$dirdeg <- ifelse(L2018dir.xyz$dir == "N", 270, L2018dir.xyz$dirdeg)
    L2018dir.xyz$dirdeg <- ifelse(L2018dir.xyz$dir == "NE", 315, L2018dir.xyz$dirdeg)
    
    L2018dir.dat <- data.frame(x=unlist(L2018dir.xyz$x), y=unlist(L2018dir.xyz$y), dirdeg=unlist(L2018dir.xyz$dirdeg))
    L2018dir.rast <- rasterFromXYZ(L2018dir.dat, crs=CRS("+proj=longlat +datum=WGS84"))
    writeRaster(L2018dir.rast, filename="L2018dir.grd", format="raster", overwrite=T)
    
    L1816.dir.xyz <- L1816.dir.xyz[-1,]
    L1816dir.xyz <- as.data.frame(L1816.dir.xyz)
    colnames(L1816dir.xyz) <- c("x","y","dir")
    L1816dir.xyz$dirdeg <- ifelse(L1816dir.xyz$dir == "E", 0, NA)
    L1816dir.xyz$dirdeg <- ifelse(L1816dir.xyz$dir == "SE", 45, L1816dir.xyz$dirdeg)
    L1816dir.xyz$dirdeg <- ifelse(L1816dir.xyz$dir == "S", 90, L1816dir.xyz$dirdeg)
    L1816dir.xyz$dirdeg <- ifelse(L1816dir.xyz$dir == "SW", 135, L1816dir.xyz$dirdeg)
    L1816dir.xyz$dirdeg <- ifelse(L1816dir.xyz$dir == "W", 180, L1816dir.xyz$dirdeg)
    L1816dir.xyz$dirdeg <- ifelse(L1816dir.xyz$dir == "NW", 235, L1816dir.xyz$dirdeg)
    L1816dir.xyz$dirdeg <- ifelse(L1816dir.xyz$dir == "N", 270, L1816dir.xyz$dirdeg)
    L1816dir.xyz$dirdeg <- ifelse(L1816dir.xyz$dir == "NE", 315, L1816dir.xyz$dirdeg)
    
    L1816dir.dat <- data.frame(x=unlist(L1816dir.xyz$x), y=unlist(L1816dir.xyz$y), dirdeg=unlist(L1816dir.xyz$dirdeg))
    L1816dir.rast <- rasterFromXYZ(L1816dir.dat, crs=CRS("+proj=longlat +datum=WGS84"))
    writeRaster(L1816dir.rast, filename="L1816dir.grd", format="raster", overwrite=T)
    
    L1614.dir.xyz <- L1614.dir.xyz[-1,]
    L1614dir.xyz <- as.data.frame(L1614.dir.xyz)
    colnames(L1614dir.xyz) <- c("x","y","dir")
    L1614dir.xyz$dirdeg <- ifelse(L1614dir.xyz$dir == "E", 0, NA)
    L1614dir.xyz$dirdeg <- ifelse(L1614dir.xyz$dir == "SE", 45, L1614dir.xyz$dirdeg)
    L1614dir.xyz$dirdeg <- ifelse(L1614dir.xyz$dir == "S", 90, L1614dir.xyz$dirdeg)
    L1614dir.xyz$dirdeg <- ifelse(L1614dir.xyz$dir == "SW", 135, L1614dir.xyz$dirdeg)
    L1614dir.xyz$dirdeg <- ifelse(L1614dir.xyz$dir == "W", 180, L1614dir.xyz$dirdeg)
    L1614dir.xyz$dirdeg <- ifelse(L1614dir.xyz$dir == "NW", 235, L1614dir.xyz$dirdeg)
    L1614dir.xyz$dirdeg <- ifelse(L1614dir.xyz$dir == "N", 270, L1614dir.xyz$dirdeg)
    L1614dir.xyz$dirdeg <- ifelse(L1614dir.xyz$dir == "NE", 315, L1614dir.xyz$dirdeg)
    
    L1614dir.dat <- data.frame(x=unlist(L1614dir.xyz$x), y=unlist(L1614dir.xyz$y), dirdeg=unlist(L1614dir.xyz$dirdeg))
    L1614dir.rast <- rasterFromXYZ(L1614dir.dat, crs=CRS("+proj=longlat +datum=WGS84"))
    writeRaster(L1614dir.rast, filename="L1614dir.grd", format="raster", overwrite=T)
    
    L1412.dir.xyz <- L1412.dir.xyz[-1,]
    L1412dir.xyz <- as.data.frame(L1412.dir.xyz)
    colnames(L1412dir.xyz) <- c("x","y","dir")
    L1412dir.xyz$dirdeg <- ifelse(L1412dir.xyz$dir == "E", 0, NA)
    L1412dir.xyz$dirdeg <- ifelse(L1412dir.xyz$dir == "SE", 45, L1412dir.xyz$dirdeg)
    L1412dir.xyz$dirdeg <- ifelse(L1412dir.xyz$dir == "S", 90, L1412dir.xyz$dirdeg)
    L1412dir.xyz$dirdeg <- ifelse(L1412dir.xyz$dir == "SW", 135, L1412dir.xyz$dirdeg)
    L1412dir.xyz$dirdeg <- ifelse(L1412dir.xyz$dir == "W", 180, L1412dir.xyz$dirdeg)
    L1412dir.xyz$dirdeg <- ifelse(L1412dir.xyz$dir == "NW", 235, L1412dir.xyz$dirdeg)
    L1412dir.xyz$dirdeg <- ifelse(L1412dir.xyz$dir == "N", 270, L1412dir.xyz$dirdeg)
    L1412dir.xyz$dirdeg <- ifelse(L1412dir.xyz$dir == "NE", 315, L1412dir.xyz$dirdeg)
    
    L1412dir.dat <- data.frame(x=unlist(L1412dir.xyz$x), y=unlist(L1412dir.xyz$y), dirdeg=unlist(L1412dir.xyz$dirdeg))
    L1412dir.rast <- rasterFromXYZ(L1412dir.dat, crs=CRS("+proj=longlat +datum=WGS84"))
    writeRaster(L1412dir.rast, filename="L1412dir.grd", format="raster", overwrite=T)
    
    L1210.dir.xyz <- L1210.dir.xyz[-1,]
    L1210dir.xyz <- as.data.frame(L1210.dir.xyz)
    colnames(L1210dir.xyz) <- c("x","y","dir")
    L1210dir.xyz$dirdeg <- ifelse(L1210dir.xyz$dir == "E", 0, NA)
    L1210dir.xyz$dirdeg <- ifelse(L1210dir.xyz$dir == "SE", 45, L1210dir.xyz$dirdeg)
    L1210dir.xyz$dirdeg <- ifelse(L1210dir.xyz$dir == "S", 90, L1210dir.xyz$dirdeg)
    L1210dir.xyz$dirdeg <- ifelse(L1210dir.xyz$dir == "SW", 135, L1210dir.xyz$dirdeg)
    L1210dir.xyz$dirdeg <- ifelse(L1210dir.xyz$dir == "W", 180, L1210dir.xyz$dirdeg)
    L1210dir.xyz$dirdeg <- ifelse(L1210dir.xyz$dir == "NW", 235, L1210dir.xyz$dirdeg)
    L1210dir.xyz$dirdeg <- ifelse(L1210dir.xyz$dir == "N", 270, L1210dir.xyz$dirdeg)
    L1210dir.xyz$dirdeg <- ifelse(L1210dir.xyz$dir == "NE", 315, L1210dir.xyz$dirdeg)
    
    L1210dir.dat <- data.frame(x=unlist(L1210dir.xyz$x), y=unlist(L1210dir.xyz$y), dirdeg=unlist(L1210dir.xyz$dirdeg))
    L1210dir.rast <- rasterFromXYZ(L1210dir.dat, crs=CRS("+proj=longlat +datum=WGS84"))
    writeRaster(L1210dir.rast, filename="L1210dir.grd", format="raster", overwrite=T)
    
    L1008.dir.xyz <- L1008.dir.xyz[-1,]
    L1008dir.xyz <- as.data.frame(L1008.dir.xyz)
    colnames(L1008dir.xyz) <- c("x","y","dir")
    L1008dir.xyz$dirdeg <- ifelse(L1008dir.xyz$dir == "E", 0, NA)
    L1008dir.xyz$dirdeg <- ifelse(L1008dir.xyz$dir == "SE", 45, L1008dir.xyz$dirdeg)
    L1008dir.xyz$dirdeg <- ifelse(L1008dir.xyz$dir == "S", 90, L1008dir.xyz$dirdeg)
    L1008dir.xyz$dirdeg <- ifelse(L1008dir.xyz$dir == "SW", 135, L1008dir.xyz$dirdeg)
    L1008dir.xyz$dirdeg <- ifelse(L1008dir.xyz$dir == "W", 180, L1008dir.xyz$dirdeg)
    L1008dir.xyz$dirdeg <- ifelse(L1008dir.xyz$dir == "NW", 235, L1008dir.xyz$dirdeg)
    L1008dir.xyz$dirdeg <- ifelse(L1008dir.xyz$dir == "N", 270, L1008dir.xyz$dirdeg)
    L1008dir.xyz$dirdeg <- ifelse(L1008dir.xyz$dir == "NE", 315, L1008dir.xyz$dirdeg)
    
    L1008dir.dat <- data.frame(x=unlist(L1008dir.xyz$x), y=unlist(L1008dir.xyz$y), dirdeg=unlist(L1008dir.xyz$dirdeg))
    L1008dir.rast <- rasterFromXYZ(L1008dir.dat, crs=CRS("+proj=longlat +datum=WGS84"))
    writeRaster(L1008dir.rast, filename="L1008dir.grd", format="raster", overwrite=T)
    
    L0806.dir.xyz <- L0806.dir.xyz[-1,]
    L0806dir.xyz <- as.data.frame(L0806.dir.xyz)
    colnames(L0806dir.xyz) <- c("x","y","dir")
    L0806dir.xyz$dirdeg <- ifelse(L0806dir.xyz$dir == "E", 0, NA)
    L0806dir.xyz$dirdeg <- ifelse(L0806dir.xyz$dir == "SE", 45, L0806dir.xyz$dirdeg)
    L0806dir.xyz$dirdeg <- ifelse(L0806dir.xyz$dir == "S", 90, L0806dir.xyz$dirdeg)
    L0806dir.xyz$dirdeg <- ifelse(L0806dir.xyz$dir == "SW", 135, L0806dir.xyz$dirdeg)
    L0806dir.xyz$dirdeg <- ifelse(L0806dir.xyz$dir == "W", 180, L0806dir.xyz$dirdeg)
    L0806dir.xyz$dirdeg <- ifelse(L0806dir.xyz$dir == "NW", 235, L0806dir.xyz$dirdeg)
    L0806dir.xyz$dirdeg <- ifelse(L0806dir.xyz$dir == "N", 270, L0806dir.xyz$dirdeg)
    L0806dir.xyz$dirdeg <- ifelse(L0806dir.xyz$dir == "NE", 315, L0806dir.xyz$dirdeg)
    
    L0806dir.dat <- data.frame(x=unlist(L0806dir.xyz$x), y=unlist(L0806dir.xyz$y), dirdeg=unlist(L0806dir.xyz$dirdeg))
    L0806dir.rast <- rasterFromXYZ(L0806dir.dat, crs=CRS("+proj=longlat +datum=WGS84"))
    writeRaster(L0806dir.rast, filename="L0806dir.grd", format="raster", overwrite=T)
    
    