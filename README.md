# Stochastic population projections in Sahul refine the human-refugia hypothesis for early Last Glacial Maximum

Directionally supervised cellular automaton to predict changes in human populations in Sahul through the Late Pleistocene into the Holocene

Data and R code associated with the paper:

<a href="http://scholar.google.com.au/citations?sortby=pubdate&hl=en&user=1sO0O3wAAAAJ&view_op=list_works">BRADSHAW, CJA</a>, <a href="https://stefanicrabtree.com/about-stefani/">SA CRABTREE</a>, <a href="https://santafe.edu/people/profile/devin-white">DA WHITE</a>, <a href="https://research.jcu.edu.au/portfolio/sean.ulm">S ULM</a>, <a href="https://research.jcu.edu.au/portfolio/michael.bird">MI BIRD</a>, <a href="https://www.emmconsulting.com.au/about/leadership-team/dr-alan-william-2/">AN WILLIAMS</a>, <a href="http://www.flinders.edu.au/people/frederik.saltre">F SALTRÉ</a>. 2022. <a href="http://doi.org/10.31219/osf.io/a45fw">Directionally supervised cellular automaton for the initial peopling of Sahul</a>. <em>OSF Preprints</em> doi:10.31219/osf.io/a45fw (in review elsewhere)

## Abstract
Reconstructing the patterns of expansion out of Africa and across the globe by modern <em>Homo sapiens</em> have been advanced using demographic and travel-cost models. However, modelled routes are ipso facto influenced by expansion rates, and vice versa. It is therefore timely to combine these two intertwined phenomena in reconstructions of the history of our species. We applied such an approach to one of the world’s earliest peopling events in Sahul covering the period from 75,000–50,000 years ago by combining recently published models producing highest-likelihood ‘superhighways’ of movement with predictions of expansion based on a demographic cellular automaton. The superhighways-supervised model approximately doubled the predicted time to continental saturation (~ 10K years; 0.4–0.5 km year<sup>-1</sup>) compared to that based on the original, directionally unsupervised model (~ 5K years; 0.7–0.9 km year<sup>-1</sup>), suggesting that rates of migration need to account for topographical constraints in addition to rate of saturation. The modification indicates a slower progression than Neolithic farmers (~ 1 km year<sup>-1</sup>) from the Near East through Europe. Additionally, the scenarios assuming two dominant entry points into Sahul exposed a previously undetected movement corridor south through the centre of Sahul early in the expansion wave. Our combined model infrastructure provides a data-driven means to examine how people initially moved through, settled, and abandoned different regions of the globe, and demonstrates how combining sophisticated continental-scale path modelling with spatial-demographic models can shed light on the complicated process of ancient peopling events. 

<br>
Prof <a href="http://scholar.google.com.au/citations?sortby=pubdate&hl=en&user=1sO0O3wAAAAJ&view_op=list_works">Corey J. A. Bradshaw</a> <br>
<a href="http://globalecologyflinders.com" target="_blank">Global Ecology</a>, <a href="http://flinders.edu.au" target="_blank">Flinders University</a>, Adelaide, Australia <br>
May 2022 <br>
<a href=mailto:corey.bradshaw@flinders.edu.au>e-mail</a> <br>
<br>
<sub><sup><sup>1</sup>Crabtree, S.A. et al. <a href="http://doi.org/10.1038/s41562-021-01106-8">Landscape rules predict optimal super-highways for the first peopling of Sahul</a>. <strong><em>Nature Human Behaviour</strong></em> 5:1303-1313, doi:10.1038/s41562-021-01106-8 (2021) (see also relevant <a href="https://github.com/dawhite/sfa">Github repository</a>)</sup></sub><br>
<sub><sup><sup>2</sup>Bradshaw, C.J.A. et al. <a href="http://doi.org/10.1038/s41467-021-21551-3">Stochastic models support rapid peopling of Late Pleistocene Sahul</a>. <strong><em>Nature Communications</strong></em> 12:2440, doi:10.1038/s41467-021-21551-3 (2021) (see also relevant <a href="https://github.com/cjabradshaw/SahulHumanSpread">Github repository</a>)</sup></sub>


## Code
The R file <a href="https://github.com/cjabradshaw/SuperhighwaysSpreadModel/blob/main/code/SHSpreadPathGithub.R"><code>SHSpreadPathGithub.R</code></a> in the <a href="https://github.com/cjabradshaw/SuperhighwaysSpreadModel/tree/main/code"><code>Code</code></a> directory produces average scenario outputs over a set number of iterations. The user can choose the particulars of the scenario (e.g., underlying <em>K</em>~NPP relationship, entry time(s), entry point(s), spatial clustering, stochastic variances, minimum viable population thresholds, etc.)

The file <a href="https://github.com/cjabradshaw/SuperhighwaysSpreadModel/blob/main/code/matrixOperators.r"><code>matrixOperators.R</code></a> includes necessary functions and is sourced directly within the R code file.

## Data
The two zipped files (<code>CSV files.zip</code> and <code>TIF files.zip</code>) in the <a href="https://github.com/cjabradshaw/SuperhighwaysSpreadModel/tree/main/data"><code>data</code></a> directory should be decompressed and their files placed in the same directory as the R code.
