# Stochastic population projections in Sahul refine the human-refugia hypothesis for early Last Glacial Maximum
<img align="right" src="www/Fig2.jpg" alt="Sahul predictions" width="400" style="margin-top: 20px">
<a href="https://doi.org/10.5281/zenodo.7263279"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.7263279.svg" alt="DOI"></a>

Directionally supervised cellular automaton to predict changes in human populations in Sahul through the Late Pleistocene into the Holocene

Data and R code associated with the paper:

<a href="http://scholar.google.com.au/citations?sortby=pubdate&hl=en&user=1sO0O3wAAAAJ&view_op=list_works">BRADSHAW, CJA</a>, <a href="https://www.emmconsulting.com.au/about/leadership-team/dr-alan-william-2/">AN WILLIAMS</a>, <a href="https://research.jcu.edu.au/portfolio/sean.ulm">S ULM</a>, <a href="https://scholar.google.com.au/citations?user=ebgyzfkAAAAJ&hl=en">H CADD</a>, <a href="https://research.jcu.edu.au/portfolio/michael.bird">MI BIRD</a>, <a href="https://stefanicrabtree.com/about-stefani/">SA CRABTREE</a>, <a href="https://santafe.edu/people/profile/devin-white">DA WHITE</a>, <a href="https://research-repository.uwa.edu.au/en/persons/pete-veth">P VETH</a>, <a href="http://www.flinders.edu.au/people/frederik.saltre">F SALTRÉ</a>. 2022. <a href="http://doi.org/10.31219/osf.io/s6hya">Stochastic population projections in Sahul refine the human-refugia hypothesis for early Last Glacial Maximum</a>. <em>OSF Preprints</em> doi:10.31219/osf.io/s6hya (in review elsewhere)

## Abstract
Pleistocene archaeology in Australia has focussed on the survival and behaviour of Indigenous populations across Sahul during the Last Glacial Maximum (28.6 ± 2.8 ka to 17.7 ± 2.2 ka). A long-standing conceptual model proposes people occupied ecological refugia while abandoning drier regions during extreme climatic conditions, with inferred patterns of subsequent recovery essential for describing the evolution of societies in the Holocene. Radiocarbon-derived population estimates partially support the conceptual model while genetic evidence does not, but how human populations were influenced by the Last Glacial Maximum remains untested. To test the refugia hypothesis, we developed a spatial-demographic model of human movement to project population patterns across Sahul from 40,000 years ago (ka) to the mid-Holocene (5 ka). The model predicts little population change in eastern Sahul or New Guinea during the Last Glacial Maximum. However, extensive movement and potential abandonment in the central-western deserts and the north-northwest coastal regions are predicted during the first half of the Last Glacial Maximum (~ 26–20 ka), with some recovery after 15 ka. The demographic implications to societies appear to have extended beyond the Last Glacial Maximum, with increasing populations not evident until the early Holocene in many regions. Our model describes a complex pattern where large areas of Sahul provided refugia that maintained populations throughout the Last Glacial Maximum, providing a possible explanation for the disparity between archaeological and genomic evidence. There was also a correlation between predicted rates of population change and those derived from radiocarbon dates, supporting the realism of applying dates to infer past demography.

<br>
Prof <a href="http://scholar.google.com.au/citations?sortby=pubdate&hl=en&user=1sO0O3wAAAAJ&view_op=list_works">Corey J. A. Bradshaw</a> <br>
<a href="http://globalecologyflinders.com" target="_blank">Global Ecology</a>, <a href="http://flinders.edu.au" target="_blank">Flinders University</a>, Adelaide, Australia <br>
July 2022 <br>
<a href=mailto:corey.bradshaw@flinders.edu.au>e-mail</a> <br>
<br>

## Code
The R file <a href="https://github.com/cjabradshaw/SahulLGMhuman/blob/main/scripts/AusHumSpreadSingleScenAvgLTNprojGithub.R"><code>AusHumSpreadSingleScenAvgLTNprojGithub.R</code></a> in the <a href="https://github.com/cjabradshaw/SahulLGMhuman/tree/main/scripts"><code>scripts</code></a> directory produces average scenario outputs over a set number of iterations.

The file <a href="https://github.com/cjabradshaw/SahulLGM/blob/main/scripts/matrixOperators.r"><code>matrixOperators.R</code></a> includes necessary functions and is sourced directly within the R code file; the file should be placed in a folder called 'scripts'

All code ran on the Flinders University <em>Deepthought</em> High-Performance Computing facility: Flinders University (2021). DeepThought (HPC). doi:<a href="https://doi.org/10.25957/FLINDERS.HPC.DEEPTHOUGHT">10.25957/FLINDERS.HPC.DEEPTHOUGHT</a>

## Data
Decompress the zipped files in the <a href="https://github.com/cjabradshaw/SahulLGMhuman/tree/main/data"><code>data</code></a> directory and place them (as well as the unzipped .csv files) into a sub-directory called 'data'.
