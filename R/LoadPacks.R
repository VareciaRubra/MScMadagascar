if(!require(plyr)) {install.packages('plyr'); library(plyr)}
if(!require(dplyr)) {install.packages('dplyr'); library(dplyr)}
if(!require(car)) {install.packages('car'); library(car)}
if(!require(magrittr)) {install.packages('magrittr'); library(magrittr)}
if(!require(lme4)) {install.packages('lme4'); library(lme4)}
if(!require(ggplot2)) {install.packages('ggplot2'); library(ggplot2)}
if(!require(GGally)) {install.packages('GGally'); library(GGally)}
if(!require(tidyr)) {install.packages('tidyr'); library(tidyr)}
if(!require(MCMCglmm)) {install.packages('MCMCglmm'); library(MCMCglmm)}
if(!require(reshape2)) {install.packages('reshape2'); library(reshape2)}
if(!require(evolqg)) {install.packages('evolqg'); library(evolqg)}
#devtools::install_github('lem-usp/evolqg'); library(evolqg)
if(!require(readr)) {devtools::install_github('hadley/readr'); library(readr)}
if(!require(doParallel)) {install.packages('doParallel'); library(doParallel)}
if(!require(phytools)) {install.packages('phytools'); library(phytools)}
if(!require(phylobase)) {install.packages('phylobase'); library(phylobase)}
if(!require(RColorBrewer)) {install.packages('RColorBrewer'); library(RColorBrewer)}
if(!require(cowplot)) {install.packages('cowplot'); library(cowplot)}
if(!require(mvoutlier)) {install.packages('mvoutlier'); library(mvoutlier)}
if(!require(MVN)) {install.packages('MVN'); library(MVN)}
if(!require(mvtnorm)) {install.packages('mvtnorm'); library(mvtnorm)}
if(!require(nortest)) {install.packages('nortest'); library(nortest)} # testar normalidade univariados por Kolmogorov-Smirnov usando a lillie.test{nortest}
if(!require(ade4)) {install.packages('ad4'); library(ad4)}
if(!require(vegan)) {install.packages('vegan'); library(vegan)}
if(!require(gdata)) {install.packages('gdata'); library(gdata)}
if(!require(psych)) {install.packages('psych'); library(psych)} # pro Bartlet test cortest.bartlett {psych}
if(!require(xtable)) {install.packages('xtable'); library(xtable)}
if(!require(geomorph)) {install.packages('geomorph'); library(geomorph)} # k.multi pra ver sinal filogenetico a la Adams
if(!require(shapes)) {install.packages('shapes'); library(shapes)} 
if(!require(mvMORPH)) {install.packages('mvMORPH'); library(mvMORPH)} # testar 
if(!require(tseries)) {install.packages('tseries'); library(tseries)} #pro MMRR
if(!require(plspm)) {install.packages('plspm'); library(plspm)} #pro MMRR





#if(!require()) {install.packages(''); library()}

################### Update do R e mais um monte de zica de pacotes e dependencias
# source("https://bioconductor.org/biocLite.R")
# biocLite()
# biocLite(pkgs = "Biostrings")
# install.packages("phangorn")
