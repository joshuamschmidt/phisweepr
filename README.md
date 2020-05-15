# phisweepr


This is an R/Rcpp implementation of Vy and Kim 2015.

In R:
```R
options(scipen=999)
require(pacman)
pacman::p_load(
devtools,
usethis,
roxygen2,
data.table,
ggplot2,
viridis,
cowplot,
Matrix,
NMOF,
parallel)
source('R/simdata_to_Robjects.R')
source('R/sfs_functions.R')
source('R/utilityFunctions.R')
source('R/MLconfig_functions.R')
source('R/MLfunctions.R')
Rcpp::sourceCpp('src/rcpp_sfs_functions.cpp')
Rcpp::sourceCpp('src/armaCubeFields.cpp')
Rcpp::sourceCpp('src/MLconfig_cpp_functions.cpp')
Rcpp::sourceCpp('src/ML_compute_cpp.cpp')

```
