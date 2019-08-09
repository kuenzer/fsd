# fsd - Functional Spatial Data in R

## Spatial Functional PCA

This package implements functional data on a grid of dimension r 
and the spatial principal components analysis.  

## Installing this package

On Windows:  
- Install Rtools if not already installed.
It can be downloaded from https://cran.r-project.org/bin/windows/Rtools/  
- To install the package, run the following lines in R

``` r
if (!require(devtools)) install.packages("devtools")
require(devtools)
if (!find_rtools()) stop("Please install Rtools first.")
install_github("kuenzer/fsd")
```