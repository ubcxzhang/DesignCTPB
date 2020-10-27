# DesignCTPB

This is the immature version of R package for designing clinical trial with potential biomarker effect. We implemented it with GPU computing and smoothing method. 

## How to install in R:

devtools::install_github("ubcxzhang/DesignCTPB")

## How to run in R:

Sys.setenv(RETICULATE_PYTHON='python_path')# Note that please specify one python version here instead of using the default r-reticulate python environment, which may encounter errors

library(DesignCTPB)

res <- plot_clinical()

res$plot_alpha # to see the 3-d rotatable plot of optimal alpha versus r2 and r3

res$plot_power # to see the 3-d rotatable plot of optimal power versus r2 and r3

res$opt_r_split

res$opt_alpha_split

res$opt_power

## R Dependencies:

R/4.0.2

reticulate(Package to interface python in R)

mnormt/fields/plotly/dply

## Python Dependencies:

Python 3.6.3 or later

numba 0.46.0 or later

scipy/numpy/pandas

## GPU and other Dependency 
### Note that we develop our package under this cuda version, while we are still testing for other versions

gcc/7.3.0

CUDA 9.2.148


## Notes:
1) This version is still under development. We are now working on testing the consistency with our original codes.
2) The documentation is also developing right now. And if readers need any documentation to test our package, please contact email: Yitao Lu: yitaolu@uvic.ca/ Belaid: bmoa@uvic.ca. But I hope you could wait untill the formal version is published, thanks for your waiting. 
3) If any suggestions and problems encountered, please contact us, thanks. 


And thanks very much for our package maintainer Belaid, he has done some great jobs. 
