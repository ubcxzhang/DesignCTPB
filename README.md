# DesignCTPB

This is the beta version of R package for designing clinical trial with potential biomarker effect. Currently we are working on the following two tasks,\
  (1) preparing documentation for this package.\
  (2) testing this package in various environments and evalueting its consistency. Our original code was developed on Compute Canada Servers with versions of dependency listed above. 
  
For fixed settings, this package can solve up to 5-dimension alpha-split problems. This can be expended to handle higher dimension problems. But in practice, we do not suggest consider too high dimensions, since considering too many subpopulation leads to too much loss in power, and not being the optimal choice.
This package also guide the choice of size of nested populations, i.e. the r-values. The function to visualize and optimize r-values only support 3-dimension. The optimization of r-values in more than 3-dimension is trivial, but visualization can be too hard.\

In our package, we restrict to user's personalized input at this stage. For example, the user can only specify the variance of drug effect's prior distribution by setting "sd" in the input value. But as developing, we will consider more flexibility for users. \

We implemented it with GPU computing and smoothing method(thin plate spline). 

## How to install in R:

devtools::install_github("ubcxzhang/DesignCTPB")

## How to run in R:

### for r-split in 3-dimension
Sys.setenv(RETICULATE_PYTHON='python_path')# Note that please specify one python version here instead of using the default r-reticulate python environment, which may encounter errors\
library(DesignCTPB)\
res <- plot_clinical()\
res$plot_alpha # to see the 3-d rotatable plot of optimal alpha versus r2 and r3\
res$plot_power # to see the 3-d rotatable plot of optimal power versus r2 and r3\
res$opt_r_split\
res$opt_alpha_split\
res$opt_power

### for one-point alpha-split
reticulate::source_python(system.file("python","power4R.py",package="DesignCTPB")) # source python4R.py into the environment\
power <- Power_sampling\
alpha_slpit(r=c(1,0.4,0.2,0.1), sd=1/sqrt(20), N1=20480, N2=10240, N3=4000, lower_bio_eff=0.1, upper_bio_effect=0.5, power=power, seed=NULL)
#### Note for selection of N3
It's still under developing for a better selcetion of N3, which should consider the proportions of each subset.

## R Dependencies:

R/4.0.2\
reticulate(Package to interface python in R)\
mnormt/fields/plotly/dply

## Python Dependencies:

Python 3.6.3 or later\
numba 0.46.0 or later\
scipy/numpy/pandas

## GPU and other Dependency 
### Note that we develop our package under this cuda version, while we are still testing for other versions

gcc/7.3.0\
CUDA 9.2.148



