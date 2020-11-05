# DesignCTPB

This is the beta version of R package for designing clinical trial with potential biomarker effect. Currently we are working on the following two tasks,\
  (1) preparing documentation for this package.\
  (2) testing this package in various environments and evalueting its consistency. Our original code was developed on Compute Canada Servers with versions of dependency listed above. 
  
For a given setting of input parameters, this package can solve up to 5-dimension alpha-split problems. This can also be expended to handle higher dimension problems. But in practice, we do not suggest consider too high dimensions, since considering too many subpopulation leads to too much loss in power, and not being the optimal choice.\
This package can also guide the choice of size of nested populations, i.e. find optimal r-values. The function visualizes and optimizes r-values, but only supports 3-dimension. The optimization of r-values in more than 3-dimension is trivial, but visualization can be too hard.

We implemented it with GPU computing and smoothing method(thin plate spline). 

## How to install in R:

devtools::install_github("ubcxzhang/DesignCTPB")

## How to run in R:

### for r-split in 3-dimension
Sys.setenv(RETICULATE_PYTHON='python_path')# Note that please specify one python version here instead of using the default r-reticulate python environment, which may encounter errors\
library(DesignCTPB)\
res <- design_ctpb(m=24, n_dim=3, N1=20480, N2=10240, N3=2000, E=NULL, SIGMA=NULL, sd_full=1/base::sqrt(20), DELTA=NULL, delta_linear_bd=c(0.2,0.8), seed=NULL)\
res$plot_alpha # to see the 3-d rotatable plot of optimal alpha versus r2 and r3\
res$plot_power # to see the 3-d rotatable plot of optimal power versus r2 and r3\
res$opt_r_split\
res$opt_alpha_split\
res$opt_power

The default inputs give the results of the strong biomarker effect in our paper. It's no practical meaning to conduct nested-population clinical trails of no biomarker effect, so the no biomarker effect example in our paper can be obtained like: (link)

In our package, the user can specify the standard deviation of each population by giving SIGMA as input, and the harzard reduction rate DELTA for each population. Just give input values to SIGMA and DELTA, but note that the inputed matrix should coincides with the matrix of r-split setting.(e.g. if m=24 and n_dim=3, which means we are going to have 276 r-split setting(like our default setting), so each row of the SIGMA(DELTA) matrix should coincides with the corresponding row of r-split setting). For obtaining the r-split setting, user can specify it personalized or follow our r_setting(m,n_dim) function. 

### for one-point alpha-split
reticulate::source_python(system.file("python","power4R.py",package="DesignCTPB")) \
Power <- Power_sampling # source python4R.py into the environment, only need to source one time if you need to run alpha_split() many times. \
alpha_slpit(r=c(1,0.5,0.3),N1=20480,N2=10240,N3=2000,E=NULL,sig=NULL,sd_full=1/base::sqrt(20),delta=NULL,delta_linear_bd = c(0.2,0.8),Power=Power,seed=NULL)
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



