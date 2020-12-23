\name{res_weak}
\alias{res_weak}
\docType{data}
\title{
The clinical trial design for weak biomarker effect
}
\description{
The dataset contains the 3-D rotable figure for the simulated nested two-subset design(weak biomarker effect), from which the optimal choice of each population's proportion and the optimal alpha split can be abtianed. 
}
\usage{data("res_weak")}
\format{
  The format is:
List of 5
 $ plot_power     :plotly object
 $ plot_alpha     :plotly object
 $ opt_r_split    : numeric vector
 $ opt_power      : numeric 
 $ opt_alpha_split: numeric vector 
}
\examples{
data(res_weak)
}
\keyword{datasets}
