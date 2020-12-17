
#' Constraint function
#' @description The constraint function of significant level for the nested-population to test the problem, typically 0.025 for two-side test and 0.05 for one-side test.
#' @param alpha_x sig.lv for the n-th(smallest) sub-population
#' @param alpha_1n vector of sig.lv for the full to (n-1)-th(except for the smallest) sub-population
#' @param r the proportion for each sub-population, r_1 is 1, r_i>r_{i+1}
#' @param  sig.lv significant level for hypothesis testing, usually 0.025 for 2-side test
#' @return the result is 0, refer to the Formula(2) of preprint paper arXiv:2005.10494 
constraint <- function(alpha_x,alpha_1n,r, sig.lv){
  rr <- base::sqrt(r)
  n_dim <- length(r)
  mat <- rr%*%(1/t(rr))
  mat[upper.tri(mat)]<- t(mat)[upper.tri(mat)]
  diag(mat) <- rep(1, n_dim)
  corr <- mat
  return(1 - mnormt::pmnorm(x  = c(stats::qnorm(1-alpha_1n),stats::qnorm(1-alpha_x)), mean = rep(0,length(r)), varcov = corr)[1]-sig.lv)
}


#' Calculate the significant level for the n-th sub-population
#' @description The function is calculating the significant level for the n-th sub-population subjects to the constraint function by giving the sig.lv of the full to the (n-1)-th sub-population
#' @param alpha_tol given the vector of significant level for the full to (n-1)-th sub-population
#' @param r the proportion for each sub-population, r1 is 1, r_i>r_{i+1}
#' @param  sig.lv significant level for hypothesis testing, usually 0.025 for 2-side test
#' @description We apply an one dimensional root finding method
#' @return the significant level for the n-th sub-population subjects to the constraint function by giving the sig.lv of the full to the (n-1)-th sub-population
alpha_kernel <- function(alpha_tol,r, sig.lv){
  aa <- try(stats::uniroot(constraint,interval = c(0,0.025),alpha_1n = alpha_tol,r=r,sig.lv = sig.lv,tol = 1e-10,maxiter = 10000, trace = 2),silent = TRUE)
  if (typeof(aa)=='list' ){
    alpha_x <- aa$root
    return(alpha_x)
  }
}

