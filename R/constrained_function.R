## Constrained function

SQRT<- function(x){
  sqrt_x <- rep(0,length(x))
  for(ii in 1:length(x)){
    sqrt_x[ii] <- base::sqrt(x[ii])
  }
  return(sqrt_x)
}




##OUTPUT
# typically 0.

#' The constraint function of significant level, typically is 0.025 for 2-side test

constraint <- function(alpha_x,alpha_1n,r, sig.lv){
  r <- SQRT(r)
  mat <- r%*%(1/t(r))
  mat[upper.tri(mat)]<- mat[lower.tri(mat)]
  corr <- mat
  if(length(r)>=4){
    return(1 - mvtnorm::pmvnorm(upper  = c(stats::qnorm(1-alpha_1n),stats::qnorm(1-alpha_x)), mean = rep(0,length(r)), corr = corr)[1]-sig.lv)
    
  }
  else{
    return(1 - mvnorm::pmnorm(upper  = c(stats::qnorm(1-alpha_1n),stats::qnorm(1-alpha_x)), mean = rep(0,length(r)), corr = corr)[1]-sig.lv)
  }
}


# OUTPUT: the satisfied alpha or NULL depends on the input

#' The function is calculating the significant level for the n-th sub-population subject to the constraint function

alpha_kernel <- function(alpha_tol,r, sig.lv){
  aa <- try(stats::uniroot(constraint,interval = c(0,0.025),alpha_1n = alpha_tol,r=r,sig.lv = sig.lv,tol = 1e-10,maxiter = 10000, trace = 2),silent = TRUE)
  if (typeof(aa)=='list' ){
    alpha_x <- aa$root
    return(alpha_x)
  }
}
