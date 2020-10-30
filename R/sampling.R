#SAMPLING FUNCTION'

#' The Alpha function is to generate valid significant level grid values

Alpha <- function(r, N3){
  n_dim <- length(r)
  if(n_dim == 3){
    a <- 1:100/4000.001
    alpha_tol <- as.matrix(expand.grid(alpha1=a,alpha2=a))
  }
  if(n_dim==4){
    a <- 1:30/1200.001
    alpha_tol <- as.matrix(expand.grid(alpha1=a,alpha2=a, alpha3=a))
  }
  
  alpha_1n <- split(alpha_tol, row(alpha_tol))
  clnum <- parallel::detectCores()
  mc <- getOption("mc.cores", clnum)
  alpha <- parallel::mclapply(alpha_1n,alpha_kernel,r, sig.lv=0.025 ,mc.cores = mc)
  alpha_tol <- cbind(alpha_tol, rep(0,(length(a))^(n_dim-1)))
  for(i in 1:(length(a))^(n_dim-1)){
    if (is.null(alpha[[i]])){
      alpha_tol[i,] <- rep(1,length(r))
      
      
    }
    else{
      alpha_tol[i,n_dim] <- alpha[[i]]
      
    }
  }
  alpha_tol <- alpha_tol[-which(alpha_tol==rep(1,length(r))),]
  if(nrow(alpha_tol)>(N3+1)){
    alpha <- alpha_tol[sample(2:nrow(alpha_tol),N3),]
  }
  return(alpha)
}


#OUTPUT:
# the drug effect given fixed proportion

#' This function if to calculate the drug effect in each sub-population given the proportion
#'


bio_effect <- function(r, lower_bio_eff, upper_bio_eff){
  return((lower_bio_eff-upper_bio_eff)*r+upper_bio_eff)
}


#Estimate the power for N3 alpha, given fixed r
#OUTPUT:
## the estimated N3 power values corresponding to fixed alpha1~alpha3, which is grid arranged.

#' @export

power_estimate_point <- function(r,sd,N1,N2,N3,lower_bio_eff, upper_bio_eff,power, seed=seed){

  rr <- SQRT(r)
  mat <- rr%*%(1/t(rr))
  mat[upper.tri(mat)]<- mat[lower.tri(mat)]
  sigma1 <- mat
  ss <- sd * (1/rr)
  sigma2 <- diag(ss)%*%sigma1%*%diag(ss)


  #overall drug effect = lower_bio_eff
  # calculate the mean for the drug effect in each subset
  mean2 <- -base::log(1-bio_effect(r, lower_bio_eff, upper_bio_eff))
  #generate random vectors for sampling
  if(is.null(seed)){
    set.seed(205851) #  for weak 76605863
  }
  else{
    set.seed(seed)
  }
  R1 <- mnormt::rmnorm(n=N1,mean=rep(0,length(r)), varcov=sigma1)
  R2 <- mnormt::rmnorm(n=N2, mean=mean2, varcov=sigma2)
  #call power in power4R.py and calculate the N power values
  alpha <- Alpha(r,N=N3)
  It <- (stats::qnorm(0.975)+stats::qnorm(0.9))^2/base::log(1-lower_bio_eff)^2

  #base::sink("py_configuration.txt")
  #print(reticulate::py_config())
  #print(reticulate::py_run_file(system.file("python","version.py",package="DesignCTBP")))
  #base::sink()
  
  pp <- power(R1,R2,r,It,alpha)
  return(list(alpha=alpha, power=pp))
}



