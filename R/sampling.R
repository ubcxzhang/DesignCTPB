#SAMPLING FUNCTION'
#' The Alpha function is to generate valid significant level grid values

Alpha <- function(r, N3){
  n_dim <- length(r)
  if(n_dim == 2){
    a <- 1:100/4000.001
    alpha_tol <- as.matrix(expand.grid(alpha1=a))
  }
  if(n_dim == 3){
    a <- 1:100/4000.001
    alpha_tol <- as.matrix(expand.grid(alpha1=a,alpha2=a))
  }
  if(n_dim==4){
    a <- 1:30/1200.001
    alpha_tol <- as.matrix(expand.grid(alpha1=a,alpha2=a, alpha3=a))
  }
  if(n_dim==5){
    a <- 1:18/720.001
    alpha_tol <- as.matrix(expand.grid(alpha1=a,alpha2=a, alpha3=a, alpha4=a))
  }
  
  alpha_1n <- split(alpha_tol, row(alpha_tol))
  clnum <- parallel::detectCores()
  mc <- getOption("mc.cores", clnum)
  alpha <- parallel::mclapply(alpha_1n,alpha_kernel,r, sig.lv=0.025 ,mc.cores = mc)
  alpha_tol <- cbind(alpha_tol, rep(0,(length(a))^(n_dim-1)))
  for(i in 1:(length(a))^(n_dim-1)){
    if (is.null(alpha[[i]])){
      alpha_tol[i,] <- rep(1,n_dim)
    }
    else{
      alpha_tol[i,n_dim] <- alpha[[i]]
      
    }
  }
  index <- which(alpha_tol==rep(1,n_dim))
  if(!setequal(which(alpha_tol==rep(1,n_dim)), integer(0))){
    alpha_tol <- alpha_tol[-index,]
  }#for 2-dimensional case there is no 1 vector
  if(nrow(alpha_tol)>=N3){
    alpha <- alpha_tol[sample(1:nrow(alpha_tol),N3),]
  }
  return(alpha)
} 




#Estimate the power for N3 alpha, given fixed r
#OUTPUT:
## the estimated N3 power values corresponding to fixed alpha1~alpha3, which is grid arranged.

power_estimator <- function(r,N1,N2,N3,E,sig,sd_full,delta,delta_linear_bd,power,seed){
  n_dim <- length(r)
  rr <- base::sqrt(r)
  mat <- rr%*%(1/t(rr))
  mat[upper.tri(mat)]<- t(mat)[upper.tri(mat)]
  diag(mat) <- rep(1, n_dim)
  sigma1 <- mat
  #verify n_dim == length(sig)==length(delta)
  if((!is.null(sig))&&(n_dim!=length(sig))){
    stop("Length of sig not coincides with the dimension!")
  }
  if((!is.null(delta))&&(n_dim!=length(delta))){
    stop("Length of delta not coincides with the dimension!")
  }
  if(is.null(sig)){
    sig <- sd_full*(1/rr)
  }
  # If user dont input the delta for each population, then the default setting is linear
  if(is.null(delta)){
    if(delta_linear_bd[2]>delta_linear_bd[1]){
      delta <- delta_linear[2]-(delta_linear_bd[2]-delta_linear_bd[1])*r
    }
    else{
      stop("Input error of upper bound and lower bound of biomarker effect!")
    }
  }
  sigma2 <- diag(sig)%*%sigma1%*%diag(sig)


  #overall drug effect = lower_bio_eff
  # calculate the mean for the drug effect in each subset
  mean2 <- -base::log(1-delta)
  #generate random vectors for sampling
  if(is.null(seed)){
    set.seed(205851) #  for weak 76605863
  }
  else{
    set.seed(seed)
  }
  R1 <- mnormt::rmnorm(n=N1,mean=rep(0,n_dim), varcov=sigma1)
  R2 <- mnormt::rmnorm(n=N2, mean=mean2, varcov=sigma2)
  #call power in power4R.py and calculate the N power values
  alpha <- Alpha(r,N=N3)
  # If user denote the number of events, then the information units in the algorithm should be E/4, 
  # else we estimate it by the following 
  if(E==NULL){
    It <- (stats::qnorm(0.975)+stats::qnorm(0.9))^2/base::log(1-delta[1])^2
  }
  else{
    It <- E/4
  }
  
  pp <- power(R1,R2,r,It,alpha)
  return(list(alpha=alpha, power=pp))
}



