# Smoothing

#' This function is to fit a smooth model given alpha and corresponding power values from Monte Carlo sampling, and in 3-dim set, we suggest thin plate splines

#' @export
alpha_split <- function(r=c(1,0.5,0.3),N1=20480,N2=10240,N3=2000,E=NULL,sig=NULL,sd_full=1/base::sqrt(20),delta=NULL,delta_linear_bd = c(0.2,0.8),Power=NULL,seed=NULL){
  n_dim <- length(r)
  if(n_dim>5){
    stop("Right now, we only support 5 dimension alpha-split!")
  }
  if(is.null(Power)){
    reticulate::source_python(system.file("python","power4R.py",package="DesignCTPB")) # source python4R.py into the environment
    Power <- Power_sampling
  }
  estimate_point <- power_estimator(r, N1, N2, N3, E, sig, sd_full, delta, delta_linear_bd ,Power,seed)
  estimate_power <- as.vector(unlist(estimate_point$power)); estimate_alpha <- as.matrix(estimate_point$alpha)
  ## Fit a thin plate splines
  Y <- estimate_power 
  eval(parse(text=paste( paste0("X",1:length(r),"="," estimate_alpha[,",1:length(r), "]",collapse = ";"), sep='')))
  #X1 <- estimate_alpha[,1]; X2 <- estimate_alpha[,2]
  estimate_model <- fields::Tps(eval(parse(text=paste("cbind(" ,paste0("X",1:(length(r)-1),collapse = ","), ")",sep=''))),Y,m = 5)
  eval(parse(text=paste( paste0("X",1:(length(r)-1),".max=","X",1:(length(r)-1),"[which.max(Y)]",collapse = ";"), sep='')))
  #X1.max = X1[which.max(Y)];X2.max = X2[which.max(Y)]
  y <- function(x){
    #new <- data.frame(X1=x[1],X2=x[2])
    x <- t(x)
    names(x) <- paste("X",1:(length(r)-1), sep='')
    new <- data.frame(x)
    p =-predict(estimate_model,new)
    return(p)
  }
  
  est <- stats::optim(eval(parse(text=paste( "c(",paste0("X",1:(length(r)-1),".max",collapse = ","), ")",sep='')))
                      ,y,lower = c(0,0),upper = c(0.025,0.025), method = "L-BFGS-B")
  alphan <- alpha_kernel(c(est$par), r=r,sig.lv = 0.025)
  res <- t(c(est$par, alphan, -est$value)); colnames(res) <- c(paste("alpha",1:length(r), sep=''),'power')
  return(res)
}


#' This function is to decide the r setting given specific density in each dimension
#' @export
r_setting <- function(m, n_dim){
  set <- seq(0.05,0.99,by=1/(m+1))
  flag_s <- n_dim -1
  flag_e <- m
  r2 <- set[flag_s:flag_e]
  r_setting <- matrix(c(rep(1, m-flag_s+1),r2), ncol=2)
  flag_s <- flag_s - 1; flag_e <- flag_e - 1
  while(flag_s >1|flag_s == 1){
    r_x <- set[flag_s:flag_e];
    len_rx <- length(r_x)
    nrow_rsetting <- nrow(r_setting)
    extend_index <- rep(1:len_rx, times = c(1:len_rx))
    r_setting.extend <- r_setting[extend_index,]
    
    rx_extend_index <- rep(0,length(extend_index))
    for(kk in 1:len_rx){
      rx_extend_index[which(extend_index == kk)] <- 1:kk
    }
    r_x.extend <- r_x[rx_extend_index]
    r_setting <- cbind(r_setting.extend, r_x.extend)
    r_setting <- r_setting[order(r_setting[,NCOL(r_setting)]), ]
    flag_s <- flag_s - 1; flag_e <- flag_e - 1
  }
  colnames(r_setting) <- c(paste("r",1:n_dim, sep=''))
  return(r_setting)
}


#' This function is to obtain the optimal results given grid points of r setting

Optim_Res<- function(m, r_set, n_dim, N1, N2, N3, E, SIGMA, sd_full, DELTA, delta_linear_bd, Power, seed){
  if(is.null(r_set)){
    if(ncol(r_set)!=n_dim){
      stop("The dimension of inputed r_set not coincides with dimension!")
    }
    r_set <- r_setting(m, n_dim)
  }
  if(!is.null(SIGMA)){
    if(ncol(SIGMA)!=n_dim){
      stop("The dimension of inputed SIGMA not coincides with dimension!")
    }
    if(nrow(SIGMA)!=nrow(r_set)){
      stop("The inputed SIGMA not coincides with r_set! Plz check with r_setting(m,n_dim) or your input r_set and decide each sig for each r setting.")
    }
  }
  if(!is.null(DELTA)){
    if(ncol(DELTA)!=n_dim){
      stop("The dimension of inputed DELTA not coincides with dimension!")
    }
    if(nrow(DELTA)!=nrow(r_set)){
      stop("The inputed DELTA not coincides with r_set! Plz check with r_setting(m,n_dim) or your input r_set and decide each delta for each r setting.")
    }
  }
  
  optim_res <- matrix(rep(0,nrow(r_set)*(n_dim+1)), nrow=nrow(r_set))
  for(ii in 1:nrow(r_set)){
    r <- r_set[ii,]
    sig <- SIGMA[ii,]# If SIGMA is null then sig is null too, else sig is the user specified value
    delta <- DELTA[ii]# If DELTA is null then delta is null too, else is the user specified value
    optim_res[ii,] <- alpha_split(r,N1,N2,N3,E,sig,sd_full,delta,delta_linear_bd,Power, seed)
  }
  Res <- cbind(r_set, optim_res); colnames(Res) <- c(paste("r", 1:n_dim, sep=""), paste("alpha", 1:n_dim, sep=""), "power")
 return(Res)
}






