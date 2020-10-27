# Smoothing

#' This function is to fit a smooth model given alpha and corresponding power values from Monte Carlo sampling, and in 3-dim set, we suggest thin plate splines

smooth_power_alpha <- function(r,sd,N1,N2,N3,lower_bio_eff, upper_bio_eff,power){
  estimate_point <- power_estimate_point(r,sd,N1,N2,N3,lower_bio_eff, upper_bio_eff,power)
  estimate_power <- as.vector(unlist(estimate_point$power)); estimate_alpha <- as.matrix(estimate_point$alpha)
  ## Fit a thin plate splines
  Y <- estimate_power ; X1 <- estimate_alpha[,1]; X2 <- estimate_alpha[,2]
  estimate_model <- fields::Tps(cbind(X1,X2),Y,m = 5)
  X1.max = X1[which.max(Y)];X2.max = X2[which.max(Y)]
  y <- function(x){
    new <- data.frame(X1 = x[1],X2 = x[2])
    p =-predict(estimate_model,new)
    return(p)
  }

  est <- stats::optim(c(X1.max, X2.max),y,lower = c(0,0),upper = c(0.025,0.025), method = "L-BFGS-B")
  alpha3 <- alpha_kernel(c(est$par), r=r,sig.lv = 0.025)
  res <- t(c(est$par, alpha3, -est$value)); colnames(res) <- c(paste("alpha",1:length(r), sep=''),'power')
  return(res)
}



#' This function is to obtain the optimal results given grid points of r setting

optim_res<- function(m,n_dim,sd,N1,N2,N3,lower_bio_eff, upper_bio_eff, power){
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
  optim_res <- matrix(rep(0,nrow(r_setting)*(n_dim+1)), nrow=nrow(r_setting))
  for(ii in 1:nrow(r_setting)){
    r <- r_setting[ii,]
    optim_res[ii,] <- smooth_power_alpha(r,sd,N1,N2,N3,lower_bio_eff, upper_bio_eff,power)
  }
  Res <- cbind(r_setting, optim_res); colnames(Res) <- c(paste("r", 1:n_dim, sep=""), paste("alpha", 1:n_dim, sep=""), "power")
 return(Res)
}






