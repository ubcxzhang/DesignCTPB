#Plot of Fig. 1
#source("R/smoothing.R")



#' This function is to obtain fig.1 and get the optimal cutoff for the clinical design
#' @export
# the default setting is our strong continuous condition in our paper
design_ctpb <- function(m=24, r_set = NULL, n_dim=3, N1=20480, N2=10240, N3=2000, E=NULL, SIGMA=NULL, sd_full=1/base::sqrt(20), DELTA=NULL, delta_linear_bd=c(0.2,0.8), seed=NULL){
  
  opt_res <- Optim_Res(m, r_set, n_dim, N1, N2, N3, E, SIGMA, sd_full, DELTA, delta_linear_bd, seed)
  r_setting <- as.matrix(opt_res[,1:n_dim]); opt_alpha <- as.matrix(opt_res[,(n_dim+1):(2*n_dim)]); opt_power <- opt_res[,NCOL(opt_res)]
  # we only develop for 3-dim right now, but we can easily extend it into higher dimensional case
  if(n_dim == 3){
    r2 <- r_setting[,2];r3 <- r_setting[,3]; opt_alpha1 <- opt_alpha[,1]; opt_alpha2 <- opt_alpha[,2];opt_alpha3 <- opt_alpha[,3]
    data1<- data.frame(alpha1 = opt_alpha1, r2 =r2, r3=r3 ); data2 <- data.frame(alpha2 = opt_alpha2, r2 =r2, r3=r3 ); data3 <- data.frame(alpha3 = opt_alpha3, r2 =r2, r3=r3 )
    model.p <- suppressWarnings(fields::Tps(cbind(r2,r3), opt_power, m=4))
    model.a1 <- suppressWarnings(fields::Tps(cbind(r2,r3), opt_alpha1, m=4))
    model.a2 <- suppressWarnings(fields::Tps(cbind(r2,r3), opt_alpha2, m=4))
    model.a3 <- suppressWarnings(fields::Tps(cbind(r2,r3), opt_alpha3, m=4))
    #power
    f.p = function(x,y){
      new=data.frame(r2=x,r3=y)
      p = stats::predict(model.p,new)
      return(p)
    }
    r2 <- seq(0,1,0.01); r3 <- seq(0,1,0.01)
    Power<- outer(r2,r3,f.p)
    for(jj in 1:101){
      for(kk in 1:101){
        if(kk<jj){ Power[jj,kk]=Power[jj,kk]}
        else{Power[jj,kk]=NA}
      }
    }
    # 3d-plot of optimal power versus r2 & r3
    try(library(dplyr), silent=TRUE)
    fig.optim.power <- plotly::plot_ly(x=r2, y=r3, z=t(Power)) %>% plotly::add_surface() %>% plotly::layout(scene=list(camera=list(eye=list(x=2, y=-1, z=0.34)),
                                                                                                                       xaxis = list(title = "r2"),
                                                                                                                       yaxis = list(title ="r3"),
                                                                                                                       zaxis = list(title = "Optimal Power ")))
    #alpha
    f1 = function(x,y){
      new=data.frame(r2=x,r3=y)
      p = stats::predict(model.a1,new)
      return(p)
    }
    
    f2 = function(x,y){
      new=data.frame(r2=x,r3=y)
      p = stats::predict(model.a2,new)
      return(p)
    }
    
    f3 = function(x,y){
      new=data.frame(r2=x,r3=y)
      p = stats::predict(model.a3,new)
      return(p)
    }
    
    r2 <- seq(0,1,0.01); r3 <- seq(0,1,0.01)
    pre_alpha1 <- outer(r2,r3,f1); pre_alpha2 <- outer(r2,r3,f2); pre_alpha3 <- outer(r2,r3,f3)
    for(jj in 1:101){
      for(kk in 1:101){
        if(kk<jj){ pre_alpha1[jj,kk]=pre_alpha1[jj,kk];pre_alpha2[jj,kk]=pre_alpha2[jj,kk];pre_alpha3[jj,kk]=pre_alpha3[jj,kk]}
        else{pre_alpha1[jj,kk]=NA; pre_alpha2[jj,kk]=NA; pre_alpha3[jj,kk]=NA}
      }
    }
    # 3d-plot of optimal alpha versus r2 & r3
    fig.alpha <- plotly::plot_ly() #showscale = FALSE)
    fig.alpha<- fig.alpha %>% plotly::add_surface(x=r2,y=r3,z=t(pre_alpha1)) %>% plotly::add_data(data1) %>% plotly::add_markers(x=~r2, y=~r3, z=~alpha1, size=2,symbol= 0,name = "alpha1")  #%>% add_text(x=0.1, y=0.21,0.0172,text = TeX("\\alpha 1")) %>% config(mathjax = "cdn")
    fig.alpha<- fig.alpha %>% plotly::add_surface(x=r2,y=r3,z= t(pre_alpha2),opacity = 0.98) %>% plotly::add_data(data2) %>% plotly::add_markers(x=~r2, y=~r3, z=~alpha2, size=2,symbol= 100,name = "alpha2") # %>% add_text(x=0.2, y=0.5, z=0.0103889,text = TeX("\\alpha 2")) %>% config(mathjax = "cdn")
    fig.alpha<- fig.alpha %>% plotly::add_surface(x=r2,y=r3,z=t(pre_alpha3),opacity = 0.98) %>% plotly::add_data(data3) %>% plotly::add_markers(x=~r2, y=~r3, z=~alpha3, size=2,symbol= 200,name = "alpha3")   %>% plotly::layout(scene=list(camera=list(eye=list(x=2, y=-1, z=0.34)),
                                                                                                                                                                                                                                             xaxis = list(title = "r2"),
                                                                                                                                                                                                                                             yaxis = list(title ="r3"),
                                                                                                                                                                                                                                             zaxis = list(title = "Optimal alpha"))) # %>% config(mathjax = "cdn")
    # obtain the optimal power at cutoff of r2 and r3 to decide whether cut or not
    y <- function(x){
      new=data.frame(r2=x[1],r3=x[2])
      p = -stats::predict(model.p,new)
      return(p)
    }
    r2.max <- r2[which.max(opt_power)]; r3.max <- r3[which.max(opt_power)]
    opt <- optim(c(r2.max, r3.max), y, upper = c(1,1),lower = c(0,0), method = "L-BFGS-B")
  }
  opt_r <- opt$par
  names(opt_r) <- c("r2",'r3')
  ct_opt_alpha1 <- f1(opt_r[1],opt_r[2])
  ct_opt_alpha2 <- f2(opt_r[1],opt_r[2])
  ct_opt_alpha3 <- f3(opt_r[1],opt_r[2])
  ct_opt_alpha <- c(ct_opt_alpha1, ct_opt_alpha2, ct_opt_alpha3)
  names(  ct_opt_alpha) <- paste("alpha",1:3, sep='')
  
  return(list(plot_power=fig.optim.power, plot_alpha = fig.alpha, opt_r_split = opt_r, opt_power = -opt$value, opt_alpha_split= ct_opt_alpha))
}

