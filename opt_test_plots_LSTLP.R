opt_test_plots_LSTLP <- function(trait,
                           TLP_e_mean,
                           TLP_e_5perc,
                           TLP_e_95perc,
                           TLP_e,
                           P50_e_mean,
                           P50_e_5perc,
                           P50_e_95perc,
                           P50_e,
                           LMA_e_mean,
                           LMA_e_5perc,
                           LMA_e_95perc,
                           LMA_e,
                           WD_e_mean,
                           WD_e_5perc,
                           WD_e_95perc,
                           WD_e,
                           slope_e_mean,
                           slope_e_5perc,
                           slope_e_95perc,
                           slope_e) {
  # Make plots which show the quality of the fit between the optimised values and the original data
  #
  # T. Pugh
  # 13.12.20
  #
  
  
  par(mfrow=c(4,4))
  par(mar=c(2,2,2,2))
  
  plot(trait$TLP,trait$P50,pch=16,xlab="TLP",ylab="P50",main="TLP vs P50")
  points(TLP_e[,1],P50_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(TLP_e[,1],P50_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(TLP_e[,1],P50_e_5perc,col="green",pch=16)
  points(TLP_e[,1],P50_e_95perc,col="green",pch=16)
  
  plot(trait$P50,trait$TLP,pch=16,xlab="P50",ylab="TLP",main="P50 vs TLP")
  points(P50_e[,1],TLP_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(P50_e[,1],TLP_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(P50_e[,1],TLP_e_5perc,col="green",pch=16)
  points(P50_e[,1],TLP_e_95perc,col="green",pch=16)
  
  plot(trait$TLP,trait$slope,pch=16,xlab="TLP",ylab="slope",main="TLP vs slope")
  points(TLP_e[,1],slope_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(TLP_e[,1],slope_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(TLP_e[,1],slope_e_5perc,col="green",pch=16)
  points(TLP_e[,1],slope_e_95perc,col="green",pch=16)
  
  plot(trait$slope,trait$TLP,pch=16,xlab="slope",ylab="TLP",main="slope vs TLP")
  points(slope_e[,1],TLP_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(slope_e[,1],TLP_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(slope_e[,1],TLP_e_5perc,col="green",pch=16)
  points(slope_e[,1],TLP_e_95perc,col="green",pch=16)
  
  plot(trait$P50,trait$slope,pch=16,xlab="P50",ylab="slope",main="P50 vs slope")
  points(P50_e[,1],slope_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(P50_e[,1],slope_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(P50_e[,1],slope_e_5perc,col="green",pch=16)
  points(P50_e[,1],slope_e_95perc,col="green",pch=16)
  
  plot(trait$slope,trait$P50,pch=16,xlab="slope",ylab="P50",main="slope vs P50")
  points(slope_e[,1],P50_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(slope_e[,1],P50_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(slope_e[,1],P50_e_5perc,col="green",pch=16)
  points(slope_e[,1],P50_e_95perc,col="green",pch=16)
  
  plot(trait$TLP,trait$WD,pch=16,xlab="TLP",ylab="WD",main="TLP vs WD")
  points(TLP_e[,1],WD_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(TLP_e[,1],WD_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(TLP_e[,1],WD_e_5perc,col="green",pch=16)
  points(TLP_e[,1],WD_e_95perc,col="green",pch=16)
  
  plot(trait$WD,trait$TLP,pch=16,xlab="WD",ylab="TLP",main="WD vs TLP")
  points(WD_e[,1],TLP_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(WD_e[,1],TLP_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(WD_e[,1],TLP_e_5perc,col="green",pch=16)
  points(WD_e[,1],TLP_e_95perc,col="green",pch=16)
  
  plot(trait$P50,trait$WD,pch=16,xlab="P50",ylab="WD",main="P50 vs WD")
  points(P50_e[,1],WD_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(P50_e[,1],WD_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(P50_e[,1],WD_e_5perc,col="green",pch=16)
  points(P50_e[,1],WD_e_95perc,col="green",pch=16)
  
  plot(trait$WD,trait$P50,pch=16,xlab="WD",ylab="P50",main="WD vs P50")
  points(WD_e[,1],P50_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(WD_e[,1],P50_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(WD_e[,1],P50_e_5perc,col="green",pch=16)
  points(WD_e[,1],P50_e_95perc,col="green",pch=16)
  
  #NOTE: From here onwards I have not made the plots in both directions
  plot(trait$TLP,trait$LMA,pch=16,xlab="TLP",ylab="LMA",main="TLP vs LMA")
  points(TLP_e[,1],LMA_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(TLP_e[,1],LMA_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(TLP_e[,1],LMA_e_5perc,col="green",pch=16)
  points(TLP_e[,1],LMA_e_95perc,col="green",pch=16)
  
  plot(trait$LS,trait$LMA,pch=16,xlab="LS",ylab="LMA",main="LS vs LMA")
  points(LS_e,LMA_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(LS_e,LMA_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(LS_e,LMA_e_5perc,col="green",pch=16)
  points(LS_e,LMA_e_95perc,col="green",pch=16)
  
  plot(trait$Ks,trait$P50,pch=16,xlab="Ks",ylab="P50",main="Ks vs P50")
  points(Ks_e,P50_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(Ks_e,P50_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(Ks_e,P50_e_5perc,col="green",pch=16)
  points(Ks_e,P50_e_95perc,col="green",pch=16)
  
  plot(trait$LS,trait$TLP,pch=16,xlab="LS",ylab="TLP",main="LS vs TLP")
  points(LS_e,TLP_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(LS_e,TLP_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(LS_e,TLP_e_5perc,col="green",pch=16)
  points(LS_e,TLP_e_95perc,col="green",pch=16)
  
  plot(trait$WD,trait$LMA,pch=16,xlab="WD",ylab="LMA",main="WD vs LMA")
  points(WD_e[,1],LMA_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(WD_e[,1],LMA_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(WD_e[,1],LMA_e_5perc,col="green",pch=16)
  points(WD_e[,1],LMA_e_95perc,col="green",pch=16)
  
  plot(trait$Ks,trait$slope,pch=16,xlab="Ks",ylab="slope",main="Ks vs slope")
  points(Ks_e,slope_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(Ks_e,slope_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(Ks_e,slope_e_5perc,col="green",pch=16)
  points(Ks_e,slope_e_95perc,col="green",pch=16)
  
  #Set back to single plot
  par(mfrow=c(1,1))
  par(mar=c(5.1,4.1,4.1,2.1))
  
}

opt_rmse <- function(trait,
                     P50_e,
                     TLP_e,
                     LMA_e,
                     WD_e,
                     slope_e,
                     ind) {
  
  nregress_opt=4 #Assume that there are 4 predictors for each of the 5 optimised traits (i.e. 5-1)
  nregress_other=5 #For the other traits (slope and WD) assume that all the optimised traits are predictors
  
  par(mfrow=c(2,3))
  par(mar=c(5,4,4,2))
  
  ind_P50=which(!is.na(trait$P50[ind]) & !is.na(P50_e[,1]))
  P50_res <- trait$P50[ind[ind_P50]]-P50_e[ind_P50,1]
  P50_e_rmse <- sqrt(mean(P50_res^2,na.rm=T))
  P50_e_ndata <- length(which(is.na(P50_res)==F))
  P50_e_R <- cor(trait$P50[ind[ind_P50]],P50_e[ind_P50,1])
  P50_e_R2 <- P50_e_R^2
  P50_e_R2adj <- 1 - ( ((1-P50_e_R2)*(P50_e_ndata-1))/(P50_e_ndata-nregress_opt-1) )
  plot(trait$P50[ind[ind_P50]],P50_e[ind_P50,1],xlab="P50",ylab="P50_opt")
  
  ind_TLP=which(!is.na(trait$TLP[ind]) & !is.na(TLP_e[,1]))
  TLP_res <- trait$TLP[ind[ind_TLP]]-TLP_e[ind_TLP,1]
  TLP_e_rmse <- sqrt(mean(TLP_res^2,na.rm=T))
  TLP_e_ndata <- length(which(is.na(TLP_res)==F))
  TLP_e_R <- cor(trait$TLP[ind[ind_TLP]],TLP_e[ind_TLP,1])
  TLP_e_R2 <- TLP_e_R^2
  TLP_e_R2adj <- 1 - ( ((1-TLP_e_R2)*(TLP_e_ndata-1))/(TLP_e_ndata-nregress_opt-1) )
  plot(trait$TLP[ind[ind_TLP]],TLP_e[ind_TLP,1],xlab="TLP",ylab="TLP_opt")
  
  ind_LMA=which(!is.na(trait$LMA[ind]) & !is.na(LMA_e[,1]))
  LMA_res <- trait$LMA[ind[ind_LMA]]-LMA_e[ind_LMA,1]
  LMA_e_rmse <- sqrt(mean(LMA_res^2,na.rm=T))
  LMA_e_ndata <- length(which(is.na(LMA_res)==F))
  LMA_e_R <- cor(trait$LMA[ind[ind_LMA]],LMA_e[ind_LMA,1])
  LMA_e_R2 <- LMA_e_R^2
  LMA_e_R2adj <- 1 - ( ((1-LMA_e_R2)*(LMA_e_ndata-1))/(LMA_e_ndata-nregress_opt-1) )
  plot(trait$LMA[ind[ind_LMA]],LMA_e[ind_LMA,1],xlab="LMA",ylab="LMA_opt")
  
  ind_WD=which(!is.na(trait$WD[ind]) & !is.na(WD_e[,1]))
  WD_res <- trait$WD[ind[ind_WD]]-WD_e[ind_WD,1]
  WD_e_rmse <- sqrt(mean(WD_res^2,na.rm=T))
  WD_e_ndata <- length(which(is.na(WD_res)==F))
  WD_e_R <- cor(trait$WD[ind[ind_WD]],WD_e[ind_WD,1])
  WD_e_R2 <- WD_e_R^2
  WD_e_R2adj <- 1 - ( ((1-WD_e_R2)*(WD_e_ndata-1))/(WD_e_ndata-nregress_other-1) )
  plot(trait$WD[ind[ind_WD]],WD_e[ind_WD,1],xlab="WD",ylab="WD_opt")
  
  ind_slope=which(!is.na(trait$slope[ind]) & !is.na(slope_e[,1]))
  slope_res <- trait$slope[ind[ind_slope]]-slope_e[ind_slope,1]
  slope_e_rmse <- sqrt(mean(slope_res^2,na.rm=T))
  slope_e_ndata <- length(which(is.na(slope_res)==F))
  slope_e_R <- cor(trait$slope[ind[ind_slope]],slope_e[ind_slope,1])
  slope_e_R2 <- slope_e_R^2
  slope_e_R2adj <- 1 - ( ((1-slope_e_R2)*(slope_e_ndata-1))/(slope_e_ndata-nregress_other-1) )
  plot(trait$slope[ind[ind_slope]],slope_e[ind_slope,1],xlab="slope",ylab="slope_opt")
  
  # Set back to single plot
  par(mfrow=c(1,1))
  par(mar=c(5.1,4.1,4.1,2.1))
  
  all_names <- c("P50","TLP","LMA","WD","Slope")
  all_e_rmse <- c(P50_e_rmse,TLP_e_rmse,LMA_e_rmse,WD_e_rmse,slope_e_rmse)
  all_e_R <- c(P50_e_R,TLP_e_R,LMA_e_R,WD_e_R,slope_e_R)
  all_e_R2 <- c(P50_e_R2,TLP_e_R2,LMA_e_R2,WD_e_R2,slope_e_R2)
  all_e_R2adj <- c(P50_e_R2adj,TLP_e_R2adj,LMA_e_R2adj,WD_e_R2adj,slope_e_R2adj)
  all_e_ndata <- c(P50_e_ndata,TLP_e_ndata,LMA_e_ndata,WD_e_ndata,slope_e_ndata)
  
  stats <- data.frame(all_names,all_e_rmse,all_e_R,all_e_R2,all_e_R2adj,all_e_ndata)
  
  return_vars <- list("stats"=stats)
  
  return(return_vars)
  
}