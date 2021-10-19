opt_test_plots <- function(trait,
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
  points(LS_e[,1],LMA_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(LS_e[,1],LMA_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(LS_e[,1],LMA_e_5perc,col="green",pch=16)
  points(LS_e[,1],LMA_e_95perc,col="green",pch=16)
  
  plot(trait$LS,trait$slope,pch=16,xlab="LS",ylab="LMA",main="LS vs slope")
  points(LS_e[,1],slope_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(LS_e[,1],slope_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(LS_e[,1],slope_e_5perc,col="green",pch=16)
  points(LS_e[,1],slope_e_95perc,col="green",pch=16)
  
  plot(trait$Ks,trait$P50,pch=16,xlab="Ks",ylab="P50",main="Ks vs P50")
  points(Ks_e[,1],P50_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(Ks_e[,1],P50_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(Ks_e[,1],P50_e_5perc,col="green",pch=16)
  points(Ks_e[,1],P50_e_95perc,col="green",pch=16)
  
  plot(trait$LS,trait$TLP,pch=16,xlab="LS",ylab="TLP",main="LS vs TLP")
  points(LS_e[,1],TLP_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(LS_e[,1],TLP_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(LS_e[,1],TLP_e_5perc,col="green",pch=16)
  points(LS_e[,1],TLP_e_95perc,col="green",pch=16)
  
  plot(trait$WD,trait$LMA,pch=16,xlab="WD",ylab="LMA",main="WD vs LMA")
  points(WD_e[,1],LMA_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(WD_e[,1],LMA_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(WD_e[,1],LMA_e_5perc,col="green",pch=16)
  points(WD_e[,1],LMA_e_95perc,col="green",pch=16)
  
  plot(trait$Ks,trait$slope,pch=16,xlab="Ks",ylab="slope",main="Ks vs slope")
  points(Ks_e[,1],slope_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(Ks_e[,1],slope_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(Ks_e[,1],slope_e_5perc,col="green",pch=16)
  points(Ks_e[,1],slope_e_95perc,col="green",pch=16)
  
  #Set back to single plot
  par(mfrow=c(1,1))
  par(mar=c(5.1,4.1,4.1,2.1))
  
}

