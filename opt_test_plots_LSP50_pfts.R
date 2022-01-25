
opt_test_plots_LSP50_pfts <- function(trait,
                                 Ks_e_mean,
                                 Ks_e_5perc,
                                 Ks_e_95perc,
                                 Ks_e,
                                 TLP_e_mean,
                                 TLP_e_5perc,
                                 TLP_e_95perc,
                                 TLP_e,
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
                                 slope_e,
                                 outs_LSP50_hv) {
  # Make plots which show the quality of the fit between the optimised values and the original data
  #
  # original version:
  # T. Pugh
  # 13.12.20
  # modified version:
  # Annemarie Eckes-Shephard May 2021 
  # Code heavily based on T. Pugh's function in opt_test_plots.R
  # Small adaptations related to LS and P50 variable and their summary stats.
  # January 2021
  # Add PFT-samples to plot, so see whether all trait spaces are adequately covered
  
  # created dataframe with strategically sampled trait-network traits (from hv-sampling)
  df <- data.frame( TLP_e =rowMeans(outs_LSP50_hv$predicted$TLP_e),Ks_e =rowMeans(outs_LSP50_hv$predicted$Ks_e),
                    LMA_e =rowMeans(outs_LSP50_hv$predicted$LMA_e),slope_e =rowMeans(outs_LSP50_hv$predicted$slope_e),
                    WD_e = rowMeans(outs_LSP50_hv$predicted$WD_e),
                    LS_start=outs_LSP50_hv$predictors$LS_e_start,P50_start=outs_LSP50_hv$predictors$P50_e_start)
  
 # df <- rbind(df,df_maxmin)
 # df <- df_maxmin
  
  par(mfrow=c(4,4))
  par(mar=c(2,2,2,2))
  
  plot(trait$P50,trait$Ks,pch=16,xlab="P50",ylab="Ks",main="P50 vs Ks")
  points(P50_e[,1],Ks_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(P50_e[,1],Ks_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(P50_e[,1],Ks_e_5perc,col="green",pch=16)
  points(P50_e[,1],Ks_e_95perc,col="green",pch=16)
  points(df$P50_start,df$Ks_e, col="purple",pch=9)
  
  plot(trait$TLP,trait$slope,pch=16,xlab="TLP",ylab="slope",main="TLP vs slope")
  points(TLP_e[,1],slope_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(TLP_e[,1],slope_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(TLP_e[,1],slope_e_5perc,col="green",pch=16)
  points(TLP_e[,1],slope_e_95perc,col="green",pch=16)
  points(df$TLP_e,df$slope_e, col="purple",pch=9)
  
  plot(trait$slope,trait$TLP,pch=16,xlab="slope",ylab="TLP",main="slope vs TLP")
  points(slope_e[,1],TLP_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(slope_e[,1],TLP_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(slope_e[,1],TLP_e_5perc,col="green",pch=16)
  points(slope_e[,1],TLP_e_95perc,col="green",pch=16)
  points(df$slope_e,df$TLP_e, col="purple",pch=16,cex=0.4)
  
  plot(trait$P50,trait$slope,pch=16,xlab="P50",ylab="slope",main="P50 vs slope")
  points(P50_e[,1],slope_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(P50_e[,1],slope_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(P50_e[,1],slope_e_5perc,col="green",pch=16)
  points(P50_e[,1],slope_e_95perc,col="green",pch=16)
  points(df$P50_start,df$slope_e, col="purple",pch=9)
  

  plot(trait$TLP,trait$WD,pch=16,xlab="TLP",ylab="WD",main="TLP vs WD")
  points(TLP_e[,1],WD_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(TLP_e[,1],WD_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(TLP_e[,1],WD_e_5perc,col="green",pch=16)
  points(TLP_e[,1],WD_e_95perc,col="green",pch=16)
  points(df$TLP_e,df$WD_e, col="purple",pch=9)
  
   
  plot(trait$WD,trait$TLP,pch=16,xlab="WD",ylab="TLP",main="WD vs TLP")
  points(WD_e[,1],TLP_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(WD_e[,1],TLP_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(WD_e[,1],TLP_e_5perc,col="green",pch=16)
  points(WD_e[,1],TLP_e_95perc,col="green",pch=16)
  points(df$WD_e,df$TLP_e, col="purple",pch=9)
  
  
  plot(trait$P50,trait$WD,pch=16,xlab="P50",ylab="WD",main="P50 vs WD")
  points(P50_e[,1],WD_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(P50_e[,1],WD_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(P50_e[,1],WD_e_5perc,col="green",pch=16)
  points(P50_e[,1],WD_e_95perc,col="green",pch=16)
  points(df$P50_start,df$WD_e, col="purple",pch=9)
  
  
  plot(trait$TLP,trait$LMA,pch=16,xlab="TLP",ylab="LMA",main="TLP vs LMA")
  points(TLP_e[,1],LMA_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(TLP_e[,1],LMA_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(TLP_e[,1],LMA_e_5perc,col="green",pch=16)
  points(TLP_e[,1],LMA_e_95perc,col="green",pch=16)
  points(df$TLP_e,df$LMA_e, col="purple",pch=9)
  
  
  plot(trait$LS,trait$LMA,pch=16,xlab="LS",ylab="LMA",main="LS vs LMA")
  points(LS_e[,1],LMA_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(LS_e[,1],LMA_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(LS_e[,1],LMA_e_5perc,col="green",pch=16)
  points(LS_e[,1],LMA_e_95perc,col="green",pch=16)
  points(df$LS_start,df$LMA_e, col="purple",pch=9)
  
  
  plot(trait$P50,trait$Ks,pch=16,xlab="P50",ylab="Ks",main="P50 vs Ks")
  points(P50_e[,1],Ks_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(P50_e[,1],Ks_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(P50_e[,1],Ks_e_5perc,col="green",pch=16)
  points(P50_e[,1],Ks_e_95perc,col="green",pch=16)
  points(df$P50_start,df$Ks_e, col="purple",pch=9)
  
  
  plot(trait$P50,trait$TLP,pch=16,xlab="P50",ylab="TLP",main="P50 vs TLP")
  points(P50_e[,1],TLP_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(P50_e[,1],TLP_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(P50_e[,1],TLP_e_5perc,col="green",pch=16)
  points(P50_e[,1],TLP_e_95perc,col="green",pch=16)
  points(df$P50_start,df$TLP_e, col="purple",pch=9)
  
  
  plot(trait$WD,trait$LMA,pch=16,xlab="WD",ylab="LMA",main="WD vs LMA")
  points(WD_e[,1],LMA_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(WD_e[,1],LMA_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(WD_e[,1],LMA_e_5perc,col="green",pch=16)
  points(WD_e[,1],LMA_e_95perc,col="green",pch=16)
  points(df$WD_e,df$LMA_e, col="purple",pch=9)
  
  
  plot(trait$Ks,trait$slope,pch=16,xlab="Ks",ylab="slope",main="Ks vs slope")
  points(Ks_e[,1],slope_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(Ks_e[,1],slope_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(Ks_e[,1],slope_e_5perc,col="green",pch=16)
  points(Ks_e[,1],slope_e_95perc,col="green",pch=16)
  points(df$Ks_e,df$slope_e, col="purple",pch=9)
  
  
  plot(trait$Ks,trait$WD,pch=16,xlab="Ks",ylab="WD",main="Ks vs WD")
  points(Ks_e[,1],WD_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(Ks_e[,1],WD_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(Ks_e[,1],WD_e_5perc,col="green",pch=16)
  points(Ks_e[,1],WD_e_95perc,col="green",pch=16)
  points(df$Ks_e,df$WD_e, col="purple",pch=9)
  
  
  plot(trait$WD,trait$slope,pch=16,xlab="WD",ylab="slope",main="WD vs slope")
  points(WD_e[,1],slope_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(WD_e[,1],slope_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(WD_e[,1],slope_e_5perc,col="green",pch=16)
  points(WD_e[,1],slope_e_95perc,col="green",pch=16)
  points(df$WD_e,df$slope_e, col="purple",pch=9)
  
  
  #Set back to single plot
  par(mfrow=c(1,1))
  par(mar=c(5.1,4.1,4.1,2.1))
  
}
