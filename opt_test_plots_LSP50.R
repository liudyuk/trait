opt_test_plots_LSP50 <- function(trait,
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
                                 slope_e) {
  # Make plots which show the quality of the fit between the optimised values and the original data
  #
  # original version:
  # T. Pugh
  # 13.12.20
  # modified version:
  # Annemarie Eckes-Shephard May 2021 
  # Code heavily based on T. Pugh's function in opt_test_plots.R
  # Small adaptations related to LS and P50 variable and their summary stats.
  
  par(mfrow=c(4,4))
  par(mar=c(2,2,2,2))
  
  plot(trait$P50,trait$LS,pch=16,xlab="P50",ylab="LS",main="P50 vs LS", col = makeTransparent('dark grey', alpha=80))
  points(P50_e[,1],LS_e[,1],col="blue",pch=16,cex=0.4) # Using central estimate coefficients
  points(P50_e[,1],LS_e[,1],col="red",pch=16,cex=0.4) # Using mean of all bootstrapped estimates 
  points(P50_e[,1],LS_e[,1],col="green",pch=16,cex=0.4)
  points(P50_e[,1],LS_e[,1],col="green",pch=16,cex=0.4)
  
  
  plot(trait$P50,trait$Ks,pch=16,xlab="P50",ylab="Ls",main="P50 vs Ks", col = makeTransparent('dark grey', alpha=80))
  points(P50_e[,1],Ks_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(P50_e[,1],Ks_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(P50_e[,1],Ks_e_5perc,col="green",pch=16)
  points(P50_e[,1],Ks_e_95perc,col="green",pch=16)
  
  plot(trait$P50,trait$TLP,pch=16,xlab="P50",ylab="TLP",main="P50 vs TLP", col = makeTransparent('dark grey', alpha=80))
  points(P50_e[,1],TLP_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(P50_e[,1],TLP_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(P50_e[,1],TLP_e_5perc,col="green",pch=16)
  points(P50_e[,1],TLP_e_95perc,col="green",pch=16)
  
  plot(trait$P50,trait$slope,pch=16,xlab="P50",ylab="slope",main="P50 vs slope", col = makeTransparent('dark grey', alpha=80))
  points(P50_e[,1],slope_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(P50_e[,1],slope_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(P50_e[,1],slope_e_5perc,col="green",pch=16)
  points(P50_e[,1],slope_e_95perc,col="green",pch=16)
  
  
  plot(trait$P50,trait$WD,pch=16,xlab="P50",ylab="WD",main="P50 vs WD", col = makeTransparent('dark grey', alpha=80))
  points(P50_e[,1],WD_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(P50_e[,1],WD_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(P50_e[,1],WD_e_5perc,col="green",pch=16)
  points(P50_e[,1],WD_e_95perc,col="green",pch=16)
  
  plot(trait$TLP,trait$slope,pch=16,xlab="TLP",ylab="slope",main="TLP vs slope", col = makeTransparent('dark grey', alpha=80))
  points(TLP_e[,1],slope_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(TLP_e[,1],slope_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(TLP_e[,1],slope_e_5perc,col="green",pch=16)
  points(TLP_e[,1],slope_e_95perc,col="green",pch=16)
  
  plot(trait$TLP,trait$WD,pch=16,xlab="TLP",ylab="WD",main="TLP vs WD", col = makeTransparent('dark grey', alpha=80))
  points(TLP_e[,1],WD_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(TLP_e[,1],WD_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(TLP_e[,1],WD_e_5perc,col="green",pch=16)
  points(TLP_e[,1],WD_e_95perc,col="green",pch=16)
  
  plot(trait$WD,trait$TLP,pch=16,xlab="WD",ylab="TLP",main="WD vs TLP", col = makeTransparent('dark grey', alpha=80))
  points(WD_e[,1],TLP_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(WD_e[,1],TLP_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(WD_e[,1],TLP_e_5perc,col="green",pch=16)
  points(WD_e[,1],TLP_e_95perc,col="green",pch=16)

  plot(trait$TLP,trait$LMA,pch=16,xlab="TLP",ylab="LMA",main="TLP vs LMA", col = makeTransparent('dark grey', alpha=80))
  points(TLP_e[,1],LMA_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(TLP_e[,1],LMA_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(TLP_e[,1],LMA_e_5perc,col="green",pch=16)
  points(TLP_e[,1],LMA_e_95perc,col="green",pch=16)
  
  plot(trait$LS,trait$LMA,pch=16,xlab="LS",ylab="LMA",main="LS vs LMA", col = makeTransparent('dark grey', alpha=80))
  points(LS_e[,1],LMA_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(LS_e[,1],LMA_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(LS_e[,1],LMA_e_5perc,col="green",pch=16)
  points(LS_e[,1],LMA_e_95perc,col="green",pch=16)
  
  plot(trait$LS,trait$WD,pch=16,xlab="LS",ylab="WD",main="LS vs WD", col = makeTransparent('dark grey', alpha=80))
  points(LS_e[,1],WD_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(LS_e[,1],WD_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(LS_e[,1],WD_e_5perc,col="green",pch=16)
  points(LS_e[,1],WD_e_95perc,col="green",pch=16)
  
  plot(trait$LS,trait$TLP,pch=16,xlab="LS",ylab="TLP",main="LS vs TLP", col = makeTransparent('dark grey', alpha=80))
  points(LS_e[,1],TLP_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(LS_e[,1],TLP_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(LS_e[,1],TLP_e_5perc,col="green",pch=16)
  points(LS_e[,1],TLP_e_95perc,col="green",pch=16)

  
  plot(trait$Ks,trait$slope,pch=16,xlab="Ks",ylab="slope",main="Ks vs slope", col = makeTransparent('dark grey', alpha=80))
  points(Ks_e[,1],slope_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(Ks_e[,1],slope_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(Ks_e[,1],slope_e_5perc,col="green",pch=16)
  points(Ks_e[,1],slope_e_95perc,col="green",pch=16)
  
  plot(trait$Ks,trait$WD,pch=16,xlab="Ks",ylab="WD",main="Ks vs WD", col = makeTransparent('dark grey', alpha=80))
  points(Ks_e[,1],WD_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(Ks_e[,1],WD_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(Ks_e[,1],WD_e_5perc,col="green",pch=16)
  points(Ks_e[,1],WD_e_95perc,col="green",pch=16)
  
  plot(trait$Ks,trait$TLP,pch=16,xlab="Ks",ylab="TLP",main="Ks vs TLP", col = makeTransparent('dark grey', alpha=80))
  points(Ks_e[,1],TLP_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(Ks_e[,1],TLP_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(Ks_e[,1],TLP_e_5perc,col="green",pch=16)
  points(Ks_e[,1],TLP_e_95perc,col="green",pch=16)
  
  plot(trait$Ks,trait$LMA,pch=16,xlab="Ks",ylab="LMA",main="Ks vs LMA", col = makeTransparent('dark grey', alpha=80))
  points(Ks_e[,1],LMA_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(Ks_e[,1],LMA_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(Ks_e[,1],LMA_e_5perc,col="green",pch=16)
  points(Ks_e[,1],LMA_e_95perc,col="green",pch=16)
  
  plot(trait$WD,trait$slope,pch=16,xlab="WD",ylab="slope",main="WD vs slope", col = makeTransparent('dark grey', alpha=80))
  points(WD_e[,1],slope_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(WD_e[,1],slope_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(WD_e[,1],slope_e_5perc,col="green",pch=16)
  points(WD_e[,1],slope_e_95perc,col="green",pch=16)
  
  plot(trait$WD,trait$LMA,pch=16,xlab="WD",ylab="LMA",main="WD vs LMA", col = makeTransparent('dark grey', alpha=80))
  points(WD_e[,1],LMA_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(WD_e[,1],LMA_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(WD_e[,1],LMA_e_5perc,col="green",pch=16)
  points(WD_e[,1],LMA_e_95perc,col="green",pch=16)
  
  plot(trait$slope,trait$TLP,pch=16,xlab="slope",ylab="TLP",main="slope vs TLP", col = makeTransparent('dark grey', alpha=80))
  points(slope_e[,1],TLP_e[,1],col="blue",pch=16) # Using central estimate coefficients
  points(slope_e[,1],TLP_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
  points(slope_e[,1],TLP_e_5perc,col="green",pch=16)
  points(slope_e[,1],TLP_e_95perc,col="green",pch=16)
  
  #Set back to single plot
  par(mfrow=c(1,1))
  par(mar=c(5.1,4.1,4.1,2.1))
  
}
