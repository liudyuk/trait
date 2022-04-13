 -m 'opt_test_plots_LSP50 <- function(trait,
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
  
  plot_scatter <- function(traitdf,trait1,trait2){
    
    if(trait2 == 'LMA' ){# plot SLA and not LMA, so do 1/LMA
      1/trait$LMA
      plot(traitdf[[trait1]],1/traitdf[[trait2]],pch=16,xlab=trait1,ylab=trait2,main=paste0(trait1,' vs SLA'), col = makeTransparent('dark grey', alpha=80),axes=FALSE)
    
      
      extended_names <- c('_e','_e_mean','_e_median','_e_5perc','_e_95perc')
      trait2_plot_names <- paste0(trait2,extended_names )
      trait1n <- get(paste0(trait1,'_e'))
      
      points(trait1n[,1], 1/get(trait2_plot_names[1])[,1], col="blue",  pch=16) # Using central estimate coefficients
      points(trait1n[,1], 1/get(trait2_plot_names[2]),     col="red",   pch=16) # Using mean of all bootstrapped estimates 
      points(trait1n[,1], 1/get(trait2_plot_names[4]),     col="green", pch=16)
      points(trait1n[,1], 1/get(trait2_plot_names[5]),     col="green", pch=16)
      abline( lm( (1/traitdf[[trait2]])             ~ traitdf[[trait1]]), col='grey', lty=2, lwd=1.2 )
      abline( lm( (1/get(trait2_plot_names[1])[,1]) ~ trait1n[,1]), col='dark green', lty=2, lwd=1.2 )
      axis(2, mgp=c(3, .6, 0), tck=-.015, labels=NA)
      axis(side = 2, lwd = 0, line = -0.8, las = 0.5)
      axis(1, mgp=c(3, .6, 0), tck=-.015, labels=NA)
      axis(side = 1, lwd = 0, line = -0.8, las = 0.5)
      box()
      mtext(side = 1,trait1, line=1.1 ,cex=0.9)
      mtext(side = 2,'SLA', line=1 ,cex=0.9)
    }else{ # plot all other traits as normal
      plot(traitdf[[trait1]],traitdf[[trait2]],pch=16,xlab=trait1,ylab=trait2,main=paste0(trait1,' vs ',trait2), col = makeTransparent('dark grey', alpha=80),axes=FALSE)
    
      
      extended_names <- c('_e','_e_mean','_e_median','_e_5perc','_e_95perc')
      trait2_plot_names <- paste0(trait2,extended_names )
      trait1n <- get(paste0(trait1,'_e'))
      
      points(trait1n[,1],get(trait2_plot_names[1])[,1],col="blue",pch=16) # Using central estimate coefficients
      points(trait1n[,1],get(trait2_plot_names[2]),col="red",pch=16) # Using mean of all bootstrapped estimates 
      points(trait1n[,1],get(trait2_plot_names[4]),col="green",pch=16)
      points(trait1n[,1],get(trait2_plot_names[5]),col="green",pch=16)
      abline(lm(traitdf[[trait2]]~traitdf[[trait1]]),col='grey',lty=2,lwd=1.2)
      abline(lm(get(trait2_plot_names[1])[,1]~trait1n[,1]),col='dark green',lty=2,lwd=1.2)
      axis(2, mgp=c(3, .6, 0), tck=-.015, labels=NA)
      axis(side = 2, lwd = 0, line = -0.8, las = 0.5)
      axis(1, mgp=c(3, .6, 0), tck=-.015, labels=NA)
      axis(side = 1, lwd = 0, line = -0.8, las = 0.5)
      box()
      mtext(side = 1,trait1, line=1.1 ,cex=0.9)
      mtext(side = 2,trait2, line=1 ,cex=0.9)
    }# end LMA test
    
  }
  
  par(mfrow=c(4,4))
  par(mar=c(2,2,2,2))
  
  #plot(trait$P50,trait$LS,pch=16,xlab="P50",ylab="LS",main="P50 vs LS", col = makeTransparent('dark grey', alpha=80))
  #abline(lm(trait$LS~trait$P50),col='black',lty=2)
  #points(P50_e[,1],LS_e[,1],col="blue",pch=16,cex=0.4) # Using central estimate coefficients
  #points(P50_e[,1],LS_e[,1],col="red",pch=16,cex=0.4) # Using mean of all bootstrapped estimates 
  #points(P50_e[,1],LS_e[,1],col="green",pch=16,cex=0.4)
  #points(P50_e[,1],LS_e[,1],col="green",pch=16,cex=0.4)
  

  plot_scatter(trait_plot,'P50','slope')
  plot_scatter(trait_plot,'P50','WD')
  plot_scatter(trait_plot,'P50','TLP')
  plot_scatter(trait_plot,'P50','Ks')
  plot_scatter(trait_plot,'P50','LMA')
 
  plot_scatter(trait_plot,'TLP','slope')
  plot_scatter(trait_plot,'TLP','WD')
  plot_scatter(trait_plot,'TLP','LMA')
  #plot_scatter(trait_plot,'TLP','P50')

  plot_scatter(trait_plot,'LS','slope')
  plot_scatter(trait_plot,'LS','WD')
  plot_scatter(trait_plot,'LS','TLP')
  plot_scatter(trait_plot,'LS','LMA')
  
  plot_scatter(trait_plot,'Ks','slope')
  plot_scatter(trait_plot,'Ks','WD')
  plot_scatter(trait_plot,'Ks','TLP')
  plot_scatter(trait_plot,'Ks','LMA')
  
  plot_scatter(trait_plot,'WD','slope')
  #plot_scatter(trait_plot,'WD','TLP')
  plot_scatter(trait_plot,'WD','LMA')
  #plot_scatter(trait_plot,'WD','P50')
  
  plot_scatter(trait_plot,'slope','WD')
  #plot_scatter(trait_plot,'slope','TLP')
  plot_scatter(trait_plot,'slope','LMA')
  #plot_scatter(trait_plot,'slope','P50')
  
  #Set back to single plot
  par(mfrow=c(1,1))
  par(mar=c(5.1,4.1,4.1,2.1))
  
}
