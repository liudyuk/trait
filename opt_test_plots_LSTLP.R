opt_test_plots_LSTLP <- function(trait,
                           Ks_e_mean,
                           Ks_e_5perc,
                           Ks_e_95perc,
                           Ks_e,
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
  # original version:
  # T. Pugh
  # 13.12.20
  # modified version:
  # Annemarie Eckes-Shephard May 2021 
  # Code heavily based on T. Pugh's function in opt_test_plots.R
  # Small adaptations related to LS and TLP variable and their summary stats.
  
  plot_scatter <- function(traitdf,trait1,trait2){
    
    # traitdf : data frame with observed traits
    if(trait2 == 'LMA' ){# plot SLA and not LMA, so do 1/LMA
      
      plot(traitdf[[trait1]],1/traitdf[[trait2]],pch=16,xlab=trait1,ylab=trait2,main=paste0(trait1,' vs SLA'), col = makeTransparent('dark grey', alpha=80),axes=FALSE)
      
      
      extended_names <- c('_e','_e_mean','_e_median','_e_5perc','_e_95perc')
      trait2_plot_names <- paste0(trait2,extended_names )
      trait1n <- get(paste0(trait1,'_e'))
      
      points(trait1n[,1], 1/get(trait2_plot_names[1])[,1], col="blue",  pch=16) # Using central estimate coefficients
      points(trait1n[,1], 1/get(trait2_plot_names[2]),     col="red",   pch=16) # Using mean of all bootstrapped estimates 
      points(trait1n[,1], 1/get(trait2_plot_names[4]),     col="green", pch=16)
      points(trait1n[,1], 1/get(trait2_plot_names[5]),     col="green", pch=16)
      ### create linear model and also 95% confidence interval
      
      traitdf1 <- traitdf[,c(trait1,trait2)]
      traitdf1[,2] <- 1/traitdf1[,2] #1/LMA to get to sla
      names(traitdf1) <- c("x","y")
      x = na.omit(traitdf[[trait1]])
      lm.out <-  lm( y ~x,data=traitdf1)
      newx = seq(min(x),max(x),by = 0.05)
      conf_interval <- predict(lm.out, newdata=data.frame(x=newx), interval="confidence",
                               level = 0.95)
      abline( lm.out, col='grey', lty=2, lwd=1.2 )
      matlines(newx, conf_interval[,2:3], col = "grey", lty=2)
      
      #testing
      #p3 <- ggplot(traitdf1, aes(x=x, y=y)) +
      #  geom_point() +
      #  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE)
      #print(p3)
      
      #estimated traits
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
      
      # observed values conf. int
      traitdf1 <- traitdf[,c(trait1,trait2)]
      names(traitdf1) <- c("x","y")
      lm.out <-  lm(y~x,data=traitdf1)
      x = na.omit(traitdf[[trait1]])
      newx = seq(min(x),max(x),by = 0.05)
      conf_interval <- predict(lm.out, newdata=data.frame(x=newx), interval="confidence",
                               level = 0.99)
      abline( lm.out, col='grey', lty=2, lwd=1.2 )
      matlines(newx, conf_interval[,2:3], col = "grey", lty=2)
      
      #estimated values conf. int
      traitdf1= data.frame(x=trait1n[,1],y=get(trait2_plot_names[1])[,1])
      lm.out_est <- lm(y~x,data = traitdf1)
      x= trait1n[,1]
      newx = seq(min(x),max(x),by = 0.05)
      abline(lm.out_est,col='dark green',lty=2,lwd=1.2)
      conf_interval <- predict(lm.out_est, newdata=data.frame(x=newx), interval="confidence",
                               level = 0.95)
      matlines(newx, conf_interval[,2:3], col = "dark green", lty=2)
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
  #plot_scatter(trait_plot,'P50','TLP')
  plot_scatter(trait_plot,'P50','Ks')
  plot_scatter(trait_plot,'P50','LMA')
  
  plot_scatter(trait_plot,'LS','slope')
  plot_scatter(trait_plot,'LS','WD')
  #plot_scatter(trait_plot,'LS','TLP')
  plot_scatter(trait_plot,'LS','LMA')
  plot_scatter(trait_plot,'LS','Ks')
  
  plot_scatter(trait_plot,'TLP','slope')
  plot_scatter(trait_plot,'TLP','WD')
  plot_scatter(trait_plot,'TLP','LMA')
  plot_scatter(trait_plot,'TLP','Ks')
  #plot_scatter(trait_plot,'TLP','P50')
  
  plot_scatter(trait_plot,'Ks','slope')
  plot_scatter(trait_plot,'Ks','WD')
  #plot_scatter(trait_plot,'Ks','TLP')
  plot_scatter(trait_plot,'Ks','LMA')
  
  plot_scatter(trait_plot,'WD','slope')
  #plot_scatter(trait_plot,'WD','TLP')
  plot_scatter(trait_plot,'WD','LMA')
  #plot_scatter(trait_plot,'WD','P50')
  
  plot_scatter(trait_plot,'slope','WD')
  #plot_scatter(trait_plot,'slope','TLP')
  plot_scatter(trait_plot,'slope','LMA')
  #plot_scatter(trait_plot,'slope','P50')
  plot.new()
  legend("center",cex=0.8,legend=c("predicted","observed","predicted","observed"),
         lty=c(NA,NA,2,2), col=c('green','grey','green','grey'),
         fill=c(NA,NA,NA,NA),pch=c(16,16,NA,NA),border=NA,box.lwd=NA)
  
  #Set back to single plot
  par(mfrow=c(1,1))
  par(mar=c(5.1,4.1,4.1,2.1))
  
}