sma_plot_stats_comp <- function(vars1,vars2,labels,nbtstrp,makeplot=F,indin=NULL) {
  # Call the sma_regress_multivar function, make plots (if required) and calculate stats
  #
  # Last column in the data.frame is the one being predicted
  #
  # Input: vars - data.frame of input variables for regression
  # Input: labels - axis labels corresponding to columns in vars
  # Input: nbtstrp - number of bootstrapped samples to take in the regression
  # Input: makeplot - make a plot of the first vs last column in vars, including the best fit line
  # Input: indin - optional input array of logicals to select rows from vars
  # Output: return_vals - output of sma_regress_multivar plus RMSE, R, R2, R2adj, estimated variable in final column and number of data points used
  #
  # T. Pugh
  # 23.07.20
  
  nvars1=ncol(vars1)
  nvars2=ncol(vars2)
  nlabels=length(labels)
  if (nvars1!=nlabels | nvars2!=nlabels) {
    stop('Number of variables and number of labels differ')
  }
  nregress=1 #Number of independent regressors
  if (is.null(indin)) {
    ind1=complete.cases(vars1)
    ind2=complete.cases(vars2)
  } else {
    ind=indin
  }
  ndata1=length(ind1[ind1])
  ndata2=length(ind2[ind2])
  
  mod1 <- sma_regress_multivar(vars1[ind1,],nbtstrp,T)
  mod2 <- sma_regress_multivar(vars2[ind2,],nbtstrp,T)
  
  var1_est_boot <- matrix(NA, nrow = ndata1, ncol = nbtstrp)
  var2_est_boot <- matrix(NA, nrow = ndata2, ncol = nbtstrp)
  if (nvars1==2) {
    var1_est <- mod1$intercept_R + mod1$slope_R.y1*vars1[ind1,1]
    var2_est <- mod2$intercept_R + mod2$slope_R.y1*vars2[ind2,1]
    for (nn in 1:nbtstrp) {
      var1_est_boot[,nn] <- mod1$boot.intercept[nn] + mod1$boot.y1[nn]*vars1[ind1,1]
      var2_est_boot[,nn] <- mod2$boot.intercept[nn] + mod2$boot.y1[nn]*vars2[ind2,1]
    }
  } else {
    stop('sma_plot_stats_comp can only handle a maximum of 2 variables')
  }
  var1_est_L95=unname(apply(var1_est_boot, 1, quantile,0.05))
  var1_est_U95=unname(apply(var1_est_boot, 1, quantile,0.95))
  var2_est_L95=unname(apply(var2_est_boot, 1, quantile,0.05))
  var2_est_U95=unname(apply(var2_est_boot, 1, quantile,0.95))
  
  #Make plot
  if (makeplot) {
    # Get axis limits
    xmax=max(max(vars1[ind1,2],na.rm=T),max(vars2[ind2,2],na.rm=T),na.rm=T)
    xmin=min(min(vars1[ind1,2],na.rm=T),min(vars2[ind2,2],na.rm=T),na.rm=T)
    ymax=max(max(vars1[ind1,1],na.rm=T),max(vars2[ind2,1],na.rm=T),na.rm=T)
    ymin=min(min(vars1[ind1,1],na.rm=T),min(vars2[ind2,1],na.rm=T),na.rm=T)
    
    plot(vars1[ind1,nvars1],vars1[ind1,1],pch=16,xlab=labels[nvars1],ylab=labels[1],main=paste(labels[nvars1]," vs ",labels[1]),xlim=c(xmin,xmax),ylim=c(ymin,ymax))
    points(vars2[ind2,nvars2],vars2[ind2,1],col="blue",pch=16)
    lines(var1_est,vars1[ind1,1],col="red")
    lines(var1_est_L95,vars1[ind1,1],col="red",lty="dotted")
    lines(var1_est_U95,vars1[ind1,1],col="red",lty="dotted")
    lines(var2_est,vars2[ind2,1],col="green",pch=16)
    lines(var2_est_L95,vars2[ind2,1],col="green",lty="dotted")
    lines(var2_est_U95,vars2[ind2,1],col="green",lty="dotted")
    #points(var_est_L95,vars[ind,1],col="green",pch=5)
    #points(var_est_U95,vars[ind,1],col="green",pch=5)
  }
  
  #Calculate RMSE
  #res <- vars[ind,nvars]-var_est
  #rmse <- sqrt(mean(res^2))
  
  #Calculate R2
  R_1 <- cor(vars1[ind1,nvars1],var1_est)
  R2_1 <- R_1^2
  R2adj_1 <- 1 - ( ((1-R2_1)*(ndata1-1))/(ndata1-nregress-1) )
  
  R_2 <- cor(vars2[ind2,nvars2],var2_est)
  R2_2 <- R_2^2
  R2adj_2 <- 1 - ( ((1-R2_2)*(ndata2-1))/(ndata2-nregress-1) )
  
  #return_vals <- list("mod"=mod,"rmse"=rmse,"R"=R,"R2"=R2,"R2adj"=R2adj,"var_est"=var_est,
  #                    "var_est_L95"=var_est_L95,"var_est_U95"=var_est_U95,"ndata"=ndata,"dataused"=ind)
  return_vals <- list("mod1"=mod1,"mod2"=mod2,"R2_1"=R2_1,"R2_2"=R2_2,"var1_est"=var1_est,"var2_est"=var2_est,
                      "var1_est_L95"=var1_est_L95,"var1_est_U95"=var1_est_U95,
                      "var2_est_L95"=var2_est_L95,"var2_est_U95"=var2_est_U95,
                      "ndata1"=ndata1,"ndata2"=ndata2,"dataused1"=ind1,"dataused2"=ind2)
  
  return(return_vals)
}

# Make some plots that compare bivariate SMA regressions for broadleaf evergreen and broadleaf deciduous.
#T. Pugh
#27.08.20

traits=read.csv("/Users/pughtam/Documents/TreeMort/Analyses/Hydraulic_modelling/Traits/woody_trait.0827.txt",sep="\t")

trait_allB<-subset(traits,group!="CC",drop = T)
trait_BE<-subset(traits,group=="BE",drop = T)
trait_BDT<-subset(traits,group=="BD" | group=="BT",drop = T)
trait_CC<-subset(traits,group=="CC",drop = T)

#source('trait_functions.R')
source('sma_multivar_regress.R')

nbtstrp=1000

#leafL_from_LMA <- sma_plot_stats_comp(data.frame(LMA,log(leafL)),c("LMA","leafL"),nbtstrp,T)

par(mfrow=c(4,4))
par(mar=c(2,2,2,2))

TLP_from_P50 <- sma_plot_stats_comp(data.frame(trait_BE$P50,trait_BE$TLP),data.frame(trait_BDT$P50,trait_BDT$TLP),c("P50","TLP"),nbtstrp,T)

TLP_from_slope <- sma_plot_stats_comp(data.frame(trait_BE$slope,trait_BE$TLP),data.frame(trait_BDT$slope,trait_BDT$TLP),c("slope","TLP"),nbtstrp,T)

P50_from_slope <- sma_plot_stats_comp(data.frame(trait_BE$slope,trait_BE$P50),data.frame(trait_BDT$slope,trait_BDT$P50),c("slope","P50"),nbtstrp,T)

TLP_from_WD <- sma_plot_stats_comp(data.frame(trait_BE$WD,trait_BE$TLP),data.frame(trait_BDT$WD,trait_BDT$TLP),c("WD","TLP"),nbtstrp,T)

P50_from_WD <- sma_plot_stats_comp(data.frame(trait_BE$WD,trait_BE$P50),data.frame(trait_BDT$WD,trait_BDT$P50),c("WD","P50"),nbtstrp,T)

TLP_from_LMA <- sma_plot_stats_comp(data.frame(trait_BE$LMA,trait_BE$TLP),data.frame(trait_BDT$LMA,trait_BDT$TLP),c("LMA","TLP"),nbtstrp,T)

LS_from_LMA <- sma_plot_stats_comp(data.frame(trait_BE$LMA,trait_BE$LS),data.frame(trait_BDT$LMA,trait_BDT$LS),c("LMA","LS"),nbtstrp,T)

Ks_from_LSHmax <- sma_plot_stats_comp(data.frame(trait_BE$LS_Hmax,trait_BE$Ks),data.frame(trait_BDT$LS_Hmax,trait_BDT$Ks),c("LS*Hmax","Ks"),nbtstrp,T)

Ks_from_P50 <- sma_plot_stats_comp(data.frame(trait_BE$P50,trait_BE$Ks),data.frame(trait_BDT$P50,trait_BDT$Ks),c("P50","Ks"),nbtstrp,T)

LS_from_TLP <- sma_plot_stats_comp(data.frame(trait_BE$TLP,trait_BE$LS),data.frame(trait_BDT$TLP,trait_BDT$LS),c("TLP","LS"),nbtstrp,T)

WD_from_LMA <- sma_plot_stats_comp(data.frame(trait_BE$LMA,trait_BE$WD),data.frame(trait_BDT$LMA,trait_BDT$WD),c("LMA","WD"),nbtstrp,T)

Ks_from_slope <- sma_plot_stats_comp(data.frame(trait_BE$slope,trait_BE$Ks),data.frame(trait_BDT$slope,trait_BDT$Ks),c("slope","Ks"),nbtstrp,T)

leafN_from_LMA <- sma_plot_stats_comp(data.frame(trait_BE$LMA,trait_BE$leafN),data.frame(trait_BDT$LMA,trait_BDT$leafN),c("LMA","leafN"),nbtstrp,T)

leafL_from_LMA <- sma_plot_stats_comp(data.frame(trait_BE$LMA,log(trait_BE$leafL)),data.frame(trait_BDT$LMA,log(trait_BDT$leafL)),c("LMA","leafL"),nbtstrp,T)

Ks_from_LMA <- sma_plot_stats_comp(data.frame(trait_BE$LMA,trait_BE$Ks),data.frame(trait_BDT$LMA,trait_BDT$Ks),c("LMA","Ks"),nbtstrp,T)

                                 