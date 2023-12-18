sma_plot_stats <- function(vars,labels,nbtstrp,makeplot=F,indin=NULL,regression_type = 'lm',makeplot_Suppl =F) {
  # Call the sma_regress_multivar lm_regress_multivar or pca_regress_multivar function, make plots (if required) and calculate stats
  #
  # Last column in the data.frame is the one being predicted
  #
  # Input: vars - data.frame of input variables for regression
  # Input: labels - axis labels corresponding to columns in vars
  # Input: nbtstrp - number of bootstrapped samples to take in the regression
  # Input: makeplot - make a plot of the first vs last column in vars, including the best fit line
  # Input: indin - optional input array of logicals to select rows from vars
  # Input: makeplot_suppl - optional, if T, creates plots of bivariate models used to predict auxiliary parameters later
  # Output: return_vals - output of sma_regress_multivar plus RMSE, R, R2, R2adj, estimated variable in final column and number of data points used
  #
  # T. Pugh
  # 23.07.20
  # AHES 06.10.2021, 07.03.2022
  # small adjustments to accommodate lm and pca regression
  
  #testing AHES
  #vars = data.frame(trait_BE$TLP,trait_BE$LS,trait_BE$LMA)
  #labels = c("TLP","LS","LMA")
  #vars =data.frame(trait_BE$TLP,trait_BE$LS,trait_BE$WD,trait_BE$LMA)
  #labels = c("TLP","LS","WD","LMA")
  ## remove above later AHES
  
  nvars=ncol(vars)
  nlabels=length(labels)
  if (nvars!=nlabels) {
    stop('Number of variables and number of labels differ')
  }
  nregress=nvars-1 #Number of independent regressors
  if (is.null(indin)) {
    ind=complete.cases(vars)
  } else {
    ind=indin
  }
  ndata=length(ind[ind])
          
  if(regression_type =='lm'){
  mod <- lm_regress_multivar(vars[ind,],nbtrstrp,T)
  mod$intercept_R
  }
  if(regression_type == 'sma'){
  mod <- sma_regress_multivar(vars[ind,],nbtstrp,T)
  }
  if(regression_type == 'pcr' | regression_type == 'plsr'){
  # scale predictor traits, assuming no NAs ( which should be true with abovestream indexing checks.) 
  vars_scaled <- vars
  vars_scaled[ind,] <- scale_traits(vars[ind,], labels, nlabels, traits_mean, traits_sd)
  mod <- pca_regress_multivar(vars_scaled[ind,], nbtstrp, T, regr_type = regression_type)
  }

  pearson_cor <- cor(vars[ind,])[1,2]

  var_est_boot <- matrix(NA, nrow = ndata, ncol = nbtstrp)
  if (nvars==2) {
    var_est <- mod$intercept_R + mod$slope_R.y1*vars[ind,1]
    for (nn in 1:nbtstrp) {
      var_est_boot[,nn] <- mod$boot.intercept[nn] + mod$boot.y1[nn]*vars[ind,1]
    }
  } else if (nvars==3) {
    var_est <- mod$intercept_R + mod$slope_R.y1*vars[ind,1] + mod$slope_R.y2*vars[ind,2]
    for (nn in 1:nbtstrp) {
      var_est_boot[,nn] <- mod$boot.intercept[nn] + mod$boot.y1[nn]*vars[ind,1] + mod$boot.y2[nn]*vars[ind,2]
    }
  } else if (nvars==4) {
    var_est <- mod$intercept_R + mod$slope_R.y1*vars[ind,1] + mod$slope_R.y2*vars[ind,2] +
      mod$slope_R.y3*vars[ind,3]
    for (nn in 1:nbtstrp) {
      var_est_boot[,nn] <- mod$boot.intercept[nn] + mod$boot.y1[nn]*vars[ind,1] + mod$boot.y2[nn]*vars[ind,2] +
        mod$boot.y3[nn]*vars[ind,3]
    }
  } else if (nvars==5) {
    var_est <- mod$intercept_R + mod$slope_R.y1*vars[ind,1] + mod$slope_R.y2*vars[ind,2] +
      mod$slope_R.y3*vars[ind,3] + mod$slope_R.y4*vars[ind,4]
    for (nn in 1:nbtstrp) {
      var_est_boot[,nn] <- mod$boot.intercept[nn] + mod$boot.y1[nn]*vars[ind,1] + mod$boot.y2[nn]*vars[ind,2] +
        mod$boot.y3[nn]*vars[ind,3] + mod$boot.y4[nn]*vars[ind,4]
    }
  } else if (nvars==6) {
    var_est <- mod$intercept_R + mod$slope_R.y1*vars[ind,1] + mod$slope_R.y2*vars[ind,2] +
      mod$slope_R.y3*vars[ind,3] + mod$slope_R.y4*vars[ind,4] + mod$slope_R.y5*vars[ind,5]
    for (nn in 1:nbtstrp) {
      var_est_boot[,nn] <- mod$boot.intercept[nn] + mod$boot.y1[nn]*vars[ind,1] + mod$boot.y2[nn]*vars[ind,2] +
        mod$boot.y3[nn]*vars[ind,3] + mod$boot.y4[nn]*vars[ind,4] + mod$boot.y5[nn]*vars[ind,5]
    }
  } else {
    stop('sma_plot_stats can only handle a maximum of 6 variables')
  }
  
  
  var_est_L95=unname(apply(var_est_boot, 1, quantile,0.05))
  var_est_U95=unname(apply(var_est_boot, 1, quantile,0.95))
  
  #Make plot against the first variable
  if (makeplot) {
    plot(vars[ind,1],vars[ind,nvars],pch=16,xlab=labels[1],ylab=labels[nvars],main=paste(labels[1]," vs ",labels[nvars]))
    points(vars[ind,1],var_est,col="red",pch=16)
    points(vars[ind,1],var_est_L95,col="green",pch=5)
    points(vars[ind,1],var_est_U95,col="green",pch=5)
  }
  if(makeplot_Suppl){ # not all variables treated here are logged (ln)
    if("leafN" %in% labels){
      trait2 = labels[2]
      ylim = c(1.5,4.5) 
    }else{
      trait2 = paste0("ln(",labels[2],")")
      ylim = c(min(vars[ind,nvars]),max(vars[ind,nvars]))
    }
    trait1 = "ln(LMA)"
    plot(vars[ind,1],vars[ind,nvars],pch=16,xlab="",ylab="", col = makeTransparent('dark grey', alpha=180),axes=FALSE, ylim = ylim)
    points(vars[ind,1],var_est,col="black",pch=16)
    points(vars[ind,1],var_est_L95,col = makeTransparent('black', alpha=80),pch=5)
    points(vars[ind,1],var_est_U95,col = makeTransparent('black', alpha=80),pch=5)
    if("leafN" %in% labels){
      mtext(side=3, line=-2,adj=0.99,paste("leafN = ",round(mod$intercept_R,3)," + ",round(mod$slope_R.y1,digits = 2), " * x"),cex=0.9)
    }else{
      mtext(side=3, line=-1,adj=0.99,paste("y = ",round(mod$intercept_R,3)," + ",round(mod$slope_R.y1,digits = 2), " * x"),cex=0.9)
    }
   
    axis(2, mgp=c(3, .6, 0), tck=-.015, labels=NA)
    axis(side = 2, lwd = 0, line = -0.8, las = 0.5)
    axis(1, mgp=c(3, .6, 0), tck=-.015, labels=NA)
    axis(side = 1, lwd = 0, line = -0.8, las = 0.5)
    box()
    mtext(side = 1,trait1, line=1.1 ,cex=0.9)
    mtext(side = 2,trait2, line=1 ,cex=0.9)
  }
  
  #Calculate RMSE
  res <- vars[ind,nvars]-var_est
  rmse <- sqrt(mean(res^2))
  
  #Calculate R2
  R <- cor(vars[ind,nvars],var_est)
  R2 <- R^2
  R2adj <- 1 - ( ((1-R2)*(ndata-1))/(ndata-nregress-1) )
  
  return_vals <- list("regression_type (lm, sma or pcr/plsr)" = regression_type,"mod"=mod,"rmse"=rmse,"R"=R,"R2"=R2,"R2adj"=R2adj,"var_est"=var_est,
                      "var_est_L95"=var_est_L95,"var_est_U95"=var_est_U95,"ndata"=ndata,"dataused"=ind,'pearson_cor'=  pearson_cor)
 
  return(return_vals)
}
