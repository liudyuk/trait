sma_plot_stats <- function(vars,labels,nbtstrp,makeplot=F,indin=NULL) {
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
  
  mod <- sma_regress_multivar(vars[ind,],nbtstrp,T)
  
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
  } else {
    stop('sma_plot_stats can only handle a maximum of 5 variables')
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
  
  #Calculate RMSE
  res <- vars[ind,nvars]-var_est
  rmse <- sqrt(mean(res^2))
  
  #Calculate R2
  R <- cor(vars[ind,nvars],var_est)
  R2 <- R^2
  R2adj <- 1 - ( ((1-R2)*(ndata-1))/(ndata-nregress-1) )
  
  return_vals <- list("mod"=mod,"rmse"=rmse,"R"=R,"R2"=R2,"R2adj"=R2adj,"var_est"=var_est,
                      "var_est_L95"=var_est_L95,"var_est_U95"=var_est_U95,"ndata"=ndata,"dataused"=ind)
  
  return(return_vals)
}

sma_plot_stats_comp <- function(vars1,vars2,labels,nbtstrp,makeplot=F,indin=NULL) {
  # Call the sma_regress_multivar function, make bivariate plots comparing two different
  # data subsets and calculate stats
  #
  # Last column in the data.frame is the one being predicted
  #
  # Input: vars1 - data.frame of input variables for regression (first subset of data, e.g. BE)
  # Input: vars2 - data.frame of input variables for regression (second subset of data, e.g. BT+BD)
  # Input: labels - axis labels corresponding to columns in vars
  # Input: nbtstrp - number of bootstrapped samples to take in the regression
  # Input: makeplot - make a plot of the first vs last column in vars, including the best fit line
  # Input: indin - optional input array of logicals to select rows from vars
  # Output: return_vals - output of sma_regress_multivar plus RMSE, R, R2, R2adj, estimated variable in final column and number of data points used
  #
  # T. Pugh
  # 27.08.20
  
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
    
    plot(vars1[ind1,1],vars1[ind1,nvars1],pch=16,xlab=labels[1],ylab=labels[nvars1],main=paste(labels[1]," vs ",labels[nvars1]),xlim=c(xmin,xmax),ylim=c(ymin,ymax))
    points(vars2[ind2,1],vars2[ind2,nvars2],col="blue",pch=16)
    lines(vars1[ind1,1],var1_est,col="red")
    lines(vars1[ind1,1],var1_est_L95,col="red",lty="dotted")
    lines(vars1[ind1,1],var1_est_U95,col="red",lty="dotted")
    lines(vars2[ind2,1],var2_est,col="green",pch=16)
    lines(vars2[ind2,1],var2_est_L95,col="green",lty="dotted")
    lines(vars2[ind2,1],var2_est_U95,col="green",lty="dotted")
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
  
  return_vals <- list("mod1"=mod1,"mod2"=mod2,"R2_1"=R2_1,"R2_2"=R2_2,"var1_est"=var1_est,"var2_est"=var2_est,
                      "var1_est_L95"=var1_est_L95,"var1_est_U95"=var1_est_U95,
                      "var2_est_L95"=var2_est_L95,"var2_est_U95"=var2_est_U95,
                      "ndata1"=ndata1,"ndata2"=ndata2,"dataused1"=ind1,"dataused2"=ind2)
  
  return(return_vals)
}


regress_limit_adjust <- function(var1,var2,var1_from_var2,thres) {
  # Function to adjust the intercept of a linear regression equation to encompass a defined selection of data below the regression line.
  #
  # Input: var1 - dependent variable
  # Input: var2 - independent variable
  # Input: var1_from_var2 - list created by sma_plot_stats describing the regression equation of interest
  # Input: thres - fraction of data to be excluded at either end of the distribution
  #
  # Output: list containing modified intercepts and the original intercept and slope
  #
  # T. Pugh
  # 03.09.20
  
  #Set a step by which to increase the intercept in each iteration
  int_step=var1_from_var2$mod$intercept_R/1000
  
  thres_upper=1-thres
  thres_lower=0+thres
  
  ind=which(is.na(var1)==F & is.na(var2)==F)
  
  #Find the upper limit
  frac_data=0.5
  cc=0
  while (frac_data<thres_upper) {
    cc=cc+1
    int=int_step*cc
    
    var1_pred_upper <- var1_from_var2$mod$intercept_R+int + var1_from_var2$mod$slope_R.y1*var2[ind] # mg N per g leaf
    
    under_thres=which(var1[ind]<var1_pred_upper)
    
    frac_data=length(under_thres)/length(ind)
  }
  intercept_upper <- var1_from_var2$mod$intercept_R+int
  
  #Find the lower limit
  frac_data=0.5
  cc=0
  while (frac_data>thres_lower) {
    cc=cc+1
    int=int_step*cc
    
    var1_pred_lower <- var1_from_var2$mod$intercept_R-int + var1_from_var2$mod$slope_R.y1*var2[ind] # mg N per g leaf
    
    under_thres=which(var1[ind]<var1_pred_lower)
    
    frac_data=length(under_thres)/length(ind)
  }
  intercept_lower <- var1_from_var2$mod$intercept_R-int
  
  return_vals <- list("intercept_upper"=intercept_upper,"intercept_lower"=intercept_lower,"intercept"=var1_from_var2$mod$intercept_R,"slope"=var1_from_var2$mod$slope_R.y1,"var1_pred_upper"=var1_pred_upper,"var1_pred_lower"=var1_pred_lower,"ind"=ind)
  
  return(return_vals)
  
}
