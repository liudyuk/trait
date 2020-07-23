sma_plot_stats <- function(vars,labels,nbtstrp,makeplot=F) {
  
  nvars=ncol(vars)
  nlabels=length(labels)
  if (nvars!=nlabels) {
    stop('Number of variables and number of labels differ')
  }
  nregress=nvars-1 #Number of independent regressors
  
  ind=complete.cases(vars)
  ndata=length(ind[ind])
  
  mod <- sma_regress_multivar(vars[ind,],nbtstrp,T)
  if (nvars==2) {
    var_est <- mod$intercept_R + mod$slope_R.y1*vars[ind,1]
  }
  else if (nvars==3) {
    var_est <- mod$intercept_R + mod$slope_R.y1*vars[ind,1] + mod$slope_R.y2*vars[ind,2]
  }
  else if (nvars==4) {
    var_est <- mod$intercept_R + mod$slope_R.y1*vars[ind,1] + mod$slope_R.y2*vars[ind,2] +
      mod$slope_R.y3*vars[ind,3]
  }
  else if (nvars==5) {
    var_est <- mod$intercept_R + mod$slope_R.y1*vars[ind,1] + mod$slope_R.y2*vars[ind,2] +
      mod$slope_R.y3*vars[ind,3] + mod$slope_R.y4*vars[ind,4]
  }
  
  #Make plot against the first variable
  if (makeplot) {
    plot(vars[ind,nvars],vars[ind,1],pch=16,xlab=labels[nvars],ylab=labels[1],main=paste(labels[nvars]," vs ",labels[1]))
    points(var_est,vars[ind,1],col="red",pch=16)
  }
  
  #Calculate RMSE
  res <- vars[ind,nvars]-var_est
  rmse <- sqrt(mean(res^2))
  
  #Calculate R2
  R <- cor(vars[ind,nvars],var_est)
  R2 <- R^2
  R2adj <- 1 - ( ((1-R2)*(ndata-1))/(ndata-nregress-1) )
  
  return_vals <- list("mod"=mod,"rmse"=rmse,"R"=R,"R2"=R2,"R2adj"=R2adj)
  
  return(return_vals)
}
