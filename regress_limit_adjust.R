regress_limit_adjust <- function(var1,var2,var1_from_var2,thres) {
  # Function to adjust the intercept of a linear regression equation to encompass a defined selection of data below the regression line (effectively a quantile regression which we can apply to the SMA regression).
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
