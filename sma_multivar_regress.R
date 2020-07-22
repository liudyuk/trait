
sma_regress_multivar <- function(yy,nbtstrp) {
  # Function to calculate a trivariate SMA regression 
  # From Richter and Stavn (2014)
  # https://journals.ametsoc.org/jtech/article/31/7/1663/4627/Determining-Functional-Relations-in-Multivariate
  #
  # Note that slopes and intercepts are returned with respect to the last column of the data frame.
  # To calculate slopes and intercepts with respect to another variable, reorder the columns in the data frame.
  # 
  # Returns intercept and slope coefficients, along with bootstrapped 95% confidence intervals for the slope.
  # Bootstrapped intervals are based on 10 000 samples.
  #
  # Input: yy - data.frame of input variables for regression
  # Input: nbtstrp - number of bootstrapped samples to take
  # Output: return_vals - array of intercept, slope, slope confidence intervals and bootstrapped samples
  #
  # T. Pugh
  # 01.07.20
  
  n=nrow(yy)
  c=ncol(yy)
  
  message("Calculating SMA regression with ",c," variables")
  
  # Compute means and standard deviations
  ybar=unname(colMeans(yy))
  sdev=unname(apply(yy, 2, sd))
  
  #Compute correlation matrix
  R=cor(yy)
  
  #Calculate eigenvectors
  PCA_R=eigen(R)
  
  #Compute intercept estimate
  if (c==3) {
    intercept_R=PCA_R$vectors[1,3]*sdev[3]/PCA_R$vectors[3,3]/sdev[1]*ybar[1]+
      PCA_R$vectors[2,3]*sdev[3]/PCA_R$vectors[3,3]/sdev[2]*ybar[2]+
      ybar[3]
  }
  else if (c==4) {
    intercept_R=PCA_R$vectors[1,4]*sdev[4]/PCA_R$vectors[4,4]/sdev[1]*ybar[1]+
      PCA_R$vectors[2,4]*sdev[4]/PCA_R$vectors[4,4]/sdev[2]*ybar[2]+
      PCA_R$vectors[3,4]*sdev[4]/PCA_R$vectors[4,4]/sdev[3]*ybar[3]+
      ybar[4]
  }
  else if (c==5) {
    intercept_R=PCA_R$vectors[1,5]*sdev[5]/PCA_R$vectors[5,5]/sdev[1]*ybar[1]+
      PCA_R$vectors[2,5]*sdev[5]/PCA_R$vectors[5,5]/sdev[2]*ybar[2]+
      PCA_R$vectors[3,5]*sdev[5]/PCA_R$vectors[5,5]/sdev[3]*ybar[3]+
      PCA_R$vectors[4,5]*sdev[5]/PCA_R$vectors[5,5]/sdev[4]*ybar[4]+
      ybar[5]
  }
  else {
    stop('sma_regress_multivar cannot accept more than 5 variables')
  }
  
  #Compute slope estimate
  slope_R.y1=-PCA_R$vectors[1,c]*sdev[c]/PCA_R$vectors[c,c]/sdev[1]
  slope_R.y2=-PCA_R$vectors[2,c]*sdev[c]/PCA_R$vectors[c,c]/sdev[2]
  if (c>3) {
    slope_R.y3=-PCA_R$vectors[3,c]*sdev[c]/PCA_R$vectors[c,c]/sdev[3]
  }
  if (c>4) {
    slope_R.y4=-PCA_R$vectors[4,c]*sdev[c]/PCA_R$vectors[c,c]/sdev[4]
  }
  
  #Generate bootstrap percentile intervals for slope estimates
  set.seed(1234)
  slope_R.y1.boot <- numeric(nbtstrp) 
  slope_R.y2.boot <- numeric(nbtstrp) 
  if (c>3) {
    slope_R.y3.boot <- numeric(nbtstrp)
  }
  if (c>4) {
    slope_R.y4.boot <- numeric(nbtstrp)
  }

  for(i in 1:nbtstrp) {
    index = 1:n
    bootindex = sample(index, n, replace=T) 
    bootsample = yy[bootindex,] 
    Rb=cor(bootsample)
    PCA_Rb=eigen(Rb) 
  
    slope_R.y1.boot[i]=-PCA_Rb$vectors[1,c]*sdev[c]/PCA_Rb$vectors[c,c]/sdev[1]
    slope_R.y2.boot[i]=-PCA_Rb$vectors[2,c]*sdev[c]/PCA_Rb$vectors[c,c]/sdev[2]
    if (c>3) {
      slope_R.y3.boot=-PCA_Rb$vectors[3,c]*sdev[c]/PCA_Rb$vectors[c,c]/sdev[3]
    }
    if (c>4) {
      slope_R.y4.boot=-PCA_Rb$vectors[4,c]*sdev[c]/PCA_Rb$vectors[c,c]/sdev[4]
    }
  }
  
  # Take the 95% confidence interval
  L95_R.y1=quantile(slope_R.y1.boot,0.025)
  U95_R.y1=quantile(slope_R.y1.boot,0.975)
  L95_R.y2=quantile(slope_R.y2.boot,0.025) 
  U95_R.y2=quantile(slope_R.y2.boot,0.975) 
  if (c>3) {
    L95_R.y3=quantile(slope_R.y3.boot,0.025)
    U95_R.y3=quantile(slope_R.y3.boot,0.975)
  }
  if (c>4) {
    L95_R.y4=quantile(slope_R.y4.boot,0.025)
    U95_R.y4=quantile(slope_R.y4.boot,0.975)
  }
  
  if (c==3) {
    return_vals <- list("intercept_R"=intercept_R,"slope_R.y1"=slope_R.y1,"slope_R.y2"=slope_R.y2,
                      "L95_R.y1"=L95_R.y1,"U95_R.y1"=U95_R.y1,"L95_R.y2"=L95_R.y2,"U95_R.y2"=U95_R.y2,
                      "boot.y1"=slope_R.y1.boot,"boot.y2"=slope_R.y2.boot)
  }
  else if (c==4) {
    return_vals <- list("intercept_R"=intercept_R,"slope_R.y1"=slope_R.y1,"slope_R.y2"=slope_R.y2,
                        "slope_R.y3"=slope_R.y3,
                        "L95_R.y1"=L95_R.y1,"U95_R.y1"=U95_R.y1,"L95_R.y2"=L95_R.y2,"U95_R.y2"=U95_R.y2,
                        "L95_R.y3"=L95_R.y3,"U95_R.y3"=U95_R.y3,
                        "boot.y1"=slope_R.y1.boot,"boot.y2"=slope_R.y2.boot,"boot.y3"=slope_R.y3.boot)
  }
  else if (c==5) {
    return_vals <- list("intercept_R"=intercept_R,"slope_R.y1"=slope_R.y1,"slope_R.y2"=slope_R.y2,
                        "slope_R.y3"=slope_R.y3,"slope_R.y4"=slope_R.y4,
                        "L95_R.y1"=L95_R.y1,"U95_R.y1"=U95_R.y1,"L95_R.y2"=L95_R.y2,"U95_R.y2"=U95_R.y2,
                        "L95_R.y3"=L95_R.y3,"U95_R.y3"=U95_R.y3,"L95_R.y4"=L95_R.y4,"U95_R.y4"=U95_R.y4,
                        "boot.y1"=slope_R.y1.boot,"boot.y2"=slope_R.y2.boot,"boot.y3"=slope_R.y3.boot,
                        "boot.y4"=slope_R.y4.boot)
  }
  else {return_vals=NA}
  
  return(return_vals)
}
