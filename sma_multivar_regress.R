
sma_regress_multivar <- function(yy,nbtstrp=10000,bootout=F) {
  # Function to calculate a multivariate SMA regression 
  # Based on Richter and Stavn (2014)
  # https://journals.ametsoc.org/jtech/article/31/7/1663/4627/Determining-Functional-Relations-in-Multivariate
  #
  # Note that slopes and intercepts are returned with respect to the last column of the data frame.
  # To calculate slopes and intercepts with respect to another variable, reorder the columns in the data frame.
  # 
  # Returns intercept and slope coefficients, along with bootstrapped 95% confidence intervals for the slope 
  # and intercept.
  #
  # Can handle between 2 and 5 input variables.
  #
  # Input: yy - data.frame of input variables for regression
  # Input: nbtstrp - number of bootstrapped samples to take (10 000 is the default)
  # Input: bootout - whether to write out the entire bootstrap sample arrays
  # Output: return_vals - array of intercept, slope, slope confidence intervals and (optionally) bootstrapped samples
  #
  # T. Pugh
  # 01.07.20
  
  nn=nrow(yy)
  cc=ncol(yy)
  
  message("Calculating SMA regression with ",cc," variables")
  
  # Compute means and standard deviations
  ybar=unname(colMeans(yy))
  sdev=unname(apply(yy, 2, sd))
  
  #Compute correlation matrix
  R=cor(yy)
  
  #Calculate eigenvectors
  PCA_R=eigen(R)
  
  #Compute intercept estimate
  if (cc==2) {
    intercept_R=PCA_R$vectors[1,2]*sdev[2]/PCA_R$vectors[2,2]/sdev[1]*ybar[1]+ybar[2]
  }
  else if (cc==3) {
    intercept_R=PCA_R$vectors[1,3]*sdev[3]/PCA_R$vectors[3,3]/sdev[1]*ybar[1]+
      PCA_R$vectors[2,3]*sdev[3]/PCA_R$vectors[3,3]/sdev[2]*ybar[2]+
      ybar[3]
  }
  else if (cc==4) {
    intercept_R=PCA_R$vectors[1,4]*sdev[4]/PCA_R$vectors[4,4]/sdev[1]*ybar[1]+
      PCA_R$vectors[2,4]*sdev[4]/PCA_R$vectors[4,4]/sdev[2]*ybar[2]+
      PCA_R$vectors[3,4]*sdev[4]/PCA_R$vectors[4,4]/sdev[3]*ybar[3]+
      ybar[4]
  }
  else if (cc==5) {
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
  slope_R.y1=-PCA_R$vectors[1,cc]*sdev[cc]/PCA_R$vectors[cc,cc]/sdev[1]
  if (cc>2) {
    slope_R.y2=-PCA_R$vectors[2,cc]*sdev[cc]/PCA_R$vectors[cc,cc]/sdev[2]
  }
  if (cc>3) {
    slope_R.y3=-PCA_R$vectors[3,cc]*sdev[cc]/PCA_R$vectors[cc,cc]/sdev[3]
  }
  if (cc>4) {
    slope_R.y4=-PCA_R$vectors[4,cc]*sdev[cc]/PCA_R$vectors[cc,cc]/sdev[4]
  }
  
  #Generate bootstrap percentile intervals for slope estimates
  set.seed(1234)
  intercept_R.boot <- numeric(nbtstrp) 
  slope_R.y1.boot <- numeric(nbtstrp) 
  if (cc>2) {
    slope_R.y2.boot <- numeric(nbtstrp) 
  }
  if (cc>3) {
    slope_R.y3.boot <- numeric(nbtstrp)
  }
  if (cc>4) {
    slope_R.y4.boot <- numeric(nbtstrp)
  }

  for(i in 1:nbtstrp) {
    index = 1:nn
    bootindex = sample(index, nn, replace=T) 
    bootsample = yy[bootindex,] 
    ybarb=unname(colMeans(bootsample))
    sdevb=unname(apply(bootsample, 2, sd)) #Calculate standard deviation for bootsample (a deviation from the code in Richter and Stavn)
    Rb=cor(bootsample)
    PCA_Rb=eigen(Rb) 
  
    #Compute intercept estimate (an addition to code included in Richter and Stavn)
    if (cc==2) {
      intercept_R.boot[i]=PCA_Rb$vectors[1,2]*sdevb[2]/PCA_Rb$vectors[2,2]/sdevb[1]*ybarb[1]+ybarb[2]
    }
    else if (cc==3) {
      intercept_R.boot[i]=PCA_Rb$vectors[1,3]*sdevb[3]/PCA_Rb$vectors[3,3]/sdevb[1]*ybarb[1]+
        PCA_Rb$vectors[2,3]*sdevb[3]/PCA_Rb$vectors[3,3]/sdevb[2]*ybarb[2]+
        ybarb[3]
    }
    else if (cc==4) {
      intercept_R.boot[i]=PCA_Rb$vectors[1,4]*sdevb[4]/PCA_Rb$vectors[4,4]/sdevb[1]*ybarb[1]+
        PCA_Rb$vectors[2,4]*sdevb[4]/PCA_Rb$vectors[4,4]/sdevb[2]*ybarb[2]+
        PCA_Rb$vectors[3,4]*sdevb[4]/PCA_Rb$vectors[4,4]/sdevb[3]*ybarb[3]+
        ybarb[4]
    }
    else if (cc==5) {
      intercept_R.boot[i]=PCA_Rb$vectors[1,5]*sdevb[5]/PCA_Rb$vectors[5,5]/sdevb[1]*ybarb[1]+
        PCA_Rb$vectors[2,5]*sdevb[5]/PCA_Rb$vectors[5,5]/sdevb[2]*ybarb[2]+
        PCA_Rb$vectors[3,5]*sdevb[5]/PCA_Rb$vectors[5,5]/sdevb[3]*ybarb[3]+
        PCA_Rb$vectors[4,5]*sdevb[5]/PCA_Rb$vectors[5,5]/sdevb[4]*ybarb[4]+
        ybarb[5]
    }
    
    slope_R.y1.boot[i]=-PCA_Rb$vectors[1,cc]*sdevb[cc]/PCA_Rb$vectors[cc,cc]/sdevb[1]
    if (cc>2) {
      slope_R.y2.boot[i]=-PCA_Rb$vectors[2,cc]*sdevb[cc]/PCA_Rb$vectors[cc,cc]/sdevb[2]
    }
    if (cc>3) {
      slope_R.y3.boot=-PCA_Rb$vectors[3,cc]*sdevb[cc]/PCA_Rb$vectors[cc,cc]/sdevb[3]
    }
    if (cc>4) {
      slope_R.y4.boot=-PCA_Rb$vectors[4,cc]*sdevb[cc]/PCA_Rb$vectors[cc,cc]/sdevb[4]
    }
  }
  
  # Take the 95% confidence interval
  L95_R.intercept=quantile(intercept_R.boot,0.025)
  U95_R.intercept=quantile(intercept_R.boot,0.975)
  
  L95_R.y1=quantile(slope_R.y1.boot,0.025)
  U95_R.y1=quantile(slope_R.y1.boot,0.975)
  if (cc>2) {
    L95_R.y2=quantile(slope_R.y2.boot,0.025) 
    U95_R.y2=quantile(slope_R.y2.boot,0.975)
  }
  if (cc>3) {
    L95_R.y3=quantile(slope_R.y3.boot,0.025)
    U95_R.y3=quantile(slope_R.y3.boot,0.975)
  }
  if (cc>4) {
    L95_R.y4=quantile(slope_R.y4.boot,0.025)
    U95_R.y4=quantile(slope_R.y4.boot,0.975)
  }
  
  # Create return values array
  r1 <- list("intercept_R"=intercept_R,"L95_R.intercept"=L95_R.intercept,"U95_R.intercept"=U95_R.intercept,
          "slope_R.y1"=slope_R.y1,"L95_R.y1"=L95_R.y1,"U95_R.y1"=U95_R.y1)
  if (cc>2) {
    r2 <- list("slope_R.y2"=slope_R.y2,"L95_R.y2"=L95_R.y2,"U95_R.y2"=U95_R.y2)
  }
  if (cc>3) {
    r3 <- list("slope_R.y3"=slope_R.y3,"L95_R.y3"=L95_R.y3,"U95_R.y3"=U95_R.y3)
  }
  if (cc>4) {
    r4 <- list("slope_R.y4"=slope_R.y4,"L95_R.y4"=L95_R.y4,"U95_R.y4"=U95_R.y4)
  }
  if (bootout) {
    b1 <- list("boot.intercept"=intercept_R.boot,"boot.y1"=slope_R.y1.boot)
    if (cc>2) {
      b1 <- list("boot.intercept"=intercept_R.boot,"boot.y1"=slope_R.y1.boot,"boot.y2"=slope_R.y2.boot)
    }
    if (cc>3) {
      b1 <- list("boot.intercept"=intercept_R.boot,"boot.y1"=slope_R.y1.boot,"boot.y2"=slope_R.y2.boot,
                 "boot.y3"=slope_R.y3.boot)
    }
    if (cc>4) {
      b1 <- list("boot.intercept"=intercept_R.boot,"boot.y1"=slope_R.y1.boot,"boot.y2"=slope_R.y2.boot,
                 "boot.y3"=slope_R.y3.boot,"boot.y4"=slope_R.y4.boot)
    }
  }

  if (cc==2) {
    if (bootout) {
      return_vals <- c(r1,b1)
    }
    else {
      return_vals <- r1
    }
  }
  else if (cc==3) {
    if (bootout) {
      return_vals <- c(r1,r2,b1)
    }
    else {
      return_vals <- c(r1,r2)
    }
  }
  else if (cc==4) {
    if (bootout) {
      return_vals <- c(r1,r2,r3,b1)
    }
    else {
      return_vals <- c(r1,r2,r3)
    }
  }
  else if (cc==5) {
    if (bootout) {
      return_vals <- c(r1,r2,r3,r4,b1)
    }
    else {
      return_vals <- c(r1,r2,r3,r4)
    }
  }
  else {return_vals=NA}
  
  return(return_vals)
}
