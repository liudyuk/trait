pca_regress_multivar <- function(yy,nbtrstrp=10000,bootout=F, regr_type='pcr'){
    # Function to flexibly calculate principal component regression (PCR)/ partial least squares regression (PLSR) regression 
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
    # Input: regr_type - pcr (takes into account variance in x) or plsr (takes into account variance in x and y)
    # Output: return_vals - array of intercept, slope, slope confidence intervals and (optionally) bootstrapped samples
    #
    # A. Eckes-Shephard
    # modified for PCA/PLSR regression 07.03.2022
    # heavily based on function sma_regress_multivar from 
    # T. Pugh
    # 01.07.20
    
    nn=nrow(yy)
    nvars = cc =ncol(yy)
    
    message("Calculating PCR/PLSR regression with ",cc," variables")
    
    #vars = vars[ind,]
    #nvars   = ncol(vars)
    if(nvars == 1){
      stop("Error: Not enough variables for model construction. You need at least 2")
      
    }else{
      
      # apply pcr/plsr model, DF2formula helps to remain flexible with input variables:
      varnames <- names(yy)
      if(regr_type =='pcr'){
        mod_tmp <- pcr(DF2formula(yy[varnames[c(nvars,1:nvars-1)]]), scale = TRUE, data = yy)
        # <- pcr(LS~., data = trait_BDT[c('P50','TLP','Ks','LS')], scale = TRUE, validation = "CV", na.action = na.omit)

      }
      if(regr_type =='plsr'){
        mod_tmp <- plsr(DF2formula(yy[varnames[c(nvars,1:nvars-1)]]),scale = TRUE, data = yy)
      }
      
      
      intercept_R <- as.numeric(coef(mod_tmp,intercept=TRUE)[1])
      
      slope_R.y1 <- as.numeric(coef(mod_tmp,intercept=TRUE)[2])
      if(nvars > 2){
        slope_R.y2 <- as.numeric(coef(mod_tmp,intercept=TRUE)[3])
      }
      if(nvars > 3){
        slope_R.y3 <- as.numeric(coef(mod_tmp,intercept=TRUE)[4])
      }
      if(nvars > 4){
        slope_R.y4 <-as.numeric(coef(mod_tmp,intercept=TRUE)[5])
      }
      if(nvars > 5){
        slope_R.y5 <-as.numeric(coef(mod_tmp,intercept=TRUE)[6])
      }
      if(nvars > 6){
        stop('sma_regress_multivar cannot accept more than 6 variables')
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
      if (cc>5) {
        slope_R.y5.boot <- numeric(nbtstrp)
      }
      
      for(i in 1:nbtstrp) {
        index = 1:nn
        bootindex = sample(index, nn, replace=T) 
        bootsample = yy[bootindex,] 
        
        #Compute intercept estimate 
        if(regr_type =='pcr'){
          mod_boot            = pcr(DF2formula(yy[varnames[c(nvars,1:nvars-1)]]), scale=TRUE,data = bootsample)
        }
        if(regr_type =='plsr'){
          mod_boot            = plsr(DF2formula(yy[varnames[c(nvars,1:nvars-1)]]), scale=TRUE,data = bootsample)
        }
        
        intercept_R.boot[i] = as.numeric(coef(mod_boot,intercept=TRUE)[1])
        
        slope_R.y1.boot[i]  = as.numeric(coef(mod_boot,intercept=TRUE)[2])
        if (cc> 2) {
          slope_R.y2.boot[i] =  as.numeric(coef(mod_boot,intercept=TRUE)[3])
        }
        if(cc> 3){
          slope_R.y3.boot[i] =  as.numeric(coef(mod_boot,intercept=TRUE)[4])
        }
        if(cc> 4){
          slope_R.y4.boot[i] =  as.numeric(coef(mod_boot,intercept=TRUE)[5])
        }
        if(cc> 5){
          slope_R.y5.boot[i] =  as.numeric(coef(mod_boot,intercept=TRUE)[6])
        }
      }#bootstrap end
      
      
      
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
      if (cc>5) {
        L95_R.y5=quantile(slope_R.y5.boot,0.025)
        U95_R.y5=quantile(slope_R.y5.boot,0.975)
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
      if (cc>5) {
        r5 <- list("slope_R.y5"=slope_R.y5,"L95_R.y5"=L95_R.y5,"U95_R.y5"=U95_R.y5)
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
        if (cc>5) {
          b1 <- list("boot.intercept"=intercept_R.boot,"boot.y1"=slope_R.y1.boot,"boot.y2"=slope_R.y2.boot,
                     "boot.y3"=slope_R.y3.boot,"boot.y4"=slope_R.y4.boot,"boot.y5"=slope_R.y5.boot)
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
      else if (cc==6) {
        if (bootout) {
          return_vals <- c(r1,r2,r3,r4,r5,b1)
        }
        else {
          return_vals <- c(r1,r2,r3,r4,r5)
        }
      }
      else {return_vals=NA}
      
      return(return_vals)
    }
    
  }

