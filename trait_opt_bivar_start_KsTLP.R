# Code heavily based on T. Pugh's function in trait_opt.R
# Code adaptations to accommodate TLP and Ks as predictors from Annemarie Eckes-Shephard, May 2021
trait_opt_bivar_start_KsTLP<- function(P50,
                                        LS,
                                        LMA,
                                        WD,
                                        slope,
                                        LMA_from_TLP,
                                        LMA_from_TLP_LS,
                                        LS_from_TLP_Ks,
                                        LS_from_LMA_TLP_Ks,
                                        P50_from_TLP_Ks,
                                        slope_from_P50_TLP_Ks,
                                        WD_from_slope_P50slope,
                                        LMA_from_LS,
                                        P50_from_Ks,
                                        LS_from_TLP,
                                        TLP_e,
                                        Ks_e,
                                        n_uncer,
                                        use_LMA_from_TLP_LS){
  # Input
  # The following are just used for calculating maximum and minimum values of each trait, beyond which the optimisation should not extend:
  # - P50,LS,LMA,WD,slope
  # The following are the relationships on which the optimisation and associated trait estimates are based:
  # (note that if the relationships used are changed, these need to be updated)
  # - LMA_from_TLP or LMA_from_TLP_LS (choice of which is used is based on use_LMA_from_TLP_LS)
  # - LS_from_TLP_Ks or LS_from_LMA_TLP_Ks (choice of which is used is based on use_LMA_from_TLP_LS)
  # - P50_from_TLP_Ks,
  # - slope_from_P50_TLP_Ks
  # - WD_from_slope_P50slope,
  # The following are all used for initial estimates before optimisation (not expected to have any effect on the results):
  # - LMA_from_LS
  # - P50_from_Ks
  # - LS_from_TLP
  # The following are the traits from which the other traits are predicted:
  # - TLP_e
  # - Ks_e
  
  # Set the tolerance level for the iteration
  tol = 0.00001
  
  #Calculate minimum and maximum values
  maxP50=max(P50,na.rm=T)
  minP50=min(P50,na.rm=T)
  maxLS=max(LS,na.rm=T)
  minLS=min(LS,na.rm=T)
  maxLMA=max(LMA,na.rm=T)
  minLMA=min(LMA,na.rm=T)
  maxWD=max(WD,na.rm=T)
  minWD=min(WD,na.rm=T)
  maxslope=max(slope,na.rm=T)
  minslope=min(slope,na.rm=T)
  
  P50_e <- matrix(NA, ncol = n_uncer) #Array now expanded to hold multiple replicate estimates based on regression coefficient uncertainty
  LMA_e <- matrix(NA, ncol = n_uncer)
  LS_e <- matrix(NA, ncol = n_uncer)
  WD_e <- matrix(NA, ncol = n_uncer)
  slope_e <- matrix(NA, ncol = n_uncer)
  
  # New outer loop which randomly samples regression coefficients from within their uncertainty bounds
  # The random sampling comes from the bootstrap sampling done in the calculations of the SMA regressions themselves. This approach has the big advantages of (a) not having to make any assumptions about the distribution of the coefficient uncertainty and (b) ensuring that the individual slope coefficients within a regression are consistent with each other.
  for (ss in 1:n_uncer) {
    # For now, make samples for the following:
    # (ensure equation choices are consistent with decisions above)
    if (ss==1) { #First pass always calculates the best estimate
      if(use_LMA_from_TLP_LS) {
        mod_LMA_intercept_sample <- LMA_from_TLP_LS$mod$intercept_R #LMA_from_TLP_LS
        mod_LMA_slope_y1_sample <- LMA_from_TLP_LS$mod$slope_R.y1
        mod_LMA_slope_y2_sample <- LMA_from_TLP_LS$mod$slope_R.y2
        
        mod_LS_intercept_sample  <- LS_from_TLP_Ks$mod$intercept_R #LS_from_TLP_Ks
        mod_LS_slope_y1_sample   <- LS_from_TLP_Ks$mod$slope_R.y1 
        mod_LS_slope_y2_sample   <- LS_from_TLP_Ks$mod$slope_R.y2 
        
      } else {
        mod_LMA_intercept_sample <- LMA_from_TLP$mod$intercept_R #LMA_from_TLP
        mod_LMA_slope_y1_sample <- LMA_from_TLP$mod$slope_R.y1
        
        mod_LS_intercept_sample  <- LS_from_LMA_TLP_Ks$mod$intercept_R #LS_from_LMA TLP_Ks
        mod_LS_slope_y1_sample   <- LS_from_LMA_TLP_Ks$mod$slope_R.y1 
        mod_LS_slope_y2_sample   <- LS_from_LMA_TLP_Ks$mod$slope_R.y2 
        mod_LS_slope_y3_sample   <- LS_from_LMA_TLP_Ks$mod$slope_R.y3 
      }

      mod_P50_intercept_sample <- P50_from_TLP_Ks$mod$intercept_R #P50_from_TLP_Ks
      mod_P50_slope_y1_sample  <- P50_from_TLP_Ks$mod$slope_R.y1
      mod_P50_slope_y2_sample  <- P50_from_TLP_Ks$mod$slope_R.y2
      mod_slope_intercept_sample <- slope_from_P50_TLP_Ks$mod$intercept_R #slope_from_P50_TLP_Ks
      mod_slope_slope_y1_sample <- slope_from_P50_TLP_Ks$mod$slope_R.y1
      mod_slope_slope_y2_sample <- slope_from_P50_TLP_Ks$mod$slope_R.y2
      mod_slope_slope_y3_sample <- slope_from_P50_TLP_Ks$mod$slope_R.y3
      mod_WD_intercept_sample <- WD_from_slope_P50slope$mod$intercept_R #WD_from_slope_P50slope
      mod_WD_slope_y1_sample <- WD_from_slope_P50slope$mod$slope_R.y1
      mod_WD_slope_y2_sample <- WD_from_slope_P50slope$mod$slope_R.y2
    } else {
      if(use_LMA_from_TLP_LS) {
        mod_LMA_intercept_sample <- LMA_from_TLP_LS$mod$boot.intercept[ss] #LMA_from_TLP_LS
        mod_LMA_slope_y1_sample  <- LMA_from_TLP_LS$mod$boot.y1[ss]
        mod_LMA_slope_y2_sample  <- LMA_from_TLP_LS$mod$boot.y2[ss]
        
        mod_LS_intercept_sample <- LS_from_TLP_Ks$mod$boot.intercept[ss] #LS_from TLP Ks
        mod_LS_slope_y1_sample <- LS_from_TLP_Ks$mod$boot.y1[ss]
        mod_LS_slope_y2_sample <- LS_from_TLP_Ks$mod$boot.y2[ss]
      } else {
        mod_LMA_intercept_sample <- LMA_from_TLP$mod$boot.intercept[ss] #LMA_from_TLP
        mod_LMA_slope_y1_sample <- LMA_from_TLP$mod$boot.y1[ss]
        
        mod_LS_intercept_sample <- LS_from_LMA_TLP_Ks$mod$boot.intercept[ss] #LS_from LMA TLP Ks
        mod_LS_slope_y1_sample <- LS_from_LMA_TLP_Ks$mod$boot.y1[ss]
        mod_LS_slope_y2_sample <- LS_from_LMA_TLP_Ks$mod$boot.y2[ss]
        mod_LS_slope_y3_sample <- LS_from_LMA_TLP_Ks$mod$boot.y3[ss]
      }

      mod_P50_intercept_sample <- P50_from_TLP_Ks$mod$boot.intercept[ss] #P50_from_TLP_Ks
      mod_P50_slope_y1_sample <- P50_from_TLP_Ks$mod$boot.y1[ss]
      mod_P50_slope_y2_sample <- P50_from_TLP_Ks$mod$boot.y2[ss]
      mod_slope_intercept_sample <- slope_from_P50_TLP_Ks$mod$boot.intercept[ss] #slope_from_P50_TLP_Ks
      mod_slope_slope_y1_sample <- slope_from_P50_TLP_Ks$mod$boot.y1[ss]
      mod_slope_slope_y2_sample <- slope_from_P50_TLP_Ks$mod$boot.y2[ss]
      mod_slope_slope_y3_sample <- slope_from_P50_TLP_Ks$mod$boot.y3[ss]
      mod_WD_intercept_sample <- WD_from_slope_P50slope$mod$boot.intercept[ss] #WD_from_slope_P50slope
      mod_WD_slope_y1_sample <- WD_from_slope_P50slope$mod$boot.y1[ss]
      mod_WD_slope_y2_sample <- WD_from_slope_P50slope$mod$boot.y2[ss]
    }
    # These regression coefficients will now be used in the optimisation calculations
    
    #LS, P50, LMA need optimising
    
    #First set some initial based on simple bivariate relationship. This is just so that the iteration has somewhere to start from. Final result should not be sensitive to these.
    LS_e_last = LS_from_TLP$mod$intercept_R + LS_from_TLP$mod$slope_R.y1*TLP_e
    P50_e_last = P50_from_Ks$mod$intercept_R + P50_from_Ks$mod$slope_R.y1*Ks_e
    LMA_e_last = LMA_from_TLP$mod$intercept_R + LMA_from_TLP$mod$slope_R.y1*LS_e_last
    
    # "diff_" variables hold the difference between the current estimate of a trait value "_e" and the previous
    # estimate "_last"
    # "diff_*_last" variables contain the differences from the last round of iteration
    # (these are compared to differences in the current round of iteration to see if changes are smaller than
    # "tol" and therefore the iteration can stop)
    # Here we initialise the "diff_*_last" variables very high
    diff_P50_last= 100
    diff_LMA_last= 100
    diff_LS_last=  100
    
    # These arrays are just for output, they store the values of every iteration for the current datapoint.
    # Useful for debugging and to check that convergence is working.
    # (only for debugging, can be commented out)
    P50_c <- matrix(NA, nrow = 100)
    LMA_c <- matrix(NA, nrow = 100)
    LS_c <- matrix(NA, nrow = 100)
    
    # Now we start the optimisation loop. Trait values are iterated until the difference between trait
    # values on successive iterations is less than "tol".
    niter=0;
    while (T) {
      niter=niter+1 # Number of iterations completed
      
      # Make estimates of trait values based on the best SMA regressions (probably multivariate in most cases)
      # The estimates of traits in each iteration are based on the estimates of their predictor traits from the previous iteration
      # When choosing TLP_e and Ks as starting values, it 'slices' through the network and determines P50 and LS. Not sure we are challenging the network's plausibility much here:
      if(use_LMA_from_TLP_LS) {
        LMA_e[ss]=mod_LMA_intercept_sample + mod_LMA_slope_y1_sample*TLP_e +
          mod_LMA_slope_y2_sample*LS_e_last
        
        LS_e[ss] = mod_LS_intercept_sample + mod_LS_slope_y1_sample*TLP_e +
          mod_LS_slope_y2_sample * Ks_e 
      } else {
        LMA_e[ss]=mod_LMA_intercept_sample + mod_LMA_slope_y1_sample*TLP_e
        
        LS_e[ss] = mod_LS_intercept_sample + mod_LS_slope_y1_sample*LMA_e +  mod_LS_slope_y2_sample*TLP_e +
          mod_LS_slope_y3_sample * Ks_e 
      }

      
      P50_e[ss]=mod_P50_intercept_sample + mod_P50_slope_y1_sample*TLP_e + 
        mod_P50_slope_y2_sample*Ks_e
      
      if (limitdataranges) {
        #Do not go beyond observed limits of data - if so, discard.
        if (P50_e[ss]>maxP50 | is.na(P50_e[ss])) {P50_e[ss]=NA; break}
        if (P50_e[ss]<minP50 | is.na(P50_e[ss])) {P50_e[ss]=NA; break}
        if (LS_e[ss]> maxLS  | is.na(LS_e[ss])) {LS_e[ss]=NA; break}
        if (LS_e[ss]< minLS  | is.na(LS_e[ss])) {LS_e[ss]=NA; break}
        if (LMA_e[ss]>maxLMA | is.na(LMA_e[ss])) {LMA_e[ss]=NA; break}
        if (LMA_e[ss]<minLMA | is.na(LMA_e[ss])) {LMA_e[ss]=NA; break}
      }
      
      # Save the values for this iteration to the output array (only for debugging, can be commented out)
      P50_c[niter] <- P50_e[ss]
      LMA_c[niter] <- LMA_e[ss]
      LS_c[niter]  <- LS_e[ss]
      
      # Calculate the difference between the current estimate of a trait value "_e" and the previous estimate "_last"
      diff_P50 = P50_e[ss]-P50_e_last
      diff_LMA = LMA_e[ss]-LMA_e_last
      diff_LS  = LS_e[ss] -LS_e_last
      
      # Now we test if the difference between trait estimates on this iteration and between trait estimates on
      # the last iteration is less than "tol" for all traits. If it is we finish the iteration.
      if (abs(diff_P50-diff_P50_last)<tol &&
          abs(diff_LMA-diff_LMA_last)<tol &&
          abs(diff_LS-diff_LS_last)<tol) {
        break
      }
      
      # Save the "diff" values ready for the next iteration
      diff_P50_last=diff_P50
      diff_LMA_last=diff_LMA
      diff_LS_last=diff_LS
      
      # Save the "_e" values ready for the next iteration
      P50_e_last=P50_e[ss]
      LMA_e_last=LMA_e[ss]
      LS_e_last =LS_e[ss]
    }
    
    # After the iteration has finished we can calculate any traits which did not need to be included in the optimisation (because they are not used in the input to calculate any other trait)
    slope_e[ss]=mod_slope_intercept_sample + mod_slope_slope_y1_sample*P50_e[ss] + 
      mod_slope_slope_y2_sample*TLP_e + mod_slope_slope_y3_sample*Ks_e
    WD_e[ss]=mod_WD_intercept_sample + mod_WD_slope_y1_sample*slope_e[ss] + 
      mod_WD_slope_y2_sample*slope_e[ss]*P50_e[ss]
    
    #if (limitdataranges) {
    #  #Do not go beyond observed limits of data
    #  if (slope_e[dd,ss]>maxslope | is.na(slope_e[dd,ss])) {slope_e[dd,ss]=NA}
    #  if (slope_e[dd,ss]<minslope | is.na(slope_e[dd,ss])) {slope_e[dd,ss]=NA}
    #  if (WD_e[dd,ss]>maxWD | is.na(WD_e[dd,ss])) {WD_e[dd,ss]=NA}
    #  if (WD_e[dd,ss]<minWD | is.na(WD_e[dd,ss])) {WD_e[dd,ss]=NA}
    #}
    
  } #Finish nbtstrp loop
  
  return_vals <- list("P50_e"=P50_e,"LS_e"=LS_e,"LMA_e"=LMA_e,"WD_e"=WD_e,"slope_e"=slope_e)
  
  return(return_vals)
}
