trait_opt_bivar_start_LSP50 <- function(P50,
                      TLP,
                      LMA,
                      WD,
                      slope,
                      LS,
                      Ks,
                      LMA_from_TLP,
                      LMA_from_TLP_LS,
                      LS_from_P50_TLP_Ks_LMA,
                      LS_from_P50_TLP_Ks,
                      Ks_from_P50_LS,
                      TLP_from_LS_LMA_P50,
                      P50_from_TLP_Ks,
                      slope_from_P50_TLP_WD,
                      WD_from_P50_Ks,
                      #WD_from_P50_LMA_Ks,
                      LMA_from_LS,
                      P50_from_Ks,
                      TLP_from_P50,
                      WD_from_LMA,
                      slope_from_WD,
                      Ks_from_LS,
                      LS_from_LMA,
                      P50_e_start,
                      LS_e_start,
                      n_uncer,
                      use_LMA_from_TLP_LS,
                      regr_type) {
  # Input
  # The following are just used for calculating maximum and minimum values of each trait, beyond which the optimisation should not extend:
  # - P50,TLP,LMA,WD,slope
  # The following are the relationships on which the optimisation and associated trait estimates are based:
  # (note that if the relationships used are changed, these need to be updated)
  # - LMA_from_TLP or LMA_from_TLP_LS (choice of which is used is based on use_LMA_from_TLP_LS)
  # - TLP_from_LS_LMA_P50
  # - P50_from_TLP_Ks,
  # - slope_from_P50_TLP_Ks
  # - WD_from_slope_P50slope,
  # The following are all used for initial estimates before optimisation (not expected to have any effect on the results):
  # - LMA_from_LS
  # - P50_from_Ks
  # - TLP_from_P50
  # The following are the traits from which the other traits are predicted:
  # - Ks_e_start
  # - LS_e_start
  
  # Set the tolerance level for the iteration
  tol=0.00001
  
  #Calculate minimum and maximum values
  maxP50=max(P50,na.rm=T)+0.1*max(P50,na.rm=T)
  minP50=min(P50,na.rm=T)-0.1*min(P50,na.rm=T)
  maxTLP=max(TLP,na.rm=T)+0.1*max(TLP,na.rm=T)
  minTLP=min(TLP,na.rm=T)-0.1*min(TLP,na.rm=T)
  maxLMA=max(LMA,na.rm=T)+0.1*max(LMA,na.rm=T)
  minLMA=min(LMA,na.rm=T)-0.1*min(LMA,na.rm=T)
  maxWD=max(WD,na.rm=T)+0.1*max(WD,na.rm=T)
  minWD=min(WD,na.rm=T)- 0.1*min(WD,na.rm=T)
  maxslope=max(slope,na.rm=T) +0.1*max(slope,na.rm=T)
  minslope=min(slope,na.rm=T) -0.1*min(slope,na.rm=T)
  maxLS=max(LS,na.rm=T)+ 0.1*max(LS,na.rm=T)
  minLS=min(LS,na.rm=T)- 0.1*min(LS,na.rm=T)
  maxKs=max(Ks,na.rm=T)+ 0.1*max(Ks,na.rm=T)
  minKs=min(Ks,na.rm=T)+ 0.1*min(Ks,na.rm=T)
  
  
  P50_e   <- matrix(NA, ncol = n_uncer) #Array now expanded to hold multiple replicate estimates based on regression coefficient uncertainty
  LMA_e   <- matrix(NA, ncol = n_uncer)
  TLP_e   <- matrix(NA, ncol = n_uncer)
  WD_e    <- matrix(NA, ncol = n_uncer)
  slope_e <- matrix(NA, ncol = n_uncer)
  LS_e    <- matrix(NA, ncol = n_uncer)
  Ks_e    <- matrix(NA, ncol = n_uncer)
  WD_e    <- matrix(NA, ncol = n_uncer)
  slope_e    <- matrix(NA, ncol = n_uncer)
  
  # New outer loop which randomly samples regression coefficients from within their uncertainty bounds
  # The random sampling comes from the bootstrap sampling done in the calculations of the SMA regressions themselves. This approach has the big advantages of 
  #(a) not having to make any assumptions about the distribution of the coefficient uncertainty and 
  #(b) ensuring that the individual slope coefficients within a regression are consistent with each other.
  for (ss in 1:n_uncer) {
    # For now, make samples for the following:
    # (ensure equation choices are consistent with decisions above)
    if (ss==1) { #First pass always calculates the best estimate
      if(use_LMA_from_TLP_LS) {  # BE:
        mod_LMA_intercept_sample <- LMA_from_TLP_LS$mod$intercept_R #LMA_from_TLP_LS
        mod_LMA_slope_y1_sample  <- LMA_from_TLP_LS$mod$slope_R.y1
        mod_LMA_slope_y2_sample  <- LMA_from_TLP_LS$mod$slope_R.y2
       
        
        mod_LS_intercept_sample  <- LS_from_P50_TLP_Ks_LMA$mod$intercept_R #LS_from_LMA TLP_Ks
        mod_LS_slope_y1_sample   <- LS_from_P50_TLP_Ks_LMA$mod$slope_R.y1 
        mod_LS_slope_y2_sample   <- LS_from_P50_TLP_Ks_LMA$mod$slope_R.y2 
        mod_LS_slope_y3_sample   <- LS_from_P50_TLP_Ks_LMA$mod$slope_R.y3 
        mod_LS_slope_y4_sample   <- LS_from_P50_TLP_Ks_LMA$mod$slope_R.y4 
        
       print(mod_LMA_intercept_sample)
       print(mod_LS_intercept_sample)
       print(regr_type)
        #no longer BE or BDT -specific: 
      #  mod_WD_intercept_sample   <- WD_from_P50_LMA_Ks$mod$intercept_R #WD_from_P50_LMA_Ks
      #  mod_WD_slope_y1_sample    <- WD_from_P50_LMA_Ks$mod$slope_R.y1
      #  mod_WD_slope_y2_sample    <- WD_from_P50_LMA_Ks$mod$slope_R.y2
      #  mod_WD_slope_y3_sample    <- WD_from_P50_LMA_Ks$mod$slope_R.y3
       
        
      } else {  # BDT:
        mod_LMA_intercept_sample <- LMA_from_TLP$mod$intercept_R #LMA_from_TLP
        mod_LMA_slope_y1_sample  <- LMA_from_TLP$mod$slope_R.y1
        
        mod_LS_intercept_sample  <- LS_from_P50_TLP_Ks$mod$intercept_R #LS_from_TLP_Ks
        mod_LS_slope_y1_sample   <- LS_from_P50_TLP_Ks$mod$slope_R.y1 
        mod_LS_slope_y2_sample   <- LS_from_P50_TLP_Ks$mod$slope_R.y2 
        
        #no longer BE or BDT -specific:
      # mod_WD_intercept_sample   <- WD_from_P50_Ks$mod$intercept_R #WD_from_P50_Ks
      #  mod_WD_slope_y1_sample    <- WD_from_P50_Ks$mod$slope_R.y1
      #  mod_WD_slope_y2_sample    <- WD_from_P50_Ks$mod$slope_R.y2
      
      }
      
      mod_WD_intercept_sample   <- WD_from_P50_Ks$mod$intercept_R #WD_from_P50_Ks
      mod_WD_slope_y1_sample    <- WD_from_P50_Ks$mod$slope_R.y1
      mod_WD_slope_y2_sample    <- WD_from_P50_Ks$mod$slope_R.y2
      
      mod_Ks_intercept_sample    <- Ks_from_P50_LS$mod$intercept #Ks_from_P50_LS
      mod_Ks_slope_y1_sample     <- Ks_from_P50_LS$mod$slope_R.y1
      mod_Ks_slope_y2_sample     <- Ks_from_P50_LS$mod$slope_R.y2

      mod_TLP_intercept_sample   <- TLP_from_LS_LMA_P50$mod$intercept_R #TLP_from_LS_LMA_P50
      mod_TLP_slope_y1_sample    <- TLP_from_LS_LMA_P50$mod$slope_R.y1
      mod_TLP_slope_y2_sample    <- TLP_from_LS_LMA_P50$mod$slope_R.y2
      mod_TLP_slope_y3_sample    <- TLP_from_LS_LMA_P50$mod$slope_R.y3
      
      mod_P50_intercept_sample   <- P50_from_TLP_Ks$mod$intercept_R #P50_from_TLP_Ks
      mod_P50_slope_y1_sample    <- P50_from_TLP_Ks$mod$slope_R.y1
      mod_P50_slope_y2_sample    <- P50_from_TLP_Ks$mod$slope_R.y2

      mod_slope_intercept_sample <-  slope_from_P50_TLP_WD$mod$intercept_R # slope_from_P50_TLP_WD
      mod_slope_slope_y1_sample <-  slope_from_P50_TLP_WD$mod$slope_R.y1
      mod_slope_slope_y2_sample <-  slope_from_P50_TLP_WD$mod$slope_R.y2
      mod_slope_slope_y3_sample <-  slope_from_P50_TLP_WD$mod$slope_R.y3
   
      
    } else {
      if(use_LMA_from_TLP_LS) {
        mod_LMA_intercept_sample <- LMA_from_TLP_LS$mod$boot.intercept[ss] #LMA_from_TLP_LS
        mod_LMA_slope_y1_sample <- LMA_from_TLP_LS$mod$boot.y1[ss]
        mod_LMA_slope_y2_sample <- LMA_from_TLP_LS$mod$boot.y2[ss]

        mod_LS_intercept_sample  <- LS_from_P50_TLP_Ks_LMA$mod$intercept_R[ss] #LS_from_LMA TLP_Ks
        mod_LS_slope_y1_sample   <- LS_from_P50_TLP_Ks_LMA$mod$slope_R.y1[ss] 
        mod_LS_slope_y2_sample   <- LS_from_P50_TLP_Ks_LMA$mod$slope_R.y2[ss] 
        mod_LS_slope_y3_sample   <- LS_from_P50_TLP_Ks_LMA$mod$slope_R.y3[ss] 
        mod_LS_slope_y4_sample   <- LS_from_P50_TLP_Ks_LMA$mod$slope_R.y4[ss] 
        
        #no longer BE or BDT -specific:
       # mod_WD_intercept_sample   <- WD_from_P50_LMA_Ks$mod$boot.intercept[ss] #WD_from_P50_LMA_Ks
      #  mod_WD_slope_y1_sample    <- WD_from_P50_LMA_Ks$mod$boot.y1[ss]
      #  mod_WD_slope_y2_sample    <- WD_from_P50_LMA_Ks$mod$boot.y2[ss]
      #  mod_WD_slope_y3_sample    <- WD_from_P50_LMA_Ks$mod$boot.y3[ss]
     
      } else {
        mod_LMA_intercept_sample <- LMA_from_TLP$mod$boot.intercept[ss] #LMA_from_TLP
        mod_LMA_slope_y1_sample  <- LMA_from_TLP$mod$boot.y1[ss]
        
        mod_LS_intercept_sample  <- LS_from_P50_TLP_Ks$mod$boot.intercept[ss] #LS_from_TLP_Ks
        mod_LS_slope_y1_sample   <- LS_from_P50_TLP_Ks$mod$boot.y1[ss] 
        mod_LS_slope_y2_sample   <- LS_from_P50_TLP_Ks$mod$boot.y2[ss] 
        
        #no longer BE or BDT -specific:
       # mod_WD_intercept_sample   <- WD_from_P50_Ks$mod$boot.intercept[ss] # WD_from_P50_slope_Ks#WD_from_slope_P50slope
      #  mod_WD_slope_y1_sample    <- WD_from_P50_Ks$mod$boot.y1[ss]
      #  mod_WD_slope_y2_sample    <- WD_from_P50_Ks$mod$boot.y2[ss]
       
      }
      
      mod_WD_intercept_sample   <- WD_from_P50_Ks$mod$boot.intercept[ss] # WD_from_P50_slope_Ks#WD_from_slope_P50slope
      mod_WD_slope_y1_sample    <- WD_from_P50_Ks$mod$boot.y1[ss]
      mod_WD_slope_y2_sample    <- WD_from_P50_Ks$mod$boot.y2[ss]
      
      mod_Ks_intercept_sample <- Ks_from_P50_LS$mod$boot.intercept[ss] #Ks_from_P50_LS
      mod_Ks_slope_y1_sample  <- Ks_from_P50_LS$mod$boot.y1[ss]
      mod_Ks_slope_y2_sample  <- Ks_from_P50_LS$mod$boot.y2[ss]
      
      mod_TLP_intercept_sample <- TLP_from_LS_LMA_P50$mod$boot.intercept[ss] #TLP_from_LS_LMA_P50
      mod_TLP_slope_y1_sample <- TLP_from_LS_LMA_P50$mod$boot.y1[ss]
      mod_TLP_slope_y2_sample <- TLP_from_LS_LMA_P50$mod$boot.y2[ss]
      mod_TLP_slope_y3_sample <- TLP_from_LS_LMA_P50$mod$boot.y3[ss]
      
      mod_P50_intercept_sample <- P50_from_TLP_Ks$mod$boot.intercept[ss] #P50_from_TLP_Ks
      mod_P50_slope_y1_sample  <- P50_from_TLP_Ks$mod$boot.y1[ss]
      mod_P50_slope_y2_sample  <- P50_from_TLP_Ks$mod$boot.y2[ss]

      mod_slope_intercept_sample <- slope_from_P50_TLP_WD$mod$boot.intercept[ss] #slope_from_P50_TLP_Ks
      mod_slope_slope_y1_sample <- slope_from_P50_TLP_WD$mod$boot.y1[ss]
      mod_slope_slope_y2_sample <- slope_from_P50_TLP_WD$mod$boot.y2[ss]
      mod_slope_slope_y3_sample <- slope_from_P50_TLP_WD$mod$boot.y3[ss]
     
      
    }
    # These regression coefficients will now be used in the optimisation calculations
    
    #TLP, P50, LMA need optimising
    
    #First set some initial based on simple bivariate (linear) relationship. This is just so that the iteration has somewhere to start from. Final result should not be sensitive to these.
    # no need to scale here for pcr and plsr, as the regression models here are derived from linear models.
     #P50_e_last = P50_from_Ks$mod$intercept_R  + P50_from_Ks$mod$slope_R.y1*Ks_e_start
    TLP_e_last = TLP_from_P50$mod$intercept_R  + TLP_from_P50$mod$slope_R.y1*P50_e_start
    LMA_e_last = LMA_from_LS$mod$intercept_R   + LMA_from_LS$mod$slope_R.y1*LS_e_start
    WD_e_last  = WD_from_LMA$mod$intercept_R   + WD_from_LMA$mod$slope_R.y1*LMA_e_last  # high pearson cor, assumed is no functional relationship
    slope_e_last  = slope_from_WD$mod$intercept_R  + slope_from_WD$mod$slope_R.y1*WD_e_last 
    Ks_e_last = Ks_from_LS$mod$intercept_R     + Ks_from_LS$mod$slope_R.y1*LS_e_start
    
    # "diff_" variables hold the difference between the current estimate of a trait value "_e" and the previous
    # estimate "_last"
    # "diff_*_last" variables contain the differences from the last round of iteration
    # (these are compared to differences in the current round of iteration to see if changes are smaller than
    # "tol" and therefore the iteration can stop)
    # Here we initialise the "diff_*_last" variables very high
    diff_P50_last=100
    diff_LMA_last=100
    diff_TLP_last=100
    diff_LS_last=100
    diff_Ks_last=100
    diff_WD_last=100
    diff_slope_last=100
    
    # These arrays are just for output, they store the values of every iteration for the current datapoint.
    # Useful for debugging and to check that convergence is working.
    # (only for debugging, can be commented out)
    P50_c <- matrix(NA, nrow = 100)
    LMA_c <- matrix(NA, nrow = 100)
    TLP_c <- matrix(NA, nrow = 100)
    LS_c <- matrix(NA, nrow = 100)
    Ks_c <- matrix(NA, nrow = 100)
    WD_c <- matrix(NA, nrow = 100)
    slope_c <- matrix(NA, nrow = 100)
    
    
    # Do the first iteration with the starting trait pair ('_start'). 
    # In the optimisation below, these traits also turn into predicted variables.
    if(use_LMA_from_TLP_LS) {
      
      ########
      #LMA_from_TLP_LS
      #scale
      if(regr_type=='pcr' || regr_type=='plsr'){
        print('scaling within the network, trait_opt_bivar_start')
      TLP_e_last <- scale_traits(TLP_e_last, "TLP_e_last", nlabels, traits_mean, traits_sd)
      LS_e_start <- scale_traits(LS_e_start, "LS_e_start", nlabels, traits_mean, traits_sd)
 
      }
      
     LMA_e_last = mod_LMA_intercept_sample + mod_LMA_slope_y1_sample*TLP_e_last +
        mod_LMA_slope_y2_sample*LS_e_start 

     #unscale
     if(regr_type=='pcr' || regr_type=='plsr'){
       TLP_e_last <- unscale_traits(TLP_e_last, "TLP_e_last", nlabels, traits_mean, traits_sd)
       LS_e_start <- unscale_traits(LS_e_start, "LS_e_start", nlabels, traits_mean, traits_sd)
     #  LMA_e_last = LMA_e_last + as.numeric(traits_mean_unscale['LMA'])
  
     }
     ########
     
     
      # LS_e_last = mod_LS_intercept_sample + mod_LS_slope_y1_sample*LMA_e_last +  mod_LS_slope_y2_sample*TLP_e_last +
      #   mod_LS_slope_y3_sample * Ks_e_start
      
      
    } else {
      
      ########
      #LMA from TLP
      #scale
      if(regr_type=='pcr' || regr_type=='plsr'){
        TLP_e_last <- scale_traits(TLP_e_last, "TLP_e_last", nlabels, traits_mean, traits_sd)
      }
      
      LMA_e_last = mod_LMA_intercept_sample + mod_LMA_slope_y1_sample*TLP_e_last
      
      #unscale
      if(regr_type=='pcr' || regr_type=='plsr'){
        TLP_e_last <- unscale_traits(TLP_e_last, "TLP_e_last", nlabels, traits_mean, traits_sd)
       # LMA_e_last = LMA_e_last + as.numeric(traits_mean_unscale['LMA'])
      }
      ########
      
      # LS_e_last = mod_LS_intercept_sample + mod_LS_slope_y1_sample*TLP_e_last +
      #    mod_LS_slope_y2_sample * Ks_e_start
      
      #WD_from_P50_Ks  moved out of loop
      #WD_e_last = mod_WD_intercept_sample  + mod_WD_slope_y1_sample*P50_e_start  +
      #  mod_WD_slope_y2_sample*Ks_e_last
    }
    
    ########
    #Ks_from_P50_LS
    #scale
    if(regr_type=='pcr' || regr_type=='plsr'){
      P50_e_start <- scale_traits(P50_e_start, "P50_e_start", nlabels, traits_mean, traits_sd)
      LS_e_start <- scale_traits(LS_e_start, "LS_e_start", nlabels, traits_mean, traits_sd)
    }
    
    Ks_e_last  = mod_Ks_intercept_sample + mod_Ks_slope_y1_sample*P50_e_start +  
      mod_Ks_slope_y2_sample* LS_e_start 
    
    #unscale
    if(regr_type=='pcr' || regr_type=='plsr'){
      P50_e_start <- unscale_traits(P50_e_start, "P50_e_start", nlabels, traits_mean, traits_sd)
      LS_e_start  <- unscale_traits(LS_e_start, "LS_e_start", nlabels, traits_mean, traits_sd)
     # Ks_e_last = Ks_e_last + as.numeric(traits_mean_unscale['Ks'])
  
    }
    ########
 
    ########
    #TLP_from_LS_LMA_P50
    #scale
    if(regr_type=='pcr' || regr_type=='plsr'){
      P50_e_start <- scale_traits(P50_e_start, "P50_e_start", nlabels, traits_mean, traits_sd)
      LS_e_start <- scale_traits(LS_e_start, "LS_e_start", nlabels, traits_mean, traits_sd)
      print('LMA_e_last')
      LMA_e_last <- scale_traits(LMA_e_last, "LMA_e_last", nlabels, traits_mean, traits_sd)
      print(LMA_e_last)
    }
    
    TLP_e_last = mod_TLP_intercept_sample + mod_TLP_slope_y1_sample*LS_e_start + 
      mod_TLP_slope_y2_sample*LMA_e_last + mod_TLP_slope_y3_sample*P50_e_start 
    
    #unscale
    if(regr_type=='pcr' || regr_type=='plsr'){
      P50_e_start <- unscale_traits(P50_e_start, "P50_e_start", nlabels, traits_mean, traits_sd)
      LS_e_start  <- unscale_traits(LS_e_start, "LS_e_start", nlabels, traits_mean, traits_sd)
      LMA_e_last  <- unscale_traits(LMA_e_last, "LMA_e_last", nlabels, traits_mean, traits_sd)
     # TLP_e_last = TLP_e_last + as.numeric(traits_mean_unscale['TLP'])
    }
    ########

   # P50_e_last = mod_P50_intercept_sample + mod_P50_slope_y1_sample*TLP_e_last +  mod_P50_slope_y2_sample*LS_e_start + 
   #   mod_P50_slope_y3_sample*Ks_e_start + mod_P50_slope_y4_sample*slope_e_last + mod_P50_slope_y5_sample*WD_e_last
    
    #slope_from_P50_TLP_WD  moved outside the loop
    #slope_e_last = mod_slope_intercept_sample + mod_slope_slope_y1_sample*P50_e_start + 
     # mod_slope_slope_y2_sample*TLP_e_last + mod_slope_slope_y3_sample*WD_e_last 
    
    # Now we start the optimisation loop. Trait values are iterated until the difference between trait
    # values on successive iterations is less than "tol".
    niter=0;
    while (T) {
      niter=niter+1 # Number of iterations completed
    
      # Make estimates of trait values based on the best SMA regressions (probably multivariate in most cases)
      # The estimates of traits in each iteration are based on the estimates of their predictor traits from the previous iteration
      if(use_LMA_from_TLP_LS) {
        
        ########
        #LMA_from_TLP_LS
        #scale
        if(regr_type=='pcr' || regr_type=='plsr'){
          TLP_e_last <- scale_traits(TLP_e_last, "TLP_e_last", nlabels, traits_mean, traits_sd)
          LS_e_start  <- scale_traits(LS_e_start, "LS_e_start", nlabels, traits_mean, traits_sd)
        
        }
        
        LMA_e[ss]= mod_LMA_intercept_sample + mod_LMA_slope_y1_sample*TLP_e_last +
          mod_LMA_slope_y2_sample*LS_e_start 
    
        #unscale
        if(regr_type=='pcr' || regr_type=='plsr'){
          TLP_e_last  <- unscale_traits(TLP_e_last, "TLP_e_last", nlabels, traits_mean, traits_sd)
          LS_e_start  <- unscale_traits(LS_e_start, "LS_e_start", nlabels, traits_mean, traits_sd)
        
         # LMA_e[ss] = LMA_e[ss] + as.numeric(traits_mean_unscale['LMA'])
        }
        #######
        
        # LS_e[ss] = mod_LS_intercept_sample + mod_LS_slope_y1_sample*LMA_e[ss] +  mod_LS_slope_y2_sample*TLP_e_last +
        #   mod_LS_slope_y3_sample * Ks_e_start
        
        #WD_from_P50_LMA_Ks  moved out of loop
       # WD_e[ss] = mod_WD_intercept_sample  + mod_WD_slope_y1_sample*P50_e_start  + 
       #   mod_WD_slope_y2_sample*LMA_e_last + mod_WD_slope_y3_sample*Ks_e_last
        
      } else {
        
        ########
        #LMA_from_TLP
        #scale
        if(regr_type=='pcr' || regr_type=='plsr'){
          TLP_e_last <- scale_traits(TLP_e_last, "TLP_e_last", nlabels, traits_mean, traits_sd)
        }
        
        LMA_e[ss] = mod_LMA_intercept_sample + mod_LMA_slope_y1_sample*TLP_e_last
        
        #unscale
        if(regr_type=='pcr' || regr_type=='plsr'){
          TLP_e_last <- unscale_traits(TLP_e_last, "TLP_e_last", nlabels, traits_mean, traits_sd)
         # LMA_e[ss] = LMA_e[ss] + as.numeric(traits_mean_unscale['LMA'])
        }
        ########
        # LS_e[ss] = mod_LS_intercept_sample + mod_LS_slope_y1_sample*TLP_e_last +
        #    mod_LS_slope_y2_sample * Ks_e_start
        
      }
      
      ########
      #TLP_from_LS_LMA_P50
      #scale
      if(regr_type=='pcr' || regr_type=='plsr'){
        P50_e_start <- scale_traits(P50_e_start, "P50_e_start", nlabels, traits_mean, traits_sd)
        LS_e_start <- scale_traits(LS_e_start, "LS_e_start", nlabels, traits_mean, traits_sd)
        LMA_e_last <- scale_traits(LMA_e_last, "LMA_e_last", nlabels, traits_mean, traits_sd)
        print('tlp scale')
      }
      
      TLP_e[ss] = mod_TLP_intercept_sample + mod_TLP_slope_y1_sample*LS_e_start + 
        mod_TLP_slope_y2_sample*LMA_e_last + mod_TLP_slope_y3_sample*P50_e_start 
      
      #unscale
      if(regr_type=='pcr' || regr_type=='plsr'){
        print('tlp unscale')
        P50_e_start <- unscale_traits(P50_e_start, "P50_e_start", nlabels, traits_mean, traits_sd)
        LS_e_start <- unscale_traits(LS_e_start, "LS_e_start", nlabels, traits_mean, traits_sd)
        LMA_e_last <- unscale_traits(LMA_e_last, "LMA_e_last", nlabels, traits_mean, traits_sd)
      #  TLP_e[ss] = TLP_e[ss] + as.numeric(traits_mean_unscale['TLP'])
      }
      ########
      
      
     # P50_e[ss] = mod_P50_intercept_sample + mod_P50_slope_y1_sample*TLP_e_last +  mod_P50_slope_y2_sample*LS_e_start + 
      #  mod_P50_slope_y3_sample*Ks_e_start + mod_P50_slope_y4_sample*slope_e_last + mod_P50_slope_y5_sample*WD_e_last
     
      ########
       #Ks_from_P50_LS_WD
      #scale
      if(regr_type=='pcr' || regr_type=='plsr'){
        P50_e_start <- scale_traits(P50_e_start, "P50_e_start", nlabels, traits_mean, traits_sd)
        LS_e_start <- scale_traits(LS_e_start, "LS_e_start", nlabels, traits_mean, traits_sd)
        print('ks unscale')
        }
     
      Ks_e[ss] = mod_Ks_intercept_sample + mod_Ks_slope_y1_sample*P50_e_start +  
          mod_Ks_slope_y2_sample* LS_e_start 
      
      #unscale
      if(regr_type=='pcr' || regr_type=='plsr'){
        print('ks unscale')
        P50_e_start <- unscale_traits(P50_e_start, "P50_e_start", nlabels, traits_mean, traits_sd)
        LS_e_start <- unscale_traits(LS_e_start, "LS_e_start", nlabels, traits_mean, traits_sd)
        #Ks_e[ss] = Ks_e[ss] + as.numeric(traits_mean_unscale['Ks'])
        }
  
        ########
      
      if (limitdataranges) {
        #Do not go beyond observed limits of data - if so, discard.
        #if (P50_e[ss]>maxP50 | is.na(P50_e[ss])) {P50_e[ss]=NA; break}
        #if (P50_e[ss]<minP50 | is.na(P50_e[ss])) {P50_e[ss]=NA; break}
        if (TLP_e[ss]>maxTLP | is.na(TLP_e[ss])) {TLP_e[ss]=NA; break}
        if (TLP_e[ss]<minTLP | is.na(TLP_e[ss])) {TLP_e[ss]=NA; break}
        if (LMA_e[ss]>maxLMA | is.na(LMA_e[ss])) {LMA_e[ss]=NA; break}
        if (LMA_e[ss]<minLMA | is.na(LMA_e[ss])) {LMA_e[ss]=NA; break}
        #  if (LS_e[ss]>maxLS | is.na(LS_e[ss])) {LS_e[ss]=NA; break}
        #  if (LS_e[ss]<minLS | is.na(LS_e[ss])) {LS_e[ss]=NA; break}
         if (Ks_e[ss]>maxKs | is.na(Ks_e[ss])) {Ks_e[ss]=NA; break}
         if (Ks_e[ss]<minKs | is.na(Ks_e[ss])) {Ks_e[ss]=NA; break}
      #  if (slope_e[ss]>maxslope | is.na(slope_e[ss])) {slope_e[ss]=NA; break}
      #  if (slope_e[ss]<minslope | is.na(slope_e[ss])) {slope_e[ss]=NA; break}
       # if (WD_e[ss]>maxWD | is.na(WD_e[ss])) {WD_e[ss]=NA; break}
       # if (WD_e[ss]<minWD | is.na(WD_e[ss])) {WD_e[ss]=NA; break}
      }
      
      # Save the values for this iteration to the output array (only for debugging, can be commented out)
      #P50_c[niter] <- P50_e[ss]
      LMA_c[niter] <- LMA_e[ss]
      TLP_c[niter] <- TLP_e[ss]
      #LS_c[niter] <- LS_e[ss]
      Ks_c[niter] <- Ks_e[ss]
      WD_c[niter] <- WD_e[ss]
      slope_c[niter] <- slope_e[ss]
      
      # Calculate the difference between the current estimate of a trait value "_e" and the previous estimate "_last"
     # diff_P50 = P50_e[ss]-P50_e_last
      diff_LMA = LMA_e[ss]-LMA_e_last
      diff_TLP = TLP_e[ss]-TLP_e_last
      #  diff_LS  = LS_e[ss] -LS_e_last
      diff_Ks  = Ks_e[ss] -Ks_e_last
    #  diff_WD  = WD_e[ss] -WD_e_last
    #  diff_slope  = slope_e[ss] -slope_e_last
      
      # Now we test if the difference between trait estimates on this iteration and between trait estimates on
      # the last iteration is less than "tol" for all traits. If it is we finish the iteration.
      if (#abs(diff_P50-diff_P50_last)<tol &&
          abs(diff_LMA-diff_LMA_last)<tol &&
          abs(diff_TLP-diff_TLP_last)<tol &&
         # abs(diff_LS-diff_LS_last)<tol &&
          abs(diff_Ks-diff_Ks_last)<tol
        #  abs(diff_WD-diff_WD_last)<tol 
        #  abs(diff_slope-diff_slope_last)<tol 
        ) {
        break
      }
      
      # Save the "diff" values ready for the next iteration
      #diff_P50_last=diff_P50
      diff_LMA_last=diff_LMA
      diff_TLP_last=diff_TLP
      # diff_LS_last=diff_LS
      diff_Ks_last=diff_Ks
    #  diff_WD_last=diff_WD
     # diff_slope_last=diff_slope
      
      # Save the "_e" values ready for the next iteration
      #P50_e_last=P50_e[ss]
      LMA_e_last=LMA_e[ss]
      TLP_e_last=TLP_e[ss]
      # LS_e_last =LS_e[ss]
      Ks_e_last =Ks_e[ss]
    #  WD_e_last =WD_e[ss]
     # slope_e_last =slope_e[ss]
    }
    
    # After the iteration has finished we can calculate any traits which did not need to be included in the optimisation (because they are not used in the input to calculate any other trait)
    
    #  slope_e[ss]=mod_slope_intercept_sample + mod_slope_slope_y1_sample*P50_e[ss] + 
    #    mod_slope_slope_y2_sample*TLP_e[ss] + mod_slope_slope_y3_sample*WD_e[ss]
    
    #  WD_e[ss]=mod_WD_intercept_sample  + mod_WD_slope_y1_sample*P50_e[ss]  + mod_WD_slope_y3_sample*slope_e[ss]+ 
    #    mod_WD_slope_y2_sample*Ks_e[ss] 
    
    
    ########
    #WD_from_P50_Ks
    #scale
    if(regr_type=='pcr' || regr_type=='plsr'){
      P50_e_start <- scale_traits(P50_e_start, "P50_e_start", nlabels, traits_mean, traits_sd)
      Ks_e[ss]    <- scale_traits(Ks_e[ss], "Ks_e", nlabels, traits_mean, traits_sd)
    }
    
    WD_e[ss] = mod_WD_intercept_sample  + mod_WD_slope_y1_sample*P50_e_start  + 
      mod_WD_slope_y2_sample*Ks_e[ss]
    
    #unscale
    if(regr_type=='pcr' || regr_type=='plsr'){
      P50_e_start <- unscale_traits(P50_e_start, "P50_e_start", nlabels, traits_mean, traits_sd)
      Ks_e[ss]    <- unscale_traits(Ks_e[ss], "Ks_e", nlabels, traits_mean, traits_sd)
      #WD_e[ss] = WD_e[ss] + as.numeric(traits_mean_unscale['WD'])
    }
    ########
    #[TODO] check all this is triggered
    ########
    #slope_from_P50_TLP_WD
    #scale
    if(regr_type=='pcr' || regr_type=='plsr'){
      P50_e_start <- scale_traits(P50_e_start, "P50_e_start", nlabels, traits_mean, traits_sd)
      TLP_e[ss]   <- scale_traits(TLP_e[ss], "TLP_e", nlabels, traits_mean, traits_sd)
      WD_e[ss]    <- scale_traits(WD_e[ss], "WD_e", nlabels, traits_mean, traits_sd)
    }
    
    slope_e[ss] = mod_slope_intercept_sample + mod_slope_slope_y1_sample*P50_e_start + 
      mod_slope_slope_y2_sample*TLP_e[ss]  + mod_slope_slope_y3_sample*WD_e[ss] 
    
    #unscale
    if(regr_type=='pcr' || regr_type=='plsr'){
      P50_e_start <- unscale_traits(P50_e_start, "P50_e_start", nlabels, traits_mean, traits_sd)
      TLP_e[ss]   <- unscale_traits(TLP_e[ss], "TLP_e", nlabels, traits_mean, traits_sd)
      WD_e[ss]    <- unscale_traits(WD_e[ss], "WD_e", nlabels, traits_mean, traits_sd)
     # slope_e[ss] = slope_e[ss] + as.numeric(traits_mean_unscale['slope'])
    }
    ########
    
    
    
    #  if( is.na(slope_e))print(dd)
    # WD_e[ss]=mod_WD_intercept_sample + mod_WD_slope_y1_sample*slope_e[ss] + 
    #  mod_WD_slope_y2_sample*slope_e[ss]*P50_e[ss]
    
    
    if (limitdataranges) {
    #  #Do not go beyond observed limits of data
      if (slope_e[ss]>maxslope | is.na(slope_e[ss])) {slope_e[ss]=NA}
      if (slope_e[ss]<minslope | is.na(slope_e[ss])) {slope_e[ss]=NA}
      if (WD_e[ss]>maxWD | is.na(WD_e[ss])) {WD_e[ss]=NA}
      if (WD_e[ss]<minWD | is.na(WD_e[ss])) {WD_e[ss]=NA}
    }
    
  } #Finish nbtstrp loop
  
  return_vals <- list("P50_e"=P50_e_start,"TLP_e"=TLP_e,"LMA_e"=LMA_e,"WD_e"=WD_e,"slope_e"=slope_e,'LS_e'= LS_e_start,'Ks_e'=Ks_e)
  
  return(return_vals)
}

