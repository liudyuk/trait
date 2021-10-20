# this function subsets the data, chooses options for the type of sampling (random or strategic subset vs. use all Ks LS combinations)
# it then carries out the optimisation. The number of outputs optimised values depends on the subset options parsed into the function
# if propagate_uncer=T, must provide nbtsrp value
trait_optim_bivar_start_KsTLP <- function(limitdataranges=T, propagate_uncer=T, nbtstrp=1000, trait_sel= T, n_trait_sel=28, spec_group_sel= 3,est_lhs){
  # ---
  if (propagate_uncer) {
    n_uncer=nbtstrp
  } else {
    n_uncer=1
  }
  
  # Get index for selected species group
  if (spec_group_sel==1) {
    ind_spec_group=which(traits$group=='BT' | traits$group=='BD')
  } else if (spec_group_sel==2) {
    ind_spec_group=which(traits$group=='BE')
  } else if (spec_group_sel==3) {
    ind_spec_group=which(traits$group=='BT')
  } else if (spec_group_sel==4) { 
    ind_spec_group=which(traits$group=='BD')
  }
  
  # Identify all combinations of Ks and LS (do this across full range of broadleaf species)
  ind=which(!is.na(traits$Ks) & !is.na(traits$TLP))
  
  Ks_comb <- traits$Ks[ind]
  TLP_comb <- traits$TLP[ind]
  
  if (trait_sel) {
    if (n_trait_sel>0) {
      # Random selection of LS and Ks values to be tested
      set.seed(1234)
      index = 1:length(ind)
      trait_samp = sample(index, n_trait_sel, replace=F) 
      Ks_e_start = Ks_comb[trait_samp] 
      TLP_e_start= TLP_comb[trait_samp] 
      
      # Plot sample against all data as a check for sampling density
      if (n_trait_sel != 4){
        plot(Ks_comb,TLP_comb)
        points(Ks_e_start,TLP_e_start,col="red")
      }
      
      #accommodate latin hypercube sampling option:
      if (n_trait_sel == 4){
        Ks_e_start  <- est_lhs$Ks_e_start
        TLP_e_start<- est_lhs$TLP_e_start 
      }
    } else {
      # Systematic sample
      
      #install.packages("hypervolume")
      library(hypervolume)
      # Fit a hypervolume (KDE at 95%)
      
      # The hypervolume algorithm uses a stochastic 'dart throwing' algorithm to determine the topography of the data. 
      # This stochasticity can sometimes lead to the hypervolume edges changing between code runs.
      # In order to reliably obtain a systematic sample 28 PFTs from the sampling space, set.seed is used
      # https://benjaminblonder.org/hypervolume_faq.html
      set.seed(1)
      # Have not rescaled trait before fitting hypervolume as the ranges of both are very similar
      hv = hypervolume(data.frame(Ks_comb,TLP_comb),method="gaussian",quantile.requested=0.95)
      # plot(hv)
      
      # Set the number of points distributed systematically across LS and TLP space to test for inclusion in the hypervolume
      sampKs=7
      sampTLP=8
      
      maxKs=max(Ks_comb,na.rm=T)
      minKs=min(Ks_comb,na.rm=T)
      maxTLP=max(TLP_comb,na.rm=T)
      minTLP=min(TLP_comb,na.rm=T)
      intKs=(maxKs-minKs)/sampKs
      intTLP=(maxTLP-minTLP)/sampTLP
      
      Ks_seq <- seq(minKs+intKs,maxKs-intKs,by=intKs)
      TLP_seq <- seq(minTLP+intTLP,maxTLP-intTLP,by=intTLP)
      Ks_TLP_seq <- expand.grid(Ks_seq,TLP_seq)
      Ks_TLP_seq2<- expand.grid(Ks_seq,TLP_seq)
      # Test those points for inclusion in the hypervolume
      in_hv <- hypervolume_inclusion_test(hv,Ks_TLP_seq,fast.or.accurate = "accurate")
      
      Ks_e_start <-  Ks_TLP_seq$Var1[in_hv]
      length( Ks_e_start )
      TLP_e_start<-  Ks_TLP_seq$Var2[in_hv]
      length( TLP_e_start)
      plot_Hypervolume(hv,Ks_e_start,TLP_e_start)
    }
  } else {
    # Go through all observed combinations of TLP and LS
    TLP_e_start=TLP_comb
    Ks_e_start=Ks_comb
  }
  
  # ---
  # Select the LMA and LS relationship to use
  if (spec_group_sel==1) {
    use_LMA_from_TLP_LS=F
    LS_from_TLP_Ks = F
  } else if (spec_group_sel==2) {
    use_LMA_from_TLP_LS=T
    LS_from_TLP_Ks = T
  } else if (spec_group_sel==3) {
    use_LMA_from_TLP_LS=F
    LS_from_TLP_Ks = F
  } else if (spec_group_sel==4) {
    use_LMA_from_TLP_LS=F
    LS_from_TLP_Ks = F
  }
  
  # ---
  # Do the optimisation
  
  ndata=length(Ks_e_start)
  
  
  P50_e <- matrix(NA, nrow= ndata, ncol = n_uncer) #Array now expanded to hold multiple replicate estimates based on regression coefficient uncertainty
  LMA_e <- matrix(NA, nrow= ndata, ncol = n_uncer)
  TLP_e <- matrix(NA, nrow= ndata, ncol = n_uncer)
  WD_e <- matrix(NA, nrow= ndata, ncol = n_uncer)
  slope_e <- matrix(NA, nrow= ndata, ncol = n_uncer)
  LS_e <- matrix(NA, nrow= ndata, ncol = n_uncer)
  Ks_e <- matrix(NA, nrow= ndata, ncol = n_uncer)
  
  
  # Loop over all the combinations of TLP and Ks
  # The new estimates of traits use the suffix "_e"
  for (dd in 1:ndata) {
    print(dd)
    
    # Carry out the optimisation - based on Ks_e_start and TLP_e_start as starting points
    opt_vals <- trait_opt_bivar_start_KsTLP(traits$P50[ind_spec_group],
                                            traits$TLP[ind_spec_group],
                                            traits$LMA[ind_spec_group],
                                            traits$WD[ind_spec_group],
                                            traits$slope[ind_spec_group],
                                            traits$LS[ind_spec_group],
                                            traits$Ks[ind_spec_group],
                                            LMA_multivar_BDT$LMA_from_TLP,
                                            LMA_multivar_BE$LMA_from_TLP_LS_WD,
                                            LS_multivar_BE$LS_from_P50_TLP_Ks_LMA,
                                            LS_multivar_BDT$LS_from_P50_TLP_Ks,
                                            Ks_multivar$Ks_from_P50_LS_WD,#Ks_from_P50_LS_slope,
                                            TLP_multivar$TLP_from_LS_LMA_P50,
                                            P50_multivar$P50_from_TLP_Ks_WD,
                                            slope_multivar$slope_from_P50_TLP_WD_Ks,#slope_from_P50_TLP_Ks,
                                            WD_multivar_BDT$WD_from_P50_slope_Ks,#WD_from_Ks_P50,#WD_from_slope_P50slope,
                                            WD_multivar_BE$WD_from_P50_slope_Ks_LMA,
                                            bivar$LMA_from_LS,
                                            bivar$P50_from_Ks,
                                            bivar$TLP_from_P50,
                                            bivar$WD_from_LMA,
                                            bivar$slope_from_WD,
                                            bivar$Ks_from_LS,
                                            bivar$LS_from_LMA,
                                            Ks_e_start[dd],
                                            TLP_e_start[dd],
                                            n_uncer,
                                            use_LMA_from_TLP_LS)
    
    
    P50_e[dd,]   <- opt_vals$P50_e
    TLP_e[dd,]   <- opt_vals$TLP_e# they are not _e estimated values, but the _start values are reported here. name used for easier downstream plotting.!! dangerous..
    LMA_e[dd,]   <- opt_vals$LMA_e
    WD_e[dd,]    <- opt_vals$WD_e
    slope_e[dd,] <- opt_vals$slope_e
    LS_e[dd,]    <- opt_vals$LS 
    Ks_e[dd,]    <- opt_vals$Ks# they are not _e estimated values, but the _start values are reported here. name used for easier downstream plotting.!! dangerous..
    
  }
  
  
  # discard values across all trait lists where observed trait limits were surpassed (labelled as NA within function trait_opt_bivar_start_LSP50)
  # this will reduce the number of predictor and predicted points in the lists
  predicted.df <- data.frame("TLP_e"=TLP_e,"Ks_e"=Ks_e,"LMA_e"=LMA_e,"WD_e"=WD_e,"slope_e"=slope_e,"LS_e" =LS_e,"P50_e" =P50_e)
  ind = complete.cases(predicted.df)
  
  #re-define list elements as matrix after sub-setting ([ind]) so that subsequent functions in analysis still work:
  predicted <- list("TLP_e"=as.matrix(TLP_e[ind,]),"Ks_e"=as.matrix(Ks_e[ind,]),"LMA_e"=as.matrix(LMA_e[ind,]),"WD_e"=as.matrix(WD_e[ind,]),"slope_e"=as.matrix(slope_e[ind,]),
                    "LS_e" =as.matrix(LS_e[ind,]),"P50_e" =as.matrix(P50_e[ind,]))
  predictors <- list('Ks_e_start' = as.matrix(Ks_e_start[ind]),'TLP_e_start' = as.matrix(TLP_e_start[ind]))
  return_vals <- list('predictors' = predictors ,'predicted' = predicted )
  
  #inform about the extend to which data was discarded during the optimisation:
  nkeep <- count(ind)[which(count(ind)$x=='TRUE'),'freq']
  ndiscard <- count(ind)[which(count(ind)$x=='FALSE'),'freq']
  
  if(length(ndiscard) != 0){
    print(paste('The optimisation started with ',dim(predicted.df)[1],'predictor values and keept ',nkeep,'values and discarded ',ndiscard))
    print('This means that predicted values could not converge, as some fell outside the range of observations')
  }
  return(return_vals)
}
