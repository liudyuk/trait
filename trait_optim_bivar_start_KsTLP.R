# this function subsets the data, chooses options for the type of sampling (random or strategic subset vs. use all Ks LS combinations)
# it then carries out the optimisation. The number of outputs optimised values depends on the subset options parsed into the function
# if propagate_uncer=T, must provide nbtsrp value
trait_optim_bivar_start_KsTLP <- function(limitdataranges=T, propagate_uncer=T, nbtstrp=1000, trait_sel= T, n_trait_sel=28, spec_group_sel= 3){
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
      Ks_e = Ks_comb[trait_samp] 
      TLP_e = TLP_comb[trait_samp] 
      # Plot sample against all data as a check for sampling density
      plot(Ks_comb,TLP_comb)
      points(Ks_e,TLP_e,col="red")
      
    } else {
      # Systematic sample
      
      #install.packages("hypervolume")
      library(hypervolume)
      # Fit a hypervolume (KDE at 95%)
      
      # The hypervolume algorithm uses a stochastic 'dart throwing' algorithm to determine the topography of the data. 
      # This stochasticity can sometimes lead to the hypervolume edges changing between code runs.
      # In order to reliably obtain a systematic sample 28 PFTs from the sampling space, set.seed is used
      # https://benjaminblonder.org/hypervolume_faq.html
      set.seed(12)
      # Have not rescaled trait before fitting hypervolume as the ranges of both are very similar
      hv = hypervolume(data.frame(Ks_comb,TLP_comb),method="gaussian",quantile.requested=0.95)
      plot(hv)
      
      # Set the number of points distributed systematically across LS and TLP space to test for inclusion in the hypervolume
      sampKs=8
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
      
      Ks_e <-  Ks_TLP_seq$Var1[in_hv]
      length( Ks_e )
      TLP_e <-  Ks_TLP_seq$Var2[in_hv]
      length( TLP_e )
    }
  } else {
    # Go through all observed combinations of TLP and LS
    TLP_e=TLP_comb
    Ks_e=Ks_comb
  }
  
  # ---
  # Select the LMA relationship to use
  if (spec_group_sel==1) {
    use_LMA_from_TLP_LS=F
  } else if (spec_group_sel==2) {
    use_LMA_from_TLP_LS=T
  } else if (spec_group_sel==3) {
    use_LMA_from_TLP_LS=F
  } else if (spec_group_sel==4) {
    use_LMA_from_TLP_LS=F
  }
  
  # ---
  # Do the optimisation
  
  ndata=length(Ks_e)
  
  
  LS_e <- matrix(NA, nrow= ndata, ncol = n_uncer) #Array now expanded to hold multiple replicate estimates based on regression coefficient uncertainty
  LMA_e <- matrix(NA, nrow= ndata, ncol = n_uncer)
  WD_e <- matrix(NA, nrow= ndata, ncol = n_uncer)
  slope_e <- matrix(NA, nrow= ndata, ncol = n_uncer)
  P50_e <- matrix(NA, nrow= ndata, ncol = n_uncer)
  
  # Loop over all the combinations of TLP and Ks
  # The new estimates of traits use the suffix "_e"
  for (dd in 1:ndata) {
    print(dd)
    
    # Carry out the optimisation - based on TLP_e and KS_e as starting points
    
    opt_vals <- trait_opt_bivar_start_KsTLP(
      traits$P50[ind_spec_group],
      traits$LS[ind_spec_group],
      traits$LMA[ind_spec_group],
      traits$WD[ind_spec_group],
      traits$slope[ind_spec_group],
      LMA_multivar_BDT$LMA_from_TLP,
      LMA_multivar_BE$LMA_from_TLP_LS,
      LS_multivar$LS_from_TLP_Ks,
      P50_multivar$P50_from_TLP_Ks,
      slope_multivar$slope_from_P50_TLP_Ks,
      WD_multivar$WD_from_slope_P50slope,
      bivar$LMA_from_TLP,
      bivar$P50_from_Ks,
      bivar$LS_from_TLP,
      TLP_e[dd],
      Ks_e[dd],
      n_uncer,
      use_LMA_from_TLP_LS)
    
    
    P50_e[dd,] <- opt_vals$P50_e
    LS_e[dd,] <- opt_vals$LS_e
    LMA_e[dd,] <- opt_vals$LMA_e
    WD_e[dd,] <- opt_vals$WD_e
    slope_e[dd,] <- opt_vals$slope_e
  }
  
  predicted <- list("P50_e"=P50_e,"LS_e"=LS_e,"LMA_e"=LMA_e,"WD_e"=WD_e,"slope_e"=slope_e)
  predictors <- list('TLP_e' = TLP_e,'Ks_e'  = Ks_e)
  return_vals <- list('predictors'=predictors ,'predicted' =predicted )
  
  return(return_vals)
}
