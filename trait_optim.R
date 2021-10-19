# Annemarie Eckes-Shephard moved this code from T. Pugh into this function ( 2021 May )
# this function subsets the data, chooses options for the type of sampling (random or strategic subset vs. use all Ks LS combinations)
# it then carries out the optimisation. The number of outputs optimised values depends on the subset options parsed into the function
# if propagate_uncer=T, must provide nbtsrp value
# This function used to be part of the main analysis script, but is now outsourced
trait_optim <- function(limitdataranges=T, propagate_uncer=T, nbtstrp=1000, trait_sel= T, n_trait_sel=28, spec_group_sel= 3,est_lhs ){
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
  
  # Identify all combinations of Ks and LS, P50 and TLP(do this across full range of broadleaf species)
  ind=which(!is.na(traits$Ks) & !is.na(traits$LS))
  
  # sample across a uniform distribution of min max of those traits
  LS_comb <- traits$LS[ind]
  Ks_comb <- traits$Ks[ind]

 
  
  if (trait_sel) {
    if (n_trait_sel>0) {
      # Random selection of LS and Ks values to be tested
      set.seed(1234)
      index = 1:length(ind)
      trait_samp = sample(index, n_trait_sel, replace=F) 
      LS_e_start = LS_comb[trait_samp] 
      Ks_e_start = Ks_comb[trait_samp] 
      # Plot sample against all data as a check for sampling density ( !=4 suppresses plot as this is now the lhcube option)
      if(n_trait_sel != 4){
        plot(LS_comb,Ks_comb)
        points(LS_e_start,Ks_e_start,col="red")
      }
     
      
      if (n_trait_sel == 4){
        LS_e_start  <- est_lhs$LS_e
        Ks_e_start <- est_lhs$Ks_e  
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
      set.seed(293)
      # Have not rescaled trait before fitting hypervolume as the ranges of both are very similar
      hv = hypervolume(data.frame(LS_comb,Ks_comb),method="gaussian",quantile.requested=0.95)
      #plot(hv)
   
     
      # Set the number of points distributed systematically across LS and Ks space to test for inclusion in the hypervolume
      sampKs=8
      sampLS=8
      
      maxLS=max(LS_comb,na.rm=T)
      minLS=min(LS_comb,na.rm=T)
      maxKs=max(Ks_comb,na.rm=T)
      minKs=min(Ks_comb,na.rm=T)
      intLS=(maxLS-minLS)/sampLS
      intKs=(maxKs-minKs)/sampKs
      
      LS_seq <- seq(minLS+intLS,maxLS-intLS,by=intLS)
      Ks_seq <- seq(minKs+intKs,maxKs-intKs,by=intKs)
      LS_Ks_seq <- expand.grid(LS_seq,Ks_seq)
      LS_Ks_seq2<- expand.grid(LS_seq,Ks_seq)
      # Test those points for inclusion in the hypervolume
      in_hv <- hypervolume_inclusion_test(hv,LS_Ks_seq,fast.or.accurate = "accurate")
      
      LS_e_start <- LS_Ks_seq$Var1[in_hv]
      length( LS_e_start )
      Ks_e_start <- LS_Ks_seq$Var2[in_hv]
      length( Ks_e_start )
  
      plot_Hypervolume(hv,LS_e_start,Ks_e_start)
    }
  } else {
    # Go through all observed combinations of KS and LS
    Ks_e_start=Ks_comb
    LS_e_start=LS_comb
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
  
  ndata=length(LS_e_start)
  
  P50_e <- matrix(NA, nrow= ndata, ncol = n_uncer) #Array now expanded to hold multiple replicate estimates based on regression coefficient uncertainty
  LMA_e <- matrix(NA, nrow= ndata, ncol = n_uncer)
  TLP_e <- matrix(NA, nrow= ndata, ncol = n_uncer)
  WD_e <- matrix(NA, nrow= ndata, ncol = n_uncer)
  slope_e <- matrix(NA, nrow= ndata, ncol = n_uncer)
  LS_e <- matrix(NA, nrow= ndata, ncol = n_uncer)
  Ks_e <- matrix(NA, nrow= ndata, ncol = n_uncer)
  
  # Loop over all the combinations of KS and LS
  # The new estimates of traits use the suffix "_e"
  for (dd in 1:ndata) {
    print(dd)
    
    # Carry out the optimisation - based on Ks_e_start and LS_e_start as starting points
    opt_vals <- trait_opt(traits$P50[ind_spec_group],
                          traits$TLP[ind_spec_group],
                          traits$LMA[ind_spec_group],
                          traits$WD[ind_spec_group],
                          traits$slope[ind_spec_group],
                          traits$LS[ind_spec_group],
                          traits$Ks[ind_spec_group],
                          LMA_multivar_BDT$LMA_from_TLP,
                          LMA_multivar_BE$LMA_from_TLP_LS_WD,
                          LS_multivar_BE$LS_from_P50_TLP_Ks, 
                          LS_multivar_BDT$LS_from_TLP_Ks,
                          Ks_multivar$Ks_from_P50_LS_slope_WD,#Ks_from_P50_LS_slope,
                          TLP_multivar$TLP_from_LS_LMA_P50_slope,
                          P50_multivar$P50_from_TLP_LS_Ks_slope_WD,
                          slope_multivar$slope_from_P50_TLP_WD_Ks,#slope_from_P50_TLP_Ks,
                          WD_multivar_BDT$WD_from_P50_slope_Ks,#WD_from_Ks_P50,#WD_from_slope_P50slope,
                          WD_multivar_BE$WD_from_P50_slope_Ks_LMA,
                          bivar$LMA_from_LS,
                          bivar$P50_from_Ks,
                          bivar$TLP_from_P50,
                          bivar$WD_from_LMA,
                          bivar$slope_from_WD,
                          bivar$Ks_from_LS,
                          bivar$LS_from_LMA, # not technically necessary here, but for consistency between functions.
                          Ks_e_start[dd],
                          LS_e_start[dd],
                          n_uncer,
                          use_LMA_from_TLP_LS)
    

    P50_e[dd,]   <- opt_vals$P50_e
    TLP_e[dd,]   <- opt_vals$TLP_e
    LMA_e[dd,]   <- opt_vals$LMA_e
    WD_e[dd,]    <- opt_vals$WD_e
    slope_e[dd,] <- opt_vals$slope_e
    LS_e[dd,]    <- opt_vals$LS # they are not _e estimated values, but the _start values are reported here. name used for easier downstream plotting.!! dangerous..
    Ks_e[dd,]    <- opt_vals$Ks
    
  }
 
  # discard values across all trait lists where observed trait limits were surpassed (labelled as NA within in function trait_opt_bivar_start_LSP50)
  # this will reduce the number of predictor and predicted points in the lists
  predicted.df <- data.frame("P50_e"=P50_e,"TLP_e"=TLP_e,"LMA_e"=LMA_e,"WD_e"=WD_e,"slope_e"=slope_e)
  ind = complete.cases(predicted.df)
  
  #re-define list elements as matrix after sub-setting ([ind]) so that subsequent functions in analysis still work:
  predicted <- list("TLP_e"= as.matrix(TLP_e[ind,]),"P50_e" = as.matrix(P50_e[ind,]),"LMA_e" = as.matrix(LMA_e[ind,]),"WD_e" = as.matrix(WD_e[ind,]),
                    "slope_e"= as.matrix(slope_e[ind,]))
  predictors <- list('Ks_e' = as.matrix(Ks_e_start[ind]),'LS_e'  = as.matrix(LS_e_start[ind]))
  return_vals <- list('predictors'=predictors ,'predicted' =predicted )
  
  #inform about the extend to which data was discarded during the optimisation:
  nkeep <- count(ind)[which(count(ind)$x=='TRUE'),'freq']
  ndiscard <- count(ind)[which(count(ind)$x=='FALSE'),'freq']
  
  if(length(ndiscard) != 0){
  print(paste('The optimisation started with ',dim(predicted.df)[1],'predictor values and keept ',nkeep,'values and discarded ',ndiscard))
  if(ndiscard>1)print('This means that predicted values could not converge, as some fell outside the range of observations')
  }
  
  return(return_vals) 
}
