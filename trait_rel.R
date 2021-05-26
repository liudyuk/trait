# Script to read in the processed trait data, make bivariate and multivariate SMA regressions and then
# carry out an optimisation procedure to unify inter-trait relationships across the whole trait dataset.
#
# Dependencies: 
# - sma_multivar_regress.R
# - trait_functions.R
# - make_bivar_plots.R
# - multivar_model_selection.R
# - whittaker_biomes_plot.R
# - trait_opt.R
# - opt_test_plots.R
# - lpjg_traits_conv.R
# - hypervolume package (if taking a systematic sample)
#
# T. Pugh
# 15.06.20
#
# Annemarie Eckes-Shephard
# May 2021
# *Minor changes to make output reproducible.
# *outsource some code bits into separate functions
# *modify existing functions to test optimisation efficacy when starting with other variable combinations (not Ks LS)


nbtstrp=1000 # Number of bootstrap samples to take in sma_multivar_regress (samples later used to calculated uncertainty in the optimisation). Was previously 10 000, using a lower number for testing, Will need to check sensitivity to this value.

#traits=read.table("/Users/liudy/trait_data/woody_trait.0625.txt")
#traits=read.csv("/Users/pughtam/Documents/TreeMort/Analyses/Hydraulic_modelling/Traits/mytrait-data/woody_trait.0803.txt",sep="\t")
#traits=read.csv("/Users/pughtam/Documents/TreeMort/Analyses/Hydraulic_modelling/Traits/woody_trait.0827.txt",sep="\t")
traits = read.csv("/Users/annemarie/Documents/1_TreeMort/2_Analysis/2_data/2_intermediate/woody_trait.0827.txt",sep="\t")

source('sma_multivar_regress.R')
source('trait_functions.R')
source('make_bivar_plots.R')
source('multivar_model_selection.R')
source("whittaker_biomes_plot.R")
source("trait_opt.R")
source("opt_test_plots.R")
source('lpjg_traits_conv.R')


source('create_uncertainty_range_stats.R')
source('trait_optim.R')
source('trait_optim_bivar_startLSTLP.R') #for testing multiple starting variables for the optimisation
source('trait_opt_bivar_start_LSTLP.R')
source('opt_test_plots_LSTLP.R')

#--- Read in the trait data ---

# Subset only the broadleaf species
trait_B<-subset(traits,group!="CC",drop = T)
#trait_B<-droplevels(trait_B)
#str(trait_B)
#attach(trait_B)

#Also create tables with just BE or BD+BT for later use
trait_BE<-subset(traits,group=="BE",drop = T)
trait_BDT<-subset(traits,group=="BD" | group=="BT",drop = T)

#--- Bivariate plots and statistics with SMA regression ---
# These are used to fill out the hypothesis framework plot with R
 
# Calculate for all broadleaf 
bivar <- make_bivar_plots(trait_B,nbtstrp)
View(bivar$all_sma_bivar)

# Calculate for all evergreen broadleaf
bivar_BE <- make_bivar_plots(trait_BE,nbtstrp)
View(bivar_BE$all_sma_bivar)

# Calculate for all deciduous broadleaf
bivar_BDT <- make_bivar_plots(trait_BDT,nbtstrp)
View(bivar_BDT$all_sma_bivar)


#--- Experiment with different plausible multivariate SMA models, based on our theory ---


# P50 fits ----------------------------------------------------------------

P50_multivar <- P50_multivar_test(trait_B)

coeffnames_P50_from_TLP_Ks <- c("Coefficient","L95","U95")
intercept_P50_from_TLP_Ks <- c(P50_multivar$P50_from_TLP_Ks$mod$intercept_R,P50_multivar$P50_from_TLP_Ks$mod$L95_R.intercept,P50_multivar$P50_from_TLP_Ks$mod$U95_R.intercept)
y1_P50_from_TLP_Ks <- c(P50_multivar$P50_from_TLP_Ks$mod$slope_R.y1,P50_multivar$P50_from_TLP_Ks$mod$L95_R.y1,P50_multivar$P50_from_TLP_Ks$mod$U95_R.y1)
y2_P50_from_TLP_Ks <- c(P50_multivar$P50_from_TLP_Ks$mod$slope_R.y2,P50_multivar$P50_from_TLP_Ks$mod$L95_R.y2,P50_multivar$P50_from_TLP_Ks$mod$U95_R.y2)

coeff_P50_from_TLP_Ks <- data.frame(coeffnames_P50_from_TLP_Ks,intercept_P50_from_TLP_Ks,y1_P50_from_TLP_Ks,y2_P50_from_TLP_Ks)
View(coeff_P50_from_TLP_Ks)

# NOTE: These coefficients are all over the place, almost certainly because we have high multicolinearity in the predictors, BUT, this is not a problem as we are not interpreting the coefficients, just using them for the prediction.


# TLP fits ----------------------------------------------------------------

TLP_multivar <- TLP_multivar_test(trait_B)


# LMA fits -----------------------------------------------------------------
# Separate for BE and BDT (BD + BT) on the basis that LMA has a very different range and set of bivariate relationships for these
# two different groups, unlike the other traits here.

LMA_multivar_BE <- LMA_multivar_test_BE(trait_BE)

LMA_multivar_BDT <- LMA_multivar_test_BDT(trait_BDT)


# WD fits -----------------------------------------------------------------

WD_multivar <- WD_multivar_test(trait_B)


# slope fits --------------------------------------------------------------

slope_multivar <- slope_multivar_test(trait_B)


# Ks fits ----------------------------------------------------------------

Ks_multivar <- Ks_multivar_test(trait_B)

# LS fits -----------------------------------------------------------------

LS_multivar <- LS_multivar_test(trait_B)


# Make plots showing quality of fits and climate coverage -----------------

MATp1 <- trait_B$MAT[P50_multivar$P50_from_TLP_Ks$dataused]
MAPp1 <- trait_B$MAP[P50_multivar$P50_from_TLP_Ks$dataused]/10
name1 <- rep("P50",length(MATp1))

MATp2 <- trait_B$MAT[TLP_multivar$TLP_from_LS_LMA_P50$dataused]
MAPp2 <- trait_B$MAP[TLP_multivar$TLP_from_LS_LMA_P50$dataused]/10
name2 <- rep("TLP",length(MATp2))

MATp3 <- trait_B$MAT[WD_multivar$WD_from_slope_P50slope$dataused]
MAPp3 <- trait_B$MAP[WD_multivar$WD_from_slope_P50slope$dataused]/10
name3 <- rep("WD",length(MATp3))

MATp4 <- trait_B$MAT[slope_multivar$slope_from_P50_TLP_Ks$dataused]
MAPp4 <- trait_B$MAP[slope_multivar$slope_from_P50_TLP_Ks$dataused]/10
name4 <- rep("Slope",length(MATp4))

MATp5 <- trait_BE$MAT[LMA_multivar_BE$LMA_from_TLP$dataused]
MAPp5 <- trait_BE$MAP[LMA_multivar_BE$LMA_from_TLP$dataused]/10
name5 <- rep("LMA (BE)",length(MATp5))

MATp6 <- trait_BDT$MAT[LMA_multivar_BDT$LMA_from_TLP$dataused]
MAPp6 <- trait_BDT$MAP[LMA_multivar_BDT$LMA_from_TLP$dataused]/10
name6 <- rep("LMA (BD)",length(MATp6))

data_MATp_MAPp <- data.frame("MATp"=c(MATp1,MATp2,MATp3,MATp4,MATp5,MATp6),
                             "MAPp"=c(MAPp1,MAPp2,MAPp3,MAPp4,MAPp5,MAPp6),
                             "name"=c(name1,name2,name3,name4,name5,name6))

whittaker_biomes_plot(data_MATp_MAPp)

# Optimisation with TLP and LS------------------------------------------------------------
# lowest bivariate relationship 0.36 when all broadleaf data is considered
# Attempt to iteratively converge on the best fit values of Ks, P50 and LMA given known TLP and LS

# Decide whether to limit the possible ranges of predicted traits to the observed values (T) or not (F)
limitdataranges=T # Currently does not converge in uncertainty propagation if not set to T

# Decide whether to run the uncertainty propagation (T) or not (F)
propagate_uncer=T

# Decide whether to run all trait combinations in the database for LS and Ks (F), or just a selection (T), T useful for generating output for LPJ-Guess
trait_sel= T 
# Number of combinations to select if trait_sel=T. Set to -1 for a systematic sample, >0 for a random sample of the size specified, we have created 28 PFTs.
n_trait_sel=28#-1

# Run for all deciduous (BT + BD) (=1), or BE (=2), or BT (=3), or BD (=4). This is used to set the maximum and minimum bounds in trait_opt().
spec_group_sel=3

outs_LSTLP <- trait_optim_bivar_startLSTLP(limitdataranges = T,propagate_uncer = T,trait_sel = F, n_trait_sel = 28, spec_group_sel = 4)
# not very elegant.. but: this is to 'release' the output from function trait_optim_bivar_startLSTLP from a list of objects into single objects
# single objects
list2env(outs_LSTLP$predictors , envir = .GlobalEnv) 
list2env(outs_LSTLP$predicted , envir = .GlobalEnv)

#Stats defining the uncertainty range for each point
create_uncertainty_range_stats(outs_LSTLP)

opt_test_plots_LSTLP(trait_plot,
                     Ks_e_mean,
                     Ks_e_5perc,
                     Ks_e_95perc,
                     Ks_e,
                     P50_e_mean,
                     P50_e_5perc,
                     P50_e_95perc,
                     P50_e,
                     LMA_e_mean,
                     LMA_e_5perc,
                     LMA_e_95perc,
                     LMA_e,
                     WD_e_mean,
                     WD_e_5perc,
                     WD_e_95perc,
                     WD_e,
                     slope_e_mean,
                     slope_e_5perc,
                     slope_e_95perc,
                     slope_e)

# Calculate the RMSE (only if running with actual values of Ks and LS [i.e. trait_sel =F in function trait_optim])
trait_sel = F

if(trait_sel==F) {
  # Identify all combinations of Ks and LS (do this across full range of broadleaf species)
  ind = which(!is.na(traits$TLP) & !is.na(traits$LS))
  
  # provide list of trait names which are the predicted traits
  trait_names = c('P50','Ks','LMA','WD','slope')
  RMSE_withTLPLS_start <- opt_rmse(traits,trait_names,ind)
}



# Optimisation with KS and TLP ------------------------------------------------------------
# 
# Attempt to iteratively converge on the best fit values of LS, P50 and LMA, given known KS and TLP

# Decide whether to limit the possible ranges of predicted traits to the observed values (T) or not (F)
limitdataranges=T # Currently does not converge in uncertainty propagation if not set to T

# Decide whether to run the uncertainty propagation (T) or not (F)
propagate_uncer=T

# Decide whether to run all trait combinations in the database for LS and Ks (F), or just a selection (T), T useful for generating output for LPJ-Guess
trait_sel= T 
# Number of combinations to select if trait_sel=T. Set to -1 for a systematic sample, >0 for a random sample of the size specified, we have created 28 PFTs.
n_trait_sel = 28#-1

# Run for all deciduous (BT + BD) (=1), or BE (=2), or BT (=3), or BD (=4). This is used to set the maximum and minimum bounds in trait_opt().
spec_group_sel=3

outs_KsTLP <- trait_optim_bivar_start_KsTLP(limitdataranges = T,propagate_uncer = T,trait_sel = F, n_trait_sel = 28, spec_group_sel = 4)
# not very elegant.. but: this is to 'release' the output from function trait_optim_bivar_startLSTLP from a list of objects into single objects
# single objects
list2env(outs_KsTLP$predictors , envir = .GlobalEnv) 
list2env(outs_KsTLP$predicted , envir = .GlobalEnv)

# Stats defining the uncertainty range for each point
create_uncertainty_range_stats(outs_KsTLP)

opt_test_plots_KsTLP(trait_plot,
                                 LS_e_mean,
                                 LS_e_5perc,
                                 LS_e_95perc,
                                 LS_e,
                                 P50_e_mean,
                                 P50_e_5perc,
                                 P50_e_95perc,
                                 P50_e,
                                 LMA_e_mean,
                                 LMA_e_5perc,
                                 LMA_e_95perc,
                                 LMA_e,
                                 WD_e_mean,
                                 WD_e_5perc,
                                 WD_e_95perc,
                                 WD_e,
                                 slope_e_mean,
                                 slope_e_5perc,
                                 slope_e_95perc,
                                 slope_e) 


# Calculate the RMSE (only if running with actual values of Ks and LS [i.e. trait_sel =F in function trait_optim])

if(trait_sel==F) {
  # Identify all combinations of Ks and LS (do this across full range of broadleaf species)
  ind = which(!is.na(traits$Ks) & !is.na(traits$TLP))
  
  # provide list of trait names which are the predicted traits
  trait_names = c('P50','LS','LMA','WD','slope')
  RMSE_withKsTLP_start <- opt_rmse(traits,trait_names,ind)
}



# Optimisation with Ks and LS------------------------------------------------------------
# Attempt to iteratively converge on the best fit values of TLP, P50 and LMA given known Ks and LS

outs_LSKs <- trait_optim(limitdataranges = T,propagate_uncer = T,trait_sel = F, n_trait_sel = 28, spec_group_sel = 4)
list2env(outs_LSKs$predictors , envir = .GlobalEnv)
list2env(outs_LSKs$predicted , envir = .GlobalEnv)

# Stats defining the uncertainty range for each point
create_uncertainty_range_stats(outs_LSKs)

# Make plots to compare with original data

if (spec_group_sel==1 | spec_group_sel==3 | spec_group_sel==4) {
  trait_plot=trait_BDT
} else if (spec_group_sel==2) {
  trait_plot=trait_BE
}

opt_test_plots(trait_plot,
               TLP_e_mean,
               TLP_e_5perc,
               TLP_e_95perc,
               TLP_e,
               P50_e_mean,
               P50_e_5perc,
               P50_e_95perc,
               P50_e,
               LMA_e_mean,
               LMA_e_5perc,
               LMA_e_95perc,
               LMA_e,
               WD_e_mean,
               WD_e_5perc,
               WD_e_95perc,
               WD_e,
               slope_e_mean,
               slope_e_5perc,
               slope_e_95perc,
               slope_e)


# Calculate the RMSE (only if running with actual values of Ks and LS [i.e. trait_sel =F in function trait_optim])

if(trait_sel==F) {
  # Identify all combinations of Ks and LS (do this across full range of broadleaf species)
  ind = which(!is.na(traits$Ks) & !is.na(traits$LS))
  
  # provide list of trait names which are the predicted traits
  trait_names = c('P50','TLP','LMA','WD','slope')
  RMSE_withKsLS_start <- opt_rmse(traits,trait_names,ind)
}


# Write the optimised trait values out to a file
traits_e_out <- data.frame(LS_e,Ks_e,
                           P50_e_mean,P50_e_5perc,P50_e_95perc,
                           TLP_e_mean,TLP_e_5perc,TLP_e_95perc,
                           LMA_e_mean,LMA_e_5perc,LMA_e_95perc,
                           WD_e_mean,WD_e_5perc,WD_e_95perc,
                           slope_e_mean,slope_e_5perc,slope_e_95perc)
write.table(format(traits_e_out, digits=3), "traits_e_out_systtraits_260820.csv", append = FALSE, sep = ",", dec = ".",row.names = F, col.names = T)





# Calculate regression of leafL from LMA ----------------------------------

leafL_from_LMA <- sma_plot_stats(data.frame(trait_B$LMA,log(trait_B$leafL)),c("LMA","leafL"),nbtstrp,T)


# Calculate limits of leafN vs LMA to allow estimate of leaf C:N ----------

leafN_from_LMA <- sma_plot_stats(data.frame(trait_B$LMA,trait_B$leafN),c("LMA","leafN"),nbtstrp,T)

leafN_from_LMA_limit <- regress_limit_adjust(trait_B$leafN,trait_B$LMA,leafN_from_LMA,0.05)

plot(trait_B$LMA,trait_B$leafN)
points(trait_B$LMA[leafN_from_LMA_limit$ind],leafN_from_LMA_limit$var1_pred_lower,col="green")
points(trait_B$LMA[leafN_from_LMA_limit$ind],leafN_from_LMA_limit$var1_pred_upper,col="red")


# Convert to the values needed in LPJ-GUESS and write out -----------------

traits_LPJG <- lpjg_traits_conv(LMA_e_mean,P50_e_mean,TLP_e_mean,slope_e_mean,
                                LS_e,WD_e_mean,Ks_e,
                                leafL_from_LMA,leafN_from_LMA,leafN_from_LMA_limit)

# Select which base PFT to use: TeBE (1), TeBS (2), IBS (3), TrBE (4) or TrBR (5)
basePFT=5

# Select output folder
#output_fol="/Users/pughtam/Documents/TreeMort/Analyses/Hydraulic_modelling/Traits/uncer_test_KsLS/revised_PFTs_141220"
output_fol="/Users/annemarie/Documents/1_TreeMort/2_Analysis/1_Inputs/"
# Set the name for the output file
if (basePFT==1) {
  LPJG_outfile <- paste(output_fol,"/LPJG_PFT_insfile_TeBE.ins",sep="")
  LPJG_summaryfile <- paste(output_fol,"/LPJG_PFT_summary_TeBE.csv",sep="")
} else if (basePFT==2) {
  LPJG_outfile <- paste(output_fol,"/LPJG_PFT_insfile_TeBS.ins",sep="")
  LPJG_summaryfile <- paste(output_fol,"/LPJG_PFT_summary_TeBS.csv",sep="")
} else if (basePFT==3) {
  LPJG_outfile <- paste(output_fol,"LPJG_PFT_insfile_IBS.ins",sep="")
  LPJG_summaryfile <- paste(output_fol,"/LPJG_PFT_summary_IBS.csv",sep="")
} else if (basePFT==4) {
  LPJG_outfile <- paste(output_fol,"/LPJG_PFT_insfile_TrBE.ins",sep="")
  LPJG_summaryfile <- paste(output_fol,"/LPJG_PFT_summary_TrBE.csv",sep="")
} else if (basePFT==5) {
  LPJG_outfile <- paste(output_fol,"/LPJG_PFT_insfile_TrBR.ins",sep="")
  LPJG_summaryfile <- paste(output_fol,"/LPJG_PFT_summary_TrBR.csv",sep="")
} else {
  stop("basePFT must be equal to 1, 2, 3, 4 or 5")
}

# Write out to LPJ-GUESS instruction file
PFTfile <- file(LPJG_outfile)
for (nn in 1:length(traits_LPJG$Ks)) {
  if (nn>1) {
    PFTfile <- file(LPJG_outfile,open="append")
  }
  
  Line1 <- paste("pft \"PFT",nn,"\" (",sep="")
  if (basePFT==1) {
    Line2 <- TeBE_header
  } else if (basePFT==2) {
    Line2 <- TeBS_header
  } else if (basePFT==3) {
    Line2 <- IBS_header
  } else if (basePFT==4) {
    Line2 <- TrBE_header
  } else if (basePFT==5) {
    Line2 <- TrBR_header
  }
  Line3 <- "\t !Hydraulics"
  Line4 <- paste("\t isohydricity ",traits_LPJG$lambda[nn],sep="")
  Line5 <- paste("\t delta_psi_max ",traits_LPJG$DeltaPsiWW[nn],sep="")
  Line6 <- paste("\t cav_slope ",traits_LPJG$polyslope[nn],sep="")
  Line7 <- paste("\t psi50_xylem ",traits_LPJG$P50[nn],sep="")
  Line8 <- paste("\t ks_max ",traits_LPJG$Ks[nn],sep="")
  Line9 <- paste("\t kr_max ",11.2e-4,sep="") # LPJ-GUESS default from Hickler et al. (2006)
  Line10 <- paste("\t kL_max ",traits_LPJG$Kleaf[nn],sep="")
  Line11 <- paste("\t wooddens ",traits_LPJG$WD[nn],sep="")
  Line12 <- paste("\t k_latosa ",traits_LPJG$LS[nn],sep="")
  Line13 <- paste("\t sla ",traits_LPJG$SLA[nn],sep="")
  Line14 <- paste("\t cton_leaf_min ",traits_LPJG$CtoNmin_LPJG[nn],sep="")
  if (basePFT==1 | basePFT==4) {
    Line15 <- paste("\t leaflong ",traits_LPJG$leaflong[nn],sep="")
    Line16 <- paste("\t turnover_leaf ",1/traits_LPJG$leaflong[nn],sep="")
  } else {
    #Use LPJ-GUESS standard values for deciduous
    Line15 <- paste("\t leaflong 0.5")
    Line16 <- paste("\t turnover_leaf 1.0")
  } 
  
  writeLines(c(Line1,Line2,Line3,Line4,Line5,Line6,Line7,Line8,Line9,Line10,Line11,Line12,Line13,Line14,Line15,Line16,"",")",""),PFTfile)
  close(PFTfile)
}

# Write out to a set of LPJ-GUESS instruction files, 1 per PFT
for (nn in 1:length(traits_LPJG$Ks)) {
  LPJG_outfile_pft <- paste(LPJG_outfile,".PFT",nn,sep="")
  file.copy("global_cf_base.ins",LPJG_outfile_pft,overwrite=T)
  PFTfile <- file(LPJG_outfile_pft,open="append")
  
  Line1 <- paste("pft \"PFT",nn,"\" (",sep="")
  if (basePFT==1) {
    Line2 <- TeBE_header
  } else if (basePFT==2) {
    Line2 <- TeBS_header
  } else if (basePFT==3) {
    Line2 <- IBS_header
  } else if (basePFT==4) {
    Line2 <- TrBE_header
  } else if (basePFT==5) {
    Line2 <- TrBR_header
  }
  Line3 <- "\t !Hydraulics"
  Line4 <- paste("\t isohydricity ",traits_LPJG$lambda[nn],sep="")
  Line5 <- paste("\t delta_psi_max ",traits_LPJG$DeltaPsiWW[nn],sep="")
  Line6 <- paste("\t cav_slope ",traits_LPJG$polyslope[nn],sep="")
  Line7 <- paste("\t psi50_xylem ",traits_LPJG$P50[nn],sep="")
  Line8 <- paste("\t ks_max ",traits_LPJG$Ks[nn],sep="")
  Line9 <- paste("\t kr_max ",11.2e-4,sep="") # LPJ-GUESS default from Hickler et al. (2006)
  Line10 <- paste("\t kL_max ",traits_LPJG$Kleaf[nn],sep="")
  Line11 <- paste("\t wooddens ",traits_LPJG$WD[nn],sep="")
  Line12 <- paste("\t k_latosa ",traits_LPJG$LS[nn],sep="")
  Line13 <- paste("\t sla ",traits_LPJG$SLA[nn],sep="")
  Line14 <- paste("\t cton_leaf_min ",traits_LPJG$CtoNmin_LPJG[nn],sep="")
  if (basePFT==1 | basePFT==4) {
    Line15 <- paste("\t leaflong ",traits_LPJG$leaflong[nn],sep="")
    Line16 <- paste("\t turnover_leaf ",1/traits_LPJG$leaflong[nn],sep="")
  } else {
    #Use LPJ-GUESS standard values for deciduous
    Line15 <- paste("\t leaflong 0.5")
    Line16 <- paste("\t turnover_leaf 1.0")
  } 
  
  writeLines(c(Line1,Line2,Line3,Line4,Line5,Line6,Line7,Line8,Line9,Line10,Line11,Line12,Line13,Line14,Line15,Line16,"",")",""),PFTfile)
  close(PFTfile)
}

# Write out a trait values table by PFT to be used for post-processing of LPJ-GUESS output
write.table(traits_LPJG,file=LPJG_summaryfile,sep=",",row.names = F)


  