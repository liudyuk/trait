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
source('trait_optim_bivar_start_KsTLP')
source('trait_opt_bivar_start_KsTLP.R')
source('opt_test_plots_KsTLP.R')
source('trait_optim_bivar_start_LSP50.R')
source('trait_opt_bivar_start_LSP50.R')
source('opt_test_plots_LSP50.R')

#---------------packages--------------------------
# PCA plotting
install_github("vqv/ggbiplot")
library('ggbiplot')
# systematic sampling for optimisation
library('hypervolume')


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
# lowest bivariate functional relationship within the network (R= 0.36) when all broadleaf data is considered
# Attempt to iteratively converge on the best fit values of Ks, P50 and LMA given known TLP and LS

# Decide whether to limit the possible ranges of predicted traits to the observed values (T) or not (F)
limitdataranges=T # Currently does not converge in uncertainty propagation if not set to T

# Decide whether to run the uncertainty propagation (T) or not (F)
propagate_uncer=T

# Decide whether to run all trait combinations in the database for LS and Ks (F), or just a selection (T), T useful for generating output for LPJ-Guess
trait_sel= T 
# Number of combinations to select if trait_sel=T. Set to -1 for a systematic sample, >0 for a random sample of the size specified, we have created 28 PFTs.
n_trait_sel=-1#-1

# Run for all deciduous (BT + BD) (=1), or BE (=2), or BT (=3), or BD (=4). This is used to set the maximum and minimum bounds in trait_opt().
spec_group_sel=3

outs_LSTLP <- trait_optim_bivar_startLSTLP(limitdataranges = limitdataranges ,propagate_uncer = propagate_uncer,trait_sel = trait_sel, n_trait_sel = n_trait_sel, spec_group_sel = spec_group_sel)
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


if(trait_sel==F) {
  # Identify all combinations of Ks and LS (do this across full range of broadleaf species)
  ind = which(!is.na(traits$TLP) & !is.na(traits$LS))
  
  # provide list of trait names which are the predicted traits
  trait_names = c('P50','Ks','LMA','WD','slope')
  RMSE_withTLPLS_start <- opt_rmse(traits,trait_names,ind)
}

# Convert to the values needed in LPJ-GUESS,PCA and write out -----------------
#KSTLP
traits_LPJG_TLPLS <- lpjg_traits_conv(LMA_e_mean,P50_e_mean,TLP_e,slope_e_mean,
                                LS_e,WD_e_mean,Ks_e_mean,
                                leafL_from_LMA,leafN_from_LMA,leafN_from_LMA_limit)

# create data frame of traits for subsequent PCA below -------
traits_TLPLS.df <- data.frame(matrix(unlist(traits_LPJG_TLPLS), ncol=length(traits_LPJG_TLPLS), byrow=FALSE))
names(traits_TLPLS.df) <- names(traits_LPJG_TLPLS)


# Optimisation with KS and TLP ------------------------------------------------------------
# no hypothesised functional link, but bivariate correlation coeff 0.3
# Attempt to iteratively converge on the best fit values of LS, P50 and LMA, given known KS and TLP


outs_KsTLP <- trait_optim_bivar_start_KsTLP(limitdataranges = limitdataranges, propagate_uncer = propagate_uncer,trait_sel = trait_sel, n_trait_sel = n_trait_sel, spec_group_sel = spec_group_sel)
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

# Convert to the values needed in LPJ-GUESS,PCA and write out -----------------
#KSTLP
traits_LPJG <- lpjg_traits_conv(LMA_e_mean,P50_e_mean,TLP_e_mean,slope_e_mean,
                                LS_e_mean,WD_e_mean,Ks_e,
                                leafL_from_LMA,leafN_from_LMA,leafN_from_LMA_limit)

# create data frame of traits for subsequent PCA below -------
traits_KSTLP.df <- data.frame(matrix(unlist(traits_LPJG), ncol=length(traits_LPJG), byrow=FALSE))
names(traits_KSTLP.df) <- names(traits_LPJG)



# Optimisation with Ks and P50 ------------------------------------------------------------
# lowest functional bivariate relationship in the network: 0.23, that is not directly used in the optimisation framweork (orange lines)
# Attempt to iteratively converge on the best fit values of Ks, TLP and LMA, given known LS and P50


outs_LSP50 <- trait_optim_bivar_start_LSP50(limitdataranges = limitdataranges ,propagate_uncer = propagate_uncer,trait_sel = trait_sel, n_trait_sel = n_trait_sel, spec_group_sel = spec_group_sel)
# not very elegant.. but: this is to 'release' the output from function trait_optim_bivar_startLSTLP from a list of objects into single objects
# single objects
list2env(outs_LSP50$predictors , envir = .GlobalEnv) 
list2env(outs_LSP50$predicted , envir = .GlobalEnv)

# Stats defining the uncertainty range for each point
create_uncertainty_range_stats(outs_LSP50)

opt_test_plots_LSP50(trait_plot,
                     Ks_e_mean,
                     Ks_e_5perc,
                     Ks_e_95perc,
                     Ks_e,
                     TLP_e_mean,
                     TLP_e_5perc,
                     TLP_e_95perc,
                     TLP_e,
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

if(trait_sel==F) {
  # Identify all combinations of Ks and LS (do this across full range of broadleaf species)
  ind = which(!is.na(traits$P50) & !is.na(traits$LS))
  
  # provide list of trait names which are the predicted traits
  trait_names = c('Ks','TLP','LMA','WD','slope')
  RMSE_withLSP50_start <- opt_rmse(traits,trait_names,ind)
}

# Convert to the values needed in LPJ-GUESS,PCA and write out -----------------
#LSP50
traits_LPJG <- lpjg_traits_conv(LMA_e_mean,P50_e,TLP_e_mean,slope_e_mean,
                                LS_e,WD_e_mean,Ks_e_mean,
                                leafL_from_LMA,leafN_from_LMA,leafN_from_LMA_limit)

# create data frame of traits for subsequent PCA below -------
traits_KSP50.df <- data.frame(matrix(unlist(traits_LPJG), ncol=length(traits_LPJG), byrow=FALSE))
names(traits_KSP50.df) <- names(traits_LPJG)

# ---------------------------------------------------------------------------------------
# Optimisation with Ks and LS------------------------------------------------------------
# No functional relationship between traits hypothesised (note: high observed correlation, R=0.4, probably via other traits), 
# so: traits are good 'outer edges' to start the optimisation from.
# Attempt to iteratively converge on the best fit values of TLP, P50 and LMA given known Ks and LS

outs_LSKs <- trait_optim(limitdataranges = limitdataranges ,propagate_uncer = propagate_uncer,trait_sel = trait_sel, n_trait_sel = n_trait_sel, spec_group_sel = spec_group_sel)
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
write.table(format(traits_e_out, digits=3), "traits_e_out_systtraits_KsLS.csv", append = FALSE, sep = ",", dec = ".",row.names = F, col.names = T)




# Calculate regression of leafL from LMA ----------------------------------

leafL_from_LMA <- sma_plot_stats(data.frame(trait_B$LMA,log(trait_B$leafL)),c("LMA","leafL"),nbtstrp,T)


# Calculate limits of leafN vs LMA to allow estimate of leaf C:N ----------

leafN_from_LMA <- sma_plot_stats(data.frame(trait_B$LMA,trait_B$leafN),c("LMA","leafN"),nbtstrp,T)

leafN_from_LMA_limit <- regress_limit_adjust(trait_B$leafN,trait_B$LMA,leafN_from_LMA,0.05)

plot(trait_B$LMA,trait_B$leafN)
points(trait_B$LMA[leafN_from_LMA_limit$ind],leafN_from_LMA_limit$var1_pred_lower,col="green")
points(trait_B$LMA[leafN_from_LMA_limit$ind],leafN_from_LMA_limit$var1_pred_upper,col="red")


# Convert to the values needed in LPJ-GUESS,PCA and write out -----------------

traits_LPJG_KSLS <- lpjg_traits_conv(LMA_e_mean,P50_e_mean,TLP_e_mean,slope_e_mean,
                                LS_e,WD_e_mean,Ks_e,
                                leafL_from_LMA,leafN_from_LMA,leafN_from_LMA_limit)

# create data frame of traits for subsequent PCA below -------
traits_LPJG_KSLS.df <- data.frame(matrix(unlist(traits_LPJG_KSLS), ncol=length(traits_LPJG_KSLS), byrow=FALSE))
names(traits_LPJG_KSLS.df) <- names(traits_LPJG_KSLS)


#------------------------------------------------------------------------------
# write .ins file for LPJGuess based on the above optimised trait combinations for 28 PFT 'variants' per PFT.

# Select which base PFT to use: TeBE (1), TeBS (2), IBS (3), TrBE (4) or TrBR (5)
basePFT=5

# Select output folder
#output_fol="/Users/pughtam/Documents/TreeMort/Analyses/Hydraulic_modelling/Traits/uncer_test_KsLS/revised_PFTs_141220"
output_fol="/Users/annemarie/Documents/1_TreeMort/2_Analysis/1_Inputs"

# create .insfiles for  LPJ-GUESS_hydro
#write_LPJG_ins.file(output_fol,basePFT = 1 ,traits_LPJG)
#write_LPJG_ins.file(output_fol,basePFT = 2 ,traits_LPJG)
#write_LPJG_ins.file(output_fol,basePFT = 3 ,traits_LPJG)
write_LPJG_ins.file(output_fol,basePFT = 4 ,traits_LPJG)
write_LPJG_ins.file(output_fol,basePFT = 5 ,traits_LPJG)


# perform PCA ------------------------------------------------------------
# as test to see whether starting from different trait combinations
# has an impact on the PFT-spread.
opt_traits <-  c('traits_LPJG_KSLS.df', 'traits_TLPLS.df', 'traits_KSTLP.df', 'traits_KSP50.df')
traits_BE  <- get(opt_traits[3])

# transform some values:
# T. Pugh
# 25.10.20
# original file: lpjg_strat_mapping_comb.m
# translated into R Annemarie Eckes-Shephard May 2021

## Convert SLA to LMA
#traits_BS.LMA=1./traits_BS.SLA;
traits_BE$LMA = 1./traits_BE$SLA

## Log traits that are non-normal
#traits_BS$P50 = log(-traits_BS$P50); 
traits_BE$P50=log(-traits_BE$P50)
#traits_BS$P88=log(-traits_BS$P88); 
traits_BE$P88=log(-traits_BE$P88)
#traits_BS$TLP=log(-traits_BS$TLP); 
traits_BE$TLP=log(-traits_BE$TLP)
#traits_BS$LS=log(traits_BS$LS); 
traits_BE$LS=log(traits_BE$LS)
#traits_BS$Ks=log(traits_BS$Ks); 
traits_BE$Ks=log(traits_BE$Ks)
#traits_BS$LMA=log(traits_BS$LMA); 
traits_BE$LMA=log(traits_BE$LMA)


#PCA on optimised trait values
opt.pca <- prcomp(traits_BE[,c(1,3,4,6,10,12,17)], center = TRUE,scale. = TRUE)

p_KsLS  <- ggbiplot(opt.pca,labels=row.names(traits_e_out))
p_TLPLS <- ggbiplot(opt.pca,labels=row.names(traits_e_out))
p_KSTLP <- ggbiplot(opt.pca,labels=row.names(traits_e_out))
p_KSP50 <- ggbiplot(opt.pca,labels=row.names(traits_e_out))

grid.arrange(p_KsLS , p_TLPLS,p_KSTLP, nrow = 1)


  