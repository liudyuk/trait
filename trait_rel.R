# Script to read in the processed trait data, make bivariate and multivariate SMA regressions and then
# carry out an optimisation procedure to unify inter-trait relationships across the whole trait dataset.
#
# [TO DO] optimisation starting with LSTLP still does not correspond to the other starting predictors' final PCA space 
# [TO DO] shorten axis labels in PCA plots
# [TO DO] Improve .ins writeout section

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
# - lhs, and QRM packages (if taking a latin hypercube sample)
#
# T. Pugh
# 15.06.20
#
# Annemarie Eckes-Shephard
# May + June 2021
# *Minor changes to make output reproducible (set.seed in hypervolume sampling)
# *outsource some code bits into separate functions
# *modify existing functions to test optimisation efficacy when starting with 
#  other variable combinations (not Ks LS)
# July 2021
# *PCA (previously a matlab script from Tom) now part of this script and applied
#  to all starting variable combinations ('predictors')
# *additional sampling method (latin hypercube sampling) implemented to exclude
#  a possible predictor range impact on final predictions 
# *reliably remove predicted trait combinations where one variable fell out 
#  the min max bounds of observations (marked with NA in trait_opt functions)
# *hard-wired (non-flexible) selection for 30 PFTs based on LS-P50 as starting
#  variables ('predictors') in the optimisation. Used for next step to check 
#  whether the emergent PFTs perform similarly to the PFTs derived through the
#  Ks-Ls trait optimisation
# October 2021
# Switch for LS predictions for BDT and BE data. 
# added linear regression model to determine trait network relations, compare against sma


nbtstrp=1000 # Number of bootstrap samples to take in sma_multivar_regress (samples later used to calculated uncertainty in the optimisation). Was previously 10 000, using a lower number for testing, Will need to check sensitivity to this value.

#traits=read.table("/Users/liudy/trait_data/woody_trait.0625.txt")
#traits=read.csv("/Users/pughtam/Documents/TreeMort/Analyses/Hydraulic_modelling/Traits/mytrait-data/woody_trait.0803.txt",sep="\t")
#traits=read.csv("/Users/pughtam/Documents/TreeMort/Analyses/Hydraulic_modelling/Traits/woody_trait.0827.txt",sep="\t")
traits = read.csv("/Users/annemarie/Documents/1_TreeMort/2_Analysis/2_data/2_intermediate/mytrait-data/woody_trait.0827.txt",sep="\t")

# related to trait network construction
source('sma_multivar_regress.R')
source('lm_regress_multivar.R')
source('trait_functions.R')
source('make_bivar_plots.R')
source('multivar_model_selection.R')
source("whittaker_biomes_plot.R")


# optimisation starting with KsLS and subsequent diagnostic plots 
source("trait_opt.R")
source('trait_optim.R')
source("opt_test_plots.R")

# for testing multiple starting variables for the optimisation
# optimisation starting with LSTLP and subsequent diagnostic plots  
source('trait_optim_bivar_startLSTLP.R') 
source('trait_opt_bivar_start_LSTLP.R')
source('opt_test_plots_LSTLP.R')

# optimisation starting with KsTLP and subsequent diagnostic plots 
source('trait_optim_bivar_start_KsTLP.R')
source('trait_opt_bivar_start_KsTLP.R')
source('opt_test_plots_KsTLP.R')

# optimisation starting with LSP50 and subsequent diagnostic plots 
source('trait_optim_bivar_start_LSP50.R')
source('trait_opt_bivar_start_LSP50.R')
source('opt_test_plots_LSP50.R')

# to create optimisation statistics
source('opt_rmse.R')
source('create_uncertainty_range_stats.R')

#functions related to latin hypercube sampling
source('estimate_trait_distr.R')
source('check_lhs_sampling.R')

# to visualise samples for optimisation
source('plot_Hypervolume.R')

# to create PFT input files for LPJGuess
source('lpjg_traits_conv.R')
source('write_LPJG_ins.file.R')


#---------------packages--------------------------
# PCA plotting
#install_github("vqv/ggbiplot") # must have devtools installed to use install_github
library('ggbiplot')
# systematic sampling for optimisation
library('hypervolume')
# efficient sampling of multiple parameter combinations using hypercube sampling
library('lhs')
# used to create samples by transforming latin hypercube sampling space
library('QRM')
#for plotting of hypercube sampling results:
library(grid)
library(ggplot2)
library("ggpubr")
# to plot PCA ggplots in grid
library(gridExtra)

# convenience package for ols normality checks, required in function lm_regress_multivar
library(olsrr)

#--- Experiment with different plausible multivariate SMA models, based on our theory ---
# multivariate lm and sma regression ppossible. change option here:
regr_type = 'lm'

#--- Read in the trait data ---

# Subset only the broadleaf species
trait_B<-subset(traits,group!="CC",drop = T)
#trait_B<-droplevels(trait_B)
#str(trait_B)
#attach(trait_B)

#BD = drought deciduous
#BT = teperature deciduous

#Also create tables with just BE or BD+BT for later use
trait_BE<-subset(traits,group=="BE",drop = T)
trait_BDT<-subset(traits,group=="BD" | group=="BT",drop = T)

#drought deciduous broadleaves only
trait_BD<-subset(traits,group=="BD",drop = T)
trait_BT<-subset(traits,group=="BT",drop = T)

#--- Bivariate plots and statistics with SMA regression ---
# These are used to fill out the hypothesis framework plot with R
 
# Calculate for all broadleaf 
bivar <- make_bivar_plots(trait_B,nbtstrp,regr_type = regr_type)
#View(bivar$all_sma_bivar)

# Calculate for all evergreen broadleaf
bivar_BE <- make_bivar_plots(trait_BE,nbtstrp,regr_type = regr_type)
# AHES: negative bivatiate correlation in data, but sma shows positive correlation 
subs <- na.omit(trait_BE[,c("P50","Ks")])
cor(subs$P50,subs$Ks)
cor(subs$Ks,subs$P50)
cor(subs$LS,subs$Ks)
cor(subs$LS,subs$P50)
#View(bivar_BE$all_sma_bivar)

# Calculate for all deciduous broadleaf
bivar_BDT <- make_bivar_plots(trait_BDT,nbtstrp, regr_type = regr_type)
# View(bivar_BDT$all_sma_bivar)
# AHES: negative bivatiate correlation in data, but sma shows positive correlation 
subs <- na.omit(trait_BD[,c("P50","LS")])
cor(subs$P50,subs$LS)
cor.test(subs$LS,subs$P50)

#----------------------------------------------------------------------------------------------------------------------
# Test single and multivariate correlation models
#----------------------------------------------------------------------------------------------------------------------



# P50 fits ----------------------------------------------------------------

P50_multivar <- P50_multivar_test(trait_B, regr_type =  regr_type)

coeffnames_P50_from_TLP_Ks <- c("Coefficient","L95","U95")
intercept_P50_from_TLP_Ks <- c(P50_multivar$P50_from_TLP_Ks$mod$intercept_R,P50_multivar$P50_from_TLP_Ks$mod$L95_R.intercept,P50_multivar$P50_from_TLP_Ks$mod$U95_R.intercept)
y1_P50_from_TLP_Ks <- c(P50_multivar$P50_from_TLP_Ks$mod$slope_R.y1,P50_multivar$P50_from_TLP_Ks$mod$L95_R.y1,P50_multivar$P50_from_TLP_Ks$mod$U95_R.y1)
y2_P50_from_TLP_Ks <- c(P50_multivar$P50_from_TLP_Ks$mod$slope_R.y2,P50_multivar$P50_from_TLP_Ks$mod$L95_R.y2,P50_multivar$P50_from_TLP_Ks$mod$U95_R.y2)

coeff_P50_from_TLP_Ks <- data.frame(coeffnames_P50_from_TLP_Ks,intercept_P50_from_TLP_Ks,y1_P50_from_TLP_Ks,y2_P50_from_TLP_Ks)
#View(coeff_P50_from_TLP_Ks)

# NOTE: These coefficients are all over the place, almost certainly because we have high multicolinearity in the predictors, BUT, this is not a problem as we are not interpreting the coefficients, just using them for the prediction.


# TLP fits ----------------------------------------------------------------

TLP_multivar <- TLP_multivar_test(trait_B, regr_type =  regr_type)


# LMA fits -----------------------------------------------------------------
# Separate for BE and BDT (BD + BT) on the basis that LMA has a very different range and set of bivariate relationships for these
# two different groups, unlike the other traits here.

LMA_multivar_BE <- LMA_multivar_test_BE(trait_BE, regr_type =  regr_type)

LMA_multivar_BDT <- LMA_multivar_test_BDT(trait_BDT, regr_type =  regr_type)



# WD fits -----------------------------------------------------------------

WD_multivar <- WD_multivar_test(trait_B,regr_type = regr_type)


# slope fits --------------------------------------------------------------

slope_multivar <- slope_multivar_test(trait_B,regr_type = regr_type)
hist(trait_B$slope[slope_multivar$slope_from_P50_TLP_Ks$dataused]-slope_multivar$slope_from_P50_TLP_Ks$var_est,xlab ="slope-slope_est", main=paste(" distribution of residuals; model: ", Ks_multivar$Ks_from_P50_L$`regression_type (lm or sma)` ))
shapiro.test(trait_B$slope[slope_multivar$slope_from_P50_TLP_Ks$dataused]-slope_multivar$slope_from_P50_TLP_Ks$var_est)
mtext()

# Ks fits ----------------------------------------------------------------

Ks_multivar     <- Ks_multivar_test(trait_B,regr_type = regr_type)
hist(trait_B$Ks[Ks_multivar$Ks_from_P50_LS_LMA$dataused]-Ks_multivar$Ks_from_P50_LS_LMA$var_est,xlab ="Ks-Ks_est", main=paste(" distribution of residuals; model: ", Ks_multivar$Ks_from_P50_L$`regression_type (lm or sma)` ))

# LS fits -----------------------------------------------------------------
# Separate for BE and BDT (BD + BT) on the basis that LS has a very different range and set of bivariate relationships for these
# two different groups, unlike the other traits here.

LS_multivar_BE <- LS_multivar_test(trait_BE, leaf_type ='BE',regr_type = regr_type) # returns LS_from_LMA_TLP_Ks


LS_multivar_BDT <- LS_multivar_test(trait_BDT, leaf_type ='BDT',regr_type = regr_type) # returns LS_from_TLP_Ks



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

MATp7 <- trait_BDT$MAT[LS_multivar_BDT$LS_from_LMA_TLP_Ks$dataused]
MAPp7 <- trait_BDT$MAP[LS_multivar_BDT$LS_from_LMA_TLP_Ks$dataused]/10
name7 <- rep("LS (BDT)",length(MATp7))

MATp8 <- trait_BE$MAT[LS_multivar_BE$LS_from_TLP_Ks$dataused]
MAPp8 <- trait_BE$MAP[LS_multivar_BE$LS_from_TLP_Ks$dataused]/10
name8 <- rep("LS (BE)",length(MATp8))

MATp9 <- trait_B$MAT[Ks_multivar$Ks_from_P50_LS$dataused]
MAPp9 <- trait_B$MAP[Ks_multivar$Ks_from_P50_LS$dataused]/10
name9 <- rep("Ks",length(MATp9))


data_MATp_MAPp <- data.frame("MATp"=c(MATp1,MATp2,MATp3,MATp4,MATp5,MATp6,MATp7,MATp8,MATp9),
                             "MAPp"=c(MAPp1,MAPp2,MAPp3,MAPp4,MAPp5,MAPp6,MAPp7,MAPp8,MAPp9),
                             "name"=c(name1,name2,name3,name4,name5,name6,name7,name8,name9))

whittaker_biomes_plot(data_MATp_MAPp)


# Derive traits that are not subject optimisation--------------------------

# Calculate regression of leafL from LMA ----------------------------------

leafL_from_LMA <- sma_plot_stats(data.frame(trait_B$LMA,log(trait_B$leafL)),c("LMA","leafL"),nbtstrp,T)


# Calculate limits of leafN vs LMA to allow estimate of leaf C:N ----------

leafN_from_LMA <- sma_plot_stats(data.frame(trait_B$LMA,trait_B$leafN),c("LMA","leafN"),nbtstrp,T)

leafN_from_LMA_limit <- regress_limit_adjust(trait_B$leafN,trait_B$LMA,leafN_from_LMA,0.05)

plot(trait_B$LMA,trait_B$leafN)
points(trait_B$LMA[leafN_from_LMA_limit$ind],leafN_from_LMA_limit$var1_pred_lower,col="green")
points(trait_B$LMA[leafN_from_LMA_limit$ind],leafN_from_LMA_limit$var1_pred_upper,col="red")

#----------------------------------------------------------------------------------------------------------------------
#Optimisation preparations: Latin Hypercube sampling
#----------------------------------------------------------------------------------------------------------------------

### latin hypercube sampling -------------------------------------------
traits <- trait_B # subset for broadleaf traits only
#traits <- trait_BE
### latin hypercube sampling of predictor traits from which to start the optimisation from.
# Identify all combinations of Ks and LS, P50 and TLP (do this across full range of broadleaf species)
ind = which(!is.na(traits$Ks) & !is.na(traits$LS) & !is.na(traits$P50) & !is.na(traits$TLP))

#create latin (hyper)cube samples for the 4 traits.
#lc_sample <- create_oalhs(28,4,bChooseLargerDesign = TRUE,bverbose = FALSE)
set.seed(1400)
# Number of samples
n = 200
#lc_sample <- optimumLHS(n,k = 4, maxSweeps = 2, eps = 0.1, verbose = FALSE)
# the drawn lc_samples from the sampling space are normally distributed and between 0 and 1

# The trait observations are not strictly normally distributed (according to shapiro- wilk test) 
# and therefore may need transformation.
# At the very least lhc- sampled value ranges need to be adapted to reflect the trait values. 
# Below are three options that were tested:
# Option 1) transform sampled marginal distributions to uniform distributions:
# that way we get stratified sample values that span the whole min-max range of the trait distribution and
# therefore the potential hydraulic strategy space.
#LS_e  <- qunif(lc_sample[,1], min(traits$LS[ind]),  max(traits$LS[ind])) 
#Ks_e  <- qunif(lc_sample[,2], min(traits$Ks[ind]),  max(traits$Ks[ind])) 
#TLP_e <- qunif(lc_sample[,3], min(traits$TLP[ind]), max(traits$TLP[ind])) 
#P50_e <- qunif(lc_sample[,4], min(traits$P50[ind]), max(traits$P50[ind])) 

# Option 2) adjust sampled marginal distributions to the variables' normal distributions ranges:
# this way we get stratified sample values that span the more commonly observed trait distribution
LS_e  <- qnorm(lc_sample[,1], mean = mean(traits$LS[ind]), sd = sd(traits$LS[ind])) 
Ks_e  <- qnorm(lc_sample[,2], mean = mean(traits$Ks[ind]),  sd = sd(traits$Ks[ind])) 
TLP_e <- qnorm(lc_sample[,3], mean = mean(traits$TLP[ind]), sd = sd(traits$TLP[ind])) 
P50_e <- qnorm(lc_sample[,4], mean = mean(traits$P50[ind]), sd = sd(traits$P50[ind])) 

# Option 3) for student t distributions
# this way we get stratified sample values that span the more commonly observed trait distribution, 
# allowing for the outer edges to have been sampled more frequently
# first, estimate the parameters for each t-distribution
#estimate_trait_distr(traits)
# then use these parameters to transform the lh_sample normal distribution into t-distr.
#LS_e  <- LaplacesDemon::qst(lc_sample[,1], mu=fitLS$par.ests[2]-1, sigma=fitLS$par.ests[3], nu=fitLS$par.ests[1], lower.tail=TRUE, log.p=FALSE)
#Ks_e  <- LaplacesDemon::qst(lc_sample[,2], mu=fitKs$par.ests[2],  sigma=fitKs$par.ests[3], nu=fitKs$par.ests[1],lower.tail=TRUE, log.p=FALSE)
#TLP_e <- LaplacesDemon::qst(lc_sample[,3], mu=fitTLP$par.ests[2], sigma=fitTLP$par.ests[3], nu=fitTLP$par.ests[1],lower.tail=TRUE, log.p=FALSE)
#P50_e <- LaplacesDemon::qst(lc_sample[,4], mu=fitP50$par.ests[2], sigma=fitP50$par.ests[3], nu=fitP50$par.ests[1],lower.tail=TRUE, log.p=FALSE)

#create objects with transformed values from above, which will act as starting values (= predictors) for the optimisation below
est_lhsLSP50 <- list(LS_e,P50_e)
names(est_lhsLSP50) <- c('LS_e','P50_e')

est_lhsLSTLP <- list(LS_e,TLP_e)
names(est_lhsLSTLP) <- c('LS_e','TLP_e')

est_lhsKsLS <- list(LS_e,Ks_e)
names(est_lhsKsLS) <- c('LS_e','Ks_e')

est_lhsKsTLP <- list(Ks_e,TLP_e)
names(est_lhsKsTLP) <- c('Ks_e','TLP_e')

# Check the distribution of the sample in relation to the observations
LSKsplot  <- check_lhs_sampling(traits$LS[ind],traits$Ks[ind],LS_e,Ks_e,'LS','Ks')
LSP50plot <- check_lhs_sampling(traits$LS[ind],traits$P50[ind],LS_e,P50_e,'LS','P50')
LSTLPplot <- check_lhs_sampling(traits$LS[ind],traits$TLP[ind],LS_e,TLP_e,'LS','TLP')

#plot the above
ggarrange(LSKsplot , LSP50plot , LSTLPplot  + rremove("x.text"), 
          labels = c("A", "B", "C"),
          ncol = 2, nrow = 2)

#----------------------------------------------------------------------------------------------------------------------
# Optimisation
#----------------------------------------------------------------------------------------------------------------------

###
# Optimisation setup
# select options that will be used in all optimisations below

# Decide whether to limit the possible ranges of predicted traits to the observed values (T) or not (F)
limitdataranges=T # Currently does not converge in uncertainty propagation if not set to T

# Decide whether to run the uncertainty propagation (T) or not (F)
propagate_uncer=F

# Decide whether to run all trait combinations in the database for LS and Ks (F), or just a selection (T), T useful for generating output for LPJ-Guess
# and useful for testing different sampling methods  ( e.g. latin hypercube vs. systematic vs. hypervolume)
trait_sel= F

# Number of combinations to select if trait_sel=T. Set to -1 for a systematic sample, >0 for a random sample of the size specified, we have created 28 PFTs.
# or set = 4 for a predefined (above) hypercube sample.
n_trait_sel= 80

# Run for all deciduous (BT + BD) (=1), or BE (=2), or BT (=3), or BD (=4). This is used to set the maximum and minimum bounds in trait_opt().
spec_group_sel = 2

#Based on the above decision, determine trait dataset to use for plotting against optimised data
if (spec_group_sel==1 | spec_group_sel==3 | spec_group_sel==4) {
  trait_plot = trait_BDT
} else if (spec_group_sel==2) {
  trait_plot = trait_BE
}



#Optimisation with TLP and LS------------------------------------------------------------
# one of the lowest bivariate functional relationship within the network (R= 0.36) 
# when all broadleaf data is considered
# Attempt to iteratively converge on the best fit values of Ks, P50 and LMA given known TLP and LS

outs_LSTLP <- trait_optim_bivar_startLSTLP(limitdataranges = limitdataranges ,propagate_uncer = propagate_uncer,trait_sel = trait_sel, n_trait_sel = n_trait_sel, spec_group_sel = spec_group_sel,est_lhs= est_lhsLSTLP)

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
  # Identify all combinations of TLP and LS (do this across full range of broadleaf species)
  ind = which(!is.na(traits$TLP) & !is.na(traits$LS))
  
  # provide list of trait names which are the predicted traits
  trait_names = c('P50','Ks','LMA','WD','slope')
  RMSE_withTLPLS_start <- opt_rmse(traits,trait_names,ind)
}

# Convert to the values needed in LPJ-GUESS-----------------
#LSTLP
traits_LPJG_LSTLP <- lpjg_traits_conv(LMA_e_mean,P50_e_mean,as.vector(TLP_e),slope_e_mean,
                                      as.vector(LS_e),WD_e_mean,Ks_e_mean,
                                leafL_from_LMA,leafN_from_LMA,leafN_from_LMA_limit)

# create data frame of traits for subsequent PCA below -------
traits_LSTLP.df <- data.frame(matrix(unlist(traits_LPJG_LSTLP), ncol=length(traits_LPJG_LSTLP), byrow=FALSE))
names(traits_LSTLP.df) <- names(traits_LPJG_LSTLP)


# Optimisation with Ks and TLP ------------------------------------------------------------
# no hypothesised functional link, bivariate correlation coeff: 0.3
# Attempt to iteratively converge on the best fit values of LS, P50 and LMA, given known Ks and TLP


outs_KsTLP <- trait_optim_bivar_start_KsTLP(limitdataranges = limitdataranges, propagate_uncer = propagate_uncer,trait_sel = trait_sel, n_trait_sel = n_trait_sel, spec_group_sel = spec_group_sel,est_lhs= est_lhsKsTLP)
# 'release' the output from function trait_optim_bivar_startKsTLP from a list of objects into single objects
# single objects into the global environment
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

# Convert to the values needed in LPJ-GUESS -----------------
traits_LPJG_KSTLP <- lpjg_traits_conv(LMA_e_mean,P50_e_mean,as.vector(TLP_e),slope_e_mean,
                                LS_e_mean,WD_e_mean,as.vector(Ks_e),
                                leafL_from_LMA,leafN_from_LMA,leafN_from_LMA_limit)

# create data frame of traits for subsequent PCA below -------
traits_KSTLP.df <- data.frame(matrix(unlist(traits_LPJG_KSTLP), ncol=length(traits_LPJG_KSTLP), byrow=FALSE))
names(traits_KSTLP.df) <- names(traits_LPJG_KSTLP)



# Optimisation with LS and P50 ------------------------------------------------------------
# lowest bivariate relationship in the trait network for evergreen subset: 0.23.
# It is not directly used in the optimisation framework (orange lines)
# as it is thought not to have a functional relationship.
# Attempt to iteratively converge on the best fit values of Ks, TLP and LMA, given known LS and P50
outs_LSP50 <- trait_optim_bivar_start_LSP50(limitdataranges = limitdataranges ,propagate_uncer = propagate_uncer,trait_sel = trait_sel, n_trait_sel = n_trait_sel, spec_group_sel = spec_group_sel,est_lhs = est_lhsLSP50)

# to 'release' the output from function trait_optim_bivar_startLSTLP from a list of objects into single objects
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

# Convert to the values needed in LPJ-GUESS -----------------
traits_LPJG_LSP50 <- lpjg_traits_conv(LMA_e_mean,as.vector(P50_e),TLP_e_mean,slope_e_mean,
                                      as.vector(LS_e),WD_e_mean,Ks_e_mean,
                                leafL_from_LMA,leafN_from_LMA,leafN_from_LMA_limit)

# create data frame of traits for subsequent PCA below -------
traits_LSP50.df <- data.frame(matrix(unlist(traits_LPJG_LSP50), ncol=length(traits_LPJG_LSP50), byrow=FALSE))
names(traits_LSP50.df) <- names(traits_LPJG_LSP50)


# Optimisation with Ks and LS------------------------------------------------------------
# No functional relationship between traits hypothesised (note: high observed correlation, R=0.4, probably via other traits), 
# so: traits are good 'outer edges' to start the optimisation from.
# Attempt to iteratively converge on the best fit values of TLP, P50 and LMA given known Ks and LS
outs_LSKs <- trait_optim(limitdataranges = limitdataranges ,propagate_uncer = propagate_uncer,trait_sel = trait_sel, n_trait_sel = n_trait_sel, spec_group_sel = spec_group_sel,est_lhs = est_lhsKsLS)

list2env(outs_LSKs$predictors , envir = .GlobalEnv)
list2env(outs_LSKs$predicted , envir = .GlobalEnv)

# Stats defining the uncertainty range for each point
create_uncertainty_range_stats(outs_LSKs)


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
#write.table(format(traits_e_out, digits=3), "traits_e_out_systtraits_KsLS_sma.csv", append = FALSE, sep = ",", dec = ".",row.names = F, col.names = T)

# Convert to the values needed in LPJ-GUESS,PCA and write out -----------------

traits_LPJG_KSLS <- lpjg_traits_conv(LMA_e_mean,P50_e_mean,TLP_e_mean,slope_e_mean,
                                as.vector(LS_e),WD_e_mean,as.vector(Ks_e),
                                leafL_from_LMA,leafN_from_LMA,leafN_from_LMA_limit)

# create data frame of traits for subsequent PCA below -------
traits_KSLS.df <- data.frame(matrix(unlist(traits_LPJG_KSLS), ncol=length(traits_LPJG_KSLS), byrow=FALSE))
names(traits_KSLS.df) <- names(traits_LPJG_KSLS)

## sanity checking parameter results plausibility:
# P50 should always be lower than TLP ?
regr_type
traits_KSLS.df$TLP > traits_KSLS.df$P50
traits_KSTLP.df$TLP > traits_KSTLP.df$P50
traits_LSP50.df$TLP > traits_LSP50.df$P50
traits_LSTLP.df$TLP > traits_LSTLP.df$P50



### collect rmse 
if(spec_group_sel==2){
rmse_be <- t(data.frame(trait_pair=c('Ks-LS','Ks-TLP','LSP-50','TLP-LS'), 
                        RMSE_BE = c(sum(RMSE_withKsLS_start$all_e_rmse),
                                 sum(RMSE_withKsTLP_start$all_e_rmse),
                                 sum(RMSE_withLSP50_start$all_e_rmse),
                                 sum(RMSE_withTLPLS_start$all_e_rmse))
                        ))
}
if(spec_group_sel==1){
  rmse_DBT <- t(data.frame(trait_pair=c('Ks-LS','Ks-TLP','LSP-50','TLP-LS'), 
                          RMSE_BDT = c(sum(RMSE_withKsLS_start$all_e_rmse),
                                      sum(RMSE_withKsTLP_start$all_e_rmse),
                                      sum(RMSE_withLSP50_start$all_e_rmse),
                                      sum(RMSE_withTLPLS_start$all_e_rmse))
  ))
}




#----------------------------------------------------------------------------------------------------------------------
# PCA
#----------------------------------------------------------------------------------------------------------------------
#trait_BE
# perform PCA ------------------------------------------------------------
# used as test to see whether starting from different trait combinations has an impact on the PFT-spread.

# save objects from analysis above, useful for quick jump to this section, if only PCA is of interest
traits_after_opt        <- c(traits_KSLS.df, traits_LSTLP.df, traits_LSP50.df, traits_KSTLP.df)
names(traits_after_opt) <- c('traits_KSLS.df', 'traits_LSTLP.df', 'traits_LSP50.df', 'traits_KSTLP.df')
#save(traits_after_opt,file='traits_after_opt_BE_sma.RData')
#load('traits_after_opt_BE_sma.RData')
#list2env(traits_after_opt , envir = .GlobalEnv)

# perform PCA with observed data:
traits_obs_BE <- transform_obs_for_PCA(trait_BE)
traits_obs_BDT <-  transform_obs_for_PCA(trait_BDT)
traits_obs_B  <-  transform_obs_for_PCA(trait_B)


par(mfrow=c(2,2))
opt_traits <-  c('traits_KSLS.df', 'traits_LSTLP.df', 'traits_KSTLP.df', 'traits_LSP50.df', 'traits_obs_BE','traits_obs_BDT','traits_obs_B')
#opt_traits_BE_sma <- 
for(o in c(1,2,3,4)){
traits_PCA  <- get(opt_traits[o])
#traits_BE <- traits_KSLS.df
#traits_BE <- traits_LSP50.df
#traits_BE <- traits_KSTLP.df
#traits_BE <- traits_LSTLP.df

# T. Pugh
# 25.10.20
# original file: lpjg_strat_mapping_comb.m
# translated into R Annemarie Eckes-Shephard May 2021

## Convert SLA to LMA
#traits_BS.LMA=1./traits_BS.SLA;
traits_PCA$LMA = 1./traits_PCA$SLA

# transform some values:
## Log traits that are non-normal
#traits_BS$P50 = log(-traits_BS$P50); 
traits_PCA$P50 = log(-traits_PCA$P50)
#traits_BS$P88 = log(-traits_BS$P88); 
traits_PCA$P88 = log(-traits_PCA$P88)
#traits_BS$TLP = log(-traits_BS$TLP); 
traits_PCA$TLP = log(-traits_PCA$TLP)
#traits_BS$LS=log(traits_BS$LS); 
traits_PCA$LS  = log(traits_PCA$LS)
#traits_BS$Ks=log(traits_BS$Ks); 
traits_PCA$Ks  = log(traits_PCA$Ks)
#traits_BS$LMA=log(traits_BS$LMA); 
traits_PCA$LMA = log(traits_PCA$LMA)


#PCA on optimised trait values and observations (0 = 5 to 7)
opt.pca <- prcomp(traits_PCA[,c(1,3,4,6,10,12,17)], center = TRUE,scale. = TRUE)


if(o==1){
  
  s <- summary(opt.pca)
  plot(opt.pca$x[,1],opt.pca$x[,2],
       xlab=paste("PCA 1 (", round(s$importance[2]*100, 1), "%)", sep = ""),
       ylab=paste("PCA 2 (", round(s$importance[5]*100, 1), "%)", sep = ""), pch=16, col='grey', cex=0.5, las=1, asp=1,ylim =c(-6,6),xlim=c(-8,8),
       main= 'Ks-LS start')
  # Benjaminbell.co.uk
  # Add grid lines
  abline(v=0, lty=2, col="grey50")
  abline(h=0, lty=2, col="grey50")
  # Add labels
  text(opt.pca$x[,1], opt.pca$x[,2], labels=row.names(opt.pca$x), pos=c(1,3,4,2), font=2)
  
  #rotation
  l.x <- opt.pca$rotation[,1]*10
  l.y <- opt.pca$rotation[,2]*10
  # Draw arrows
  arrows(x0=0, x1=l.x, y0=0, y1=l.y, col="red", length=0.15, lwd=1.5)
  
  # Label position
  l.pos <- l.y # Create a vector of y axis coordinates
  lo <- which(l.y < 0) # Get the variables on the bottom half of the plot
  hi <- which(l.y > 0) # Get variables on the top half
  # Replace values in the vector
  l.pos <- replace(l.pos, lo, "1")
  l.pos <- replace(l.pos, hi, "3")
  # Variable labels
  text(l.x, l.y, labels=row.names(opt.pca$rotation), col="red", pos=l.pos)
p_KsLS  <- ggbiplot(opt.pca,labels=1:dim(traits_PCA)[1],labels.size = 2) +
  #labs(y= "y axis name", x = "x axis name") +
  ggtitle('KsLS')
}
if(o==2){
  
  s <- summary(opt.pca)
  plot(opt.pca$x[,1],opt.pca$x[,2],
       xlab=paste("PCA 1 (", round(s$importance[2]*100, 1), "%)", sep = ""),
       ylab=paste("PCA 2 (", round(s$importance[5]*100, 1), "%)", sep = ""), pch=16, col='grey', cex=0.5, las=1, asp=1,ylim =c(-6,6),xlim=c(-8,8),
       main= 'LS - TLP start')
  # Benjaminbell.co.uk
  # Add grid lines
  abline(v=0, lty=2, col="grey50")
  abline(h=0, lty=2, col="grey50")
  # Add labels
  text(opt.pca$x[,1], opt.pca$x[,2], labels=row.names(opt.pca$x), pos=c(1,3,4,2), font=2)
  
  #rotation
  l.x <- opt.pca$rotation[,1]*10
  l.y <- opt.pca$rotation[,2]*10
  # Draw arrows
  arrows(x0=0, x1=l.x, y0=0, y1=l.y, col="red", length=0.15, lwd=1.5)

  
  # Label position
  l.pos <- l.y # Create a vector of y axis coordinates
  lo <- which(l.y < 0) # Get the variables on the bottom half of the plot
  hi <- which(l.y > 0) # Get variables on the top half
  # Replace values in the vector
  l.pos <- replace(l.pos, lo, "1")
  l.pos <- replace(l.pos, hi, "3")
  # Variable labels
  text(l.x, l.y, labels=row.names(opt.pca$rotation), col="red", pos=l.pos)
p_LSTLP <- ggbiplot(opt.pca,labels=1:dim(traits_PCA)[1],labels.size = 2) +
  #labs(y= "y axis name", x = "x axis name") +
  ggtitle('LSTLP')
}
if(o==3){
  
  s <- summary(opt.pca)
  plot(opt.pca$x[,1],opt.pca$x[,2],
       xlab=paste("PCA 1 (", round(s$importance[2]*100, 1), "%)", sep = ""),
       ylab=paste("PCA 2 (", round(s$importance[5]*100, 1), "%)", sep = ""), pch=16, col='grey', cex=0.5, las=1, asp=1,ylim =c(-6,6),xlim=c(-8,8),
       main= 'Ks - TLP start')
  # Benjaminbell.co.uk
  # Add grid lines
  abline(v=0, lty=2, col="grey50")
  abline(h=0, lty=2, col="grey50")
  # Add labels
  text(opt.pca$x[,1], opt.pca$x[,2], labels=row.names(opt.pca$x), pos=c(1,3,4,2), font=2)
  
  #rotation
  l.x <- opt.pca$rotation[,1]*10
  l.y <- opt.pca$rotation[,2]*10
  # Draw arrows
  arrows(x0=0, x1=l.x, y0=0, y1=l.y, col="red", length=0.15, lwd=1.5)
  
  # Label position
  l.pos <- l.y # Create a vector of y axis coordinates
  lo <- which(l.y < 0) # Get the variables on the bottom half of the plot
  hi <- which(l.y > 0) # Get variables on the top half
  # Replace values in the vector
  l.pos <- replace(l.pos, lo, "1")
  l.pos <- replace(l.pos, hi, "3")
  # Variable labels
  text(l.x, l.y, labels=row.names(opt.pca$rotation), col="red", pos=l.pos)
p_KSTLP <- ggbiplot(opt.pca,labels=1:dim(traits_PCA)[1],labels.size = 2) +
  #labs(y= "y axis name", x = "x axis name") +
  ggtitle('KsTLP')
}
if(o==4){
  
  s <- summary(opt.pca)
  plot(opt.pca$x[,1],opt.pca$x[,2],
       xlab=paste("PCA 1 (", round(s$importance[2]*100, 1), "%)", sep = ""),
       ylab=paste("PCA 2 (", round(s$importance[5]*100, 1), "%)", sep = ""), pch=16, col='grey', cex=0.5, las=1, asp=1,ylim =c(-6,6),xlim=c(-8,8),
       main= 'LS - P50 start')
  # Benjaminbell.co.uk
  # Add grid lines
  abline(v=0, lty=2, col="grey50")
  abline(h=0, lty=2, col="grey50")
  # Add labels
  text(opt.pca$x[,1], opt.pca$x[,2], labels=row.names(opt.pca$x), pos=c(1,3,4,2), font=2)
  
  #rotation
  l.x <- opt.pca$rotation[,1]*10
  l.y <- opt.pca$rotation[,2]*10
  # Draw arrows
  arrows(x0=0, x1=l.x, y0=0, y1=l.y, col="red", length=0.15, lwd=1.5)
  
  # Label position
  l.pos <- l.y # Create a vector of y axis coordinates
  lo <- which(l.y < 0) # Get the variables on the bottom half of the plot
  hi <- which(l.y > 0) # Get variables on the top half
  # Replace values in the vector
  l.pos <- replace(l.pos, lo, "1")
  l.pos <- replace(l.pos, hi, "3")
  # Variable labels
  text(l.x, l.y, labels=row.names(opt.pca$rotation), col="red", pos=l.pos)
  
#p_LSP50 <- ggbiplot(opt.pca,labels=1:dim(traits_PCA)[1],labels.size = 2) +
  #labs(y= "y axis name", x = "x axis name") +
 # ggtitle('LSP50')
}
if(o==5){
  p_obs_BE <- ggbiplot(opt.pca,labels=1:dim(traits_PCA)[1],labels.size = 2) +
    #labs(y= "y axis name", x = "x axis name") +
    ggtitle('obs_BE')
}
if(o==6){
  p_obs_BDT <- ggbiplot(opt.pca,labels=1:dim(traits_PCA)[1],labels.size = 2) +
    #labs(y= "y axis name", x = "x axis name") +
    ggtitle('obs_BDT')
}
if(o==7){
  p_obs_B <- ggbiplot(opt.pca,labels=1:dim(traits_PCA)[1],labels.size = 2) +
    #labs(y= "y axis name", x = "x axis name") +
    ggtitle('obs_B')
}
}
if(spec_group_sel==1){
  mtext(outer=TRUE, 'Broadleaf deciduous', line = -1.5)
}
if(spec_group_sel ==2){
  mtext(outer=TRUE, 'Broadleaf evergreen', line = -1.5)
}
# plot standardised, scaled PCA outputs
dev.new()
grid.arrange(p_KsLS, p_LSTLP, p_KSTLP, p_LSP50, nrow = 2)
#grid.arrange(p_LSTLP_old,p_LSP50_old, nrow = 1)
grid.arrange(p_KsLS, p_LSP50, nrow = 1)


ggplot(opt.pca,aes(x=PC1,y=PC2 ))
plot(opt.pca$x[,1],opt.pca$x[,2],pch=16,cex=0.5,col='grey',ylim =c(-6,6),xlim=c(-8,8))

s <- summary(opt.pca)
plot(opt.pca$x[,1],opt.pca$x[,2],
     xlab=paste("PCA 1 (", round(s$importance[2]*100, 1), "%)", sep = ""),
     ylab=paste("PCA 2 (", round(s$importance[5]*100, 1), "%)", sep = ""), pch=16, col='grey', cex=0.5, las=1, asp=1,ylim =c(-6,6),xlim=c(-8,8),
     main= opt_traits[o])
# Benjaminbell.co.uk
# Add grid lines
abline(v=0, lty=2, col="grey50")
abline(h=0, lty=2, col="grey50")
# Add labels
text(opt.pca$x[,1], opt.pca$x[,2], labels=row.names(opt.pca$x), pos=c(1,3,4,2), font=2)

#rotation
l.x <- opt.pca$rotation[,1]*10
l.y <- opt.pca$rotation[,2]*10
# Draw arrows
arrows(x0=0, x1=l.x, y0=0, y1=l.y, col="red", length=0.15, lwd=1.5)
lines(opt.pca)

# Label position
l.pos <- l.y # Create a vector of y axis coordinates
lo <- which(l.y < 0) # Get the variables on the bottom half of the plot
hi <- which(l.y > 0) # Get variables on the top half
# Replace values in the vector
l.pos <- replace(l.pos, lo, "1")
l.pos <- replace(l.pos, hi, "3")
# Variable labels
text(l.x, l.y, labels=row.names(opt.pca$rotation), col="red", pos=l.pos)

# screeplot(opt.pca)

# PCA on observed data
# not enough data points
# dev.new()
# grid.arrange(p_obs_BE, p_obs_BDT,  p_obs_B ,ncol = 3)
# not enough data points

#----------------------------------------------------------------------------------------------------------------------
# Write out LPJGuess .ins files
#----------------------------------------------------------------------------------------------------------------------

#Create LPJGuess insfiles-----------------------------------------------------------------------------
###old: write .ins file for LPJGuess based on the above optimised trait combinations for 28 PFT 'variants' per PFT.
# write .ins file for LPJGuess based on the above optimised trait combinations for 28 PFT 'variants' per PFT.

# Select which base PFT to use: TeBE (1), TeBS (2), IBS (3), TrBE (4) or TrBR (5)
#basePFT = 4

#NOT IN USE (yet) Select trait combination from which to write output using their location in the list traits_LPJG:
# traits_LPJG_KSLS (1) , traits_LPJG_LSTLP (2), etc
#traits_LPJG = c(traits_LPJG_KSLS, traits_LPJG_LSTLP, traits_LPJG_KSTLP, traits_LPJG_LSP50)
#tc = 4

# Selected traits from latin hypercube sampling with n=200 starting values in the predictors.
# each optimisation (KSLSstart or LSP50start) removed some predicted values. The PCA plots were used to manually select
# a PFT-subset for each predictor starting combination.
# The locations for the loadings in the PCA can shift depending on the computer and R version used,
# so this selection may not reflect the local selection. The subset selected here is saved as .RData 
# I selected traits relative to the loading-variables' locations and not based on their position along the PC1 and PC2
# axis. Not sure that was correct..


if(spec_group_sel==2){# evergreen
# for traits_LPJG_KSLS
ind = c(1,115,79,33,118,151,101,34,153,63,23,22,179,145,192,131,148,185,132,25,123,119,15,43,71,177,50,190,26,157)
traits_LPJG_KSLS_BE <- lapply(seq_along(traits_LPJG_KSLS), function(x) traits_LPJG_KSLS[[x]][ind])
names(traits_LPJG_KSLS_BE) <- names(traits_LPJG_KSLS)

# for traits_LPJG_LSP50
ind =  c(1,115,119,31,123,104,182,34,106,102,23,143,153,47,85,185,29,59,2,184,177,81,186,187,126,171,63,192,125,49)
traits_LPJG_LSP50_BE <- lapply(seq_along(traits_LPJG_LSP50), function(x) traits_LPJG_LSP50[[x]][ind])
names(traits_LPJG_LSP50_BE) <- names(traits_LPJG_LSP50)

#save new PFT subset or load existing one: 
#save(traits_LPJG_KSLS_BE, traits_LPJG_LSP50_BE, file = 'LPJGuessPFTS_BE.RData')
#load('LPJGuessPFTS_BE.RData')
}

if(spec_group_sel==4){# deciduous
  # for traits_LPJG_KSLS
  ind = c(150,52,116,1,48,11,99,93,122,46,80,164,161,65,97,166,128,67,149,123,119,88,89,127,96,136,105,17,23,16,19)
  traits_LPJG_KSLS_BDT <- lapply(seq_along(traits_LPJG_KSLS), function(x) traits_LPJG_KSLS[[x]][ind])
  names(traits_LPJG_KSLS_BDT) <- names(traits_LPJG_KSLS)
  
  # for traits_LPJG_LSP50
  ind =  c(117,69,140,176,53,55,113,5,115,26,15,146,149,45,111,17,97,67,159,27,75,89,163,38,132,126,152,41,123,184,19)
  traits_LPJG_LSP50_BDT <- lapply(seq_along(traits_LPJG_LSP50), function(x) traits_LPJG_LSP50[[x]][ind])
  names(traits_LPJG_LSP50_BDT) <- names(traits_LPJG_LSP50)
  
  #save new PFT subset or load existing one: 
  save(traits_LPJG_KSLS_BDT, traits_LPJG_LSP50_BDT, file = 'LPJGuessPFTS_BDT.RData')
  #load('LPJGuessPFTS_BE.RData')
}



# Select output folder
#output_fol="/Users/pughtam/Documents/TreeMort/Analyses/Hydraulic_modelling/Traits/uncer_test_KsLS/revised_PFTs_141220"
#AHES commented out for now during testing
#output_fol="/Users/annemarie/Documents/1_TreeMort/2_Analysis/1_Inputs"
output_fol="/Users/annemarie/Desktop"
if(spec_group_sel==2){# broadleaf evergreen in TRY
  # Select which base PFT to use: TeBE (1), TeBS (2), IBS (3), TrBE (4) or TrBR (5)
  basePFT = 4  # tropical broadleaf evergreen PFT for LPJGuess
  # create .ins files for  LPJ-GUESS_hydro
  #started with KSLS
  write_LPJG_ins.file(output_fol,basePFT = basePFT ,traits_LPJG = traits_LPJG_KSLS_BE)
  #started with LSP50
  write_LPJG_ins.file(output_fol,basePFT = basePFT ,traits_LPJG = traits_LPJG_LSP50_BE)
}
if(spec_group_sel==4){# 4 = TBD tropical broadleaf deciduous in TRY
  basePFT= 5 # tropical raingreen PFT for LPJGuess
  # create .ins files for  LPJ-GUESS_hydro
  #started with KSLS
  write_LPJG_ins.file(output_fol,basePFT = basePFT ,traits_LPJG = traits_LPJG_KSLS_BDT)
  #started with LSP50
  write_LPJG_ins.file(output_fol,basePFT = basePFT ,traits_LPJG = traits_LPJG_LSP50_BDT)
}

#write_LPJG_ins.file(output_fol,basePFT = 2 ,traits_LPJG)
#write_LPJG_ins.file(output_fol,basePFT = 3 ,traits_LPJG)
#write_LPJG_ins.file(output_fol,basePFT = 4 ,traits_LPJG)
#write_LPJG_ins.file(output_fol,basePFT = 5 ,traits_LPJG)

