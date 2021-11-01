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

# for sanity checking optimisation results
source('test_cor_signs.R')

# pretty plotting of PCA results
source('pca_with_pretty_biplot.R')
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

#--- Experiment with different plausible multivariate linear models, based on our hypotheses ---
# multivariate lm and sma regression possible. change option here:
regr_type = 'lm'

#--- Read in the trait data ---

######
##temporary testing of TLP P50 things:
#trait_CC<-subset(traits,group=="CC",drop = T)
#plot(abs(exp(traits$P50)) -abs(traits$TLP) )
#points( abs(exp(trait_CC$P50)) -abs(trait_CC$TLP), col = 'purple')
#points(abs(exp(trait_BE$P50)) -abs(trait_BE$TLP), col = 'green')
#points(abs(exp(trait_BDT$P50)) -abs(trait_BDT$TLP), col = 'brown')

## replace 36 potentially problematic values of P50 with NA. Rather arbitraty cutoff value of 0.5.
## Conifers don't show these low values, and P50 measurements in them are usually more accurate because Conifers have shorter tracheids.
## idx <- traits[which(abs(exp(traits$P50)) - abs(traits$TLP) < 0.5),]
######


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
subs <- na.omit(trait_BDT[,c("WD","slope")])

cor(subs$WD,subs$slope)
cor.test(subs$LMA,subs$Ks)
plot(-exp(subs$P50),subs$slope)
plot(subs$Ks,subs$WD)
#View(bivar_BE$all_sma_bivar)

# Calculate for all deciduous broadleaf
bivar_BDT <- make_bivar_plots(trait_BDT,nbtstrp, regr_type = regr_type)
# View(bivar_BDT$all_sma_bivar)
# AHES: negative bivatiate correlation in data, but sma shows positive correlation 
subs <- na.omit(trait_BD[,c("P50","LS")])
cor(subs$P50,subs$LS)
cor.test(subs$WD,subs$P50)

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

WD_multivar_BDT <- WD_multivar_test(trait_BDT,leaf_type='BDT',regr_type = regr_type)

WD_multivar_BE <- WD_multivar_test(trait_BE,leaf_type='BE',regr_type = regr_type)

# slope fits --------------------------------------------------------------

slope_multivar <- slope_multivar_test(trait_B,regr_type = regr_type)
#hist(trait_B$slope[slope_multivar$slope_from_P50_TLP_Ks$dataused]-slope_multivar$slope_from_P50_TLP_Ks$var_est,xlab ="slope-slope_est", main=paste(" distribution of residuals; model: ", Ks_multivar$Ks_from_P50_L$`regression_type (lm or sma)` ))
#shapiro.test(trait_B$slope[slope_multivar$slope_from_P50_TLP_Ks$dataused]-slope_multivar$slope_from_P50_TLP_Ks$var_est)


# Ks fits ----------------------------------------------------------------

Ks_multivar     <- Ks_multivar_test(trait_B,regr_type = regr_type)

# LS fits -----------------------------------------------------------------
# Separate for BE and BDT (BD + BT) on the basis that LS has a very different range and set of bivariate relationships for these
# two different groups, unlike the other traits here.

LS_multivar_BE <- LS_multivar_test(trait_BE, leaf_type ='BE',regr_type = regr_type) # returns LS_from_LMA_TLP_Ks

LS_multivar_BDT <- LS_multivar_test(trait_BDT, leaf_type ='BDT',regr_type = regr_type) # returns LS_from_TLP_Ks



# Make plots showing quality of fits and climate coverage -----------------


# Derive traits that are not subject to optimisation--------------------------

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
lc_sample <- optimumLHS(n,k = 4, maxSweeps = 2, eps = 0.1, verbose = FALSE)
# the drawn lc_samples from the sampling space are normally distributed and between 0 and 1

# The trait observations are not strictly normally distributed (according to shapiro- wilk test) 
# and therefore may need transformation.
# At the very least lhc- sampled value ranges need to be adapted to reflect the trait values. 
# Below are three options that were tested:
# Option 1) transform sampled marginal distributions to uniform distributions:
# that way we get stratified sample values that span the whole min-max range of the trait distribution and
# therefore the potential hydraulic strategy space.
#LS_e_start  <- qunif(lc_sample[,1], min(traits$LS[ind]),  max(traits$LS[ind])) 
#Ks_e_start  <- qunif(lc_sample[,2], min(traits$Ks[ind]),  max(traits$Ks[ind])) 
#TLP_e_start <- qunif(lc_sample[,3], min(traits$TLP[ind]), max(traits$TLP[ind])) 
#P50_e_start <- qunif(lc_sample[,4], min(traits$P50[ind]), max(traits$P50[ind])) 

# Option 2) adjust sampled marginal distributions to the variables' normal distributions ranges:
# this way we get stratified sample values that span the more commonly observed trait distribution
LS_e_start  <- qnorm(lc_sample[,1], mean = mean(traits$LS[ind]), sd = sd(traits$LS[ind])) 
Ks_e_start  <- qnorm(lc_sample[,2], mean = mean(traits$Ks[ind]),  sd = sd(traits$Ks[ind])) 
TLP_e_start <- qnorm(lc_sample[,3], mean = mean(traits$TLP[ind]), sd = sd(traits$TLP[ind])) 
P50_e_start <- qnorm(lc_sample[,4], mean = mean(traits$P50[ind]), sd = sd(traits$P50[ind])) 

# Option 3) for student t distributions
# this way we get stratified sample values that span the more commonly observed trait distribution, 
# allowing for the outer edges to have been sampled more frequently
# first, estimate the parameters for each t-distribution
#estimate_trait_distr(traits)
# then use these parameters to transform the lh_sample normal distribution into t-distr.
#LS_e_start  <- LaplacesDemon::qst(lc_sample[,1], mu=fitLS$par.ests[2]-1, sigma=fitLS$par.ests[3], nu=fitLS$par.ests[1], lower.tail=TRUE, log.p=FALSE)
#Ks_e_start  <- LaplacesDemon::qst(lc_sample[,2], mu=fitKs$par.ests[2],  sigma=fitKs$par.ests[3], nu=fitKs$par.ests[1],lower.tail=TRUE, log.p=FALSE)
#TLP_e_start <- LaplacesDemon::qst(lc_sample[,3], mu=fitTLP$par.ests[2], sigma=fitTLP$par.ests[3], nu=fitTLP$par.ests[1],lower.tail=TRUE, log.p=FALSE)
#P50_e_start <- LaplacesDemon::qst(lc_sample[,4], mu=fitP50$par.ests[2], sigma=fitP50$par.ests[3], nu=fitP50$par.ests[1],lower.tail=TRUE, log.p=FALSE)

#create objects with transformed values from above, which will act as starting values (= predictors) for the optimisation below
est_lhsLSP50 <- list(LS_e_start,P50_e_start)
names(est_lhsLSP50) <- c('LS_e_start','P50_e_start')

est_lhsLSTLP <- list(LS_e_start,TLP_e_start)
names(est_lhsLSTLP) <- c('LS_e_start','TLP_e_start')

est_lhsKsLS <- list(LS_e_start,Ks_e_start)
names(est_lhsKsLS) <- c('LS_e_start','Ks_e_start')

est_lhsKsTLP <- list(Ks_e_start,TLP_e_start)
names(est_lhsKsTLP) <- c('Ks_e_start','TLP_e_start')

# Check the distribution of the sample in relation to the observations
LSKsplot  <- check_lhs_sampling(traits$LS[ind],traits$Ks[ind],LS_e_start,Ks_e_start,'LS','Ks')
LSP50plot <- check_lhs_sampling(traits$LS[ind],traits$P50[ind],LS_e_start,P50_e_start,'LS','P50')
LSTLPplot <- check_lhs_sampling(traits$LS[ind],traits$TLP[ind],LS_e_start,TLP_e_start,'LS','TLP')

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
propagate_uncer=T

# Decide whether to run all trait combinations in the database for LS and Ks (F), or just a selection (T), T useful for generating output for LPJ-Guess
# and useful for testing different sampling methods  ( e.g. latin hypercube vs. systematic vs. hypervolume)
trait_sel= T

# Number of combinations to select if trait_sel=T. Set to -1 for a systematic sample, >0 for a random sample of the size specified, we have created 28 PFTs.
# or set = 4 for a predefined (above) hypercube sample.
n_trait_sel= -1

# Run for all deciduous (BT + BD) (=1), or BE (=2), or BT (=3), or BD (=4). This is used to set the maximum and minimum bounds in trait_opt().
spec_group_sel = 1

#Based on the above decision, determine trait dataset to use for plotting against optimised data
if (spec_group_sel==1 | spec_group_sel==3 | spec_group_sel==4) {
  trait_plot = trait_BDT
} else if (spec_group_sel==2) {
  trait_plot = trait_BE
}



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
traits_LPJG_KSTLP <- lpjg_traits_conv(LMA_e_mean,P50_e_mean,as.vector(TLP_e[,1]),slope_e_mean,
                                LS_e_mean,WD_e_mean,as.vector(Ks_e[,1]),
                                leafL_from_LMA,leafN_from_LMA,leafN_from_LMA_limit)

# create data frame of traits for subsequent PCA below -------
traits_KSTLP.df <- data.frame(matrix(unlist(traits_LPJG_KSTLP), ncol=length(traits_LPJG_KSTLP), byrow=FALSE))
names(traits_KSTLP.df) <- names(traits_LPJG_KSTLP)


# # create correlation table to report and test for sign in correlation being correct --------
tt <- test_cor_signs(trait_BDT,data.frame(LMA=LMA_e_mean,P50=P50_e_mean,TLP=as.vector(TLP_e[,1]),slope=slope_e_mean,LS = LS_e_mean,WD =WD_e_mean,Ks =as.vector(Ks_e[,1])))
tt
# clean up to not contaminate the below workflow
rm(list=(ls(pattern="_e_")))
rm(Ks_e,TLP_e)

# Optimisation with LS and P50 ------------------------------------------------------------
# lowest bivariate relationship in the trait network for evergreen subset: parson cor = 0.23.
# it is thought to have a functional relationship.
# Attempt to iteratively converge on the best fit values of Ks, TLP ,slope, WD and LMA, given known LS and P50
outs_LSP50 <- trait_optim_bivar_start_LSP50(limitdataranges = limitdataranges ,propagate_uncer = propagate_uncer,trait_sel = trait_sel, n_trait_sel = n_trait_sel, spec_group_sel = spec_group_sel,est_lhs = est_lhsLSP50)

# to 'release' the output from function trait_optim_bivar_startLSP50 from a list of objects into single objects
# single objects, 
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
traits_LPJG_LSP50 <- lpjg_traits_conv(LMA_e_mean,as.vector(P50_e[,1]),TLP_e_mean,slope_e_mean,
                                      as.vector(LS_e[,1]),WD_e_mean,Ks_e_mean,
                                leafL_from_LMA,leafN_from_LMA,leafN_from_LMA_limit)

# create data frame of traits for subsequent PCA below -------
traits_LSP50.df <- data.frame(matrix(unlist(traits_LPJG_LSP50), ncol=length(traits_LPJG_LSP50), byrow=FALSE))
names(traits_LSP50.df) <- names(traits_LPJG_LSP50)

# # create correlation table to report and test for sign in correlation being correct --------
tt <- test_cor_signs(trait_BDT,data.frame(LMA=LMA_e_mean,P50=as.vector(P50_e[,1]),TLP=TLP_e_mean,slope=slope_e_mean,LS = as.vector(LS_e[,1]),WD =WD_e_mean,Ks = Ks_e_mean))
tt
# clean up to not contaminate the below workflow
rm(list=(ls(pattern="_e_")))
rm(LS_e,P50_e)

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

## create correlation table to report and test for sign in correlation being correct --------
tt <- test_cor_signs(trait_BDT,data.frame(LMA = LMA_e_mean, P50 = P50_e_mean, TLP = TLP_e_mean, slope = slope_e_mean, LS = as.vector(LS_e), WD = WD_e_mean, Ks = as.vector(Ks_e)))
tt
# clean up
rm(list=(ls(pattern="_e_")))
rm(Ks_e,LS_e)


## sanity checking parameter results plausibility:
# P50 should always be lower than TLP ?

traits_KSLS.df$TLP > traits_KSLS.df$P50
traits_KSTLP.df$TLP > traits_KSTLP.df$P50
traits_LSP50.df$TLP > traits_LSP50.df$P50
traits_LSTLP.df$TLP > traits_LSTLP.df$P50



### collect rmse 
if(spec_group_sel==2){
rmse_BE <- t(data.frame(trait_pair=c('Ks-LS','Ks-TLP','LSP-50'), 
                        RMSE_BE = c(sum(RMSE_withKsLS_start$all_e_rmse),
                                 sum(RMSE_withKsTLP_start$all_e_rmse),
                                 sum(RMSE_withLSP50_start$all_e_rmse),
                                # sum(RMSE_withTLPLS_start$all_e_rmse)
                                 )
                        ))
}
if(spec_group_sel==1){
  rmse_DBT <- t(data.frame(trait_pair=c('Ks-LS','Ks-TLP','LSP-50'), 
                          RMSE_BDT = c(sum(RMSE_withKsLS_start$all_e_rmse),
                                      sum(RMSE_withKsTLP_start$all_e_rmse),
                                      sum(RMSE_withLSP50_start$all_e_rmse),
                                     # sum(RMSE_withTLPLS_start$all_e_rmse)
                                     )
  ))
}




#----------------------------------------------------------------------------------------------------------------------
# PCA
#----------------------------------------------------------------------------------------------------------------------
#trait_BE
# perform PCA ------------------------------------------------------------
# used as test to see whether starting from different trait combinations has an impact on the PFT-spread.

# save objects from analysis above, useful for quick jump to this section, if only PCA is of interest
traits_after_opt        <- c(traits_KSLS.df, traits_KSTLP.df, traits_LSP50.df)
names(traits_after_opt) <- c('traits_KSLS.df', 'traits_KSTLP.df', 'traits_LSP50.df')
#save(traits_after_opt,file='traits_after_opt_BDT_lm.RData')

#load('traits_after_opt_BE_sma.RData')
#list2env(traits_after_opt , envir = .GlobalEnv)

# perform PCA with observed data:
#traits_obs_BE  <- transform_obs_for_PCA(trait_BE)
#traits_obs_BDT <- transform_obs_for_PCA(trait_BDT)
#traits_obs_B   <- transform_obs_for_PCA(trait_B)


par(mfrow=c(1,3))
#par(mfrow=c(1,1))
opt_traits <-  c('traits_KSLS.df', 'traits_KSTLP.df', 'traits_LSP50.df')#
for(o in 1:3){
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
  
  p_KsLS  <- ggbiplot(opt.pca,labels=1:dim(traits_PCA)[1],labels.size = 2,circle=TRUE) +
  #labs(y= "y axis name", x = "x axis name") +
  ggtitle('KsLS')

pca_with_pretty_biplot(traits_PCA[,c(1,3,4,6,10,12,17)])
##sign of correlations correct?
test_cor_signs(trait_BE, traits_PCA)

}


if(o==2){
  pca_with_pretty_biplot(traits_PCA[,c(1,3,4,6,10,12,17)])
  
 
p_LSTLP <- ggbiplot(opt.pca,labels=1:dim(traits_PCA)[1],labels.size = 2) +
  #labs(y= "y axis name", x = "x axis name") +
  ggtitle('KsTLP')

##sign of correlations correct?
test_cor_signs(trait_BE, traits_PCA)

}
if(o==3){
  pca_with_pretty_biplot(traits_PCA[,c(1,3,4,6,10,12,17)])
  
p_KSTLP <- ggbiplot(opt.pca,labels=1:dim(traits_PCA)[1],labels.size = 2) +
  #labs(y= "y axis name", x = "x axis name") +
  ggtitle('LSP50')
##sign of correlations correct?
test_cor_signs(trait_BE, traits_PCA)
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
  # for traits_LPJG_LSP50
  traits_LPJG_LSP50_BE <- traits_LPJG_LSP50
  
  #save new PFT subset or load existing one: 
  save(traits_LPJG_LSP50_BE, file = 'LPJGuessPFTS_BE.RData')
  #load('LPJGuessPFTS_BE.RData')
}

if(spec_group_sel==1){# deciduous
  # for traits_LPJG_LSP50
  traits_LPJG_LSP50_BDT <- traits_LPJG_LSP50
  
  #save new PFT subset or load existing one: 
  save(traits_LPJG_LSP50_BDT, file = 'LPJGuessPFTS_BDT.RData')
  #load('LPJGuessPFTS_BE.RData')
}



# Select output folder
#output_fol="/Users/pughtam/Documents/TreeMort/Analyses/Hydraulic_modelling/Traits/uncer_test_KsLS/revised_PFTs_141220"
#AHES commented out for now during testing
#output_fol="/Users/annemarie/Documents/1_TreeMort/2_Analysis/1_Inputs"

if(spec_group_sel==2){# broadleaf evergreen in TRY
  output_fol="/Users/annemarie/Desktop/TrBE"
  # Select which base PFT to use: TeBE (1), TeBS (2), IBS (3), TrBE (4) or TrBR (5)
  basePFT = 4  # tropical broadleaf evergreen PFT for LPJGuess
  # create .ins files for  LPJ-GUESS_hydro
  #started with KSLS
  #write_LPJG_ins.file(output_fol,basePFT = basePFT ,traits_LPJG = traits_LPJG_KSLS_BE)
  #started with LSP50
  write_LPJG_ins.file(output_fol,basePFT = basePFT ,traits_LPJG = traits_LPJG_LSP50_BE)
}
if(spec_group_sel==1){# 1 = TBD tropical and temperate broadleaf deciduous in TRY
  output_fol="/Users/annemarie/Desktop/TrBR"
  basePFT= 5 # tropical raingreen PFT for LPJGuess
  # create .ins files for  LPJ-GUESS_hydro
  #started with KSLS
  #write_LPJG_ins.file(output_fol,basePFT = basePFT ,traits_LPJG = traits_LPJG_KSLS_BDT)
  #started with LSP50
  write_LPJG_ins.file(output_fol,basePFT = basePFT ,traits_LPJG = traits_LPJG_LSP50_BDT)
}

#write_LPJG_ins.file(output_fol,basePFT = 2 ,traits_LPJG)
#write_LPJG_ins.file(output_fol,basePFT = 3 ,traits_LPJG)
#write_LPJG_ins.file(output_fol,basePFT = 4 ,traits_LPJG)
#write_LPJG_ins.file(output_fol,basePFT = 5 ,traits_LPJG)

