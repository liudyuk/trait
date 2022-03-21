# moved other starting-pair network optimisations here
# Script to read in the processed trait data, make bivariate and multivariate SMA regressions and then
# carry out an optimisation procedure to unify inter-trait relationships across the whole trait dataset.
#
# [TO DO] shorten axis labels in PCA plots

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
# *Switch for LS predictions for BDT and BE data. 
# *added linear regression model to determine trait network relations, compare against sma
# *include WD into trait-network predictions
# December 2021
# *new ins-files after bugfix
# January:
# outsourced trait-pair starting optimisations for network and the PCA viarual and comparisons into trait_rel_supplementaries.R
# latin hypercube sampling is also in that script


#clear working space
rm(list=ls())

nbtstrp=1000 # Number of bootstrap samples to take in sma_multivar_regress (samples later used to calculated uncertainty in the optimisation). Was previously 10 000, using a lower number for testing, Will need to check sensitivity to this value.

#traits=read.table("/Users/liudy/trait_data/woody_trait.0625.txt")
#traits=read.csv("/Users/pughtam/Documents/TreeMort/Analyses/Hydraulic_modelling/Traits/mytrait-data/woody_trait.0803.txt",sep="\t")
#traits=read.csv("/Users/pughtam/Documents/TreeMort/Analyses/Hydraulic_modelling/Traits/woody_trait.0827.txt",sep="\t")
traits = read.csv("/Users/annemarie/Documents/1_TreeMort/2_Analysis/2_data/2_intermediate/mytrait-data/woody_trait.0827.txt",sep="\t")

# related to trait network construction
source('sma_multivar_regress.R')
source('lm_regress_multivar.R')
source('pca_regress_multivar.R')
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

# show PFT-variant location within observed trait space and sample space
source('opt_test_plots_LSP50_pfts.R')

# pretty plotting of PCA results
source('pca_with_pretty_biplot.R')
source('pca_with_pretty_plot_PFTs.R')

#---------------packages--------------------------
# PCA plotting
#install_github("vqv/ggbiplot") # must have devtools installed to use install_github
library('ggbiplot')

#other aesthetics:
# to distinguish between PFTvariants
library(RColorBrewer) # to plot the PFTs in different colours
cols <- c(brewer.pal(9, "Set1") ,brewer.pal(8, "Set2"),brewer.pal(11, "Set3"),brewer.pal(3,'Paired'))

#for plotting of hypercube sampling results:
library(grid)
library(ggplot2)
library("ggpubr")
# to plot PCA ggplots in grid
library(gridExtra)
# convenience package for ols normality checks, required in function lm_regress_multivar
library(olsrr)
# for PCR regression models, now used in the network
library(pls)



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
trait_BE  <- subset(traits,group=="BE",drop = T)
trait_BDT <- subset(traits,group=="BD" | group=="BT",drop = T)

#drought deciduous broadleaves only
trait_BD <- subset(traits,group=="BD",drop = T)
trait_BT <- subset(traits,group=="BT",drop = T)

#----
# stretch the TLP values

#--- Bivariate plots and statistics with SMA regression ---
# These are used to fill out the hypothesis framework plot with R

# Calculate for all broadleaf 
bivar <- make_bivar_plots(trait_B,nbtstrp,regr_type = 'lm')
#View(bivar$all_sma_bivar)

# Calculate for all evergreen broadleaf
bivar_BE <- make_bivar_plots(trait_BE,nbtstrp,regr_type = 'lm')
#View(bivar_BE$all_sma_bivar)

# Calculate for all deciduous broadleaf
bivar_BDT <- make_bivar_plots(trait_BDT,nbtstrp, regr_type = 'lm')
# View(bivar_BDT$all_sma_bivar)


#--- Experiment with different plausible multivariate linear models, based on our hypotheses ---
# multivariate lm, sma, pc(r) or pls(r) regression possible. change option here:
regr_type = 'plsr'

if (regr_type == 'pcr' || regr_type == 'plsr'){
 
  traits_sd   = as.data.frame(t(apply(traits[13:44],2,sd,na.rm=TRUE)))
  traits_mean_unscale = as.data.frame(t(apply(traits[13:44],2,mean,na.rm=TRUE)))
  # now that centering is se to TRUE, and it is applied to all values, including Y, do not standardise using the mean, but only sd only.
  #Test that this is the case, setting mean to 0 for now, so that 0 will be subtracted from every X
   traits_mean = traits[1,13:44]
   traits_mean[1,] <- 0  # legacy- set traits_mean to 0, as the centering is now done by the pcr() function itself.
}

# for now here, in to its own file later:

scale_traits <- function(varst, labels, nlabels, traits_mean, traits_sd){ 
  # input:
  # varst dataframe of which the first nlables-1 variables will be scaled
  # labels =names of all traits in the dataframe
  # nlabels number of traits in the dataframe
  # traits_mean =  mean trait value from all traits in the dataset. will be subsetted for the correct one using 'labels'.
  # straits_sd = standard deviation from all traits in the dataset. will be subsetted for the correct one using 'labels'.
  # scaling is only performed on the explanatory variables X, because scale=TRUE only scales X variables. Therefore I do nlables-1
  # as index nlabels contains Y.
  
  # determine context in which the function is used ( e.g. in generating the model or predictions)
  # for scaling in generating the models, "varst" is multidimensional:
  if(!is.null(dim(varst))){ 
    ns <- names(varst)[1:(nlabels-1)]
    print(paste('scaling',ns))
    for (n in 1:length(ns)){ # normalise predictor traits only
      varst[,ns[n]] <- (varst[,ns[n]] - traits_mean[,labels[n]] ) / traits_sd[,labels[n]] 
    }
  }else{  # for scaling of single values during the network predictions, varst is a single value ( could be changed to a list, but done like this for now)
    if(labels == "TLP_e_last"  || labels == "TLP_e_start" || labels == "TLP_e" ){
      varst <- (varst - traits_mean[,"TLP"] ) / traits_sd[,"TLP"] 
    }
    if(labels == "LS_e_last"   || labels == "LS_e_start" || labels == "LS_e" ){
      varst <- (varst - traits_mean[,"LS"] ) / traits_sd[,"LS"] 
    }
    if(labels == "LMA_e_last"   || labels == "LMA_e_start" || labels == "LMA_e" ){
      varst <- (varst - traits_mean[,"LMA"] ) / traits_sd[,"LMA"] 
    }
    if(labels == "P50_e_last"   || labels == "P50_e_start" || labels == "P50_e" ){
      varst <- (varst - traits_mean[,"P50"] ) / traits_sd[,"P50"] 
    }
    if(labels == "Ks_e_last"   || labels == "Ks_e_start" || labels == "Ks_e" ){
      varst <- (varst - traits_mean[,"Ks"] ) / traits_sd[,"Ks"] 
    }
    if(labels == "WD_e_last"   || labels == "WD_e_start" || labels == "WD_e" ){
      varst <- (varst - traits_mean[,"WD"] ) / traits_sd[,"WD"] 
    }
    if(labels == "slope_e_last"   || labels == "slope_e_start" || labels == "slope_e" ){
      varst <- (varst - traits_mean[,"slope"] ) / traits_sd[,"slope"] 
    }
  }
  return(varst)
}
unscale_traits <- function(varst, labels, nlabels, traits_mean, traits_sd){ 
# for unscaling of single values during the network predictions, varst is a single value ( could be changed to a list, but done like this for now)
    if(labels == "TLP_e_last"  || labels == "TLP_e_start" || labels == "TLP_e" ){
      varst <- varst * traits_sd[,"TLP"] + traits_mean[,"TLP"] 
    }
    if(labels == "LS_e_last"   || labels == "LS_e_start" || labels == "LS_e" ){
      varst <- varst * traits_sd[,"LS"]  + traits_mean[,"LS"] 
    }
    if(labels == "LMA_e_last"   || labels == "LMA_e_start" || labels == "LMA_e" ){
      varst <- varst * traits_sd[,"LMA"] + traits_mean[,"LMA"] 
    }
    if(labels == "P50_e_last"   || labels == "P50_e_start" || labels == "P50_e" ){
      varst <- varst * traits_sd[,"P50"] + traits_mean[,"P50"] 
    }
    if(labels == "Ks_e_last"   || labels == "Ks_e_start" || labels == "Ks_e" ){
      varst <- varst * traits_sd[,"Ks"] + traits_mean[,"Ks"] 
    }
    if(labels == "WD_e_last"   || labels == "WD_e_start" || labels == "WD_e" ){
      varst <- varst * traits_sd[,"WD"] + traits_mean[,"WD"] 
    }
    if(labels == "slope_e_last"   || labels == "slope_e_start" || labels == "slope_e" ){
      varst <- varst * traits_sd[,"slope"] + - traits_mean[,"slope"] 
    }
  
  return(varst)
}


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

LMA_multivar_BE  <- LMA_multivar_test_BE(trait_BE, regr_type =  regr_type)

LMA_multivar_BDT <- LMA_multivar_test_BDT(trait_BDT, regr_type =  regr_type)



# WD fits -----------------------------------------------------------------

#WD_multivar_BDT <- WD_multivar_test(trait_BE,leaf_type='BDT',regr_type = regr_type)

#WD_multivar_BE <- WD_multivar_test(trait_BE,leaf_type='BE',regr_type = regr_type)

WD_multivar <- WD_multivar_test(trait_B,leaf_type=NULL,regr_type = regr_type)

# slope fits --------------------------------------------------------------

slope_multivar <- slope_multivar_test(trait_B,regr_type = regr_type)
#AHES hist(trait_B$slope[slope_multivar$slope_from_P50_TLP_Ks$dataused]-slope_multivar$slope_from_P50_TLP_Ks$var_est,xlab ="slope-slope_est", main=paste(" distribution of residuals; model: ", Ks_multivar$Ks_from_P50_L$`regression_type (lm or sma)` ))
#AHES shapiro.test(trait_B$slope[slope_multivar$slope_from_P50_TLP_Ks$dataused]-slope_multivar$slope_from_P50_TLP_Ks$var_est)


# Ks fits ----------------------------------------------------------------

Ks_multivar     <- Ks_multivar_test(trait_B,regr_type = regr_type)

# LS fits -----------------------------------------------------------------
# Separate for BE and BDT (BD + BT) on the basis that LS has a very different range and set of bivariate relationships for these
# two different groups, unlike the other traits here.

LS_multivar_BE  <- LS_multivar_test(trait_BE, leaf_type ='BE',regr_type = regr_type) # returns LS_from_LMA_TLP_Ks

LS_multivar_BDT <- LS_multivar_test(trait_BDT, leaf_type ='BDT',regr_type = regr_type) # returns LS_from_TLP_Ks


# Make plots showing quality of fits and climate coverage -----------------

#[TO DO] where has this code gone? retrieve from previous version


# Derive traits that are not subject to optimisation--------------------------

# Calculate regression of leafL from LMA ----------------------------------
# sticking with linear regression here:

leafL_from_LMA <- sma_plot_stats(data.frame(trait_B$LMA,log(trait_B$leafL)),c("LMA","leafL"),nbtstrp,T)


# Calculate limits of leafN vs LMA to allow estimate of leaf C:N ----------
# sticking with linear regression here:
leafN_from_LMA <- sma_plot_stats(data.frame(trait_B$LMA,trait_B$leafN),c("LMA","leafN"),nbtstrp,T)

leafN_from_LMA_limit <- regress_limit_adjust(trait_B$leafN,trait_B$LMA,leafN_from_LMA,0.05)

plot(trait_B$LMA,trait_B$leafN)
points(trait_B$LMA[leafN_from_LMA_limit$ind],leafN_from_LMA_limit$var1_pred_lower,col="green")
points(trait_B$LMA[leafN_from_LMA_limit$ind],leafN_from_LMA_limit$var1_pred_upper,col="red")


# Create PCA regression-models (later called PCR models) on what we have already determined to be the relationships between the traits 
# through functional knowledge on that they should interact, and checked above using 
# correlation analysis to double-check that they indeed also interact in the dataset ( otherwise we won't be able to predict anything)
# - avoid multicollinearity
# - reduce dimensions, but keep most variability in predictors
# - reduce risk of overfitting

testing=FALSE
if(testing==TRUE){
  set.seed (1234)
  
  #TEST new function:
  
  #[TODO] scale + rescale manually
  sdtraits  <- as.data.frame(t(apply(na.omit(trait_BDT[c("TLP","Ks","WD","P50")])[1:50,],2,sd,na.rm=TRUE)))
  mtraits   <- as.data.frame(t(apply(na.omit(trait_BDT[c("TLP","Ks","WD","P50")])[1:50,],2,mean,na.rm=TRUE)))
  
  #training dataset:
  train <- na.omit(trait_BDT[c("TLP","Ks","WD","P50")])[1:50,]
  
  # model construction and coefficient exctraction
  pcr_model <- plsr(P50~., data = train[c("TLP","Ks","WD","P50")], scale = TRUE,center=TRUE, validation = "LOO",
                     na.action = na.omit,x=TRUE,y=TRUE)
  #opls_model = oscorespls.fit(X=as.matrix(train[c("TLP","Ks","WD")]) ,Y = as.matrix(train[c("P50")]),ncomp=3)
  l_model <- lm(P50~., data = train[c("TLP","Ks","WD","P50")], na.action = na.omit)
  l_coef   <- coef(l_model)
  pcr_coef <- coef(pcr_model,intercept=TRUE,ncomp=3)
  
  #testing dataset:
  test <- na.omit(trait_BDT[c("TLP","Ks","WD","P50")])[50:133,]
  
  #apply model for prediction
  # Why does scaling TLP, Ks and WD with their sd not give better results? I thought the model was constructed with them scaled?
  P50pred_pcr <- pcr_coef[1] + pcr_coef[2]*test$TLP + pcr_coef[3]*test$Ks + pcr_coef[4]*test$WD
  P50pred_l   <-  l_coef[1] + l_coef[2]*test$TLP +  l_coef[3]*test$Ks +  l_coef[4]*test$WD
  P50pred_pcr_auto <- predict(pcr_model,newdata = test,ncomp=1 )
 
  # plot
  plot(test$P50,P50pred_pcr_auto,pch=10,ylim=c(-0.5,2),xlim=c(-0.5,2))
  # centering = TRUE is the default in pcr, applied to X and Ys.  # therefore,if center?TRUE I thought I must 'uncenter' Y ( by doing + mtraits$P50)
  # or use centering= FALSE. I don't understand why 'uncentering' doesn't give the same result as the above though.
  points(test$P50,P50pred_pcr,col='orange',pch=5) 
  points(test$P50,P50pred_pcr+mtraits$P50,col='orange',pch=5) 
  points(test$P50,P50pred_l,pch=11,col='blue',cex=0.8) 
  legend("topleft",legend=c('lm','pcr manual','pcr+mean manual', 'pcr auto'),pch=c(11,5,5,10),col=c("black","orange","orange","blue"))
  
  #
  plot(RMSEP(pcr_model)) # shows that 1 PC achieves a Root Mean Squared Error of Prediction of < 53%
  #test individual points:
  points(test$P50[1], predict(pcr_model,newdata = test[1,],ncomp=1),col='red' )
  
  #test impact of ncomps
  points(test$P50, predict(pcr_model,newdata = test,ncomp=1), col ="grey" )
  points(test$P50, predict(pcr_model,newdata = test,ncomp=2), col ="dark grey" )
  points(test$P50, predict(pcr_model,newdata = test,ncomp=3), col ="yellow" )
  
  
  # the same as LM. Why?
  plot(test$P50,P50pred_l,pch=11,col='blue',cex=0.8) 
  points(test$P50, predict(pcr_model,newdata = test,ncomp=3), col ="yellow" )
  
  
  # what does predict() do that I don't understand
  # why is predict() using the pcr model giving the same result as the linear model? I thought the pcr model would be somehow different.
  # I thought the coefficients in the pcr model would consider the interactions between all Xes, and in the
  
  
  
  mod_boot <- pcr(DF2formula(yy[varnames[c(nvars,1:nvars-1)]]),scale=FALSE,center=TRUE, data = bootsample,x=TRUE,y=TRUE)
  
  DF2formula(yy[varnames[c(nvars,1:nvars-1)]])
  
  coef(pcr_model,intercept = TRUE)
  
  
  trait = trait_BDT
  P50_from_TLP_Ks_WD_plsr <- sma_plot_stats(data.frame(trait$TLP,trait$Ks,trait$WD,trait$P50),c("TLP","Ks","WD","P50"),nbtstrp,T, regression_type= 'plsr')
  P50_from_TLP_Ks_WD_pca  <- sma_plot_stats(data.frame(trait$TLP,trait$Ks,trait$WD,trait$P50),c("TLP","Ks","WD","P50"),nbtstrp,T, regression_type= 'pcr')
  P50_from_TLP_Ks_WD_lm   <- sma_plot_stats(data.frame(trait$TLP,trait$Ks,trait$WD,trait$P50),c("TLP","Ks","WD","P50"),nbtstrp,T, regression_type= 'lm')
  #P50_from_TLP_Ks_WD_sma  <- sma_plot_stats(data.frame(trait$TLP,trait$Ks,trait$WD,trait$P50),c("TLP","Ks","WD","P50"),nbtstrp,T, regression_type= 'sma')
 
  P50_from_TLP_Ks_WD_pca$mod$intercept_R ==  P50_from_TLP_Ks_WD_lm$mod$intercept_R
  P50_from_TLP_Ks_WD_pca$mod$slope_R.y1  ==  P50_from_TLP_Ks_WD_lm$mod$slope_R.y1
  P50_from_TLP_Ks_WD_pca$mod$slope_R.y2  ==  P50_from_TLP_Ks_WD_lm$mod$slope_R.y2
  P50_from_TLP_Ks_WD_pca$mod$slope_R.y3  ==  P50_from_TLP_Ks_WD_lm$mod$slope_R.y3
  
  P50_from_TLP_Ks_WD_lm$mod$intercept_R
  P50_from_TLP_Ks_WD_lm$mod$slope_R.y1
  P50_from_TLP_Ks_WD_lm$mod$slope_R.y2
  P50_from_TLP_Ks_WD_lm$mod$slope_R.y3
  
  # 490
  [1] 4.599416
  [1] -0.324917
  [1] "plsr"
  
  LS_multivar_BDT <- LS_multivar_test(LS~.,trait_BDT[c('P50','TLP','Ks','LS')], leaf_type ='BDT',regr_type = regr_type) # returns LS_from_TLP_Ks
  
  pcr_pred <- predict.pcr(pcr_model,trait_BDT[108,c('P50','TLP','Ks','LS')], ncomp = pcr_model$ncomp)
  predplot(pcr_model)
  
  pcr_model2 <- plsr(LS~., data = trait_BDT[c('P50','TLP','Ks','LS')], validation = "LOO", jackknife=T)
  plot(pcr_model2,col='green')
  validationplot(pcr_model2, val.type = "MSEP")
  validationplot()
  
  
  coef(pcr_model2,intercept=TRUE)
  
  jack.test(pcr_model2, ncomp = 3)
  pcr_model$Ymeans  # is this the intercept?
  # make sure newdata is a data frame when doing predictions or be sure the newdata matrix contains only predictors. From help:
  #It is also possible to supply a matrix instead of a data frame as newdata, which is then assumed to be the X data matrix. 
  # Note that the usual checks for the type of the data are then omitted. 
  #Also note that this is only possible with predict; it will not work in functions like predplot, 
  # RMSEP or R2, because they also need the response variable of the new data.
  
  regression_type
  leafN_from_LMA      <- sma_plot_stats(data.frame(trait_B$LMA,trait_B$leafN),c("LMA","leafN"),nbtstrp,T)
  leafL_from_LMA_sma  <- sma_plot_stats(data.frame(trait_B$LMA,log(trait_B$leafL)),c("LMA","leafL"),nbtstrp,T,regression_type = 'sma')
  leafL_from_LMA_lm   <- sma_plot_stats(data.frame(trait_B$LMA,log(trait_B$leafL)),c("LMA","leafL"),nbtstrp,T,regression_type = 'lm')
  leafL_from_LMA_pca  <- sma_plot_stats(data.frame(trait_B$LMA,log(trait_B$leafL)),c("LMA","leafL"),nbtstrp,T,regression_type = 'pcr')
  leafL_from_LMA_plsr <- sma_plot_stats(data.frame(trait_B$LMA,log(trait_B$leafL)),c("LMA","leafL"),nbtstrp,T,regression_type = 'plsr')
  
  trait = trait_B
  P50_from_TLP_Ks_WD_plsr <- sma_plot_stats(data.frame(trait$TLP,trait$Ks,trait$WD,trait$P50),c("TLP","Ks","WD","P50"),nbtstrp,T, regression_type= 'plsr')
  P50_from_TLP_Ks_WD_pca  <- sma_plot_stats(data.frame(trait$TLP,trait$Ks,trait$WD,trait$P50),c("TLP","Ks","WD","P50"),nbtstrp,T, regression_type= 'pcr')
  P50_from_TLP_Ks_WD_lm   <- sma_plot_stats(data.frame(trait$TLP,trait$Ks,trait$WD,trait$P50),c("TLP","Ks","WD","P50"),nbtstrp,T, regression_type= 'lm')
  P50_from_TLP_Ks_WD_sma  <- sma_plot_stats(data.frame(trait$TLP,trait$Ks,trait$WD,trait$P50),c("TLP","Ks","WD","P50"),nbtstrp,T, regression_type= 'sma')
  
  P50_from_TLP_Ks_WD_plsr$dataused == P50_from_TLP_Ks_WD_lm$dataused
  P50_from_TLP_Ks_WD_plsr$R == P50_from_TLP_Ks_WD_lm$R
  P50_from_TLP_Ks_WD_plsr$mod$intercept_R == P50_from_TLP_Ks_WD_lm$mod$intercept_R
  P50_from_TLP_Ks_WD_plsr$mod$slope_R.y1 == P50_from_TLP_Ks_WD_lm$mod$slope_R.y1
  P50_from_TLP_Ks_WD_sma$mod$slope_R.y1
  P50_from_TLP_Ks_WD_plsr$mod$slope_R.y1
  P50_from_TLP_Ks_WD_lm$mod$slope_R.y1
  P50_from_TLP_Ks_WD_pca$mod$slope_R.y1
  
  varnames =c("TLP","Ks","WD","P50")
  nvars=length(varnames)
  yy= vars=data.frame(trait$TLP,trait$Ks,trait$WD,trait$P50)
  pcr_model       <- pcr(DF2formula(yy[varnames[c(nvars,1:nvars-1)]]), data = yy, scale = TRUE, validation = "CV", na.action = na.omit)
  coef(pcr_model,intercept = TRUE)
  P50_from_TLP_Ks_WD_plsr <- sma_plot_stats(data.frame(trait$TLP,trait$Ks,trait$WD,trait$P50),c("TLP","Ks","WD","P50"),nbtstrp,T, regression_type= 'plsr')
  P50_from_TLP_Ks_WD_pca <- sma_plot_stats(data.frame(trait$TLP,trait$Ks,trait$WD,trait$P50),c("TLP","Ks","WD","P50"),nbtstrp,T, regression_type= 'pcr')
  
  P50_from_TLP_Ks_WD_pca$mod$intercept_R
  P50_from_TLP_Ks_WD_pca$mod$slope_R.y1
  P50_from_TLP_Ks_WD_pca$mod$slope_R.y2
  P50_from_TLP_Ks_WD_pca$mod$slope_R.y3
  
  P50_from_TLP_Ks_WD_lm  <- sma_plot_stats(data.frame(trait$TLP,trait$Ks,trait$WD,trait$P50),c("TLP","Ks","WD","P50"),nbtstrp,T, regression_type= 'lm')
  
  P50_from_TLP_Ks_WD_lm$mod$intercept_R
  P50_from_TLP_Ks_WD_lm$mod$slope_R.y1
  P50_from_TLP_Ks_WD_lm$mod$slope_R.y2
  P50_from_TLP_Ks_WD_lm$mod$slope_R.y3
  
  
}




#----------------------------------------------------------------------------------------------------------------------
# Optimisation
#----------------------------------------------------------------------------------------------------------------------
# do analysis only on Broadleaf data:
#traits <- trait_B
###
# Optimisation setup
# select options that will be used in all optimisations below

# Decide whether to limit the possible ranges of predicted traits to the observed values (T) or not (F)
limitdataranges=T # Currently does not converge in uncertainty propagation if not set to T

# Decide whether to run the uncertainty propagation (T) or not (F)
propagate_uncer= T

# Decide whether to run all trait combinations in the database for LS and Ks (F), or just a selection (T), T useful for generating output for LPJ-Guess
# and useful for testing different sampling methods  ( e.g. latin hypercube vs. systematic vs. hypervolume)
trait_sel= F

# Number of combinations to select if trait_sel=T. Set to -1 for a systematic sample, >0 for a random sample of the size specified, we have created 28 PFTs.
# or set = 4 for a predefined (above) hypercube sample.
n_trait_sel= 28

# Run for all deciduous (BT + BD) (=1), or BE (=2), or BT (=3), or BD (=4). This is used to set the maximum and minimum bounds in trait_opt().
spec_group_sel = 2

#Based on the above decision, determine trait dataset to use for plotting against optimised data
if (spec_group_sel==1 | spec_group_sel==3 | spec_group_sel==4) {
  trait_plot = trait_BDT
} else if (spec_group_sel==2) {
  trait_plot = trait_BE
}


# Optimisation with LS and P50 ------------------------------------------------------------
# lowest bivariate relationship in the trait network for evergreen subset: pearson cor = 0.23.
# it is thought to have no functional relationship.
# Attempt to iteratively converge on the best fit values of Ks, TLP ,slope, WD and LMA, given known LS and P50
outs_LSP50     <- trait_optim_bivar_start_LSP50(limitdataranges = limitdataranges ,propagate_uncer = propagate_uncer,trait_sel = trait_sel, n_trait_sel = n_trait_sel, spec_group_sel = spec_group_sel,est_lhs = est_lhsLSP50,regr_type = regr_type)
outs_LSP50_hv  <- trait_optim_bivar_start_LSP50(limitdataranges = limitdataranges ,propagate_uncer = propagate_uncer,trait_sel = T, n_trait_sel = -1, spec_group_sel = spec_group_sel,est_lhs = est_lhsLSP50,regr_type = regr_type)
outs_LSP50_lhc <- trait_optim_bivar_start_LSP50(limitdataranges = limitdataranges ,propagate_uncer = propagate_uncer,trait_sel = T, n_trait_sel = 4, spec_group_sel = spec_group_sel,est_lhs = est_lhsLSP50,regr_type = regr_type)
#save(outs_LSP50 , file= 'outs_LSP50.RData')
# to 'release' the output from function trait_optim_bivar_startLSP50 from a list of objects into single objects
# single objects, 
list2env(outs_LSP50$predictors , envir = .GlobalEnv) 
list2env(outs_LSP50$predicted , envir = .GlobalEnv)

# Stats defining the uncertainty range for each point
create_uncertainty_range_stats(outs_LSP50)

#trait_BE = trait_BE_save

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
tt <- test_cor_signs(trait_plot,data.frame(LMA=LMA_e_mean,P50=as.vector(P50_e[,1]),TLP=TLP_e_mean,slope=slope_e_mean,LS = as.vector(LS_e[,1]),WD =WD_e_mean,Ks = Ks_e_mean))
tt


#----------------------------------------------------------------------------------------------------------------------
# PCA
#----------------------------------------------------------------------------------------------------------------------
# perform PCA ------------------------------------------------------------
# used as test to see whether starting from different trait combinations has an impact on the PFT-spread.

traits_PCA  <- traits_LSP50.df

# T. Pugh
# 25.10.20
# original file: lpjg_strat_mapping_comb.m
# translated into R Annemarie Eckes-Shephard May 2021

## Convert SLA to LMA
traits_PCA$LMA = 1./traits_PCA$SLA

# transform some values:
## Log traits that are non-normal
traits_PCA$P50 = log(-traits_PCA$P50)
traits_PCA$P88 = log(-traits_PCA$P88)
traits_PCA$TLP = log(-traits_PCA$TLP)
traits_PCA$LS  = log(traits_PCA$LS)
traits_PCA$Ks  = log(traits_PCA$Ks)
traits_PCA$LMA = log(traits_PCA$LMA)

#PCA on optimised trait values (0 = 5 to 7)
pca_with_pretty_biplot(traits_PCA[,c(1,3,4,6,10,12,17)])

if(spec_group_sel==1){
  mtext(outer=TRUE, 'Broadleaf deciduous')
}
if(spec_group_sel ==2){
  mtext(outer=TRUE, 'Broadleaf evergreen')
}

##DECISION: take P50LS as trait pair to generate parametersets using the network


#----------------------------------------------------------------------------------------------------------------------
# Trait sampling for .insfile
#----------------------------------------------------------------------------------------------------------------------
# hyper-volume sampled PFTs, alongside extreme-value outer edges:
outs_LSP50_hv  <- trait_optim_bivar_start_LSP50(limitdataranges = limitdataranges ,propagate_uncer = propagate_uncer,trait_sel = T, n_trait_sel = -1, spec_group_sel = spec_group_sel,est_lhs = est_lhsLSP50,regr_type = regr_type)
save(outs_LSP50_hv,file='outs_LSP50_hv.RData')
#display trait values that will be selected for PFTs(purple), show that their spread is across a wide range of values 
opt_test_plots_LSP50_pfts(traits,#trait_B,#trait_plot,
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
                          slope_e,
                          outs_LSP50_hv = outs_LSP50_hv,
                          cols=cols) 

# to 'release' the output from function trait_optim_bivar_startLSP50 from a list of objects into single objects
# single objects, 
list2env(outs_LSP50_hv$predictors , envir = .GlobalEnv) 
list2env(outs_LSP50_hv$predicted , envir = .GlobalEnv)

# Stats defining the uncertainty range for each point
create_uncertainty_range_stats(outs_LSP50_hv)

traits_LPJG_LSP50_pft <- lpjg_traits_conv(LMA_e_mean,as.vector(P50_e[,1]),TLP_e_mean,slope_e_mean,
                                      as.vector(LS_e[,1]),WD_e_mean,Ks_e_mean,
                                      leafL_from_LMA,leafN_from_LMA,leafN_from_LMA_limit)

# create data frame of traits for another PCA only on the PFT subvariants below -------
traits_LSP50.df <- data.frame(matrix(unlist(traits_LPJG_LSP50_pft), ncol=length(traits_LPJG_LSP50_pft), byrow=FALSE))
names(traits_LSP50.df) <- names(traits_LPJG_LSP50_pft)

traits_PCA  <- traits_LSP50.df

# T. Pugh
# 25.10.20
# original file: lpjg_strat_mapping_comb.m
# translated into R Annemarie Eckes-Shephard May 2021

## Convert SLA to LMA
traits_PCA$LMA = 1./traits_PCA$SLA

# transform some values:
## Log traits that are non-normal
traits_PCA$P50 = log(-traits_PCA$P50)
traits_PCA$P88 = log(-traits_PCA$P88)
traits_PCA$TLP = log(-traits_PCA$TLP)
traits_PCA$LS  = log(traits_PCA$LS)
traits_PCA$Ks  = log(traits_PCA$Ks)
traits_PCA$LMA = log(traits_PCA$LMA)

#PCA on optimised trait values (0 = 5 to 7)
pca_with_pretty_biplot_Pfts(traits_PCA[,c(1,3,4,6,10,12,17)])

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
  traits_LPJG_LSP50_BE <- traits_LPJG_LSP50_pft
  
  #save new PFT subset or load existing one: 
  save(traits_LPJG_LSP50_BE, file = 'LPJGuessPFTS_BE21-03-2022plsr.RData')
  #load('LPJGuessPFTS_BE.RData')
}

if(spec_group_sel==1){# deciduous
  # for traits_LPJG_LSP50
  traits_LPJG_LSP50_BDT <- traits_LPJG_LSP50_pft
  
  #save new PFT subset or load existing one: 
  save(traits_LPJG_LSP50_BDT, file = 'LPJGuessPFTS_BDT14-03-2022lm.RData')
  #load('LPJGuessPFTS_BE.RData')
}



# Select output folder
#output_fol="/Users/pughtam/Documents/TreeMort/Analyses/Hydraulic_modelling/Traits/uncer_test_KsLS/revised_PFTs_141220"
#AHES commented out for now during testing
#output_fol="/Users/annemarie/Documents/1_TreeMort/2_Analysis/1_Inputs"
output_fol="/Users/annemarie/Desktop/"

if(spec_group_sel==2){# broadleaf evergreen in TRY
  output_fol="/Users/annemarie/Desktop/TrBE/plsr/"
  # Select which base PFT to use: TeBE (1), TeBS (2), IBS (3), TrBE (4) or TrBR (5)
  basePFT = 4  # tropical broadleaf evergreen PFT for LPJGuess
  # create .ins files for  LPJ-GUESS_hydro
  #started with KSLS
  #write_LPJG_ins.file(output_fol,basePFT = basePFT ,traits_LPJG = traits_LPJG_KSLS_BE)
  #started with LSP50
  write_LPJG_ins.file(output_fol,basePFT = basePFT ,traits_LPJG = traits_LPJG_LSP50_BE)
}
if(spec_group_sel==1){# 1 = TBD tropical and temperate broadleaf deciduous in TRY
  output_fol="/Users/annemarie/Desktop/TrBR/plsr/"
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

##Compare range of traits
pcr_TrBE_new <- read.csv(file='~/Desktop/TrBE/pcr_new/LPJG_PFT_summary_TrBE.csv',header = TRUE)
plsr_TrBE<- read.csv(file='~/Desktop/TrBE/plsr_new/LPJG_PFT_summary_TrBE.csv',header = TRUE)
#pcr_TrBE <- read.csv(file='~/Desktop/TrBE/pcr/LPJG_PFT_summary_TrBE.csv',header = TRUE)
lm_TrBE <- read.csv(file='~/Desktop/TrBE/lm/LPJG_PFT_summary_TrBE.csv',header = TRUE)
previous_TrBE <- read.csv(file='~/Documents/1_TreeMort/2_Analysis/1_Inputs/ins_files/LPJG_PFT_summary_TrBE_old.csv',header = TRUE)

par(mfrow=c(4,4))

if(testing==TRUE){
for (n in names(pcr_TrBE)){
  if(n %in% names(trait_plot)){
    if(n=='WD'){
      df <- data.frame(Obs=c(na.omit(trait_plot[n]))
      boxplot(as.matrix(data.frame(Obs=c(max(na.omit(trait_plot[n])),min(na.omit(trait_plot[n]))),PFTs=c(NA,NA))), boxfill = NA, border = NA)
      boxplot(at="Obs",t(trait_plot[n]),xlab='',ylab='',col='purple',xlim=c(0,6),main=n) 
    }else if(n=='P50'){
      boxplot(add=TRUE, rep(1,length(t(trait_plot[n]))),t(trait_plot[n]),xlab='',ylab='',col='purple',xlim=c(0,6),main=n) 
    }else if(n=='Ks'){
      plot(rep(1,length(t(trait_plot[n]))),exp(t(trait_plot[n])),xlab='',ylab='',col='purple',xlim=c(0,6),main=n)  
    }else if(n=='TLP'){
      plot(rep(1,length(t(trait_plot[n]))),-exp(t(trait_plot[n])),xlab='',ylab='',col='purple',xlim=c(0,6),main=n)  
    }else if (n=='LS'){
      plot(rep(1,length(t(trait_plot[n]))),exp(t(trait_plot[n]))*10000,xlab='',ylab='',col='purple',xlim=c(0,6),main=n)  
    }else if (n=='SLA'){# (1/LMA_e_mean_unlogged)*1000*2
      plot(rep(1,length(t(trait_plot['LMA']))),(1/exp(t(trait_plot['LMA'])))*1000*2,xlab='',ylab='',col='purple',xlim=c(0,6),main=n)  
    }else if(n=='slope'){
      plot(rep(1,length(t(trait_plot[n]))),exp(t(trait_plot[n])),xlab='',ylab='',col='purple',xlim=c(0,6),main=n)  
    }else{
      plot(rep(1,length(t(trait_plot[n]))),t(trait_plot[n]),xlab='',ylab='',col='purple',xlim=c(0,6),main=n)  
    }
    
    
    boxplot(add=TRUE,2,t(previous_TrBE[n]))
    }else{
      plot(rep(2,28),t(previous_TrBE[n]),xlab='', ylim= c(min(c(t(previous_TrBE[n]),t(pcr_TrBE[n]),t(lm_TrBE[n]))),max(c(t(previous_TrBE[n]),t(pcr_TrBE[n]),t(lm_TrBE[n])))),main = n,xlim=c(0,6),ylab='')
      }
 # points(rep(3,30),t(pcr_TrBE[n]), col='blue', pch=2)
  boxplot(add=TRUE,4,t(lm_TrBE[n]), col='orange', pch=5)
  boxplot(add=TRUE,(5,31),t(plsr_TrBE_new[n]), col='green', pch=2)
}
}

par(mfrow=c(3,2))
for (n in names(pcr_TrBE)){
  if(n %in% names(trait_plot)){
    if(n=='WD'){
      df <- data.frame(Obs      = c(na.omit(trait_plot[[n]])),
                       lm   = c(t(lm_TrBE[n])/1000*2, rep(NA,length(c(na.omit(trait_plot[[n]])))-length(c(t(lm_TrBE[n]))))), 
                       plsr = c(t(plsr_TrBE[n])/1000*2, rep(NA,length(c(na.omit(trait_plot[[n]])))-length(c(t(plsr_TrBE[n]))))))
                       
      boxplot(df,xlab='',ylab='',col=2:4,xlim=c(0,4),main=n) 
      points(rep(3,length(df$plsr)),df$plsr)
      points(rep(2,length(df$lm)),df$lm)
    }else if(n=='P50'){
      df <- data.frame(Obs      = c(-exp(na.omit(trait_plot[[n]]))),
                       lm   = c((t(lm_TrBE[n])), rep(NA,length(c(na.omit(trait_plot[[n]])))-length(c(t(lm_TrBE[n]))))), 
                       plsr = c((t(plsr_TrBE[n])), rep(NA,length(c(na.omit(trait_plot[[n]])))-length(c(t(plsr_TrBE[n]))))))
      
      boxplot(df,xlab='',ylab='',col=2:4,xlim=c(0,4),main=n) 
      points(rep(3,length(df$plsr)),df$plsr)
      points(rep(2,length(df$lm)),df$lm)
    }else if(n=='Ks'){
      df <- data.frame(Obs      = c((na.omit(trait_plot[[n]]))),
                       lm   = c(log((t(lm_TrBE[n]))), rep(NA,length(c(na.omit(trait_plot[[n]])))-length(c(t(lm_TrBE[n]))))), 
                       plsr = c(log((t(plsr_TrBE[n]))), rep(NA,length(c(na.omit(trait_plot[[n]])))-length(c(t(plsr_TrBE[n]))))))
      
      boxplot(df,xlab='',ylab='',col=2:4,xlim=c(0,4),main=paste0('log(',n,')')) 
      points(rep(3,length(df$plsr)),df$plsr)
      points(rep(2,length(df$lm)),df$lm) 
    }else if(n=='TLP'){
      df <- data.frame(Obs      = c((na.omit(trait_plot[[n]]))),
                       lm   = c(log(-(t(lm_TrBE[n]))), rep(NA,length(c(na.omit(trait_plot[[n]])))-length(c(t(lm_TrBE[n]))))), 
                       plsr = c(log(-(t(plsr_TrBE[n]))), rep(NA,length(c(na.omit(trait_plot[[n]])))-length(c(t(plsr_TrBE[n]))))))
      
      boxplot(df,xlab='',ylab='',col=2:4,xlim=c(0,4),main=paste0('log(-',n,')')) 
      points(rep(3,length(df$plsr)),df$plsr)
      points(rep(2,length(df$lm)),df$lm) 
    }else if (n=='LS'){
      df <- data.frame(Obs      = c((na.omit(trait_plot[[n]]))),
                       lm   = c(log((t(lm_TrBE[n])/10000)), rep(NA,length(c(na.omit(trait_plot[[n]])))-length(c(t(lm_TrBE[n]))))), 
                       plsr = c(log((t(plsr_TrBE[n])/10000)), rep(NA,length(c(na.omit(trait_plot[[n]])))-length(c(t(plsr_TrBE[n]))))))
      
      boxplot(df,xlab='',ylab='',col=2:4,xlim=c(0,4),main=paste0('log(',n,')')) # m2 (leaf) cm-2 (sap)
      points(rep(3,length(df$plsr)),df$plsr)
      points(rep(2,length(df$lm)),df$lm) 
    }else if (n=='SLA'){
      df <- data.frame(Obs      = c( (1/exp(na.omit(trait_plot[['LMA']]))*1000*2) ),
                       lm   = c( (t(lm_TrBE[n]) ), rep(NA,length(c(na.omit(trait_plot[['LMA']])))-length(c(t(lm_TrBE[n]))))), 
                       plsr = c( (t(plsr_TrBE[n]) ), rep(NA,length(c(na.omit(trait_plot[['LMA']])))-length(c(t(plsr_TrBE[n]))))))
      
      boxplot(df,xlab='',ylab='',col=2:4,xlim=c(0,4),main=paste0(n)) # m2 kgC-1
      points(rep(3,length(df$plsr)),df$plsr)
      points(rep(2,length(df$lm)),df$lm) 
    }else if(n=='slope'){
      df <- data.frame(Obs      = c(na.omit(trait_plot[[n]]) ),
                       lm   = c( log((t(lm_TrBE[n]) )), rep(NA,length(c(na.omit(trait_plot[[n]])))-length(c(t(lm_TrBE[n]))))), 
                       plsr = c( log((t(plsr_TrBE[n]) )), rep(NA,length(c(na.omit(trait_plot[[n]])))-length(c(t(plsr_TrBE[n]))))))
      
      boxplot(df,xlab='',ylab='',col=2:4,xlim=c(0,4),main=paste0('log(',n,')')) # m2 kgC-1
      points(rep(3,length(df$plsr)),df$plsr)
      points(rep(2,length(df$lm)),df$lm) 
      
    }
  }
}

plot(rep(1,length(t(trait_plot['WD']))), -0.5571 + (2.9748*(t(trait_plot['WD']))),xlab='',ylab='',col='purple',xlim=c(0,6),main="DeltaPsiWW") 
points(rep(2,28),t(previous_TrBE["DeltaPsiWW"]))
points(rep(3,30),t(pcr_TrBE["DeltaPsiWW"]), col='blue', pch=2)
points(rep(4,31),t(lm_TrBE["DeltaPsiWW"]), col='orange', pch=4)
points(rep(4,31),t(lm_TrBE["DeltaPsiWW"]), col='orange', pch=5)
points(rep(5,31),t(plsr_TrBE_new["DeltaPsiWW"]), col='green', pch=2)

plot(rep(1,length(t(trait_plot['TLP']))), -0.188+(-0.3*(-exp(t(trait_plot['TLP'])))),xlab='',ylab='',col='purple',xlim=c(0,6),main="lambda") 
points(rep(2,28),t(previous_TrBE["lambda"]))
points(rep(3,30),t(pcr_TrBE["lambda"]), col='blue', pch=2)
points(rep(4,31),t(lm_TrBE["lambda"]), col='orange', pch=5)
points(rep(5,31),t(plsr_TrBE_new["lambda"]), col='green', pch=2)


plot.new()
legend('center',legend=c('obs','old','lm','pls','pls PCoptim'),col=c('purple','black','orange','blue','green'),pch=c(1,1,2,5),cex=0.65)


