# Script to read in the processed trait data, make bivariate and multivariate SMA regressions and then
# carry out an optimisation procedure to unify inter-trait relationships across the whole trait dataset.
#
# Dependencies: sma_multivar_regress.R
#
# T. Pugh
# 15.06.20

nbtstrp=1000 # Number of bootstrap samples to take in sma_multivar_regress (samples later used to calculated uncertainty in the optimisation). Was previously 10 000, using a lower number for testing, Will need to check sensitivity to this value.

#traits=read.table("/Users/liudy/trait_data/woody_trait.0625.txt")
#traits=read.csv("/Users/pughtam/Documents/TreeMort/Analyses/Hydraulic_modelling/Traits/mytrait-data/woody_trait.0803.txt",sep="\t")
traits=read.csv("/Users/pughtam/Documents/TreeMort/Analyses/Hydraulic_modelling/Traits/woody_trait.0827.txt",sep="\t")

source('sma_multivar_regress.R')
source('trait_functions.R')
source('make_bivar_plots.R')
source('multivar_model_selection.R')

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
intercept_P50_from_TLP_Ks <- c(P50_from_TLP_Ks$mod$intercept_R,P50_from_TLP_Ks$mod$L95_R.intercept,P50_from_TLP_Ks$mod$U95_R.intercept)
y1_P50_from_TLP_Ks <- c(P50_from_TLP_Ks$mod$slope_R.y1,P50_from_TLP_Ks$mod$L95_R.y1,P50_from_TLP_Ks$mod$U95_R.y1)
y2_P50_from_TLP_Ks <- c(P50_from_TLP_Ks$mod$slope_R.y2,P50_from_TLP_Ks$mod$L95_R.y2,P50_from_TLP_Ks$mod$U95_R.y2)

coeff_P50_from_TLP_Ks <- data.frame(coeffnames_P50_from_TLP_Ks,intercept_P50_from_TLP_Ks,y1_P50_from_TLP_Ks,y2_P50_from_TLP_Ks)
View(coeff_P50_from_TLP_Ks)

# NOTE: These coefficients are all over the place, almost certainly because we have high multicolinearity in the predictors, BUT, this is not a problem as we are not interpreting the coefficients, just using them for the prediction.


# TLP fits ----------------------------------------------------------------

TLP_multivar <- TLP_multivar_test(trait_B)


# LMA fits -----------------------------------------------------------------

LMA_multivar <- LMA_multivar_test(trait_B,trait_BE,trait_BDT)


# WD fits -----------------------------------------------------------------

WD_multivar <- WD_multivar_test(trait_B)


# slope fits --------------------------------------------------------------

slope_multivar <- slope_multivar_test(trait_B)


# Optimisation ------------------------------------------------------------
# Attempt to iteratively converge on the best fit values of TLP, P50 and LMA given known Ks and LS

# Decide whether to limit the possible ranges of predicted traits to the observed values (T) or not (F)
limitdataranges=T # Currently does not converge in uncertainty propagation if not set to T

# Decide whether to run the uncertainty propagation (T) or not (F)
propagate_uncer=T

# Decide whether to run all trait combinations in the database for LS and Hmax (F), or just a selection (T)
trait_sel=T
# Number of combinations to select if trait_sel=T. Set to -1 for a systematic sample, >0 for a random sample of the size specified
n_trait_sel=-1

# Run for all Broadleaved (i.e. BE + BT + BD) (=1), or all deciduous (BT + BD) (=2), or BE (=3), or BT (=4), or BD (=5). This is used to set the maximum and minimum bounds in trait_opt().
spec_group_sel=3

# ---
if (propagate_uncer) {
  n_uncer=nbtstrp
} else {
  n_uncer=1
}

# Get index for selected species group
if (spec_group_sel==1) {
  ind_spec_group=1:length(traits$group) # Take everything (assumes that conifers were scopes out at the top of the file)
} else if (spec_group_sel==2) {
  ind_spec_group=which(traits$group=='BT' | traits$group=='BD')
} else if (spec_group_sel==3) {
  ind_spec_group=which(traits$group=='BE')
} else if (spec_group_sel==4) {
  ind_spec_group=which(traits$group=='BT')
} else if (spec_group_sel==5) {
  ind_spec_group=which(traits$group=='BD')
}

# Identify all combinations of Ks and LS (do this across full range of broadleaf species)
ind=which(!is.na(Ks) & !is.na(LS))

LS_comb <- LS[ind]
Ks_comb <- Ks[ind]

if (trait_sel) {
  if (n_trait_sel>0) {
    # Random selection of LS and Hmax values to be tested
    set.seed(1234)
    index = 1:length(ind)
    trait_samp = sample(index, n_trait_sel, replace=F) 
    LS_e = LS_comb[trait_samp] 
    Ks_e = Ks_comb[trait_samp] 
    # Plot sample against all data as a check for sampling density
    plot(LS_comb,Ks_comb)
    points(LS_e,Ks_e,col="red")
    
  } else {
    # Systematic sample
    
    library(hypervolume)
    # Fit a hypervolume (KDE at 95%)
    # Have not rescaled trait before fitting hypervolume as the ranges of both are very similar
    hv = hypervolume(data.frame(LS_comb,Ks_comb),method="gaussian",quantile.requested=0.95)
    plot(hv)
    
    # Set the number of points distributed systematically across LS and Hmax space to test for inclusion in the hypervolume
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
    
    # Test those points for inclusion in the hypervolume
    in_hv <- hypervolume_inclusion_test(hv,LS_Ks_seq,fast.or.accurate = "accurate")
    
    LS_e <- LS_Ks_seq$Var1[in_hv]
    Ks_e <- LS_Ks_seq$Var2[in_hv]
  }
} else {
  # Go through all observed combinations of Hmax and LS
  Ks_e=Ks_comb
  LS_e=LS_comb
}

# ---
# Select the LMA relationship to use
if (spec_group_sel==1) {
  LMA_from_TLP_select=LMA_from_TLP
} else if (spec_group_sel==2) {
  LMA_from_TLP_select=LMA_from_TLP_BDT
} else if (spec_group_sel==3) {
  LMA_from_TLP_select=LMA_from_TLP_BE
} else if (spec_group_sel==4) {
  LMA_from_TLP_select=LMA_from_TLP_BDT
} else if (spec_group_sel==5) {
  LMA_from_TLP_select=LMA_from_TLP_BDT
}

# ---
# Actually do the optimisation

ndata=length(LS_e)

P50_e <- matrix(NA, nrow= ndata, ncol = n_uncer) #Array now expanded to hold multiple replicate estimates based on regression coefficient uncertainty
LMA_e <- matrix(NA, nrow= ndata, ncol = n_uncer)
TLP_e <- matrix(NA, nrow= ndata, ncol = n_uncer)
WD_e <- matrix(NA, nrow= ndata, ncol = n_uncer)
slope_e <- matrix(NA, nrow= ndata, ncol = n_uncer)

# Loop over all the combinations of Hmax and LS
# The new estimates of traits use the suffix "_e"
for (dd in 1:ndata) {
  print(dd)
  
  # Carry out the optimisation
  opt_vals <- trait_opt(P50[ind_spec_group],TLP[ind_spec_group],LMA[ind_spec_group],WD[ind_spec_group],slope[ind_spec_group],LMA_from_TLP_select,TLP_from_LS_LMA_P50,P50_from_TLP_Ks,slope_from_P50_TLP_Ks,WD_from_slope_P50slope,LMA_from_LS,P50_from_Ks,TLP_from_P50,Ks_e[dd],LS_e[dd],n_uncer)
  
  P50_e[dd,] <- opt_vals$P50_e
  TLP_e[dd,] <- opt_vals$TLP_e
  LMA_e[dd,] <- opt_vals$LMA_e
  WD_e[dd,] <- opt_vals$WD_e
  slope_e[dd,] <- opt_vals$slope_e
}
  
#Stats defining the uncertainty range for each point
TLP_e_mean=unname(apply(TLP_e, 1, mean,na.rm=T))
TLP_e_median=unname(apply(TLP_e, 1, median,na.rm=T))
TLP_e_5perc=unname(apply(TLP_e, 1, quantile,0.05,na.rm=T))
TLP_e_95perc=unname(apply(TLP_e, 1, quantile,0.95,na.rm=T))

P50_e_mean=unname(apply(P50_e, 1, mean,na.rm=T))
P50_e_median=unname(apply(P50_e, 1, median,na.rm=T))
P50_e_5perc=unname(apply(P50_e, 1, quantile,0.05,na.rm=T))
P50_e_95perc=unname(apply(P50_e, 1, quantile,0.95,na.rm=T))

LMA_e_mean=unname(apply(LMA_e, 1, mean,na.rm=T))
LMA_e_median=unname(apply(LMA_e, 1, median,na.rm=T))
LMA_e_5perc=unname(apply(LMA_e, 1, quantile,0.05,na.rm=T))
LMA_e_95perc=unname(apply(LMA_e, 1, quantile,0.95,na.rm=T))

WD_e_mean=unname(apply(WD_e, 1, mean,na.rm=T))
WD_e_median=unname(apply(WD_e, 1, median,na.rm=T))
WD_e_5perc=unname(apply(WD_e, 1, quantile,0.05,na.rm=T))
WD_e_95perc=unname(apply(WD_e, 1, quantile,0.95,na.rm=T))

slope_e_mean=unname(apply(slope_e, 1, mean,na.rm=T))
slope_e_median=unname(apply(slope_e, 1, median,na.rm=T))
slope_e_5perc=unname(apply(slope_e, 1, quantile,0.05,na.rm=T))
slope_e_95perc=unname(apply(slope_e, 1, quantile,0.95,na.rm=T))


#Make plots to compare with original data
par(mfrow=c(4,4))
par(mar=c(2,2,2,2))

plot(TLP,P50,pch=16,xlab="TLP",ylab="P50",main="TLP vs P50")
points(TLP_e[,1],P50_e[,1],col="blue",pch=16) # Using central estimate coefficients
points(TLP_e[,1],P50_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
points(TLP_e[,1],P50_e_5perc,col="green",pch=16)
points(TLP_e[,1],P50_e_95perc,col="green",pch=16)

plot(P50,TLP,pch=16,xlab="P50",ylab="TLP",main="P50 vs TLP")
points(P50_e[,1],TLP_e[,1],col="blue",pch=16) # Using central estimate coefficients
points(P50_e[,1],TLP_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
points(P50_e[,1],TLP_e_5perc,col="green",pch=16)
points(P50_e[,1],TLP_e_95perc,col="green",pch=16)

plot(TLP,slope,pch=16,xlab="TLP",ylab="slope",main="TLP vs slope")
points(TLP_e[,1],slope_e[,1],col="blue",pch=16) # Using central estimate coefficients
points(TLP_e[,1],slope_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
points(TLP_e[,1],slope_e_5perc,col="green",pch=16)
points(TLP_e[,1],slope_e_95perc,col="green",pch=16)

plot(slope,TLP,pch=16,xlab="slope",ylab="TLP",main="slope vs TLP")
points(slope_e[,1],TLP_e[,1],col="blue",pch=16) # Using central estimate coefficients
points(slope_e[,1],TLP_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
points(slope_e[,1],TLP_e_5perc,col="green",pch=16)
points(slope_e[,1],TLP_e_95perc,col="green",pch=16)

plot(P50,slope,pch=16,xlab="P50",ylab="slope",main="P50 vs slope")
points(P50_e[,1],slope_e[,1],col="blue",pch=16) # Using central estimate coefficients
points(P50_e[,1],slope_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
points(P50_e[,1],slope_e_5perc,col="green",pch=16)
points(P50_e[,1],slope_e_95perc,col="green",pch=16)

plot(slope,P50,pch=16,xlab="slope",ylab="P50",main="slope vs P50")
points(slope_e[,1],P50_e[,1],col="blue",pch=16) # Using central estimate coefficients
points(slope_e[,1],P50_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
points(slope_e[,1],P50_e_5perc,col="green",pch=16)
points(slope_e[,1],P50_e_95perc,col="green",pch=16)

plot(TLP,WD,pch=16,xlab="TLP",ylab="WD",main="TLP vs WD")
points(TLP_e[,1],WD_e[,1],col="blue",pch=16) # Using central estimate coefficients
points(TLP_e[,1],WD_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
points(TLP_e[,1],WD_e_5perc,col="green",pch=16)
points(TLP_e[,1],WD_e_95perc,col="green",pch=16)

plot(WD,TLP,pch=16,xlab="WD",ylab="TLP",main="WD vs TLP")
points(WD_e[,1],TLP_e[,1],col="blue",pch=16) # Using central estimate coefficients
points(WD_e[,1],TLP_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
points(WD_e[,1],TLP_e_5perc,col="green",pch=16)
points(WD_e[,1],TLP_e_95perc,col="green",pch=16)

plot(P50,WD,pch=16,xlab="P50",ylab="WD",main="P50 vs WD")
points(P50_e[,1],WD_e[,1],col="blue",pch=16) # Using central estimate coefficients
points(P50_e[,1],WD_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
points(P50_e[,1],WD_e_5perc,col="green",pch=16)
points(P50_e[,1],WD_e_95perc,col="green",pch=16)

plot(WD,P50,pch=16,xlab="WD",ylab="P50",main="WD vs P50")
points(WD_e[,1],P50_e[,1],col="blue",pch=16) # Using central estimate coefficients
points(WD_e[,1],P50_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
points(WD_e[,1],P50_e_5perc,col="green",pch=16)
points(WD_e[,1],P50_e_95perc,col="green",pch=16)

#NOTE: From here onwards I have not made the plots in both directions
plot(TLP,LMA,pch=16,xlab="TLP",ylab="LMA",main="TLP vs LMA")
points(TLP_e[,1],LMA_e[,1],col="blue",pch=16) # Using central estimate coefficients
points(TLP_e[,1],LMA_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
points(TLP_e[,1],LMA_e_5perc,col="green",pch=16)
points(TLP_e[,1],LMA_e_95perc,col="green",pch=16)

plot(LS,LMA,pch=16,xlab="LS",ylab="LMA",main="LS vs LMA")
points(LS_e,LMA_e[,1],col="blue",pch=16) # Using central estimate coefficients
points(LS_e,LMA_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
points(LS_e,LMA_e_5perc,col="green",pch=16)
points(LS_e,LMA_e_95perc,col="green",pch=16)

plot(Ks,P50,pch=16,xlab="Ks",ylab="P50",main="Ks vs P50")
points(Ks_e,P50_e[,1],col="blue",pch=16) # Using central estimate coefficients
points(Ks_e,P50_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
points(Ks_e,P50_e_5perc,col="green",pch=16)
points(Ks_e,P50_e_95perc,col="green",pch=16)

plot(LS,TLP,pch=16,xlab="LS",ylab="TLP",main="LS vs TLP")
points(LS_e,TLP_e[,1],col="blue",pch=16) # Using central estimate coefficients
points(LS_e,TLP_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
points(LS_e,TLP_e_5perc,col="green",pch=16)
points(LS_e,TLP_e_95perc,col="green",pch=16)

plot(WD,LMA,pch=16,xlab="WD",ylab="LMA",main="WD vs LMA")
points(WD_e[,1],LMA_e[,1],col="blue",pch=16) # Using central estimate coefficients
points(WD_e[,1],LMA_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
points(WD_e[,1],LMA_e_5perc,col="green",pch=16)
points(WD_e[,1],LMA_e_95perc,col="green",pch=16)

plot(Ks,slope,pch=16,xlab="Ks",ylab="slope",main="Ks vs slope")
points(Ks_e,slope_e[,1],col="blue",pch=16) # Using central estimate coefficients
points(Ks_e,slope_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
points(Ks_e,slope_e_5perc,col="green",pch=16)
points(Ks_e,slope_e_95perc,col="green",pch=16)

#Use up remaining unallocated plots and set back to single plot
par(mfrow=c(1,1))
par(mar=c(5.1,4.1,4.1,2.1))


# Calculate the RMSE

#P50_res <- P50[ind]-P50_e_mean
P50_res <- P50[ind]-P50_e[,1]
P50_e_rmse <- sqrt(mean(P50_res^2,na.rm=T))
P50_e_ndata <- length(which(is.na(P50_res)==F))

#TLP_res <- TLP[ind]-TLP_e_mean
TLP_res <- TLP[ind]-TLP_e[,1]
TLP_e_rmse <- sqrt(mean(TLP_res^2,na.rm=T))
TLP_e_ndata <- length(which(is.na(TLP_res)==F))

#LMA_res <- LMA[ind]-LMA_e_mean
LMA_res <- LMA[ind]-LMA_e[,1]
LMA_e_rmse <- sqrt(mean(LMA_res^2,na.rm=T))
LMA_e_ndata <- length(which(is.na(LMA_res)==F))

#WD_res <- WD[ind]-WD_e_mean
WD_res <- WD[ind]-WD_e[,1]
WD_e_rmse <- sqrt(mean(WD_res^2,na.rm=T))
WD_e_ndata <- length(which(is.na(WD_res)==F))

#slope_res <- slope[ind]-slope_e_mean
slope_res <- slope[ind]-slope_e[,1]
slope_e_rmse <- sqrt(mean(slope_res^2,na.rm=T))
slope_e_ndata <- length(which(is.na(slope_res)==F))

# Compare to RMSE from best multivariate regression

all_e_rmse <- c(P50_e_rmse,TLP_e_rmse,LMA_e_rmse,WD_e_rmse,slope_e_rmse)
all_e_ndata <- c(P50_e_ndata,TLP_e_ndata,LMA_e_ndata,WD_e_ndata,slope_e_ndata)
all_multivar_rmse <- c(P50_from_TLP_Ks$rmse,TLP_from_LS_LMA_P50$rmse,LMA_from_TLP$rmse,WD_from_slope_P50slope$rmse,slope_from_P50_TLP_Ks$rmse)
all_multivar_ndata <- c(P50_from_TLP_Ks$ndata,TLP_from_LS_LMA_P50$ndata,LMA_from_TLP$ndata,WD_from_slope_P50slope$ndata,slope_from_P50_TLP_Ks$ndata)

all_rmse_comp <- data.frame(all_e_rmse,all_e_ndata,all_multivar_rmse,all_multivar_ndata)
View(all_rmse_comp)


# Write the optimised trait values out to a file
traits_e_out <- data.frame(LS_e,Ks_e,
                           P50_e_mean,P50_e_5perc,P50_e_95perc,
                           TLP_e_mean,TLP_e_5perc,TLP_e_95perc,
                           LMA_e_mean,LMA_e_5perc,LMA_e_95perc,
                           WD_e_mean,WD_e_5perc,WD_e_95perc,
                           slope_e_mean,slope_e_5perc,slope_e_95perc)
write.table(format(traits_e_out, digits=3), "traits_e_out_systtraits_260820.csv", append = FALSE, sep = ",", dec = ".",row.names = F, col.names = T)



# Calculate regression of leafL from LMA ----------------------------------

leafL_from_LMA <- sma_plot_stats(data.frame(LMA,log(leafL)),c("LMA","leafL"),nbtstrp,T)


# Calculate limits of leafN vs LMA to allow estimate of leaf C:N ----------

leafN_from_LMA <- sma_plot_stats(data.frame(LMA,leafN),c("LMA","leafN"),nbtstrp,T)

leafN_from_LMA_limit <- regress_limit_adjust(leafN,LMA,leafN_from_LMA,0.05)

plot(LMA,leafN)
points(LMA[leafN_from_LMA_limit$ind],leafN_from_LMA_limit$var1_pred_lower,col="green")
points(LMA[leafN_from_LMA_limit$ind],leafN_from_LMA_limit$var1_pred_upper,col="red")


# Convert to the values needed in LPJ-GUESS and write out -----------------

source('lpjg_traits_conv.R')

traits_LPJG <- lpjg_traits_conv(LMA_e_mean,P50_e_mean,TLP_e_mean,slope_e_mean,
                                LS_e,WD_e_mean,Ks_e,
                                leafL_from_LMA,leafN_from_LMA,leafN_from_LMA_limit)

# Select which base PFT to use: TeBE (1), TeBS (2), IBS (3), TrBE (4) or TrBR (5)
basePFT=5

# Select output folder
output_fol="/Users/pughtam/Documents/TreeMort/Analyses/Hydraulic_modelling/Traits/uncer_test_KsLS/TrBR_with_BE_traits_test"

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


  