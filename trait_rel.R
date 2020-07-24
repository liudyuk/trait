# Script to read in the processed trait data, make bivariate and multivariate SMA regressions and then
# carry out an optimisation procedure to unify inter-trait relationships across the whole trait dataset.
#
# Dependencies: sma_multivar_regress.R
#
# T. Pugh
# 15.06.20

nbtstrp=1000 # Number of bootstrap samples to take in sma_multivar_regress (samples later used to calculated uncertainty in the optimisation). Was previously 10 000, using a lower number for testing, Will need to check sensitivity to this value.

#traits=read.table("/Users/liudy/trait_data/woody_trait.0625.txt")
traits=read.table("/Users/pughtam/Documents/TreeMort/Analyses/Hydraulic_modelling/Traits/mytrait-data/woody_trait.0625.txt")

source('sma_multivar_regress.R')
source('trait_functions.R')

#--- Read in the trait data ---
#traits=read.table(traitfile)
#attach(traits)
#detach(traits)

traitb<-subset(traits,group!="CC",drop = T)
#this is to delete the level that should not exist after subsetting
traitb<-droplevels(traitb)
str(traitb)
attach(traitb)

#--- Bivariate plots with SMA regression ---

par(mfrow=c(4,4))
par(mar=c(2,2,2,2))

TLP_from_P50 <- sma_plot_stats(data.frame(P50,TLP),c("P50","TLP"),nbtstrp,T)

TLP_from_slope <- sma_plot_stats(data.frame(slope,TLP),c("slope","TLP"),nbtstrp,T)

P50_from_slope <- sma_plot_stats(data.frame(slope,P50),c("slope","P50"),nbtstrp,T)

TLP_from_WD <- sma_plot_stats(data.frame(WD,TLP),c("WD","TLP"),nbtstrp,T)

P50_from_WD <- sma_plot_stats(data.frame(WD,P50),c("WD","P50"),nbtstrp,T)

TLP_from_LMA <- sma_plot_stats(data.frame(LMA,TLP),c("LMA","TLP"),nbtstrp,T)

LS_from_LMA <- sma_plot_stats(data.frame(LMA,LS),c("LMA","LS"),nbtstrp,T)

Ks_from_LSHmax <- sma_plot_stats(data.frame(LS_Hmax,Ks),c("LS*Hmax","Ks"),nbtstrp,T)

Ks_from_P50 <- sma_plot_stats(data.frame(P50,Ks),c("P50","Ks"),nbtstrp,T)

LS_from_TLP <- sma_plot_stats(data.frame(TLP,LS),c("TLP","LS"),nbtstrp,T)

WD_from_LMA <- sma_plot_stats(data.frame(LMA,WD),c("LMA","WD"),nbtstrp,T)

Ks_from_slope <- sma_plot_stats(data.frame(slope,Ks),c("slope","Ks"),nbtstrp,T)

Ks_from_Hmax <- sma_plot_stats(data.frame(Hmax,Ks),c("Hmax","Ks"),nbtstrp,T)

leafN_from_LMA <- sma_plot_stats(data.frame(LMA,leafN),c("LMA","leafN"),nbtstrp,T)

LS_from_KsHmax <- sma_plot_stats(data.frame(Ks_Hmax,LS),c("Ks/Hmax","LS"),nbtstrp,T)

Ks_from_LMA <- sma_plot_stats(data.frame(LMA,Ks),c("LMA","Ks"),nbtstrp,T)

# Make a data frame summarising the fits of the regressions
all_label1 <- c("TLP","slope","slope","WD","WD","LMA","LMA","LS*Hmax","P50","TLP","LMA","slope","Hmax","LMA","Ks/Hmax","LMA")
all_label2 <- c("P50","TLP","P50","TLP","P50","TLP","LS","Ks","Ks","LS","WD","Ks","Ks","leafN","LS","Ks")
all_R2 <- c(sma_TLP_P50$R2,sma_TLP_slope$R2,sma_P50_slope$R2,sma_TLP_WD$R2,sma_P50_WD$R2,
               sma_TLP_LMA$R2,sma_LS_LMA$R2,sma_Ks_LSHmax$R2,sma_Ks_P50$R2,sma_LS_TLP$R2,
               sma_WD_LMA$R2,sma_Ks_slope$R2,sma_Ks_Hmax$R2,sma_leafN_LMA$R2,sma_LS_KsHmax$R2,
               sma_Ks_LMA$R2)
all_R2adj <- c(sma_TLP_P50$R2adj,sma_TLP_slope$R2adj,sma_P50_slope$R2adj,sma_TLP_WD$R2adj,sma_P50_WD$R2adj,
            sma_TLP_LMA$R2adj,sma_LS_LMA$R2adj,sma_Ks_LSHmax$R2adj,sma_Ks_P50$R2adj,sma_LS_TLP$R2adj,
            sma_WD_LMA$R2adj,sma_Ks_slope$R2adj,sma_Ks_Hmax$R2adj,sma_leafN_LMA$R2adj,sma_LS_KsHmax$R2adj,
            sma_Ks_LMA$R2adj)
all_rmse <- c(sma_TLP_P50$rmse,sma_TLP_slope$rmse,sma_P50_slope$rmse,sma_TLP_WD$rmse,sma_P50_WD$rmse,
            sma_TLP_LMA$rmse,sma_LS_LMA$rmse,sma_Ks_LSHmax$rmse,sma_Ks_P50$rmse,sma_LS_TLP$rmse,
            sma_WD_LMA$rmse,sma_Ks_slope$rmse,sma_Ks_Hmax$rmse,sma_leafN_LMA$rmse,sma_LS_KsHmax$rmse,
            sma_Ks_LMA$rmse)

all_sma_bivar <- data.frame(all_label1,all_label2,all_R2,all_R2adj,all_rmse)

#Use up remaining unallocated plots and set back to single plot
#plot.new()
#plot.new()
#plot.new()
par(mfrow=c(1,1))
par(mar=c(5.1,4.1,4.1,2.1))

#--- Experiment with different plausible multivariate SMA models, based on our theory ---


# P50 fits ----------------------------------------------------------------

# P50 from TLP, Ks and WD
P50_from_TLP_Ks_WD <- sma_plot_stats(data.frame(TLP,Ks,WD,P50),c("TLP","Ks","WD","P50"),nbtstrp)
plot(P50[P50_from_TLP_Ks_WD$dataused],P50_from_TLP_Ks_WD$var_est,pch=16,xlab="P50",ylab="P50_est",main="P50 vs P50_est")
#NOTE: This actually reduces R2 compared to the next model which has one variable less. That is strange behaviour (and cannot be attributed to the reduction in species number, as testing that does not explain the difference)

# P50 from TLP and Ks
P50_from_TLP_Ks <- sma_plot_stats(data.frame(TLP,Ks,P50),c("TLP","Ks","P50"),nbtstrp)
plot(P50[P50_from_TLP_Ks$dataused],P50_from_TLP_Ks$var_est,pch=16,xlab="P50",ylab="P50_est",main="P50 vs P50_est")

# P50 from TLP
P50_from_TLP <- sma_plot_stats(data.frame(TLP,P50),c("TLP","P50"),nbtstrp)

# P50 from Ks
P50_from_Ks <- sma_plot_stats(data.frame(Ks,P50),c("Ks","P50"),nbtstrp)

# P50 from TLP (same species as for TLP and Ks)
P50_from_TLP_limitspec <- sma_plot_stats(data.frame(TLP,P50),c("TLP","P50"),nbtstrp,F,P50_from_TLP_Ks$dataused)

# P50 from Ks (same species as for TLP and Ks)
P50_from_Ks_limitspec <- sma_plot_stats(data.frame(Ks,P50),c("Ks","P50"),nbtstrp,F,P50_from_TLP_Ks$dataused)

# Summarise statistics
all_testnames_P50 <- c("P50_from_TLP_Ks_WD","P50_from_TLP_Ks","P50_from_TLP","P50_from_Ks","P50_from_TLP_limitspec","P50_from_Ks_limitspec")
all_R2_P50 <- c(P50_from_TLP_Ks_WD$R2,P50_from_TLP_Ks$R2,P50_from_TLP$R2,P50_from_Ks$R2,P50_from_TLP_limitspec$R2,P50_from_Ks_limitspec$R2)
all_R2adj_P50 <- c(P50_from_TLP_Ks_WD$R2adj,P50_from_TLP_Ks$R2adj,P50_from_TLP$R2adj,P50_from_Ks$R2adj,P50_from_TLP_limitspec$R2adj,P50_from_Ks_limitspec$R2adj)
all_rmse_P50 <- c(P50_from_TLP_Ks_WD$rmse,P50_from_TLP_Ks$rmse,P50_from_TLP$rmse,P50_from_Ks$rmse,P50_from_TLP_limitspec$rmse,P50_from_Ks_limitspec$rmse)
all_ndata_P50 <- c(P50_from_TLP_Ks_WD$ndata,P50_from_TLP_Ks$ndata,P50_from_TLP$ndata,P50_from_Ks$ndata,P50_from_TLP_limitspec$ndata,P50_from_Ks_limitspec$ndata)

all_P50 <- data.frame(all_testnames_P50,all_R2_P50,all_R2adj_P50,all_rmse_P50,all_ndata_P50)
View(all_P50)

# BEST MODEL: P50_from_TLP_Ks
# Test MAT and PPT coverage of species for best model
plot(MAT[P50_from_TLP_Ks$dataused],MAP[P50_from_TLP_Ks$dataused])
# WIDE CLIMATE COVERAGE

# DECISION: P50_from_TLP_Ks


# TLP fits ----------------------------------------------------------------

# TLP from LS, LMA, P50 and WD
TLP_from_LS_LMA_P50_WD <- sma_plot_stats(data.frame(LS,LMA,P50,WD,TLP),c("LS","LMA","P50","WD","TLP"),nbtstrp)
plot(TLP[TLP_from_LS_LMA_P50_WD$dataused],TLP_from_LS_LMA_P50_WD$var_est,pch=16,xlab="TLP",ylab="TLP_est",main="TLP vs TLP_est")

# TLP from LS, LMA, P50 and slope
TLP_from_LS_LMA_P50_slope <- sma_plot_stats(data.frame(LS,LMA,P50,slope,TLP),c("LS","LMA","P50","slope","TLP"),nbtstrp)
plot(TLP[TLP_from_LS_LMA_P50_slope$dataused],TLP_from_LS_LMA_P50_slope$var_est,pch=16,xlab="TLP",ylab="TLP_est",main="TLP vs TLP_est")

# TLP from LS, LMA and P50
TLP_from_LS_LMA_P50 <- sma_plot_stats(data.frame(LS,LMA,P50,TLP),c("LS","LMA","P50","TLP"),nbtstrp)
plot(TLP[TLP_from_LS_LMA_P50$dataused],TLP_from_LS_LMA_P50$var_est,pch=16,xlab="TLP",ylab="TLP_est",main="TLP vs TLP_est")

# TLP from LS and LMA
TLP_from_LS_LMA <- sma_plot_stats(data.frame(LS,LMA,TLP),c("LS","LMA","TLP"),nbtstrp)
plot(TLP[TLP_from_LS_LMA$dataused],TLP_from_LS_LMA$var_est,pch=16,xlab="TLP",ylab="TLP_est",main="TLP vs TLP_est")

# TLP from P50 and LMA
TLP_from_P50_LMA <- sma_plot_stats(data.frame(P50,LMA,TLP),c("P50","LMA","TLP"),nbtstrp)
plot(TLP[TLP_from_P50_LMA$dataused],TLP_from_P50_LMA$var_est,pch=16,xlab="TLP",ylab="TLP_est",main="TLP vs TLP_est")

# TLP from P50 and LMA (same species as for LS, LMA and P50)
TLP_from_P50_LMA_limitspec <- sma_plot_stats(data.frame(P50,LMA,TLP),c("P50","LMA","TLP"),nbtstrp,F,TLP_from_LS_LMA_P50$dataused)

# TLP from P50 and LS
TLP_from_P50_LS <- sma_plot_stats(data.frame(P50,LS,TLP),c("P50","LS","TLP"),nbtstrp)
plot(TLP[TLP_from_P50_LS$dataused],TLP_from_P50_LS$var_est,pch=16,xlab="TLP",ylab="TLP_est",main="TLP vs TLP_est")

# TLP from LS
TLP_from_LS <- sma_plot_stats(data.frame(LS,TLP),c("LS","TLP"),nbtstrp)

# TLP from P50 (same species as for LS, LMA and P50)
TLP_from_P50_limitspec <- sma_plot_stats(data.frame(P50,TLP),c("P50","TLP"),nbtstrp,F,TLP_from_LS_LMA_P50$dataused)

# Summarise statistics
all_testnames_TLP <- c("TLP_from_LS_LMA_P50_WD","TLP_from_LS_LMA_P50_slope","TLP_from_LS_LMA_P50","TLP_from_LS_LMA","TLP_from_P50_LMA","TLP_from_P50_LMA_limitspec","TLP_from_P50_LS","TLP_from_LS","TLP_from_P50","TLP_from_LMA","TLP_from_P50_limitspec")
all_R2_TLP <- c(TLP_from_LS_LMA_P50_WD$R2,TLP_from_LS_LMA_P50_slope$R2,TLP_from_LS_LMA_P50$R2,TLP_from_LS_LMA$R2,TLP_from_P50_LMA$R2,TLP_from_P50_LMA_limitspec$R2,TLP_from_P50_LS$R2,TLP_from_LS$R2,TLP_from_P50$R2,TLP_from_LMA$R2,TLP_from_P50_limitspec$R2)
all_R2adj_TLP <- c(TLP_from_LS_LMA_P50_WD$R2adj,TLP_from_LS_LMA_P50_slope$R2adj,TLP_from_LS_LMA_P50$R2adj,TLP_from_LS_LMA$R2adj,TLP_from_P50_LMA$R2adj,TLP_from_P50_LMA_limitspec$R2adj,TLP_from_P50_LS$R2adj,TLP_from_LS$R2adj,TLP_from_P50$R2adj,TLP_from_LMA$R2adj,TLP_from_P50_limitspec$R2adj)
all_rmse_TLP <- c(TLP_from_LS_LMA_P50_WD$rmse,TLP_from_LS_LMA_P50_slope$rmse,TLP_from_LS_LMA_P50$rmse,TLP_from_LS_LMA$rmse,TLP_from_P50_LMA$rmse,TLP_from_P50_LMA_limitspec$rmse,TLP_from_P50_LS$rmse,TLP_from_LS$rmse,TLP_from_P50$rmse,TLP_from_LMA$rmse,TLP_from_P50_limitspec$rmse)
all_ndata_TLP <- c(TLP_from_LS_LMA_P50_WD$ndata,TLP_from_LS_LMA_P50_slope$ndata,TLP_from_LS_LMA_P50$ndata,TLP_from_LS_LMA$ndata,TLP_from_P50_LMA$ndata,TLP_from_P50_LMA_limitspec$ndata,TLP_from_P50_LS$ndata,TLP_from_LS$ndata,TLP_from_P50$ndata,TLP_from_LMA$ndata,TLP_from_P50_limitspec$ndata)

all_TLP <- data.frame(all_testnames_TLP,all_R2_TLP,all_R2adj_TLP,all_rmse_TLP,all_ndata_TLP)
View(all_TLP)

# CHOICE: Although WD improves the fit, do not use it as our hypothesis framework does not posit a direct link between TLP and WD (only indirect via P50 and slope)

# BEST MODEL: TLP_from_LS_LMA_P50
# NO EVIDENCE that the smaller selection of species for 3 variables is the reason behind the fit
# Test MAT and PPT coverage of species for best model
plot(MAT[TLP_from_LS_LMA_P50$dataused],MAP[TLP_from_LS_LMA_P50$dataused])
# Test MAT and PPT coverage of species from TLP_from_P50_LMA
plot(MAT[TLP_from_P50_LMA$dataused],MAP[TLP_from_P50_LMA$dataused])
# Test MAT and PPT coverage of species from TLP_from_LMA
plot(MAT[TLP_from_LMA$dataused],MAP[TLP_from_LMA$dataused])
# WIDE CLIMATE COVERAGE for TLP_from_LS_LMA_P50 EXCEPT temperate rainforest (which is captured by TLP_from_P50_LMA)

# Test if relationships are consistent in character despite regardless of climate zone differences
plot(TLP[TLP_from_LS_LMA_P50$dataused],TLP_from_LS_LMA_P50$var_est,pch=16,xlab="TLP",ylab="TLP_est",main="TLP vs TLP_est")
points(TLP[TLP_from_P50_LMA$dataused],TLP_from_P50_LMA$var_est,pch=16,col="red")

# DECISION: TLP_from_LS_LMA_P50


# LMA fits -----------------------------------------------------------------

# LMA from TLP, LS and WD
LMA_from_TLP_LS_WD <- sma_plot_stats(data.frame(TLP,LS,WD,LMA),c("TLP","LS","WD","LMA"),nbtstrp)
plot(LMA[LMA_from_TLP_LS_WD$dataused],LMA_from_TLP_LS_WD$var_est,pch=16,xlab="LMA",ylab="LMA_est",main="LMA vs LMA_est")

# LMA from TLP and LS
LMA_from_TLP_LS <- sma_plot_stats(data.frame(TLP,LS,LMA),c("TLP","LS","LMA"),nbtstrp)
plot(LMA[LMA_from_TLP_LS$dataused],LMA_from_TLP_LS$var_est,pch=16,xlab="LMA",ylab="LMA_est",main="LMA vs LMA_est")

# LMA from TLP and WD
LMA_from_TLP_WD <- sma_plot_stats(data.frame(TLP,WD,LMA),c("TLP","WD","LMA"),nbtstrp)
plot(LMA[LMA_from_TLP_WD$dataused],LMA_from_TLP_WD$var_est,pch=16,xlab="LMA",ylab="LMA_est",main="LMA vs LMA_est")

# LMA from LS and WD
LMA_from_LS_WD <- sma_plot_stats(data.frame(LS,WD,LMA),c("LS","WD","LMA"),nbtstrp)
plot(LMA[LMA_from_LS_WD$dataused],LMA_from_LS_WD$var_est,pch=16,xlab="LMA",ylab="LMA_est",main="LMA vs LMA_est")

# LMA from TLP
LMA_from_TLP <- sma_plot_stats(data.frame(TLP,LMA),c("TLP","LMA"),nbtstrp)
plot(LMA[LMA_from_TLP$dataused],LMA_from_TLP$var_est,pch=16,xlab="LMA",ylab="LMA_est",main="LMA vs LMA_est")

# LMA from LS
LMA_from_LS <- sma_plot_stats(data.frame(LS,LMA),c("LS","LMA"),nbtstrp)

# LMA from WD
LMA_from_WD <- sma_plot_stats(data.frame(WD,LMA),c("WD","LMA"),nbtstrp)

# Summarise statistics
all_testnames_LMA <- c("LMA_from_TLP_LS_WD","LMA_from_TLP_LS","LMA_from_TLP_WD","LMA_from_LS_WD","LMA_from_TLP","LMA_from_LS","LMA_from_WD")
all_R2_LMA <- c(LMA_from_TLP_LS_WD$R2,LMA_from_TLP_LS$R2,LMA_from_TLP_WD$R2,LMA_from_LS_WD$R2,LMA_from_TLP$R2,LMA_from_LS$R2,LMA_from_WD$R2)
all_R2adj_LMA <- c(LMA_from_TLP_LS_WD$R2adj,LMA_from_TLP_LS$R2adj,LMA_from_TLP_WD$R2adj,LMA_from_LS_WD$R2adj,LMA_from_TLP$R2adj,LMA_from_LS$R2adj,LMA_from_WD$R2adj)
all_rmse_LMA <- c(LMA_from_TLP_LS_WD$rmse,LMA_from_TLP_LS$rmse,LMA_from_TLP_WD$rmse,LMA_from_LS_WD$rmse,LMA_from_TLP$rmse,LMA_from_LS$rmse,LMA_from_WD$rmse)
all_ndata_LMA <- c(LMA_from_TLP_LS_WD$ndata,LMA_from_TLP_LS$ndata,LMA_from_TLP_WD$ndata,LMA_from_LS_WD$ndata,LMA_from_TLP$ndata,LMA_from_LS$ndata,LMA_from_WD$ndata)

all_LMA <- data.frame(all_testnames_LMA,all_R2_LMA,all_R2adj_LMA,all_rmse_LMA,all_ndata_LMA)
View(all_LMA)

# CHOICE: best model in R2 terms is LMA_from_LS_WD, but in RMSE is (marginally) LMA_from_TLP
# Prefer to go with LMA_from_TLP on the basis that it better fits our hypothesis framework, there is also a clearer conceptual link between TLP and LMA than LS and LMA (and WD link is only expected for evergreen species)

# Test MAT and PPT coverage of species for chosen model
plot(MAT[LMA_from_TLP$dataused],MAP[LMA_from_TLP$dataused])
# WIDE CLIMATE COVERAGE

# DECISION: LMA_from_TLP


#2)Estimate LMA from LS and TLP
#note: LS is not good and consider to not use
ind_LMA=which(!is.na(LMA) & !is.na(LS) & !is.na(TLP))
mod_LMA <- sma_regress_multivar(data.frame(TLP[ind_LMA],LS[ind_LMA],LMA[ind_LMA]),nbtstrp)
LMA_est <- mod_LMA$intercept_R + mod_LMA$slope_R.y1*TLP[ind_LMA] + mod_LMA$slope_R.y2*LS[ind_LMA]
#Cross-checks to make sure that result is sensible
plot(mod_LS_LMA,"SMA",pch=16,xlab="LS",ylab="LMA",main="LS vs LMA")
points(LS[ind_LMA],LMA_est,col="red",pch=16)
plot(mod_TLP_LMA,"SMA",pch=16,xlab="TLP",ylab="LMA",main="TLP vs LMA")
points(TLP[ind_LMA],LMA_est,col="red",pch=16)
#here TLP is a goood predictor for LMA_est, but not LS.

#3)Estimate WD from P50 and LMA
ind_WD=which(!is.na(WD) & !is.na(P50) & !is.na(LMA))
mod_WD <- sma_regress_multivar(data.frame(P50[ind_WD],LMA[ind_WD],P50[ind_WD]),nbtstrp)
WD_est <- mod_WD$intercept_R + mod_WD$slope_R.y1*P50[ind_WD] + mod_WD$slope_R.y2*LMA[ind_WD]
#Cross-checks to make sure that result is sensible
plot(mod_P50_WD,"SMA",pch=16,xlab="P50",ylab="WD",main="P50 vs WD")
points(P50[ind_WD],WD_est,col="red",pch=16)
plot(mod_WD_LMA,"SMA",pch=16,xlab="WD",ylab="LMA",main="WD vs LMA")
points(WD_est,LMA[ind_WD],col="red",pch=16)
#why the regression for P50 and WD is so good? sel-correlated?
#test P50 TLP LMA
ind_WD=which(!is.na(WD) & !is.na(P50) & !is.na(LMA) & !is.na(TLP))
mod_WD <- sma_regress_multivar(data.frame(P50[ind_WD],LMA[ind_WD],TLP[ind_WD],WD[ind_WD]),nbtstrp)
WD_est <- mod_WD$intercept_R + mod_WD$slope_R.y1*P50[ind_WD] + mod_WD$slope_R.y2*LMA[ind_WD]+ mod_WD$slope_R.y3*TLP[ind_WD]
#Cross-checks to make sure that result is sensible
plot(mod_P50_WD,"SMA",pch=16,xlab="P50",ylab="WD",main="P50 vs WD")
points(P50[ind_WD],WD_est,col="red",pch=16)
plot(mod_WD_LMA,"SMA",pch=16,xlab="WD",ylab="LMA",main="WD vs LMA")
points(WD_est,LMA[ind_WD],col="red",pch=16)
plot(mod_TLP_WD,"SMA",pch=16,xlab="TLP",ylab="WD",main="TLP vs WD")
points(WD_est,TLP[ind_WD],col="red",pch=16)
#test P50 TLP 
ind_WD=which(!is.na(WD) & !is.na(P50)  & !is.na(TLP))
mod_WD <- sma_regress_multivar(data.frame(P50[ind_WD],TLP[ind_WD],P50[ind_WD]),nbtstrp)
WD_est <- mod_WD$intercept_R + mod_WD$slope_R.y1*P50[ind_WD] + mod_WD$slope_R.y2*TLP[ind_WD]
#Cross-checks to make sure that result is sensible
plot(mod_P50_WD,"SMA",pch=16,xlab="P50",ylab="WD",main="P50 vs WD")
points(P50[ind_WD],WD_est,col="red",pch=16)
plot(mod_TLP_WD,"SMA",pch=16,xlab="TLP",ylab="WD",main="TLP vs WD")
points(WD_est,TLP[ind_WD],col="red",pch=16)

#test TLP LMA
ind_WD=which(!is.na(WD) & !is.na(LMA) & !is.na(TLP))
mod_WD <- sma_regress_multivar(data.frame(LMA[ind_WD],TLP[ind_WD],WD[ind_WD]),nbtstrp)
WD_est <- mod_WD$intercept_R + mod_WD$slope_R.y1*LMA[ind_WD]+ mod_WD$slope_R.y2*TLP[ind_WD]
#Cross-checks to make sure that result is sensible
plot(mod_TLP_WD,"SMA",pch=16,xlab="TLP",ylab="WD",main="TLP vs WD")
points(TLP[ind_WD],WD_est,col="red",pch=16)
plot(mod_WD_LMA,"SMA",pch=16,xlab="WD",ylab="LMA",main="WD vs LMA")
points(WD_est,LMA[ind_WD],col="red",pch=16)

library(smatr) #NOTE TO DAIJUN: I can't use this package (it's for an older version of R). Why not use the lmodel2 package as at the start of this script?
#test LMA
ind_WD=which(!is.na(WD) & !is.na(LMA))
mod_WD <- sma(data.frame(WD[ind_WD],LMA[ind_WD]),method=c("SMA"))
summary(mod_WD)
WD_est <- mod_WD$coef[[1]][1,1] + mod_WD$coef[[1]][2,1]*LMA[ind_WD]
plot(mod_WD_LMA,"SMA",pch=16,xlab="WD",ylab="LMA",main="WD vs LMA")
points(WD_est,LMA[ind_WD],col="red",pch=16)
#test TLP
ind_WD=which(!is.na(WD) & !is.na(TLP))
mod_WD <- sma(data.frame(WD[ind_WD],TLP[ind_WD]),method=c("SMA"))
summary(mod_WD)
WD_est <- mod_WD$coef[[1]][1,1] + mod_WD$coef[[1]][2,1]*TLP[ind_WD]

plot(mod_TLP_WD,"SMA",pch=16,xlab="TLP",ylab="WD",main="TLP vs WD")
points(TLP[ind_WD],WD_est,col="red",pch=16)
#P50
ind_WD=which(!is.na(WD) & !is.na(P50))
mod_WD <- sma(data.frame(WD[ind_WD],P50[ind_WD]),method=c("SMA"))
summary(mod_WD)
WD_est <- mod_WD$coef[[1]][1,1] + mod_WD$coef[[1]][2,1]*P50[ind_WD]
#WD_est <- 0.4031998 + 0.2291030*P50[ind_WD]
plot(mod_P50_WD,"SMA",pch=16,xlab="P50",ylab="WD",main="P50 vs WD")
points(P50[ind_WD],WD_est,col="red",pch=16)

#something weird so strong correlation




#4)Estimate TLP from P50,LS,LMA and WD (223 sp)
ind_TLP=which(!is.na(P50) & !is.na(TLP) & !is.na(LMA) & !is.na(LS)& !is.na(WD))
mod_TLP <- sma_regress_multivar(data.frame(P50[ind_TLP],LS[ind_TLP],LMA[ind_TLP],WD[ind_TLP],TLP[ind_TLP]),nbtstrp)
TLP_est <- mod_TLP$intercept_R + mod_TLP$slope_R.y1*P50[ind_TLP] + mod_TLP$slope_R.y2*LS[ind_TLP] + mod_TLP$slope_R.y3*LMA[ind_TLP]+mod_TLP$slope_R.y4*WD[ind_TLP]
#Cross-checks to make sure that result is sensible
plot(mod_TLP_P50,"SMA",pch=16,xlab="TLP",ylab="P50",main="TLP vs P50")
points(TLP_est,P50[ind_TLP],col="red",pch=16)
plot(mod_LS_TLP,"SMA",pch=16,xlab="LS",ylab="TLP",main="LS vs TLP")
points(LS[ind_TLP],TLP_est,col="red",pch=16)
plot(mod_TLP_LMA,"SMA",pch=16,xlab="TLP",ylab="LMA",main="TLP vs LMA")
points(TLP_est,LMA[ind_TLP],col="red",pch=16)
plot(mod_TLP_WD,"SMA",pch=16,xlab="TLP",ylab="WD",main="TLP vs WD")
points(TLP_est,WD[ind_TLP],col="red",pch=16)

##checcked total SP by no WD is 226, no LS is 352, no LMA is 279,noP50 is 264
#4_1)Estimate TLP from P50,LS,LMA,WD 
ind_TLP=which(!is.na(P50) & !is.na(TLP) & !is.na(LMA) & !is.na(WD))
mod_TLP <- sma_regress_multivar(data.frame(P50[ind_TLP],LMA[ind_TLP],WD[ind_TLP],TLP[ind_TLP]),nbtstrp)
TLP_est <- mod_TLP$intercept_R + mod_TLP$slope_R.y1*P50[ind_TLP] + mod_TLP$slope_R.y2*LMA[ind_TLP] + mod_TLP$slope_R.y3*WD[ind_TLP]
#Cross-checks to make sure that result is sensible
plot(mod_TLP_P50,"SMA",pch=16,xlab="TLP",ylab="P50",main="TLP vs P50")
points(TLP_est,P50[ind_TLP],col="red",pch=16)
plot(mod_TLP_LMA,"SMA",pch=16,xlab="TLP",ylab="LMA",main="TLP vs LMA")
points(TLP_est,LMA[ind_TLP],col="red",pch=16)
plot(mod_TLP_WD,"SMA",pch=16,xlab="TLP",ylab="WD",main="TLP vs WD")
points(TLP_est,WD[ind_TLP],col="red",pch=16)
#if do not incclud LS, the total sp is 352. So it has large number of SP. And the regrressionss are ok even though for WD and TLP iis not strong


#5)Estimate slope from TLP, P50 and kstem
ind_slope=which(!is.na(P50) & !is.na(TLP) & !is.na(Ks) & !is.na(slope))
mod_slope <- sma_regress_multivar(data.frame(P50[ind_slope],TLP[ind_slope],Ks[ind_slope],slope[ind_slope]),nbtstrp)
slope_est <- mod_slope$intercept_R + mod_slope$slope_R.y1*P50[ind_slope] + mod_slope$slope_R.y2*TLP[ind_slope] + mod_slope$slope_R.y3*Ks[ind_slope]
#Cross-checks to make sure that result is sensible
plot(mod_Ks_slope,"SMA",pch=16,xlab="kstem",ylab="slope",main="kstem vs slope")
points(Ks[ind_slope],slope_est,col="red",pch=16)
plot(mod_TLP_slope,"SMA",pch=16,xlab="TLP",ylab="slope",main="TLP vs slope")
points(TLP[ind_slope],slope_est,col="red",pch=16)
plot(mod_P50_slope,"SMA",pch=16,xlab="P50",ylab="slope",main="P50 vs slope")
points(P50[ind_slope],slope_est,col="red",pch=16)


#6) Estimate Ks from LMA and LS/Hmax.
ind_Ks=which(!is.na(P50) & !is.na(LS_Hmax)  & !is.na(Ks))
mod_Ks <- sma_regress_multivar(data.frame(P50[ind_Ks],LS_Hmax[ind_Ks],Ks[ind_Ks]),nbtstrp)
Ks_est <- mod_Ks$intercept_R + mod_Ks$slope_R.y1*P50[ind_Ks] + mod_Ks$slope_R.y2*LS_Hmax[ind_Ks]
#Cross-checks to make sure that result is sensible
plot(mod_Ks_LS_Hmax,"SMA",pch=16,xlab="LS_Hmax",ylab="Ks",main="Ksvs LS_Hmax")
points(LS_Hmax[ind_Ks],Ks_est,col="red",pch=16)
plot(mod_Ks_P50,"SMA",pch=16,xlab="Ks",ylab="P50",main="Ks vs P50")
points(P50[ind_Ks],Ks_est,col="red",pch=16)

#--- Attempt to iteratively converge on the best fit values of Ks, TLP, P50 and LMA given known Hmax and LS ---

#Calculate minimum and maximum values
maxP50=max(P50,na.rm=T)
minP50=min(P50,na.rm=T)
maxTLP=max(TLP,na.rm=T)
minTLP=min(TLP,na.rm=T)
maxKs=max(Ks,na.rm=T)
minKs=min(Ks,na.rm=T)
maxLS=max(LS,na.rm=T)
minLS=min(LS,na.rm=T)
maxLMA=max(LMA,na.rm=T)
minLMA=min(LMA,na.rm=T)
maxWD=max(WD,na.rm=T)
minWD=min(WD,na.rm=T)
maxslope=max(slope,na.rm=T)
minslope=min(slope,na.rm=T)

# Set the tolerance level for the iteration
tol=0.00001

#Go through all observed combinations of Hmax and LS

ind=which(!is.na(Hmax) & !is.na(LS))

Hmax_e=Hmax[ind]
LS_e=LS[ind]

ndata=length(LS_e)

P50_e <- matrix(NA, nrow = ndata, ncol = n_uncer) #Array now expanded to hold multiple replicate estimates based on regression coefficient uncertainty
LMA_e <- matrix(NA, nrow = ndata, ncol = n_uncer)
TLP_e <- matrix(NA, nrow = ndata, ncol = n_uncer)
WD_e <- matrix(NA, nrow = ndata, ncol = n_uncer)
slope_e <- matrix(NA, nrow = ndata, ncol = n_uncer)

# Decide whether to limit the possible ranges of predicted traits to the observed values (T) or not (F)
limitdataranges=F

# New outer loop which randomly samples regression coefficients from within their uncertainty bounds
# The random sampling comes from the bootstrap sampling done in the calculations of the SMA regressions themselves. This approach has the big advantages of (a) not having to make any assumptions about the distribution of the coefficient uncertainty and (b) ensuring that the individual slope coefficients within a regression are consistent with each other.
for (ss in 1:nbtstrp) {
# For now, make samples for the following:
mod_LMA_slope_y1_sample <- mod_LMA$boot.y1[ss]
mod_LMA_slope_y2_sample <- mod_LMA$boot.y2[ss]
mod_TLP_slope_y1_sample <- mod_TLP$boot.y1[ss]
mod_TLP_slope_y2_sample <- mod_TLP$boot.y2[ss]
mod_TLP_slope_y3_sample <- mod_TLP$boot.y3[ss]
mod_P50_slope_y1_sample <- mod_P50$boot.y1[ss]
mod_P50_slope_y2_sample <- mod_P50$boot.y2[ss]
# Not implemented for WD yet, as I'm not sure which regression equation to use
# These regression coefficients will now be used in the optimisation calculations


# Outer loop over all the combinations of Hmax and LS
# The new estimates of traits use the suffix "_e"
for (dd in 1:ndata) {

  #Calculate the value for LS*Hmax, which comes direct from the two input traits
  LS_Hmax_e=log(exp(LS_e[dd])*Hmax_e[dd],base=exp(1))
  
  #Calculate the values of Ks based only on the bivariate relationships with LS*Hmax
  Ks_e[dd]= mod_Ks_LS_Hmax$regression.results$Slope[3]*LS_Hmax_e +
    mod_Ks_LS_Hmax$regression.results$Intercept[3]
    
  if (limitdataranges) {
    #Do not go beyond observed limits of data
    if (Ks_e[dd]>maxKs) {Ks_e[dd]=maxKs}
    if (Ks_e[dd]<minKs) {Ks_e[dd]=minKs}
  }

  #TLP, P50, LMA, WD need optimising

  #First set some initial based on simple bivariate relationship. This is just so that the iteration has somewhere to start from. Final result should not be sensitive to these.
  LMA_e_last = mod_LS_LMA$regression.results$Slope[3]*LS_e[dd] +
    mod_LS_LMA$regression.results$Intercept[3]
  TLP_e_last = mod_LS_TLP$regression.results$Slope[3]*LS_e[dd] +
    mod_LS_TLP$regression.results$Intercept[3]
  P50_e_last = mod_TLP_P50$regression.results$Slope[3]*TLP_e_last +
    mod_TLP_P50$regression.results$Intercept[3]
  WD_e_last = mod_TLP_WD$regression.results$Slope[3]*TLP_e_last +
    mod_TLP_WD$regression.results$Intercept[3] #Not sure if TLP vs WD is the best choice, but as only for initialisation shouldn't be too important.

  # "diff_" variables hold the difference between the current estimate of a trait value "_e" and the previous
  # estimate "_last"
  # "diff_*_last" variables contain the differences from the last round of iteration
  # (these are compared to differences in the current round of iteration to see if changes are smaller than
  # "tol" and therefore the iteration can stop)
  # Here we initialise the "diff_*_last" variables very high (why set 100)
  diff_P50_last=100
  diff_LMA_last=100
  diff_TLP_last=100
  diff_WD_last=100

  # These arrays are just for output, they store the values of every iteration for the current datapoint.
  # Useful for debugging and to check that convergence is working.
  # (only for debugging, can be commented out)
  P50_c <- matrix(NA, nrow = 100)
  LMA_c <- matrix(NA, nrow = 100)
  TLP_c <- matrix(NA, nrow = 100)
  WD_c <- matrix(NA, nrow = 100)
  
  # Now we start the optimisation loop. Trait values are iterated until the difference between trait
  # values on successive iterations is less than "tol".
  niter=0;
  while (T) {
    niter=niter+1 # Number of iterations completed

    # Make estimates of trait values based on the best SMA regressions (probably multivariate in most cases)
    # The estimates of traits in each iteration are based on the estimates of their predictor traits from the previous iteration
    LMA_e[dd,ss]=mod_LMA$intercept_R + mod_LMA_slope_y1_sample*TLP_e_last + mod_LMA_slope_y2_sample*LS_e[dd]
    TLP_e[dd,ss]=mod_TLP$intercept_R + mod_TLP_slope_y1_sample*P50_e_last +  mod_TLP_slope_y2_sample*LMA_e_last+mod_TLP_slope_y3_sample*WD_e_last
    P50_e[dd,ss]=mod_P50$intercept_R + mod_P50_slope_y1_sample*TLP_e_last + mod_P50_slope_y2_sample*Ks_e[dd]
    #WD_e[dd,ss]=mod_WD$intercept_R + mod_WD$slope_R.y1*TLP_e_last + mod_WD$slope_R.y2*P50_e_last + mod_WD$slope_R.y3*LMA_e_last
    # Test taking WD from a simple bivariate relationship (because not converging with the above multivariate one - need to sort out the multivariate fit!)
    WD_e[dd,ss] = mod_TLP_WD$regression.results$Slope[3]*TLP_e_last + mod_TLP_WD$regression.results$Intercept[3]
    
    if (limitdataranges) {
      #Do not go beyond observed limits of data
      if (P50_e[dd,ss]>maxP50 | is.na(P50_e[dd,ss])) {P50_e[dd,ss]=NA; break}
      if (P50_e[dd,ss]<minP50 | is.na(P50_e[dd,ss])) {P50_e[dd,ss]=NA; break}
      if (TLP_e[dd,ss]>maxTLP | is.na(TLP_e[dd,ss])) {TLP_e[dd,ss]=NA; break}
      if (TLP_e[dd,ss]<minTLP | is.na(TLP_e[dd,ss])) {TLP_e[dd,ss]=NA; break}
      if (LMA_e[dd,ss]>maxLMA | is.na(LMA_e[dd,ss])) {LMA_e[dd,ss]=NA; break}
      if (LMA_e[dd,ss]<minLMA | is.na(LMA_e[dd,ss])) {LMA_e[dd,ss]=NA; break}
      if (WD_e[dd,ss]>maxWD | is.na(WD_e[dd,ss])) {WD_e[dd,ss]=NA; break}
      if (WD_e[dd,ss]<minWD | is.na(WD_e[dd,ss])) {WD_e[dd,ss]=NA; break}
    }
    
    # Save the values for this iteration to the output array (only for debugging, can be commented out)
    P50_c[niter] <- P50_e[dd,ss]
    LMA_c[niter] <- LMA_e[dd,ss]
    TLP_c[niter] <- TLP_e[dd,ss]
    WD_c[niter] <- WD_e[dd,ss]

    # Calculate the difference between the current estimate of a trait value "_e" and the previous estimate "_last"
    diff_P50 = P50_e[dd,ss]-P50_e_last
    diff_LMA = LMA_e[dd,ss]-LMA_e_last
    diff_TLP = TLP_e[dd,ss]-TLP_e_last
    diff_WD = WD_e[dd,ss]-WD_e_last

    # Now we test if the difference between trait estimates on this iteration and between trait estimates on
    # the last iteration is less than "tol" for all traits. If it is we finish the iteration.
    if (abs(diff_P50-diff_P50_last)<tol &&
      abs(diff_LMA-diff_LMA_last)<tol &&
      abs(diff_WD-diff_WD_last)<tol &&
      abs(diff_TLP-diff_TLP_last)<tol) {
      break
    }

    # Save the "diff" values ready for the next iteration
    diff_P50_last=diff_P50
    diff_LMA_last=diff_LMA
    diff_TLP_last=diff_TLP
    diff_WD_last=diff_WD

    # Save the "_e" values ready for the next iteration
    P50_e_last=P50_e[dd,ss]
    LMA_e_last=LMA_e[dd,ss]
    TLP_e_last=TLP_e[dd,ss]
    WD_e_last=WD_e[dd,ss]
  }

  # After the iteration has finished we can calculate any traits which did not need to be included in the optimisation (because they are not used in the input to calculate any other trait)
  #WD_e[dd]=mod_WD$intercept_R + mod_WD$slope_R.y1*TLP_e[dd] + mod_WD$slope_R.y2*P50_e[dd] + mod_WD$slope_R.y3*LMA_e[dd]
  slope_e[dd,ss]=mod_slope$intercept_R + mod_slope$slope_R.y1*P50_e[dd,ss] + mod_slope$slope_R.y2*TLP_e[dd,ss] + mod_slope$slope_R.y3*Ks_e[dd]
  
  if (limitdataranges) {
    #Do not go beyond observed limits of data
    if (slope_e[dd,ss]>maxslope | is.na(slope_e[dd,ss])) {slope_e[dd,ss]=NA}
    if (slope_e[dd,ss]<minslope | is.na(slope_e[dd,ss])) {slope_e[dd,ss]=NA}
  }

}

} #Finish nbtstrp loop

#Optionally limit to ranges of observed traits
WD_e[WD_e>maxWD]=maxWD
WD_e[WD_e<minWD]=minWD
LMA_e[LMA_e>maxLMA]=maxLMA
LMA_e[LMA_e<minLMA]=minLMA
LS_e[LS_e>maxLS]=maxLS
LS_e[LS_e<minLS]=minLS
TLP_e[TLP_e>maxTLP]=maxTLP
TLP_e[TLP_e<minTLP]=minTLP
P50_e[P50_e>maxP50]=maxP50
P50_e[P50_e<minP50]=minP50
slope_e[slope_e>maxslope]=maxslope
slope_e[slope_e<minslope]=minslope

Hmax_e
exp(Ks_e)
WD_e
exp(LMA_e)
exp(LS_e)
exp(slope_e)
-exp(P50_e)
-exp(TLP_e)

#Make plots to compare with original data
par(mfrow=c(4,4))
par(mar=c(2,2,2,2))

plot(mod_TLP_P50,"SMA",pch=16,xlab="TLP",ylab="P50",main="TLP vs P50")
points(TLP_e,P50_e,col="red",pch=16)

plot(mod_TLP_slope,"SMA",pch=16,xlab="TLP",ylab="slope",main="TLP vs slope")
points(TLP_e,slope_e,col="red",pch=16)

plot(mod_P50_slope,"SMA",pch=16,xlab="P50",ylab="slope",main="P50 vs slope")
points(P50_e,slope_e,col="red",pch=16)

plot(mod_TLP_WD,"SMA",pch=16,xlab="TLP",ylab="WD",main="TLP vs WD")
points(TLP_e,WD_e,col="red",pch=16)

plot(mod_P50_WD,"SMA",pch=16,xlab="P50",ylab="WD",main="P50 vs WD")
points(P50_e,WD_e,col="red",pch=16)

plot(mod_TLP_LMA,"SMA",pch=16,xlab="TLP",ylab="LMA",main="TLP vs LMA")
points(TLP_e,LMA_e,col="red",pch=16)

plot(mod_LS_LMA,"SMA",pch=16,xlab="LS",ylab="LMA",main="LS vs LMA")
points(LS_e,LMA_e,col="red",pch=16)

plot(mod_Ks_LS_Hmax,"SMA",pch=16,xlab="LS_Hmax",ylab="LS",main="LS_Hmax vs Ks")
points(log(exp(LS_e)*Hmax_e),LS_e,col="red",pch=16)
#plot(log(exp(kstem_e)/Hmax_e),LS_e,col="red",pch=16)

plot(mod_Ks_P50,"SMA",pch=16,xlab="Ks",ylab="P50",main="Ks vs P50")
points(Ks_e,P50_e,col="red",pch=16)

plot(mod_LS_TLP,"SMA",pch=16,xlab="LS",ylab="TLP",main="LS vs TLP")
points(LS_e,TLP_e,col="red",pch=16)

plot(mod_WD_LMA,"SMA",pch=16,xlab="WD",ylab="LMA",main="WD vs LMA")
points(WD_e,LMA_e,col="red",pch=16) 

plot(mod_Ks_slope,"SMA",pch=16,xlab="Ks",ylab="slope",main="Ks vs slope")
points(Ks_e,slope_e,col="red",pch=16) 

plot(mod_Ks_Hmax,"SMA",pch=16,xlab="Ks",ylab="Hmax",main="Ks vs Hmax")
points(Ks_e,Hmax_e,col="red",pch=16) 

#Use up remaining unallocated plots and set back to single plot
par(mfrow=c(1,1))
par(mar=c(5.1,4.1,4.1,2.1))

#--- Calculate regressions that characterise the optimised relationship ---

par(mfrow=c(4,4))
par(mar=c(2,2,2,2))

mod_TLP_P50_opt <- lmodel2(P50_e ~ TLP_e, data.frame(TLP_e,P50_e))
plot(mod_TLP_P50_opt,"SMA",pch=16,xlab="TLP",ylab="P50",main="TLP vs P50")

mod_TLP_slope_opt <- lmodel2(slope_e ~ TLP_e, data.frame(TLP_e,slope_e))
plot(mod_TLP_slope_opt,"SMA",pch=16,xlab="TLP",ylab="slope",main="TLP vs slope")

mod_P50_slope_opt <- lmodel2(slope_e ~ P50_e, data.frame(P50_e,slope_e))
plot(mod_P50_slope_opt,"SMA",pch=16,xlab="P50",ylab="slope",main="P50 vs slope")

mod_TLP_WD_opt <- lmodel2(WD_e ~ TLP_e, data.frame(TLP_e,WD_e))
plot(mod_TLP_WD_opt,"SMA",pch=16,xlab="TLP",ylab="WD",main="TLP vs WD")

mod_P50_WD_opt <- lmodel2(WD_e ~ P50_e, data.frame(P50_e,WD_e))
plot(mod_P50_WD_opt,"SMA",pch=16,xlab="P50",ylab="WD",main="P50 vs WD")

mod_TLP_LMA_opt <- lmodel2(LMA_e ~ TLP_e, data.frame(TLP_e,LMA_e))
plot(mod_TLP_LMA_opt,"SMA",pch=16,xlab="TLP",ylab="LMA",main="TLP vs LMA")

mod_LS_LMA_opt <- lmodel2(LMA_e ~ LS_e, data.frame(LS_e,LMA_e))
plot(mod_LS_LMA_opt,"SMA",pch=16,xlab="LS",ylab="LMA",main="LS vs LMA")

#mod_kstemHmax_LS_opt <- lmodel2(LS_e ~ kstemHmax_e, data.frame(kstemHmax_e,LS_e))
#plot(mod_kstemHmax_LS_opt,"SMA",pch=16,xlab="kstemHmax",ylab="LS",main="kstem/Hmax vs LS")

mod_Ks_P50_opt <- lmodel2(P50_e ~ Ks_e, data.frame(Ks_e,P50_e))
plot(mod_Ks_P50_opt,"SMA",pch=16,xlab="Ks",ylab="P50",main="Ks vs P50")

mod_LS_TLP_opt <- lmodel2(TLP_e ~ LS_e, data.frame(LS_e,TLP_e))
plot(mod_LS_TLP_opt,"SMA",pch=16,xlab="LS",ylab="TLP",main="LS vs TLP")

mod_WD_LMA_opt <- lmodel2(LMA_e ~ WD_e, data.frame(WD_e,LMA_e))
plot(mod_WD_LMA_opt,"SMA",pch=16,xlab="WD",ylab="LMA",main="WD vs LMA")

mod_Ks_slope_opt <- lmodel2(slope_e ~ Ks_e, data.frame(Ks_e,slope_e))

plot(mod_Ks_slope_opt,"SMA",pch=16,xlab="Ks",ylab="slope",main="Ks vs slope")

mod_Ks_Hmax <- lmodel2(Hmax ~ Ks, data.frame(Ks,Hmax))
plot(mod_Ks_Hmax,"SMA",pch=16,xlab="Ks",ylab="Hmax",main="Ks vs Hmax")

#Use up remaining unallocated plots and set back to single plot
#plot.new()
#plot.new()
#plot.new()
par(mfrow=c(1,1))
par(mar=c(5.1,4.1,4.1,2.1))
