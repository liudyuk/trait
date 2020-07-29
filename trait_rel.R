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
coeffnames_P50_from_TLP_Ks <- c("Coefficient","L95","U95")
intercept_P50_from_TLP_Ks <- c(P50_from_TLP_Ks$mod$intercept_R,P50_from_TLP_Ks$mod$L95_R.intercept,P50_from_TLP_Ks$mod$U95_R.intercept)
y1_P50_from_TLP_Ks <- c(P50_from_TLP_Ks$mod$slope_R.y1,P50_from_TLP_Ks$mod$L95_R.y1,P50_from_TLP_Ks$mod$U95_R.y1)
y2_P50_from_TLP_Ks <- c(P50_from_TLP_Ks$mod$slope_R.y2,P50_from_TLP_Ks$mod$L95_R.y2,P50_from_TLP_Ks$mod$U95_R.y2)

coeff_P50_from_TLP_Ks <- data.frame(coeffnames_P50_from_TLP_Ks,intercept_P50_from_TLP_Ks,y1_P50_from_TLP_Ks,y2_P50_from_TLP_Ks)
View(coeff_P50_from_TLP_Ks)
# NOTE: These coefficients are all over the place, almost certainly because we have high multicolinearity in the predictors, BUT, this is not a problem as we are not interpreting the coefficients, just using them for the prediction.


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


# Ks fits -----------------------------------------------------------------

# Ks from LS_Hmax, P50 and slope
Ks_from_LSHmax_P50_slope <- sma_plot_stats(data.frame(LS_Hmax,P50,slope,Ks),c("LS*Hmax","P50","slope","Ks"),nbtstrp)
plot(Ks[Ks_from_LSHmax_P50_slope$dataused],Ks_from_LSHmax_P50_slope$var_est,pch=16,xlab="Ks",ylab="Ks_est",main="Ks vs Ks_est")

# Ks from LS_Hmax and P50
Ks_from_LSHmax_P50 <- sma_plot_stats(data.frame(LS_Hmax,P50,Ks),c("LS*Hmax","P50","Ks"),nbtstrp)
plot(Ks[Ks_from_LSHmax_P50$dataused],Ks_from_LSHmax_P50$var_est,pch=16,xlab="Ks",ylab="Ks_est",main="Ks vs Ks_est")

# Ks from LS_Hmax and slope
Ks_from_LSHmax_slope <- sma_plot_stats(data.frame(LS_Hmax,slope,Ks),c("LS*Hmax","slope","Ks"),nbtstrp)
plot(Ks[Ks_from_LSHmax_slope$dataused],Ks_from_LSHmax_slope$var_est,pch=16,xlab="Ks",ylab="Ks_est",main="Ks vs Ks_est")

# Summarise statistics
all_testnames_Ks <- c("Ks_from_LSHmax_P50_slope","Ks_from_LSHmax_P50","Ks_from_LSHmax_slope","Ks_from_LSHmax","Ks_from_P50","Ks_from_slope")
all_R2_Ks <- c(Ks_from_LSHmax_P50_slope$R2,Ks_from_LSHmax_P50$R2,Ks_from_LSHmax_slope$R2,Ks_from_LSHmax$R2,Ks_from_P50$R2,Ks_from_slope$R2)
all_R2adj_Ks <- c(Ks_from_LSHmax_P50_slope$R2adj,Ks_from_LSHmax_P50$R2adj,Ks_from_LSHmax_slope$R2adj,Ks_from_LSHmax$R2adj,Ks_from_P50$R2adj,Ks_from_slope$R2adj)
all_rmse_Ks <- c(Ks_from_LSHmax_P50_slope$rmse,Ks_from_LSHmax_P50$rmse,Ks_from_LSHmax_slope$rmse,Ks_from_LSHmax$rmse,Ks_from_P50$rmse,Ks_from_slope$rmse)
all_ndata_Ks <- c(Ks_from_LSHmax_P50_slope$ndata,Ks_from_LSHmax_P50$ndata,Ks_from_LSHmax_slope$ndata,Ks_from_LSHmax$ndata,Ks_from_P50$ndata,Ks_from_slope$ndata)

all_Ks <- data.frame(all_testnames_Ks,all_R2_Ks,all_R2adj_Ks,all_rmse_Ks,all_ndata_Ks)
View(all_Ks)

# CHOICE: Ks_from_LSHmax_P50 has the best combination of R2adj and RMSE

# Test MAT and PPT coverage of species for chosen model
plot(MAT[Ks_from_LSHmax_P50$dataused],MAP[Ks_from_LSHmax_P50$dataused])
# WIDE CLIMATE COVERAGE

# DECISION: Ks_from_LSHmax_P50


# WD fits -----------------------------------------------------------------

# NOTE: Ignoring fits with TLP and LMA because they do not fit conceptually with the hypothesis framework

# WD from P50 and slope
WD_from_P50_slope <- sma_plot_stats(data.frame(P50,slope,WD),c("P50","slope","WD"),nbtstrp)
plot(WD[WD_from_P50_slope$dataused],WD_from_P50_slope$var_est,pch=16,xlab="WD",ylab="WD_est",main="WD vs WD_est")

# WD from slope and P50*slope (interaction term) - NOTE: Testing the interactions because P50 and slope are so highly correlated with each other
WD_from_slope_P50slope <- sma_plot_stats(data.frame(slope,P50*slope,WD),c("slope","P50*slope","WD"),nbtstrp)
plot(WD[WD_from_slope_P50slope$dataused],WD_from_slope_P50slope$var_est,pch=16,xlab="WD",ylab="WD_est",main="WD vs WD_est")

# WD from P50 and P50*slope (interaction term)
WD_from_P50_P50slope <- sma_plot_stats(data.frame(P50,P50*slope,WD),c("P50","P50*slope","WD"),nbtstrp)
plot(WD[WD_from_P50_P50slope$dataused],WD_from_P50_P50slope$var_est,pch=16,xlab="WD",ylab="WD_est",main="WD vs WD_est")

# WD from P50 and slope
WD_from_P50 <- sma_plot_stats(data.frame(P50,WD),c("P50","WD"),nbtstrp)
plot(WD[WD_from_P50$dataused],WD_from_P50$var_est,pch=16,xlab="WD",ylab="WD_est",main="WD vs WD_est")

# WD from slope
WD_from_slope <- sma_plot_stats(data.frame(slope,WD),c("slope","WD"),nbtstrp)
plot(WD[WD_from_slope$dataused],WD_from_slope$var_est,pch=16,xlab="WD",ylab="WD_est",main="WD vs WD_est")

# Summarise statistics
all_testnames_WD <- c("WD_from_P50_slope","WD_from_slope_P50slope","WD_from_P50_P50slope","WD_from_P50","WD_from_slope")
all_R2_WD <- c(WD_from_P50_slope$R2,WD_from_slope_P50slope$R2,WD_from_P50_P50slope$R2,WD_from_P50$R2,WD_from_slope$R2)
all_R2adj_WD <- c(WD_from_P50_slope$R2adj,WD_from_slope_P50slope$R2adj,WD_from_P50_P50slope$R2adj,WD_from_P50$R2adj,WD_from_slope$R2adj)
all_rmse_WD <- c(WD_from_P50_slope$rmse,WD_from_slope_P50slope$rmse,WD_from_P50_P50slope$rmse,WD_from_P50$rmse,WD_from_slope$rmse)
all_ndata_WD <- c(WD_from_P50_slope$ndata,WD_from_slope_P50slope$ndata,WD_from_P50_P50slope$ndata,WD_from_P50$ndata,WD_from_slope$ndata)

all_WD <- data.frame(all_testnames_WD,all_R2_WD,all_R2adj_WD,all_rmse_WD,all_ndata_WD)
View(all_WD)

# CHOICE: WD_from_slope_P50slope has the best combination of R2adj and RMSE

# Test MAT and PPT coverage of species for chosen model
plot(MAT[WD_from_slope_P50slope$dataused],MAP[WD_from_slope_P50slope$dataused])
# WIDE CLIMATE COVERAGE

# DECISION: WD_from_slope_P50slope

# slope fits --------------------------------------------------------------

# slope from P50, TLP, WD and Ks
slope_from_P50_TLP_WD_Ks <- sma_plot_stats(data.frame(P50,TLP,WD,Ks,slope),c("P50","TLP","WD","Ks","slope"),nbtstrp)
plot(slope[slope_from_P50_TLP_WD_Ks$dataused],slope_from_P50_TLP_WD_Ks$var_est,pch=16,xlab="slope",ylab="slope_est",main="slope vs slope_est")

# slope from P50, TLP and WD
slope_from_P50_TLP_WD <- sma_plot_stats(data.frame(P50,TLP,WD,slope),c("P50","TLP","WD","slope"),nbtstrp)
plot(slope[slope_from_P50_TLP_WD$dataused],slope_from_P50_TLP_WD$var_est,pch=16,xlab="slope",ylab="slope_est",main="slope vs slope_est")

# slope from P50, TLP and Ks
slope_from_P50_TLP_Ks <- sma_plot_stats(data.frame(P50,TLP,Ks,slope),c("P50","TLP","Ks","slope"),nbtstrp)
plot(slope[slope_from_P50_TLP_Ks$dataused],slope_from_P50_TLP_Ks$var_est,pch=16,xlab="slope",ylab="slope_est",main="slope vs slope_est")

# slope from P50, WD and Ks
slope_from_P50_WD_Ks <- sma_plot_stats(data.frame(P50,WD,Ks,slope),c("P50","WD","Ks","slope"),nbtstrp)
plot(slope[slope_from_P50_WD_Ks$dataused],slope_from_P50_WD_Ks$var_est,pch=16,xlab="slope",ylab="slope_est",main="slope vs slope_est")

# slope from TLP, WD and Ks
slope_from_TLP_WD_Ks <- sma_plot_stats(data.frame(TLP,WD,Ks,slope),c("TLP","WD","Ks","slope"),nbtstrp)
plot(slope[slope_from_TLP_WD_Ks$dataused],slope_from_TLP_WD_Ks$var_est,pch=16,xlab="slope",ylab="slope_est",main="slope vs slope_est")

# slope from P50 and TLP
slope_from_P50_TLP <- sma_plot_stats(data.frame(P50,TLP,slope),c("P50","TLP","slope"),nbtstrp)
plot(slope[slope_from_P50_TLP$dataused],slope_from_P50_TLP$var_est,pch=16,xlab="slope",ylab="slope_est",main="slope vs slope_est")

# slope from P50 and TLPP50 (interaction)
slope_from_P50_TLPP50 <- sma_plot_stats(data.frame(P50,TLP*P50,slope),c("P50","TLP*P50","slope"),nbtstrp)
plot(slope[slope_from_P50_TLPP50$dataused],slope_from_P50_TLPP50$var_est,pch=16,xlab="slope",ylab="slope_est",main="slope vs slope_est")

# slope from TLP and TLPP50 (interaction)
slope_from_TLP_TLPP50 <- sma_plot_stats(data.frame(TLP,TLP*P50,slope),c("TLP","TLP*P50","slope"),nbtstrp)
plot(slope[slope_from_TLP_TLPP50$dataused],slope_from_TLP_TLPP50$var_est,pch=16,xlab="slope",ylab="slope_est",main="slope vs slope_est")

# slope from TLP and Ks
slope_from_TLP_Ks <- sma_plot_stats(data.frame(TLP,Ks,slope),c("TLP","Ks","slope"),nbtstrp)
plot(slope[slope_from_TLP_Ks$dataused],slope_from_TLP_Ks$var_est,pch=16,xlab="slope",ylab="slope_est",main="slope vs slope_est")

# slope from TLP and WD
slope_from_TLP_WD <- sma_plot_stats(data.frame(TLP,WD,slope),c("TLP","WD","slope"),nbtstrp)
plot(slope[slope_from_TLP_WD$dataused],slope_from_TLP_WD$var_est,pch=16,xlab="slope",ylab="slope_est",main="slope vs slope_est")

# slope from TLP
slope_from_TLP <- sma_plot_stats(data.frame(TLP,slope),c("TLP","slope"),nbtstrp)
plot(slope[slope_from_TLP$dataused],slope_from_TLP$var_est,pch=16,xlab="slope",ylab="slope_est",main="slope vs slope_est")

# slope from P50
slope_from_P50 <- sma_plot_stats(data.frame(P50,slope),c("P50","slope"),nbtstrp)
plot(slope[slope_from_P50$dataused],slope_from_P50$var_est,pch=16,xlab="slope",ylab="slope_est",main="slope vs slope_est")

# slope from WD
slope_from_WD <- sma_plot_stats(data.frame(WD,slope),c("WD","slope"),nbtstrp)
plot(slope[slope_from_WD$dataused],slope_from_WD$var_est,pch=16,xlab="slope",ylab="slope_est",main="slope vs slope_est")

# slope from Ks
slope_from_Ks <- sma_plot_stats(data.frame(Ks,slope),c("Ks","slope"),nbtstrp)
plot(slope[slope_from_Ks$dataused],slope_from_Ks$var_est,pch=16,xlab="slope",ylab="slope_est",main="slope vs slope_est")

# slope from P50, WD and Ks (same species as slope_from_P50_TLP_WD_Ks) (NOTE: testing because the fit for slope_from_P50_WD_Ks is only marginally worse than for slope_from_P50_TLP_Ks, but it has many more species)
slope_from_P50_WD_Ks_limitspec <- sma_plot_stats(data.frame(P50,WD,Ks,slope),c("P50","WD","Ks","slope"),nbtstrp,F,slope_from_P50_TLP_WD_Ks$dataused)
plot(slope[slope_from_P50_WD_Ks_limitspec$dataused],slope_from_P50_WD_Ks_limitspec$var_est,pch=16,xlab="slope",ylab="slope_est",main="slope vs slope_est")

# slope from P50, TLP and Ks (same species as slope_from_P50_TLP_WD_Ks) (NOTE: testing because the fit for slope_from_P50_WD_Ks is only marginally worse than for slope_from_P50_TLP_Ks, but it has many more species)
slope_from_P50_TLP_Ks_limitspec <- sma_plot_stats(data.frame(P50,TLP,Ks,slope),c("P50","TLP","Ks","slope"),nbtstrp,F,slope_from_P50_TLP_WD_Ks$dataused)
plot(slope[slope_from_P50_TLP_Ks_limitspec$dataused],slope_from_P50_TLP_Ks_limitspec$var_est,pch=16,xlab="slope",ylab="slope_est",main="slope vs slope_est")

# slope from P50 (same species as slope_from_P50_TLP_WD_Ks)
slope_from_P50_limitspec <- sma_plot_stats(data.frame(P50,slope),c("P50","slope"),nbtstrp,F,slope_from_P50_TLP_WD_Ks$dataused)
plot(slope[slope_from_P50_limitspec$dataused],slope_from_P50_limitspec$var_est,pch=16,xlab="slope",ylab="slope_est",main="slope vs slope_est")

# Summarise statistics
all_testnames_slope <- c("slope_from_P50_TLP_WD_Ks","slope_from_P50_TLP_WD","slope_from_P50_TLP_Ks","slope_from_P50_WD_Ks","slope_from_TLP_WD_Ks","slope_from_P50_TLP","slope_from_P50_TLPP50","slope_from_TLP_TLPP50","slope_from_TLP_Ks","slope_from_TLP_WD","slope_from_TLP","slope_from_P50","slope_from_WD","slope_from_Ks","slope_from_P50_WD_Ks_limitspec","slope_from_P50_TLP_Ks_limitspec","slope_from_P50_limitspec")
all_R2_slope <- c(slope_from_P50_TLP_WD_Ks$R2,slope_from_P50_TLP_WD$R2,slope_from_P50_TLP_Ks$R2,slope_from_P50_WD_Ks$R2,slope_from_TLP_WD_Ks$R2,slope_from_P50_TLP$R2,slope_from_P50_TLPP50$R2,slope_from_TLP_TLPP50$R2,slope_from_TLP_Ks$R2,slope_from_TLP_WD$R2,slope_from_TLP$R2,slope_from_P50$R2,slope_from_WD$R2,slope_from_Ks$R2,slope_from_P50_WD_Ks_limitspec$R2,slope_from_P50_TLP_Ks_limitspec$R2,slope_from_P50_limitspec$R2)
all_R2adj_slope <- c(slope_from_P50_TLP_WD_Ks$R2adj,slope_from_P50_TLP_WD$R2adj,slope_from_P50_TLP_Ks$R2adj,slope_from_P50_WD_Ks$R2adj,slope_from_TLP_WD_Ks$R2adj,slope_from_P50_TLP$R2adj,slope_from_P50_TLPP50$R2adj,slope_from_TLP_TLPP50$R2adj,slope_from_TLP_Ks$R2adj,slope_from_TLP_WD$R2adj,slope_from_TLP$R2adj,slope_from_P50$R2adj,slope_from_WD$R2adj,slope_from_Ks$R2adj,slope_from_P50_WD_Ks_limitspec$R2adj,slope_from_P50_TLP_Ks_limitspec$R2adj,slope_from_P50_limitspec$R2adj)
all_rmse_slope <- c(slope_from_P50_TLP_WD_Ks$rmse,slope_from_P50_TLP_WD$rmse,slope_from_P50_TLP_Ks$rmse,slope_from_P50_WD_Ks$rmse,slope_from_TLP_WD_Ks$rmse,slope_from_P50_TLP$rmse,slope_from_P50_TLPP50$rmse,slope_from_TLP_TLPP50$rmse,slope_from_TLP_Ks$rmse,slope_from_TLP_WD$rmse,slope_from_TLP$rmse,slope_from_P50$rmse,slope_from_WD$rmse,slope_from_Ks$rmse,slope_from_P50_WD_Ks_limitspec$rmse,slope_from_P50_TLP_Ks_limitspec$rmse,slope_from_P50_limitspec$rmse)
all_ndata_slope <- c(slope_from_P50_TLP_WD_Ks$ndata,slope_from_P50_TLP_WD$ndata,slope_from_P50_TLP_Ks$ndata,slope_from_P50_WD_Ks$ndata,slope_from_TLP_WD_Ks$ndata,slope_from_P50_TLP$ndata,slope_from_P50_TLPP50$ndata,slope_from_TLP_TLPP50$ndata,slope_from_TLP_Ks$ndata,slope_from_TLP_WD$ndata,slope_from_TLP$ndata,slope_from_P50$ndata,slope_from_WD$ndata,slope_from_Ks$ndata,slope_from_P50_WD_Ks_limitspec$ndata,slope_from_P50_TLP_Ks_limitspec$ndata,slope_from_P50_limitspec$ndata)

all_slope <- data.frame(all_testnames_slope,all_R2_slope,all_R2adj_slope,all_rmse_slope,all_ndata_slope)
View(all_slope)

# CHOICE: slope_from_P50_TLP_Ks

# Test MAT and PPT coverage of species for chosen model
plot(MAT[slope_from_P50_TLP_Ks$dataused],MAP[slope_from_P50_TLP_Ks$dataused])
# WIDE CLIMATE COVERAGE

# DECISION: slope_from_P50_TLP_Ks


# Optimisation ------------------------------------------------------------
# Attempt to iteratively converge on the best fit values of Ks, TLP, P50 and LMA given known Hmax and LS

# Decide whether to limit the possible ranges of predicted traits to the observed values (T) or not (F)
limitdataranges=T # Currently does not converge in uncertainty propagation if not set to T

# Decide whether to run the uncertainty propagation (T) or not (F)
propagate_uncer=T

# Decide whether to run all trait combinations in the database for LS and Hmax (F), or just a selection (T)
trait_sel=T
# Number of combinations to select if trait_sel=T
n_trait_sel=20

# ---
if (propagate_uncer) {
  n_uncer=nbtstrp
} else {
  n_uncer=1
}

# Identify all combinations of Hmax and LS
ind=which(!is.na(Hmax) & !is.na(LS))

LS_comb <- LS[ind]
Hmax_comb <- Hmax[ind]
LS_Hmax_comb <- LS_Hmax[ind]

if (trait_sel) {
  # Random selection of LS and Hmax values to be tested
  set.seed(1234)
  index = 1:length(ind)
  trait_samp = sample(index, n_trait_sel, replace=F) 
  LS_e = LS_comb[trait_samp] 
  Hmax_e = Hmax_comb[trait_samp] 
  LS_Hmax_e = LS_Hmax_comb[trait_samp] 
  # Plot sample against all data as a check for sampling density
  plot(LS_comb,Hmax_comb)
  points(LS_e,Hmax_e,col="red")
} else {
  # Go through all observed combinations of Hmax and LS
  Hmax_e=Hmax_comb
  LS_e=LS_comb
  LS_Hmax_e=LS_Hmax_comb
}

ndata=length(LS_e)

P50_e <- matrix(NA, nrow= ndata, ncol = n_uncer) #Array now expanded to hold multiple replicate estimates based on regression coefficient uncertainty
LMA_e <- matrix(NA, nrow= ndata, ncol = n_uncer)
TLP_e <- matrix(NA, nrow= ndata, ncol = n_uncer)
Ks_e <- matrix(NA, nrow= ndata, ncol = n_uncer)
WD_e <- matrix(NA, nrow= ndata, ncol = n_uncer)
slope_e <- matrix(NA, nrow= ndata, ncol = n_uncer)

# Loop over all the combinations of Hmax and LS
# The new estimates of traits use the suffix "_e"
for (dd in 1:ndata) {
  print(dd)
  
  # Carry out the optimisation
  opt_vals <- trait_opt(P50,TLP,Ks,LS,LMA,WD,slope,LMA_from_TLP,TLP_from_LS_LMA_P50,P50_from_TLP_Ks,Ks_from_LSHmax_P50,slope_from_P50_TLP_Ks,WD_from_slope_P50slope,Hmax_e[dd],LS_e[dd],LS_Hmax_e[dd],n_uncer)
  
  P50_e[dd,] <- opt_vals$P50_e
  TLP_e[dd,] <- opt_vals$TLP_e
  LMA_e[dd,] <- opt_vals$LMA_e
  Ks_e[dd,] <- opt_vals$Ks_e
  WD_e[dd,] <- opt_vals$WD_e
  slope_e[dd,] <- opt_vals$slope_e
}
  
#Stats defining the uncertainty range for each point
Ks_e_mean=unname(apply(Ks_e, 1, mean))
Ks_e_median=unname(apply(Ks_e, 1, median))
Ks_e_5perc=unname(apply(Ks_e, 1, quantile,0.05))
Ks_e_95perc=unname(apply(Ks_e, 1, quantile,0.95))

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

plot(P50,WD,,pch=16,xlab="P50",ylab="WD",main="P50 vs WD")
points(P50_e[,1],WD_e[,1],col="blue",pch=16) # Using central estimate coefficients
points(P50_e[,1],WD_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
points(P50_e[,1],WD_e_5perc,col="green",pch=16)
points(P50_e[,1],WD_e_95perc,col="green",pch=16)

plot(WD,P50,,pch=16,xlab="WD",ylab="P50",main="WD vs P50")
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

plot(LS_Hmax,Ks,pch=16,xlab="LS*Hmax",ylab="Ks",main="LS_Hmax vs Ks")
points(log(exp(LS_e)*Hmax_e),Ks_e[,1],col="blue",pch=16) # Using central estimate coefficients
points(log(exp(LS_e)*Hmax_e),Ks_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
points(log(exp(LS_e)*Hmax_e),Ks_e_5perc,col="green",pch=16)
points(log(exp(LS_e)*Hmax_e),Ks_e_95perc,col="green",pch=16)

plot(Ks,P50,pch=16,xlab="Ks",ylab="P50",main="Ks vs P50")
points(Ks_e[,1],P50_e[,1],col="blue",pch=16) # Using central estimate coefficients
points(Ks_e[,1],P50_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
points(Ks_e[,1],P50_e_5perc,col="green",pch=16)
points(Ks_e[,1],P50_e_95perc,col="green",pch=16)

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

#plot(Ks,slope,pch=16,xlab="Ks",ylab="slope",main="Ks vs slope")
#points(Ks_e[,1],slope_e[,1],col="blue",pch=16) # Using central estimate coefficients
#points(Ks_e_mean,slope_e_mean,col="red",pch=16) # Using mean of all bootstrapped estimates 
#points(Ks_e_5perc,slope_e_5perc,col="green",pch=16)
#points(Ks_e_95perc,slope_e_95perc,col="green",pch=16)

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

#Ks_res <- Ks[ind]-Ks_e_mean
Ks_res <- Ks[ind]-Ks_e[,1]
Ks_e_rmse <- sqrt(mean(Ks_res^2,na.rm=T))
Ks_e_ndata <- length(which(is.na(Ks_res)==F))

#WD_res <- WD[ind]-WD_e_mean
WD_res <- WD[ind]-WD_e[,1]
WD_e_rmse <- sqrt(mean(WD_res^2,na.rm=T))
WD_e_ndata <- length(which(is.na(WD_res)==F))

#slope_res <- slope[ind]-slope_e_mean
slope_res <- slope[ind]-slope_e[,1]
slope_e_rmse <- sqrt(mean(slope_res^2,na.rm=T))
slope_e_ndata <- length(which(is.na(slope_res)==F))

# Compare to RMSE from best multivarirate regression

all_e_rmse <- c(P50_e_rmse,TLP_e_rmse,LMA_e_rmse,Ks_e_rmse,WD_e_rmse,slope_e_rmse)
all_e_ndata <- c(P50_e_ndata,TLP_e_ndata,LMA_e_ndata,Ks_e_ndata,WD_e_ndata,slope_e_ndata)
all_multivar_rmse <- c(P50_from_TLP_Ks$rmse,TLP_from_LS_LMA_P50$rmse,LMA_from_TLP$rmse,Ks_from_LSHmax_P50$rmse,WD_from_slope_P50slope$rmse,slope_from_P50_TLP_Ks$rmse)
all_multivar_ndata <- c(P50_from_TLP_Ks$ndata,TLP_from_LS_LMA_P50$ndata,LMA_from_TLP$ndata,Ks_from_LSHmax_P50$ndata,WD_from_slope_P50slope$ndata,slope_from_P50_TLP_Ks$ndata)

all_rmse_comp <- data.frame(all_e_rmse,all_e_ndata,all_multivar_rmse,all_multivar_ndata)
View(all_rmse_comp)
