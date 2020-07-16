# Script to read in the processed trait data, make bivariate and multivariate SMA regressions and then
# carry out an optimisation procedure to unify inter-trait relationships across the whole trait dataset.
#
# Dependencies: sma_multivar_regress.R
#
# T. Pugh
# 15.06.20

#setwd('/Users/pughtam/Documents/TreeMort/Analyses/Hydraulic_modelling/Traits/Daijun_processed_data_240520')

traitfile='woody_trait.0625.txt'

source('sma_multivar_regress.R')

#--- Read in the trait data ---

traits=read.table(traitfile)

attach(traits)
detach(traits)
detach(traits[traits$group=="BE",])
detach(traits[traits$group!="CC",])

traitb<-subset(traits,group!="CC",drop = T)
#this is to delete the level that should not exist after subsetting
traitb<-droplevels(traitb)
str(traitb)
attach(traitb)

P50
#TOM CODE
#--- Select a subset of the trait variables to work with ---
#treesonly=T
#groupsel='BE'
#
#group=traits$group
#kstem=traits$Ks[group==groupsel]
#Hmax=traits$Hmax[group==groupsel]
#WD=traits$WD[group==groupsel]
#LS=traits$LS[group==groupsel]
#LMA=traits$LMA[group==groupsel]
#P50=traits$P50[group==groupsel]
#TLP=traits$TLP[group==groupsel]
#slope=traits$slope[group==groupsel]

#if (treesonly) {
#  kstem[Hmax<5 | is.nan(Hmax)]=NA
#  P50[Hmax<5 | is.nan(Hmax)]=NA
#  TLP[Hmax<5 | is.nan(Hmax)]=NA
#  slope[Hmax<5 | is.nan(Hmax)]=NA
#  LS[Hmax<5 | is.nan(Hmax)]=NA
#  Hmax[Hmax<5 | is.nan(Hmax)]=NA
#  LMA[Hmax<5 | is.nan(Hmax)]=NA
#  WD[Hmax<5 | is.nan(Hmax)]=NA
#}
#
# Calculate transformed traits
#P50log=log(-P50,base=exp(1))
#P50log[-P50<0]=NA
#TLPlog=log(-TLP,base=exp(1))
#TLPlog[-TLP<0]=NA
#
#kstemHmax=log(exp(kstem)/Hmax,base=exp(1))
#
##Remove high LMA values (where are they coming from - need a good justification to remove but they are skewing the analysis)
#LMA[LMA>log(600,base=exp(1))]=NA

#--- Bivariate plots with SMA regression ---

library(lmodel2)

par(mfrow=c(4,4))
par(mar=c(2,2,2,2))

mod_TLP_P50 <- lmodel2(P50 ~ TLP)
plot(mod_TLP_P50,"SMA",pch=16,xlab="TLP",ylab="P50",main="TLP vs P50")

mod_TLP_slope <- lmodel2(slope ~ TLP)
plot(mod_TLP_slope,"SMA",pch=16,xlab="TLP",ylab="slope",main="TLP vs slope")

mod_P50_slope <- lmodel2(slope ~ P50)
plot(mod_P50_slope,"SMA",pch=16,xlab="P50",ylab="slope",main="P50 vs slope")

mod_TLP_WD <- lmodel2(WD ~ TLP)
plot(mod_TLP_WD,"SMA",pch=16,xlab="TLP",ylab="WD",main="TLP vs WD")

mod_P50_WD <- lmodel2(WD ~ P50)
plot(mod_P50_WD,"SMA",pch=16,xlab="P50",ylab="WD",main="P50 vs WD")

mod_TLP_LMA <- lmodel2(LMA ~ TLP)
plot(mod_TLP_LMA,"SMA",pch=16,xlab="TLP",ylab="LMA",main="TLP vs LMA")

mod_LS_LMA <- lmodel2(LMA ~ LS)
plot(mod_LS_LMA,"SMA",pch=16,xlab="LS",ylab="LMA",main="LS vs LMA")

mod_Ks_LS_Hmax <- lmodel2(Ks ~ LS_Hmax)
plot(mod_Ks_LS_Hmax,"SMA",pch=16,xlab="LS_Hmax",ylab="Ks",main="LS_Hmax vs Ks")

mod_Ks_P50 <- lmodel2(P50 ~ Ks)
plot(mod_Ks_P50,"SMA",pch=16,xlab="Ks",ylab="P50",main="Ks vs P50")

mod_LS_TLP <- lmodel2(TLP ~ LS)
plot(mod_LS_TLP,"SMA",pch=16,xlab="LS",ylab="TLP",main="LS vs TLP")

mod_WD_LMA <- lmodel2(LMA ~ WD)
plot(mod_WD_LMA,"SMA",pch=16,xlab="WD",ylab="LMA",main="WD vs LMA")

mod_Ks_slope <- lmodel2(slope ~ Ks)
plot(mod_Ks_slope,"SMA",pch=16,xlab="Ks",ylab="slope",main="Ks vs slope")

mod_Ks_Hmax <- lmodel2(Hmax ~ Ks)
plot(mod_Ks_Hmax,"SMA",pch=16,xlab="Ks",ylab="Hmax",main="Ks vs Hmax")

mod_leafN_LMA <- lmodel2(LMA ~ leafN)
plot(mod_leafN_LMA,"SMA",pch=16,xlab="leafN",ylab="LMA",main="leafN vs LMA")

mod_Ks_Hmax_LS <- lmodel2(LS ~ Ks_Hmax)
plot(mod_Ks_Hmax_LS,"SMA",pch=16,xlab="Ks_Hmax",ylab="LS",main="Ks_Hmax vs LS")


#Use up remaining unallocated plots and set back to single plot
#plot.new()
#plot.new()
#plot.new()
par(mfrow=c(1,1))
par(mar=c(5.1,4.1,4.1,2.1))

#--- Partial correlations for linked traits ---

library(ppcor)

#LMA, LS, WD, TLP
ind=which(!is.na(LMA) & !is.na(LS) & !is.na(WD) & !is.na(TLP))
pcor_LMA <- pcor(data.frame(TLP[ind],LS[ind],WD[ind],LMA[ind]))

#TLP, LS, WD, LMA, P50 (do not include slope due to lack of data)
ind=which(!is.na(LMA) & !is.na(LS) & !is.na(WD) & !is.na(TLP) & !is.na(P50))
pcor_TLP <- pcor(data.frame(LMA[ind],LS[ind],WD[ind],P50[ind],TLP[ind]))

#P50, TLP, kstem (do not include slope due to lack of data)
ind=which(!is.na(P50) & !is.na(TLP) & !is.na(Ks)) 
pcor_P50 <- pcor(data.frame(P50[ind],TLP[ind],Ks[ind]))

#slope, P50, TLP, kstem
ind=which(!is.na(P50) & !is.na(TLP) & !is.na(slope) & !is.na(Ks)) 
pcor_slope <- pcor(data.frame(P50[ind],TLP[ind],Ks[ind],slope[ind]))

#WD, LMA, P50
ind=which(!is.na(P50) & !is.na(LMA) & !is.na(TLP)& !is.na(WD)) 
pcor_WD <- pcor(data.frame(P50[ind],LMA[ind],TLP[ind],WD[ind]))

# How many species do we have all traits for?
#Including slope
length(which(!is.na(LMA) & !is.na(LS) & !is.na(WD) & !is.na(TLP) & !is.na(P50)& !is.na(Ks)& !is.na(Hmax) & !is.na(slope)))
#Excluding slope
length(which(!is.na(LMA) & !is.na(LS) & !is.na(WD) & !is.na(TLP) & !is.na(P50)& !is.na(Ks)& !is.na(Hmax)))


#--- Multivariate SMA models ---

#Estimate P50 from TLP and kstem
ind=which(!is.na(P50) & !is.na(TLP) & !is.na(Ks))
mod_P50 <- sma_regress_multivar(data.frame(TLP[ind],Ks[ind],P50[ind]))
P50_est <- mod_P50$intercept_R + mod_P50$slope_R.y1*TLP[ind] + mod_P50$slope_R.y2*Ks[ind]
#Cross-checks to make sure that result is sensible
plot(mod_TLP_P50,"SMA",pch=16,xlab="TLP",ylab="P50",main="TLP vs P50")
points(TLP[ind],P50_est,col="red",pch=16)
plot(mod_Ks_P50,"SMA",pch=16,xlab="Ks",ylab="P50",main="Ks vs P50")
points(Ks[ind],P50_est,col="red",pch=16)

#Estimate LMA from LS and TLP
#note: My test is good frr the LS and consider to not use
ind=which(!is.na(LS) & !is.na(TLP) & !is.na(LMA))
mod_LMA <- sma_regress_multivar(data.frame(TLP[ind],LS[ind],LMA[ind]))
LMA_est <- mod_LMA$intercept_R + mod_LMA$slope_R.y1*TLP[ind] + mod_LMA$slope_R.y2*LS[ind]
#Cross-checks to make sure that result is sensible
plot(mod_LS_LMA,"SMA",pch=16,xlab="LS",ylab="LMA",main="LS vs LMA")
points(LS[ind],LMA_est,col="red",pch=16)
plot(mod_TLP_LMA,"SMA",pch=16,xlab="TLP",ylab="LMA",main="TLP vs LMA")
points(TLP[ind],LMA_est,col="red",pch=16)

#Estimate WD from P50, TLP and LMA
ind=which(!is.na(P50) & !is.na(TLP) & !is.na(LMA) & !is.na(WD))
mod_WD <- sma_regress_multivar(data.frame(TLP[ind],P50[ind],LMA[ind],WD[ind]))
WD_est <- mod_WD$intercept_R + mod_WD$slope_R.y1*TLP[ind] + mod_WD$slope_R.y2*P50[ind] + mod_WD$slope_R.y3*LMA[ind]
#Cross-checks to make sure that result is sensible
plot(mod_TLP_WD,"SMA",pch=16,xlab="TLP",ylab="WD",main="TLP vs WD")
points(TLP[ind],WD_est,col="red",pch=16)
plot(mod_P50_WD,"SMA",pch=16,xlab="P50",ylab="WD",main="P50 vs WD")
points(P50[ind],WD_est,col="red",pch=16)
plot(mod_WD_LMA,"SMA",pch=16,xlab="WD",ylab="LMA",main="WD vs LMA")
points(WD_est,LMA[ind],col="red",pch=16)

#Estimate TLP from P50,LS,LMA and WD
ind=which(!is.na(P50) & !is.na(TLP) & !is.na(LMA) & !is.na(LS)& !is.na(WD))
mod_TLP <- sma_regress_multivar(data.frame(P50[ind],LS[ind],LMA[ind],TLP[ind],TLP[ind]))
TLP_est <- mod_TLP$intercept_R + mod_TLP$slope_R.y1*P50[ind] + mod_TLP$slope_R.y2*LS[ind] + mod_TLP$slope_R.y3*LMA[ind]+mod_TLP$slope_R.y4*WD[ind]
#Cross-checks to make sure that result is sensible
plot(mod_TLP_P50,"SMA",pch=16,xlab="TLP",ylab="P50",main="TLP vs P50")
points(TLP_est,P50[ind],col="red",pch=16)
plot(mod_LS_TLP,"SMA",pch=16,xlab="LS",ylab="TLP",main="LS vs TLP")
points(LS[ind],TLP_est,col="red",pch=16)
plot(mod_TLP_LMA,"SMA",pch=16,xlab="TLP",ylab="LMA",main="TLP vs LMA")
points(TLP_est,LMA[ind],col="red",pch=16)
plot(mod_TLP_WD,"SMA",pch=16,xlab="TLP",ylab="WD",main="TLP vs WD")
points(TLP_est,WD[ind],col="red",pch=16)
#why the regression forr TLP and WD is so good? sel-correlated?

#Estimate slope from TLP, P50 and kstem
ind=which(!is.na(P50) & !is.na(TLP) & !is.na(Ks) & !is.na(slope))
mod_slope <- sma_regress_multivar(data.frame(P50[ind],TLP[ind],Ks[ind],slope[ind]))
slope_est <- mod_slope$intercept_R + mod_slope$slope_R.y1*P50[ind] + mod_slope$slope_R.y2*TLP[ind] + mod_slope$slope_R.y3*Ks[ind]
#Cross-checks to make sure that result is sensible
plot(mod_Ks_slope,"SMA",pch=16,xlab="kstem",ylab="slope",main="kstem vs slope")
points(Ks[ind],slope_est,col="red",pch=16)
plot(mod_TLP_slope,"SMA",pch=16,xlab="TLP",ylab="slope",main="TLP vs slope")
points(TLP[ind],slope_est,col="red",pch=16)
plot(mod_P50_slope,"SMA",pch=16,xlab="P50",ylab="slope",main="P50 vs slope")
points(P50[ind],slope_est,col="red",pch=16)


#Estimate LS from LMA and kstem/Hmax.
ind=which(!is.na(LMA) & !is.na(Ks_Hmax) & !is.na(LS))
mod_LS <- sma_regress_multivar(data.frame(LMA[ind],Ks_Hmax[ind],LS[ind]))
LS_est <- mod_LS$intercept_R + mod_LS$slope_R.y1*LMA[ind] + mod_LS$slope_R.y2*Ks_Hmax[ind]
#Cross-checks to make sure that result is sensible
plot(mod_Ks_Hmax_LS,"SMA",pch=16,xlab="Ks_Hmax",ylab="LS",main="Ks_Hmax vs LS")
points(Ks_Hmax[ind],LS_est,col="red",pch=16)
plot(mod_LS_LMA,"SMA",pch=16,xlab="LS",ylab="LMA",main="LS vs LMA")
points(LS_est,LMA[ind],col="red",pch=16)





#--- Attempt to iteratively converge on the best fit values of LS, TLP, P50 and LMA given known Hmax and kstem ---

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

#Go through all observed combinations of Hmax and kstem

ind=which(!is.na(Hmax) & !is.na(Ks))

tol=0.00001

Hmax_e=Hmax[ind]
Ks_e=Ks[ind]
#Hmax_e=10
#kstem_e=-1.5

ndata=length(Ks_e)

LS_e <- matrix(NA, nrow = ndata)
P50_e <- matrix(NA, nrow = ndata)
LMA_e <- matrix(NA, nrow = ndata)
TLP_e <- matrix(NA, nrow = ndata)
WD_e <- matrix(NA, nrow = ndata)
slope_e <- matrix(NA, nrow = ndata)

optLS=F
limitdataranges=F

for (dd in 1:ndata) {

  Ks_Hmax_e=log(exp(Ks_e[dd])/Hmax_e[dd],base=exp(1))
  
  if (!optLS) {
    LS_e[dd]= mod_Ks_Hmax_LS$regression.results$Slope[3]*Ks_Hmax_e +
    mod_Ks_Hmax_LS$regression.results$Intercept[3]
    
    if (limitdataranges) {
      #Do not go beyond observed limits of data
      if (LS_e[dd]>maxLS) {LS_e[dd]=maxLS}
      if (LS_e[dd]<minLS) {LS_e[dd]=minLS}
    }
  }

  #TLP, P50, LMA (and possibly LS) need optimising

  #Initialise based on simple bivariate relationship
  if (optLS) {
    LS_e_last = mod_Ks_LS_Hmax$regression.results$Slope[3]*Ks_Hmax_e +
      mod_Ks_LS_Hmax$regression.results$Intercept[3]
    LMA_e_last = mod_LS_LMA$regression.results$Slope[3]*LS_e_last +
      mod_LS_LMA$regression.results$Intercept[3]
    TLP_e_last = mod_LS_TLP$regression.results$Slope[3]*LS_e_last +
      mod_LS_TLP$regression.results$Intercept[3]
  }
  else {
    LMA_e_last = mod_LS_LMA$regression.results$Slope[3]*LS_e[dd] +
      mod_LS_LMA$regression.results$Intercept[3]
    TLP_e_last = mod_LS_TLP$regression.results$Slope[3]*LS_e[dd] +
      mod_LS_TLP$regression.results$Intercept[3]
  }
  P50_e_last = mod_TLP_P50$regression.results$Slope[3]*TLP_e_last +
    mod_TLP_P50$regression.results$Intercept[3]

  #Set diff array very high initially
  if (optLS) {
    diff_LS_last=100
  }
  diff_P50_last=100
  diff_LMA_last=100
  diff_TLP_last=100

  #LS_c <- matrix(NA, nrow = 100)
  #P50_c <- matrix(NA, nrow = 100)
  #LMA_c <- matrix(NA, nrow = 100)
  #TLP_c <- matrix(NA, nrow = 100)
  
  niter=0;
  while (T) {
    niter=niter+1

    if (optLS) {
      LS_e[dd]=mod_LS$intercept_R + mod_LS$slope_R.y1*LMA_e_last + mod_LS$slope_R.y2*LS_Hmax_e
      LMA_e[dd]=mod_LMA$intercept_R + mod_LMA$slope_R.y1*TLP_e_last + mod_LMA$slope_R.y2*LS_e_last
      TLP_e[dd]=mod_TLP$intercept_R + mod_TLP$slope_R.y1*P50_e_last + mod_TLP$slope_R.y2*LS_e_last + mod_TLP$slope_R.y3*LMA_e_last+mod_TLP$slope_R.y4*WD_e_last
    }
    else {
      LMA_e[dd]=mod_LMA$intercept_R + mod_LMA$slope_R.y1*TLP_e_last + mod_LMA$slope_R.y2*LS_e[dd]
      TLP_e[dd]=mod_TLP$intercept_R + mod_TLP$slope_R.y1*P50_e_last + mod_TLP$slope_R.y2*LS_e[dd] + mod_TLP$slope_R.y3*LMA_e_last+mod_TLP$slope_R.y4*WD_e_last
    }
    P50_e[dd]=mod_P50$intercept_R + mod_P50$slope_R.y1*TLP_e_last + mod_P50$slope_R.y2*Ks_e[dd]
    
    
    if (limitdataranges) {
      #Do not go beyond observed limits of data
      if (P50_e[dd]>maxP50 | is.na(P50_e[dd])) {P50_e[dd]=NA; break}
      if (P50_e[dd]<minP50 | is.na(P50_e[dd])) {P50_e[dd]=NA; break}
      if (TLP_e[dd]>maxTLP | is.na(TLP_e[dd])) {TLP_e[dd]=NA; break}
      if (TLP_e[dd]<minTLP | is.na(TLP_e[dd])) {TLP_e[dd]=NA; break}
      if (LMA_e[dd]>maxLMA | is.na(LMA_e[dd])) {LMA_e[dd]=NA; break}
      if (LMA_e[dd]<minLMA | is.na(LMA_e[dd])) {LMA_e[dd]=NA; break}
      if (optLS) {
        if (LS_e[dd]>maxLS | is.na(LS_e[dd])) {LS_e[dd]=NA; break}
        if (LS_e[dd]<minLS | is.na(LS_e[dd])) {LS_e[dd]=NA; break}
      }
    }
    
    #P50_c[niter] <- P50_e[dd]
    #LMA_c[niter] <- LMA_e[dd]
    #TLP_c[niter] <- TLP_e[dd]
    #LS_c[niter] <- LS_e[dd]

    if (optLS) {
      diff_LS <- LS_e[dd]-LS_e_last
    }
    diff_P50 = P50_e[dd]-P50_e_last
    diff_LMA = LMA_e[dd]-LMA_e_last
    diff_TLP = TLP_e[dd]-TLP_e_last

    if (optLS) {
      if (abs(diff_P50-diff_P50_last)<tol &&
        abs(diff_LMA-diff_LMA_last)<tol &&
        abs(diff_TLP-diff_TLP_last)<tol &&
        abs(diff_LS-diff_LS_last)<tol) {
        break
      }
    }
    else {
      if (abs(diff_P50-diff_P50_last)<tol &&
        abs(diff_LMA-diff_LMA_last)<tol &&
        abs(diff_TLP-diff_TLP_last)<tol) {
        break
      }
    }

    if (optLS) {
      diff_LS_last=diff_LS
    }
    diff_P50_last=diff_P50
    diff_LMA_last=diff_LMA
    diff_TLP_last=diff_TLP

    if (optLS) {
      LS_e_last=LS_e[dd]
    }
    P50_e_last=P50_e[dd]
    LMA_e_last=LMA_e[dd]
    TLP_e_last=TLP_e[dd]
  }

  WD_e[dd]=mod_WD$intercept_R + mod_WD$slope_R.y1*TLP_e[dd] + mod_WD$slope_R.y2*P50_e[dd] + mod_WD$slope_R.y3*LMA_e[dd]
  slope_e[dd]=mod_slope$intercept_R + mod_slope$slope_R.y1*P50_e[dd] + mod_slope$slope_R.y2*TLP_e[dd] + mod_slope$slope_R.y3*Ks[dd]
  
  if (limitdataranges) {
    #Do not go beyond observed limits of data
    if (WD_e[dd]>maxWD | is.na(WD_e[dd])) {WD_e[dd]=NA}
    if (WD_e[dd]<minWD | is.na(WD_e[dd])) {WD_e[dd]=NA}
    if (slope_e[dd]>maxslope | is.na(slope_e[dd])) {slope_e[dd]=NA}
    if (slope_e[dd]<minslope | is.na(slope_e[dd])) {slope_e[dd]=NA}
  }

}

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
