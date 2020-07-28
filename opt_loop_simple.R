# Extraction of just the optimisation loop from trait_rel.R. The code here is designed for playing with
# offline and only takes 1 value each of Hmax and LS as inputs.
#
# Note that the relevant correlations have to be calculated using trait_rel.R before running this code.
#
# Attempt to iteratively converge on the best fit values of Ks, TLP, P50 and LMA given known Hmax and LS

# Set Hmax and LS values
Hmax_e=40 #Not logged
LS_e=0 #Logged

# Decide whether to limit the possible ranges of predicted traits to the observed values (T) or not (F)
limitdataranges=F

#--- Calculations start here ---

if (limitdataranges) {
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
}

# Set the tolerance level for the iteration
tol=0.00001

# The new estimates of traits use the suffix "_e"
  
#Calculate the value for LS*Hmax, which comes direct from the two input traits
LS_Hmax_e=log(exp(LS_e)*Hmax_e,base=exp(1))
  
#Calculate the values of Ks based only on the bivariate relationships with LS*Hmax
Ks_e= mod_Ks_LS_Hmax$regression.results$Slope[3]*LS_Hmax_e +
    mod_Ks_LS_Hmax$regression.results$Intercept[3]
  
if (limitdataranges) {
  #Do not go beyond observed limits of data
  if (Ks_e[dd]>maxKs) {Ks_e[dd]=maxKs}
  if (Ks_e[dd]<minKs) {Ks_e[dd]=minKs}
}
  
#TLP, P50, LMA, WD need optimising
  
#First set some initial based on simple bivariate relationship. This is just so that the iteration has somewhere to start from. Final result should not be sensitive to these.
LMA_e_last = mod_LS_LMA$regression.results$Slope[3]*LS_e +
  mod_LS_LMA$regression.results$Intercept[3]
TLP_e_last = mod_LS_TLP$regression.results$Slope[3]*LS_e +
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
# Here we initialise the "diff_*_last" variables very high, could be any large value.
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
  LMA_e=mod_LMA$intercept_R + mod_LMA$slope_R.y1*TLP_e_last + mod_LMA$slope_R.y2*LS_e
  TLP_e=mod_TLP$intercept_R + mod_TLP$slope_R.y1*P50_e_last + mod_TLP$slope_R.y2*LMA_e_last+mod_TLP$slope_R.y3*WD_e_last
  P50_e=mod_P50$intercept_R + mod_P50$slope_R.y1*TLP_e_last + mod_P50$slope_R.y2*Ks_e
  #WD_e=mod_WD$intercept_R + mod_WD$slope_R.y1*TLP_e_last + mod_WD$slope_R.y2*P50_e_last + mod_WD$slope_R.y3*LMA_e_last
  # Test taking WD from a simple bivariate relationship (because not converging with the above multivariate one - need to sort out the multivariate fit!)
  WD_e = mod_TLP_WD$regression.results$Slope[3]*TLP_e_last + mod_TLP_WD$regression.results$Intercept[3]
    
  if (limitdataranges) {
    #Do not go beyond observed limits of data
    if (P50_e>maxP50 | is.na(P50_e)) {P50_e=NA; break}
    if (P50_e<minP50 | is.na(P50_e)) {P50_e=NA; break}
    if (TLP_e>maxTLP | is.na(TLP_e)) {TLP_e=NA; break}
    if (TLP_e<minTLP | is.na(TLP_e)) {TLP_e=NA; break}
    if (LMA_e>maxLMA | is.na(LMA_e)) {LMA_e=NA; break}
    if (LMA_e<minLMA | is.na(LMA_e)) {LMA_e=NA; break}
    if (WD_e>maxWD | is.na(WD_e)) {WD_e=NA; break}
    if (WD_e<minWD | is.na(WD_e)) {WD_e=NA; break}
  }
    
  # Save the values for this iteration to the output array (only for debugging, can be commented out)
  P50_c[niter] <- P50_e
  LMA_c[niter] <- LMA_e
  TLP_c[niter] <- TLP_e
  WD_c[niter] <- WD_e
    
  # Calculate the difference between the current estimate of a trait value "_e" and the previous estimate "_last"
  diff_P50 = P50_e-P50_e_last
  diff_LMA = LMA_e-LMA_e_last
  diff_TLP = TLP_e-TLP_e_last
  diff_WD = WD_e-WD_e_last
    
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
  P50_e_last=P50_e
  LMA_e_last=LMA_e
  TLP_e_last=TLP_e
  WD_e_last=WD_e
}
  
# After the iteration has finished we can calculate any traits which did not need to be included in the optimisation (because they are not used in the input to calculate any other trait)
slope_e=mod_slope$intercept_R + mod_slope$slope_R.y1*P50_e + mod_slope$slope_R.y2*TLP_e + mod_slope$slope_R.y3*Ks_e
  
if (limitdataranges) {
  #Do not go beyond observed limits of data
  if (slope_e>maxslope | is.na(slope_e)) {slope_e=NA}
  if (slope_e<minslope | is.na(slope_e)) {slope_e=NA}
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

#We can plot the values of the traits on each iteration to see how the optimisation converges
plot(P50_c)
plot(TLP_c)
plot(LMA_c)
plot(WD_c)
