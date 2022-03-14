# Convert trait database and estimated values to the values/variables needed in LPJ-GUESS
#
# T. Pugh
# 04.08.20

SLA_to_leaflong <- function(SLA_LPJG) {
  # Calculate leaf longevity based on equation for SLA from leaf longevity in LPJ-GUESS v4.0
  leaflong_exponent <- (log10(SLA_LPJG/0.2)-2.41)/-0.38
  leaflong_raw <- (10^leaflong_exponent)/12
  # Limit values to a maximum of 10 years
  leaflong <- matrix(NA,ncol=length(leaflong_raw))
  for (nn in 1:length(leaflong_raw)) {
    leaflong[nn] <- max(min(leaflong_raw[nn],10),0.5)
  }
  leaflong <- t(leaflong) # Rotate array for consistency with other arrays
  return(leaflong)
}

lpjg_traits_conv <- function(LMA_e_mean,P50_e_mean,TLP_e_mean,slope_e_mean,
                             LS_e,WD_e_mean,Ks_e,
                             leafL_from_LMA=NULL,
                             leafN_from_LMA=NULL,
                             leafN_from_LMA_limit=NULL) {
  # Unlog traits
  LMA_e_mean_unlogged <- exp(LMA_e_mean)
  P50_e_mean_unlogged <- -exp(P50_e_mean)
  TLP_e_mean_unlogged <- -exp(TLP_e_mean)
  slope_e_mean_unlogged <- exp(slope_e_mean)
  LS_e_mean_unlogged <- exp(LS_e)
  Ks_e_mean_unlogged <- exp(Ks_e)
  
  # Unit conversions or analytical transformations
  WD_LPJG <- (WD_e_mean*1000)/2 # From g cm-3 to kgC m-3
  SLA_LPJG <- (1/LMA_e_mean_unlogged)*1000*2 # From log g m-2 to m2 kgC-1
  P88_LPJG <- P50_e_mean_unlogged - (38/slope_e_mean_unlogged)
  polyslope_LPJG <- 2/(log(P50_e_mean_unlogged/P88_LPJG)) # Hill/polynomial slope for input into LPJ-GUESS using equation from Phillip Papastefanou
  LS_LPJG <- LS_e_mean_unlogged*10000 # From m2 (leaf) cm-2 (sap) to m2 (leaf) m-2 (sap)
  
  # Empirical transformations
  if (is.null(leafL_from_LMA)) {
    #Calculate leaf longevity using equation taken from LPJ-GUESS
    leaflong_LPJG <- SLA_to_leaflong(SLA_LPJG)
  } else {
    #Base the calculation on dataset herein
    leaflong_LPJG <- exp(leafL_from_LMA$mod$intercept_R + leafL_from_LMA$mod$slope_R.y1*LMA_e_mean)
    #Apply a minimum of 0.5 for LPJ-GUESS
    #leaflong_LPJG <- pmax(leaflong_LPJG,0.5)
  }
  
  leafN_LPJG <- exp(leafN_from_LMA$mod$intercept_R + leafN_from_LMA$mod$slope_R.y1*LMA_e_mean) # mg N per g leaf
  CtoN_LPJG <- 1000/(leafN_LPJG*2) #Assume 0.5 g C per g leaf
  
  leafNmax_LPJG <- exp(leafN_from_LMA_limit$intercept_upper + leafN_from_LMA_limit$slope*LMA_e_mean) # mg N per g leaf
  CtoNmin_LPJG <- 1000/(leafNmax_LPJG*2) #Assume 0.5 g C per g leaf
  
  lambda_LPJG <- -0.188+(-0.3*TLP_e_mean_unlogged) # Based on https://docs.google.com/spreadsheets/d/1KvB97vt6OVdy3GPUdOdlxd8c9TdjEk9z4qIl2mcfPHk/edit#gid=51069844
  DeltaPsiWW_LPJG <- -0.5571 + (2.9748*WD_e_mean) # Based on https://docs.google.com/spreadsheets/d/1KvB97vt6OVdy3GPUdOdlxd8c9TdjEk9z4qIl2mcfPHk/edit#gid=51069844
  Kleaf_LPJG <- 3.3682 + (-1.21*TLP_e_mean_unlogged) #### [TODO] Check relationship between TLP and Kleaf using other regression methods
  
  #Limit ranges
  lambda_LPJG[lambda_LPJG  < -0.3] = -0.3 # From Papastefanou et al. (2020, Front. Plant Sci.)
  lambda_LPJG[lambda_LPJG  > 1.0]  = 1.0  # From Papastefanou et al. (2020, Front. Plant Sci.)
  DeltaPsiWW_LPJG[DeltaPsiWW_LPJG < 0.3]  = 0.3 # From Papastefanou et al. (2020, Front. Plant Sci.)
  DeltaPsiWW_LPJG[DeltaPsiWW_LPJG > 10.0] = 10.0
  
  # No transformation needed
  Ks_LPJG <- Ks_e_mean_unlogged
  P50_LPJG <- P50_e_mean_unlogged
  
  traits_LPJG <- list("WD"=round(WD_LPJG,digits=4),"SLA"=round(SLA_LPJG,digits=4),"P50"=round(P50_LPJG,digits=4),
                      "P88"=round(P88_LPJG,digits=4),"polyslope"=round(polyslope_LPJG,digits=4),"LS"=round(LS_LPJG,digits=4),
                      "leaflong"=round(leaflong_LPJG,digits=4),"lambda"=round(lambda_LPJG,digits=4),
                      "DeltaPsiWW"=round(DeltaPsiWW_LPJG,digits=4),"Ks"=round(Ks_LPJG,digits=4),
                      "Kleaf"=round(Kleaf_LPJG,digits=4),"TLP"=round(TLP_e_mean_unlogged,digits=4),
                      "slope"=round(slope_e_mean_unlogged,digits=4),"leafN_LPJG"=round(leafN_LPJG,digits=4),
                      "CtoN_LPJG"=round(CtoN_LPJG,digits=4),"CtoNmin_LPJG"=round(CtoNmin_LPJG,digits=4))
  return(traits_LPJG)
}


# LPJ-GUESS instruction file information ----------------------------------

TeBE_header <- c("\t ! Temperate broadleaved evergreen tree",
                 "",
                 "\t include 1",
                 "\t tree",
                 "\t broadleaved",
                 "\t shade_tolerant",
                 "\t evergreen",
                 "\t temperate",
                 "\t !leaflong 3",
                 "\t !turnover_leaf 0.33",
                 "\t tcmin_surv -1",
                 "\t tcmin_est 0",
                 "\t tcmax_est 18.8",
                 "\t twmin_est 5",
                 "\t gdd5min_est 2000",
                 "\t longevity 300 !from TH 2010-04-07 was 350 AA",
                 "\t fireresist 0.3",
                 "\t eps_iso 24.0",
                 "\t seas_iso 0",
                 "\t eps_mon 1.6",
                 "\t storfrac_mon 0.",
                 "")

TeBS_header <- c(  "\t ! Shade-tolerant temperate broadleaved summergreen tree",
                   "",
                   "\t include 1",
                   "\t tree",
                   "\t broadleaved",
                   "\t shade_tolerant",
                   "\t summergreen",
                   "\t temperate",
                   "\t tcmin_surv -14",
                   "\t tcmin_est -13",
                   "\t tcmax_est 6",
                   "\t twmin_est 5",
                   "\t gdd5min_est 1100",
                   "\t longevity 400",
                   "\t fireresist 0.1",
                   "\t eps_iso 45.0",
                   "\t seas_iso 1",
                   "\t eps_mon 1.6",
                   "\t storfrac_mon 0.",
                   "")

IBS_header <- c(  "\t ! Shade-intolerant broadleaved summergreen tree",
                  "",
                  "\t include 1",
                  "\t tree",
                  "\t broadleaved",
                  "\t shade_intolerant",
                  "\t summergreen",
                  "\t boreal",
                  "\t tcmin_surv -30",
                  "tcmin_est -30",
                  "\t tcmax_est 7",
                  "\t twmin_est -1000  !no limit",
                  "\t gdd5min_est 350 !from TH 2010-03-10 AA",
                  "\t longevity 300 !from TH 2010-04-07 was 300 AA",
                  "\t fireresist 0.1",
                  "\t eps_iso 45.0",
                  "\t seas_iso 1",
                  "\t eps_mon 1.6",
                  "\t storfrac_mon 0.",
                  "")

TrBE_header <- c(  "\t ! Tropical broadleaved evergreen tree",
                   "",
                   "\t include 1",
                   "\t tree",
                   "\t broadleaved",
                   "\t shade_tolerant",
                   "\t evergreen",
                   "\t tropical",
                   "\t !leaflong 2",
                   "\t !turnover_leaf 0.5",
                   "\t longevity 500   !from Thomas H 2010-03-30 new 500 instead of 600 2010-04-07",
                   "\t fireresist 0.1",
                   "\t eps_iso 24.0",
                   "\t seas_iso 0",
                   "\t eps_mon 0.8",
                   "\t storfrac_mon 0.",
                   "")

TrBR_header <- c(  "\t ! Tropical broadleaved raingreen tree",
                   "",
                   "\t include 1",
                   "\t tree",
                   "\t broadleaved",
                   "\t shade_intolerant",
                   "\t tropical",
                   "\t phenology \"raingreen\"",
                   "\t fnstorage 0.15",
                   "\t !leaflong 0.5",
                   "\t !turnover_leaf 1",
                   "\t longevity 400    ! from Thomas h 2010-03-30",
                   "\t fireresist 0.3",
                   "\t eps_iso 45.0",
                   "\t seas_iso 0",
                   "\t eps_mon 2.4",
                   "\t storfrac_mon 0.",
                   "")
