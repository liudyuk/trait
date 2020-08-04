# Convert trait database values to the values/variables needed in LPJ-GUESS
#
# T. Pugh
# 04.08.20

SLA_to_leaflong <- function(SLA_LPJG) {
  # Calculate leaf longevity based on equation for SLA from leaf longevsity in LPJ-GUESS v4.0
  leaflong_exponent <- (log10(SLA_LPJG/0.2)-2.41)/-0.38
  leaflong_raw <- (10^leaflong_exponent)/12
  # Limit values to a maximum of 10 years
  leaflong <- matrix(NA,ncol=length(leaflong_raw))
  for (nn in 1:length(leaflong_raw)) {
    leaflong[nn] <- max(min(leaflong_raw[nn],10),0.5)
  }
  return(leaflong)
}

lpjg_traits_conv <- function(LMA_e_mean,P50_e_mean,TLP_e_mean,slope_e_mean,
                             LS_e,WD_e_mean,Ks_e_mean) {
  # Unlog traits
  LMA_e_mean_unlogged <- exp(LMA_e_mean)
  P50_e_mean_unlogged <- -exp(P50_e_mean)
  TLP_e_mean_unlogged <- -exp(TLP_e_mean)
  slope_e_mean_unlogged <- exp(slope_e_mean)
  LS_e_mean_unlogged <- exp(LS_e)
  Ks_e_mean_unlogged <- exp(Ks_e_mean)
  
  # Unit conversions or analytical transformations
  WD_LPJG <- (WD_e_mean*1000)/2 # From g cm-3 to kgC m-3
  SLA_LPJG <- (1/LMA_e_mean_unlogged)*1000/2 # From log g m-2 to m2 kgC-1
  P88_LPJG <- P50_e_mean_unlogged - (38/slope_e_mean_unlogged)
  polyslope_LPJG <- 2/(log(P50_e_mean_unlogged/P88_LPJG)) # Hill/polynomial slope for input into LPJ-GUESS using equation from Phillip Papastefanou
  LS_LPJG <- LS_e_mean_unlogged*10000 # From m2 (leaf) cm-2 (sap) to m2 (leaf) m-2 (sap)
  
  # Empirical transformations
  leaflong_LPJG <- SLA_to_leaflong(SLA_LPJG)
  lambda_LPJG <- -0.188+(-0.3*TLP_e_mean_unlogged) # Based on https://docs.google.com/spreadsheets/d/1KvB97vt6OVdy3GPUdOdlxd8c9TdjEk9z4qIl2mcfPHk/edit#gid=51069844
  DeltaPsiWW_LPJG <- -0.5571 + (2.9748*WD_e_mean) # Based on https://docs.google.com/spreadsheets/d/1KvB97vt6OVdy3GPUdOdlxd8c9TdjEk9z4qIl2mcfPHk/edit#gid=51069844
  Kleaf_LPJG <- 3.3682 + (-1.21*TLP_e_mean_unlogged)
  
  # No transformation needed
  Ks_LPJG <- Ks_e_mean_unlogged
  P50_LPJG <- P50_e_mean_unlogged
  
  traits_LPJG <- list("WD"=WD_LPJG,"SLA"=SLA_LPJG,"P50"=P50_LPJG,"P88"=P88_LPJG,"polyslope"=polyslope_LPJG,"LS"=LS_LPJG,"leaflong"=leaflong_LPJG,"lambda"=lambda_LPJG,"DeltaPsiWW"=DeltaPsiWW_LPJG,"Ks"=Ks_LPJG,"Kleaf"=Kleaf_LPJG)
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
                 "\t turnover_leaf 0.33",
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

