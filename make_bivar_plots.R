make_bivar_plots <- function(trait,nbtstrp) {
  # Make the bivariate SMA regressions necessary to populate the hypothesis framework
  #
  # Dependencies: trait_functions.R
  #
  # T. Pugh
  # 12.12.20
  
  par(mfrow=c(4,4))
  par(mar=c(2,2,2,2))
  
  TLP_from_P50 <- sma_plot_stats(data.frame(trait$P50,trait$TLP),c("P50","TLP"),nbtstrp,T)
  
  slope_from_TLP <- sma_plot_stats(data.frame(trait$TLP,trait$slope),c("TLP","slope"),nbtstrp,T)
  
  slope_from_P50 <- sma_plot_stats(data.frame(trait$P50,trait$slope),c("P50","slope"),nbtstrp,T)
  
  TLP_from_WD <- sma_plot_stats(data.frame(trait$WD,trait$TLP),c("WD","TLP"),nbtstrp,T)
  
  WD_from_P50 <- sma_plot_stats(data.frame(trait$P50,trait$WD),c("P50","WD"),nbtstrp,T)
  
  WD_from_slope <- sma_plot_stats(data.frame(trait$slope,trait$WD),c("slope","WD"),nbtstrp,T)
  
  LMA_from_TLP <- sma_plot_stats(data.frame(trait$TLP,trait$LMA),c("TLP","LMA"),nbtstrp,T)
  
  LS_from_LMA <- sma_plot_stats(data.frame(trait$LMA,trait$LS),c("LMA","LS"),nbtstrp,T)
  
  P50_from_Ks <- sma_plot_stats(data.frame(trait$Ks,trait$P50),c("Ks","P50"),nbtstrp,T)
  
  P50_from_LS <- sma_plot_stats(data.frame(trait$LS,trait$P50),c("LS","P50"),nbtstrp,T)
  
  LS_from_TLP <- sma_plot_stats(data.frame(trait$TLP,trait$LS),c("TLP","LS"),nbtstrp,T)
  
  WD_from_LMA <- sma_plot_stats(data.frame(trait$LMA,trait$WD),c("LMA","WD"),nbtstrp,T)
  
  slope_from_Ks <- sma_plot_stats(data.frame(trait$Ks,trait$slope),c("Ks","slope"),nbtstrp,T)
  
  # Set back to single plot
  par(mfrow=c(1,1))
  par(mar=c(5.1,4.1,4.1,2.1))
  
  # Make a data frame summarising the fits of the regressions
  all_label1 <- c("P50","TLP",   "P50", "WD", "P50","slope","LMA","LMA","Ks", "LS", "TLP","LMA","Ks")
  all_label2 <- c("TLP","slope","slope","TLP","WD", "WD",   "TLP","LS", "P50","P50","LS", "WD", "slope")
  all_R <- c(TLP_from_P50$R,
              slope_from_TLP$R,
              slope_from_P50$R,
              TLP_from_WD$R,
              WD_from_P50$R,
              WD_from_slope$R,
              LMA_from_TLP$R,
              LS_from_LMA$R,
              P50_from_Ks$R,
              P50_from_LS$R,
              LS_from_TLP$R,
              WD_from_LMA$R,
              slope_from_Ks$R)
  all_R2 <- c(TLP_from_P50$R2,
              slope_from_TLP$R2,
              slope_from_P50$R2,
              TLP_from_WD$R2,
              WD_from_P50$R2,
              WD_from_slope$R2,
              LMA_from_TLP$R2,
              LS_from_LMA$R2,
              P50_from_Ks$R2,
              P50_from_LS$R2,
              LS_from_TLP$R2,
              WD_from_LMA$R2,
              slope_from_Ks$R2)
  all_R2adj <- c(TLP_from_P50$R2adj,
                 slope_from_TLP$R2adj,
                 slope_from_P50$R2adj,
                 TLP_from_WD$R2adj,
                 WD_from_P50$R2adj,
                 WD_from_slope$R2adj,
                 LMA_from_TLP$R2adj,
                 LS_from_LMA$R2adj,
                 P50_from_Ks$R2adj,
                 P50_from_LS$R2adj,
                 LS_from_TLP$R2adj,
                 WD_from_LMA$R2adj,
                 slope_from_Ks$R2adj)
  all_rmse <- c(TLP_from_P50$rmse,
                slope_from_TLP$rmse,
                slope_from_P50$rmse,
                TLP_from_WD$rmse,
                WD_from_P50$rmse,
                WD_from_slope$rmse,
                LMA_from_TLP$rmse,
                LS_from_LMA$rmse,
                P50_from_Ks$rmse,
                P50_from_LS$rmse,
                LS_from_TLP$rmse,
                WD_from_LMA$rmse,
                slope_from_Ks$rmse)
  
  all_sma_bivar <- data.frame(all_label1,all_label2,all_R,all_R2,all_R2adj,all_rmse)
  
  return_vals <- list("TLP_from_P50"=TLP_from_P50,
                      "slope_from_TLP"=slope_from_TLP,
                      "slope_from_P50"=slope_from_P50,
                      "TLP_from_WD"=TLP_from_WD,
                      "WD_from_P50"=WD_from_P50,
                      "WD_from_slope"=WD_from_slope,
                      "LMA_from_TLP"=LMA_from_TLP,
                      "LS_from_LMA"=LS_from_LMA,
                      "P50_from_Ks"=P50_from_Ks,
                      "P50_from_LS"=P50_from_LS,
                      "LS_from_TLP"=LS_from_TLP,
                      "WD_from_LMA"=WD_from_LMA,
                      "slope_from_Ks"=slope_from_Ks,
                      "all_sma_bivar"=all_sma_bivar)
  
  return(return_vals)
}