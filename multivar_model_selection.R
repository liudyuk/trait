# Functions in this file go through the logical process of selecting the best performing multivariate models for the different traits following links in the hypothesis framework.
#
# T. Pugh
# 12.12.20


P50_multivar_test <- function(trait,view_stats=FALSE, regr_type = 'lm' )  {
  
  if(regr_type=='lm'){#lm allows for mor than 5 variables 
    # P50 from TLP, LS , Ks slope and WD
    P50_from_TLP_LS_Ks_slope_WD <- sma_plot_stats(data.frame(trait$TLP,trait$LS,trait$Ks,trait$slope,trait$WD,trait$P50),c("TLP","LS","Ks","slope","WD","P50"),nbtstrp,regression_type = regr_type)
    plot(trait$P50[P50_from_TLP_LS_Ks_slope_WD$dataused],P50_from_TLP_LS_Ks_slope_WD$var_est,pch=16,xlab="P50",ylab="P50_est",main="P50 vs P50_est")
  }else{ # sma
    P50_from_TLP_LS_Ks_slope_WD <- list(R=NA, R2adj=NA,rmse=NA, ndata=NA)
  }
    # P50 from TLP, LS , Ks and slope
    P50_from_TLP_LS_Ks_slope <- sma_plot_stats(data.frame(trait$TLP,trait$LS,trait$Ks,trait$slope,trait$P50),c("TLP","LS","Ks","slope","P50"),nbtstrp,regression_type = regr_type)
    plot(trait$P50[P50_from_TLP_LS_Ks_slope$dataused],P50_from_TLP_LS_Ks_slope$var_est,pch=16,xlab="P50",ylab="P50_est",main="P50 vs P50_est")
  
    # P50 from TLP, Ks and WD
    P50_from_TLP_Ks_WD <- sma_plot_stats(data.frame(trait$TLP,trait$Ks,trait$WD,trait$P50),c("TLP","Ks","WD","P50"),nbtstrp, regression_type = regr_type )
    plot(trait$P50[P50_from_TLP_Ks_WD$dataused],P50_from_TLP_Ks_WD$var_est,pch=16,xlab="P50",ylab="P50_est",main="P50 vs P50_est")
    
    # P50 from TLP, LS and Ks
    P50_from_TLP_LS_Ks <- sma_plot_stats(data.frame(trait$TLP,trait$LS,trait$Ks,trait$P50),c("TLP","LS","Ks","P50"),nbtstrp,regression_type = regr_type)
    plot(trait$P50[P50_from_TLP_LS_Ks$dataused],P50_from_TLP_LS_Ks$var_est,pch=16,xlab="P50",ylab="P50_est",main="P50 vs P50_est")
    
    # P50 from TLP and Ks
    P50_from_TLP_Ks <- sma_plot_stats(data.frame(trait$TLP,trait$Ks,trait$P50),c("TLP","Ks","P50"),nbtstrp,regression_type = regr_type)
    plot(trait$P50[P50_from_TLP_Ks$dataused],P50_from_TLP_Ks$var_est,pch=16,xlab="P50",ylab="P50_est",main="P50 vs P50_est")
    
    # P50 from TLP and LS
    P50_from_TLP_LS <- sma_plot_stats(data.frame(trait$TLP,trait$LS,trait$P50),c("TLP","LS","P50"),nbtstrp,regression_type = regr_type)
    plot(trait$P50[P50_from_TLP_LS$dataused],P50_from_TLP_LS$var_est,pch=16,xlab="P50",ylab="P50_est",main="P50 vs P50_est")
    
    # P50 from Ks and LS
    P50_from_Ks_LS <- sma_plot_stats(data.frame(trait$Ks,trait$LS,trait$P50),c("Ks","LS","P50"),nbtstrp,regression_type = regr_type)
    plot(trait$P50[P50_from_Ks_LS$dataused],P50_from_Ks_LS$var_est,pch=16,xlab="P50",ylab="P50_est",main="P50 vs P50_est")
    
    # P50 from TLP
    P50_from_TLP <- sma_plot_stats(data.frame(trait$TLP,trait$P50),c("TLP","P50"),nbtstrp,regression_type = regr_type)
    plot(trait$P50[P50_from_TLP$dataused],P50_from_TLP$var_est,pch=16,xlab="P50",ylab="P50_est",main="P50 vs P50_est")
    
    # P50 from Ks
    P50_from_Ks <- sma_plot_stats(data.frame(trait$Ks,trait$P50),c("Ks","P50"),nbtstrp,regression_type = regr_type)
    
    # P50 from LS
    P50_from_LS <- sma_plot_stats(data.frame(trait$LS,trait$P50),c("LS","P50"),nbtstrp,regression_type = regr_type)
    
    # Some of the above combinations have relatively high R2adj, but very reduced data points. Therefore test if the worse-performing combinations are only worse performing because they include a greater diversity of data.
    
    # P50 from TLP (same species as for TLP and LS)
    P50_from_TLP_limitTLPLS <- sma_plot_stats(data.frame(trait$TLP,trait$P50),c("TLP","P50"),nbtstrp,F,indin = P50_from_TLP_LS$dataused,regression_type = regr_type)
    
    # P50 from Ks (same species as for TLP and Ks)
    P50_from_Ks_limitTLPKs <- sma_plot_stats(data.frame(trait$Ks,trait$P50),c("Ks","P50"),nbtstrp,F,indin = P50_from_TLP_Ks$dataused,regression_type = regr_type)
    
    # P50 from LS (same species as for TLP and LS)
    P50_from_LS_limitTLPLS <- sma_plot_stats(data.frame(trait$LS,trait$P50),c("LS","P50"),nbtstrp,F,indin = P50_from_TLP_LS$dataused,regression_type = regr_type)
    
    # Summarise statistics
    all_testnames_P50 <- c("P50_from_TLP_LS_Ks_slope_WD","P50_from_TLP_LS_Ks_slope",  "P50_from_TLP_Ks_WD","P50_from_TLP_Ks","P50_from_TLP","P50_from_Ks","P50_from_TLP_limitTLPLS","P50_from_Ks_limitTLPKs","P50_from_TLP_LS_Ks","P50_from_TLP_LS","P50_from_Ks_LS","P50_from_LS_limitTLPLS","P50_from_LS")
    all_R2_P50 <- c(P50_from_TLP_LS_Ks_slope_WD$R,       P50_from_TLP_LS_Ks_slope$R,   P50_from_TLP_Ks_WD$R2,P50_from_TLP_Ks$R2,P50_from_TLP$R2,P50_from_Ks$R2,P50_from_TLP_limitTLPLS$R2,P50_from_Ks_limitTLPKs$R2,P50_from_TLP_LS_Ks$R2,P50_from_TLP_LS$R2,P50_from_Ks_LS$R2,P50_from_LS_limitTLPLS$R2,P50_from_LS$R2)
    all_R2adj_P50 <- c(P50_from_TLP_LS_Ks_slope_WD$R2adj,P50_from_TLP_LS_Ks_slope$R2adj,P50_from_TLP_Ks_WD$R2adj,P50_from_TLP_Ks$R2adj,P50_from_TLP$R2adj,P50_from_Ks$R2adj,P50_from_TLP_limitTLPLS$R2adj,P50_from_Ks_limitTLPKs$R2adj,P50_from_TLP_LS_Ks$R2adj,P50_from_TLP_LS$R2adj,P50_from_Ks_LS$R2adj,P50_from_LS_limitTLPLS$R2adj,P50_from_LS$R2adj)
    all_rmse_P50 <- c(P50_from_TLP_LS_Ks_slope_WD$rmse, P50_from_TLP_LS_Ks_slope$rmse,  P50_from_TLP_Ks_WD$rmse,P50_from_TLP_Ks$rmse,P50_from_TLP$rmse,P50_from_Ks$rmse,P50_from_TLP_limitTLPLS$rmse,P50_from_Ks_limitTLPKs$rmse,P50_from_TLP_LS_Ks$rmse,P50_from_TLP_LS$rmse,P50_from_Ks_LS$rmse,P50_from_LS_limitTLPLS$rmse,P50_from_LS$rmse)
    all_ndata_P50 <- c(P50_from_TLP_LS_Ks_slope_WD$ndata,P50_from_TLP_LS_Ks_slope$ndata,P50_from_TLP_Ks_WD$ndata,P50_from_TLP_Ks$ndata,P50_from_TLP$ndata,P50_from_Ks$ndata,P50_from_TLP_limitTLPLS$ndata,P50_from_Ks_limitTLPKs$ndata,P50_from_TLP_LS_Ks$ndata,P50_from_TLP_LS$ndata,P50_from_Ks_LS$ndata,P50_from_LS_limitTLPLS$ndata,P50_from_LS$ndata)
   
    all_P50 <- data.frame(all_testnames_P50,all_R2_P50,all_R2adj_P50,all_rmse_P50,all_ndata_P50)
    if(view_stats==TRUE)View(all_P50)
    
    # BEST MODEL: P50_from_TLP_Ks
    # Test MAT and PPT coverage of species for best model
    plot(trait$MAT[P50_from_TLP_Ks$dataused],trait$MAP[P50_from_TLP_Ks$dataused])
    #plot(trait$MAT[P50_from_TLP_LS_Ks_slope$dataused],trait$MAP[P50_from_TLP_LS_Ks_slope$dataused])
    #P50_from_TLP_LS_Ks_slope_WD
    plot(trait$MAT[P50_from_TLP_Ks_WD$dataused],trait$MAP[P50_from_TLP_Ks_WD$dataused])
    # RESULT: WIDE CLIMATE COVERAGE
    
    # DECISION: P50_from_TLP_Ks_WD
    
    return_vals <- list("all_P50"=all_P50,
                        "P50_from_TLP_Ks_WD" = P50_from_TLP_Ks_WD )
 
  return(return_vals)
  
}

P50_multivar_test_save <- function(trait,view_stats=FALSE, regr_type = 'lm' ,test=T)  {
  
  if(test==T){
  # P50 from TLP, Ks and WD
  P50_from_TLP_Ks_WD <- sma_plot_stats(data.frame(trait$TLP,trait$Ks,trait$WD,trait$P50),c("TLP","Ks","WD","P50"),nbtstrp, regression_type = regr_type )
  plot(trait$P50[P50_from_TLP_Ks_WD$dataused],P50_from_TLP_Ks_WD$var_est,pch=16,xlab="P50",ylab="P50_est",main="P50 vs P50_est")
  
  # P50 from TLP, LS and Ks
  P50_from_TLP_LS_Ks <- sma_plot_stats(data.frame(trait$TLP,trait$LS,trait$Ks,trait$P50),c("TLP","LS","Ks","P50"),nbtstrp,regression_type = regr_type)
  plot(trait$P50[P50_from_TLP_LS_Ks$dataused],P50_from_TLP_LS_Ks$var_est,pch=16,xlab="P50",ylab="P50_est",main="P50 vs P50_est")
  
  # P50 from TLP and Ks
  P50_from_TLP_Ks <- sma_plot_stats(data.frame(trait$TLP,trait$Ks,trait$P50),c("TLP","Ks","P50"),nbtstrp,regression_type = regr_type)
  plot(trait$P50[P50_from_TLP_Ks$dataused],P50_from_TLP_Ks$var_est,pch=16,xlab="P50",ylab="P50_est",main="P50 vs P50_est")
  
  # P50 from TLP and LS
  P50_from_TLP_LS <- sma_plot_stats(data.frame(trait$TLP,trait$LS,trait$P50),c("TLP","LS","P50"),nbtstrp,regression_type = regr_type)
  plot(trait$P50[P50_from_TLP_LS$dataused],P50_from_TLP_LS$var_est,pch=16,xlab="P50",ylab="P50_est",main="P50 vs P50_est")
  
  # P50 from Ks and LS
  P50_from_Ks_LS <- sma_plot_stats(data.frame(trait$Ks,trait$LS,trait$P50),c("Ks","LS","P50"),nbtstrp,regression_type = regr_type)
  plot(trait$P50[P50_from_Ks_LS$dataused],P50_from_Ks_LS$var_est,pch=16,xlab="P50",ylab="P50_est",main="P50 vs P50_est")
  
  # P50 from TLP
  P50_from_TLP <- sma_plot_stats(data.frame(trait$TLP,trait$P50),c("TLP","P50"),nbtstrp,regression_type = regr_type)
  plot(trait$P50[P50_from_TLP$dataused],P50_from_TLP$var_est,pch=16,xlab="P50",ylab="P50_est",main="P50 vs P50_est")
  
  # P50 from Ks
  P50_from_Ks <- sma_plot_stats(data.frame(trait$Ks,trait$P50),c("Ks","P50"),nbtstrp,regression_type = regr_type)
  
  # P50 from LS
  P50_from_LS <- sma_plot_stats(data.frame(trait$LS,trait$P50),c("LS","P50"),nbtstrp,regression_type = regr_type)
  
  # Some of the above combinations have relatively high R2adj, but very reduced data points. Therefore test if the worse-performing combinations are only worse performing because they include a greater diversity of data.
  
  # P50 from TLP (same species as for TLP and LS)
  P50_from_TLP_limitTLPLS <- sma_plot_stats(data.frame(trait$TLP,trait$P50),c("TLP","P50"),nbtstrp,F,indin = P50_from_TLP_LS$dataused,regression_type = regr_type)
  
  # P50 from Ks (same species as for TLP and Ks)
  P50_from_Ks_limitTLPKs <- sma_plot_stats(data.frame(trait$Ks,trait$P50),c("Ks","P50"),nbtstrp,F,indin = P50_from_TLP_Ks$dataused,regression_type = regr_type)
  
  # P50 from LS (same species as for TLP and LS)
  P50_from_LS_limitTLPLS <- sma_plot_stats(data.frame(trait$LS,trait$P50),c("LS","P50"),nbtstrp,F,indin = P50_from_TLP_LS$dataused,regression_type = regr_type)
  
  # Summarise statistics
  all_testnames_P50 <- c("P50_from_TLP_Ks_WD","P50_from_TLP_Ks","P50_from_TLP","P50_from_Ks","P50_from_TLP_limitTLPLS","P50_from_Ks_limitTLPKs","P50_from_TLP_LS_Ks","P50_from_TLP_LS","P50_from_Ks_LS","P50_from_LS_limitTLPLS","P50_from_LS")
  all_R2_P50 <- c(P50_from_TLP_Ks_WD$R2,P50_from_TLP_Ks$R2,P50_from_TLP$R2,P50_from_Ks$R2,P50_from_TLP_limitTLPLS$R2,P50_from_Ks_limitTLPKs$R2,P50_from_TLP_LS_Ks$R2,P50_from_TLP_LS$R2,P50_from_Ks_LS$R2,P50_from_LS_limitTLPLS$R2,P50_from_LS$R2)
  all_R2adj_P50 <- c(P50_from_TLP_Ks_WD$R2adj,P50_from_TLP_Ks$R2adj,P50_from_TLP$R2adj,P50_from_Ks$R2adj,P50_from_TLP_limitTLPLS$R2adj,P50_from_Ks_limitTLPKs$R2adj,P50_from_TLP_LS_Ks$R2adj,P50_from_TLP_LS$R2adj,P50_from_Ks_LS$R2adj,P50_from_LS_limitTLPLS$R2adj,P50_from_LS$R2adj)
  all_rmse_P50 <- c(P50_from_TLP_Ks_WD$rmse,P50_from_TLP_Ks$rmse,P50_from_TLP$rmse,P50_from_Ks$rmse,P50_from_TLP_limitTLPLS$rmse,P50_from_Ks_limitTLPKs$rmse,P50_from_TLP_LS_Ks$rmse,P50_from_TLP_LS$rmse,P50_from_Ks_LS$rmse,P50_from_LS_limitTLPLS$rmse,P50_from_LS$rmse)
  all_ndata_P50 <- c(P50_from_TLP_Ks_WD$ndata,P50_from_TLP_Ks$ndata,P50_from_TLP$ndata,P50_from_Ks$ndata,P50_from_TLP_limitTLPLS$ndata,P50_from_Ks_limitTLPKs$ndata,P50_from_TLP_LS_Ks$ndata,P50_from_TLP_LS$ndata,P50_from_Ks_LS$ndata,P50_from_LS_limitTLPLS$ndata,P50_from_LS$ndata)
  
  all_P50 <- data.frame(all_testnames_P50,all_R2_P50,all_R2adj_P50,all_rmse_P50,all_ndata_P50)
   if(view_stats==TRUE)View(all_P50)
  
  # BEST MODEL: P50_from_TLP_Ks
  # Test MAT and PPT coverage of species for best model
  plot(trait$MAT[P50_from_TLP_Ks$dataused],trait$MAP[P50_from_TLP_Ks$dataused])
  # RESULT: WIDE CLIMATE COVERAGE
  
  # DECISION: P50_from_TLP_Ks
  
   return_vals <- list("all_P50"=all_P50,
                      "P50_from_TLP_Ks"=P50_from_TLP_Ks)
   save(all_P50, file= 'data/all_P50.RData')
  }
  else{
    load('data/all_P50.RData')
    
    # P50 from TLP and Ks
    P50_from_TLP_Ks <- sma_plot_stats(data.frame(trait$TLP,trait$Ks,trait$P50),c("TLP","Ks","P50"),nbtstrp,regression_type = regr_type)
    plot(trait$P50[P50_from_TLP_Ks$dataused],P50_from_TLP_Ks$var_est,pch=16,xlab="P50",ylab="P50_est",main="P50 vs P50_est")
    
    
    return_vals <- list("all_P50"=all_P50,
                        "P50_from_TLP_Ks"=P50_from_TLP_Ks)
                        
  }
  return(return_vals)
   
}

TLP_multivar_test <- function(trait,view_stats=FALSE, regr_type = 'lm' )  {
  
  # TLP from LS, LMA, P50 and WD
  TLP_from_LS_LMA_P50_WD <- sma_plot_stats(data.frame(trait$LS,trait$LMA,trait$P50,trait$WD,trait$TLP),c("LS","LMA","P50","WD","TLP"),nbtstrp,regression_type = regr_type)
  plot(trait$TLP[TLP_from_LS_LMA_P50_WD$dataused],TLP_from_LS_LMA_P50_WD$var_est,pch=16,xlab="TLP",ylab="TLP_est",main="TLP vs TLP_est")
  
  # TLP from LS, LMA, P50 and slope
  TLP_from_LS_LMA_P50_slope <- sma_plot_stats(data.frame(trait$LS,trait$LMA,trait$P50,trait$slope,trait$TLP),c("LS","LMA","P50","slope","TLP"),nbtstrp,regression_type = regr_type)
  plot(trait$TLP[TLP_from_LS_LMA_P50_slope$dataused],TLP_from_LS_LMA_P50_slope$var_est,pch=16,xlab="TLP",ylab="TLP_est",main="TLP vs TLP_est")
  
  # TLP from LS, LMA and P50
  TLP_from_LS_LMA_P50 <- sma_plot_stats(data.frame(trait$LS,trait$LMA,trait$P50,trait$TLP),c("LS","LMA","P50","TLP"),nbtstrp,regression_type = regr_type)
  plot(trait$TLP[TLP_from_LS_LMA_P50$dataused],TLP_from_LS_LMA_P50$var_est,pch=16,xlab="TLP",ylab="TLP_est",main="TLP vs TLP_est")
  
  # TLP from LS and LMA
  TLP_from_LS_LMA <- sma_plot_stats(data.frame(trait$LS,trait$LMA,trait$TLP),c("LS","LMA","TLP"),nbtstrp,regression_type = regr_type)
  plot(trait$TLP[TLP_from_LS_LMA$dataused],TLP_from_LS_LMA$var_est,pch=16,xlab="TLP",ylab="TLP_est",main="TLP vs TLP_est")
  
  # TLP from P50 and LMA
  TLP_from_P50_LMA <- sma_plot_stats(data.frame(trait$P50,trait$LMA,trait$TLP),c("P50","LMA","TLP"),nbtstrp,regression_type = regr_type)
  plot(trait$TLP[TLP_from_P50_LMA$dataused],TLP_from_P50_LMA$var_est,pch=16,xlab="TLP",ylab="TLP_est",main="TLP vs TLP_est")
  
  # TLP from P50 and LS
  TLP_from_P50_LS <- sma_plot_stats(data.frame(trait$P50,trait$LS,trait$TLP),c("P50","LS","TLP"),nbtstrp,regression_type = regr_type)
  plot(trait$TLP[TLP_from_P50_LS$dataused],TLP_from_P50_LS$var_est,pch=16,xlab="TLP",ylab="TLP_est",main="TLP vs TLP_est")
  
  # TLP from LS
  TLP_from_LS <- sma_plot_stats(data.frame(trait$LS,trait$TLP),c("LS","TLP"),nbtstrp,regression_type = regr_type)
  
  # TLP from P50
  TLP_from_P50 <- sma_plot_stats(data.frame(trait$P50,trait$TLP),c("P50","TLP"),nbtstrp,regression_type = regr_type)
  
  # TLP from LMA
  TLP_from_LMA <- sma_plot_stats(data.frame(trait$LMA,trait$TLP),c("LMA","TLP"),nbtstrp,regression_type = regr_type)
  
  # Some of the above combinations have relatively high R2adj, but very reduced data points. Therefore test if the worse-performing combinations are only worse performing because they include a greater diversity of data.
  
  # TLP from P50 (same species as for LS, LMA and P50)
  TLP_from_P50_limitLSLMAP50 <- sma_plot_stats(data.frame(trait$P50,trait$TLP),c("P50","TLP"),nbtstrp,F,indin = TLP_from_LS_LMA_P50$dataused,regression_type = regr_type)
  
  # TLP from P50 and LMA (same species as for LS, LMA and P50)
  TLP_from_P50_LMA_limitLALMAP50 <- sma_plot_stats(data.frame(trait$P50,trait$LMA,trait$TLP),c("P50","LMA","TLP"),nbtstrp,F,indin = TLP_from_LS_LMA_P50$dataused,regression_type = regr_type)
  
  # Summarise statistics
  all_testnames_TLP <- c("TLP_from_LS_LMA_P50_WD","TLP_from_LS_LMA_P50_slope","TLP_from_LS_LMA_P50","TLP_from_LS_LMA","TLP_from_P50_LMA","TLP_from_P50_LMA_limitLALMAP50","TLP_from_P50_LS","TLP_from_LS","TLP_from_P50","TLP_from_LMA","TLP_from_P50_limitLSLMAP50")
  all_R2_TLP <- c(TLP_from_LS_LMA_P50_WD$R2,TLP_from_LS_LMA_P50_slope$R2,TLP_from_LS_LMA_P50$R2,TLP_from_LS_LMA$R2,TLP_from_P50_LMA$R2,TLP_from_P50_LMA_limitLALMAP50$R2,TLP_from_P50_LS$R2,TLP_from_LS$R2,TLP_from_P50$R2,TLP_from_LMA$R2,TLP_from_P50_limitLSLMAP50$R2)
  all_R2adj_TLP <- c(TLP_from_LS_LMA_P50_WD$R2adj,TLP_from_LS_LMA_P50_slope$R2adj,TLP_from_LS_LMA_P50$R2adj,TLP_from_LS_LMA$R2adj,TLP_from_P50_LMA$R2adj,TLP_from_P50_LMA_limitLALMAP50$R2adj,TLP_from_P50_LS$R2adj,TLP_from_LS$R2adj,TLP_from_P50$R2adj,TLP_from_LMA$R2adj,TLP_from_P50_limitLSLMAP50$R2adj)
  all_rmse_TLP <- c(TLP_from_LS_LMA_P50_WD$rmse,TLP_from_LS_LMA_P50_slope$rmse,TLP_from_LS_LMA_P50$rmse,TLP_from_LS_LMA$rmse,TLP_from_P50_LMA$rmse,TLP_from_P50_LMA_limitLALMAP50$rmse,TLP_from_P50_LS$rmse,TLP_from_LS$rmse,TLP_from_P50$rmse,TLP_from_LMA$rmse,TLP_from_P50_limitLSLMAP50$rmse)
  all_ndata_TLP <- c(TLP_from_LS_LMA_P50_WD$ndata,TLP_from_LS_LMA_P50_slope$ndata,TLP_from_LS_LMA_P50$ndata,TLP_from_LS_LMA$ndata,TLP_from_P50_LMA$ndata,TLP_from_P50_LMA_limitLALMAP50$ndata,TLP_from_P50_LS$ndata,TLP_from_LS$ndata,TLP_from_P50$ndata,TLP_from_LMA$ndata,TLP_from_P50_limitLSLMAP50$ndata)
  
  all_TLP <- data.frame(all_testnames_TLP,all_R2_TLP,all_R2adj_TLP,all_rmse_TLP,all_ndata_TLP)
   if(view_stats==TRUE)View(all_TLP)
  
  # CHOICE: Although WD improves the fit, do not use it as our hypothesis framework does not posit a direct link between TLP and WD (only indirect via P50 and slope)
  
  # BEST MODEL: TLP_from_LS_LMA_P50
  # NO EVIDENCE that the smaller selection of species for 3 variables is the reason behind the fit
  # Test MAT and PPT coverage of species for best model
  plot(trait$MAT[TLP_from_LS_LMA_P50$dataused],trait$MAP[TLP_from_LS_LMA_P50$dataused])
  # Test MAT and PPT coverage of species from TLP_from_P50_LMA
  plot(trait$MAT[TLP_from_P50_LMA$dataused],trait$MAP[TLP_from_P50_LMA$dataused])
  # Test MAT and PPT coverage of species from TLP_from_LMA
  plot(trait$MAT[TLP_from_LMA$dataused],trait$MAP[TLP_from_LMA$dataused])
  # RESULT: WIDE CLIMATE COVERAGE for TLP_from_LS_LMA_P50 EXCEPT temperate rainforest (which is captured by TLP_from_P50_LMA)
  
  # Test if relationships are consistent in character despite regardless of climate zone differences
  plot(trait$TLP[TLP_from_LS_LMA_P50$dataused],TLP_from_LS_LMA_P50$var_est,pch=16,xlab="TLP",ylab="TLP_est",main="TLP vs TLP_est")
  points(trait$TLP[TLP_from_P50_LMA$dataused],TLP_from_P50_LMA$var_est,pch=16,col="red")
  
  # DECISION: TLP_from_LS_LMA_P50
  
  return_vals <- list("all_TLP"=all_TLP,
                      "TLP_from_LS_LMA_P50"=TLP_from_LS_LMA_P50)
  
  return(return_vals)
}

LMA_multivar_test_BDT <- function(trait_BDT,view_stats=FALSE, regr_type = 'lm' )  {
  # This version for BD and BT only.
  
  # LMA from TLP, LS and WD
  LMA_from_TLP_LS_WD <- sma_plot_stats(data.frame(trait_BDT$TLP,trait_BDT$LS,trait_BDT$WD,trait_BDT$LMA),c("TLP","LS","WD","LMA"),nbtstrp,regression_type = regr_type)
  plot(trait_BDT$LMA[LMA_from_TLP_LS_WD$dataused],LMA_from_TLP_LS_WD$var_est,pch=16,xlab="LMA",ylab="LMA_est",main="LMA vs LMA_est")
  
  # LMA from TLP and LS
  LMA_from_TLP_LS <- sma_plot_stats(data.frame(trait_BDT$TLP,trait_BDT$LS,trait_BDT$LMA),c("TLP","LS","LMA"),nbtstrp,regression_type = regr_type)
  plot(trait_BDT$LMA[LMA_from_TLP_LS$dataused],LMA_from_TLP_LS$var_est,pch=16,xlab="LMA",ylab="LMA_est",main="LMA vs LMA_est")
  
  # LMA from TLP and WD
  LMA_from_TLP_WD <- sma_plot_stats(data.frame(trait_BDT$TLP,trait_BDT$WD,trait_BDT$LMA),c("TLP","WD","LMA"),nbtstrp,regression_type = regr_type)
  plot(trait_BDT$LMA[LMA_from_TLP_WD$dataused],LMA_from_TLP_WD$var_est,pch=16,xlab="LMA",ylab="LMA_est",main="LMA vs LMA_est")
  
  # LMA from LS and WD
  LMA_from_LS_WD <- sma_plot_stats(data.frame(trait_BDT$LS,trait_BDT$WD,trait_BDT$LMA),c("LS","WD","LMA"),nbtstrp,regression_type = regr_type)
  plot(trait_BDT$LMA[LMA_from_LS_WD$dataused],LMA_from_LS_WD$var_est,pch=16,xlab="LMA",ylab="LMA_est",main="LMA vs LMA_est")
  
  # LMA from TLP
  LMA_from_TLP <- sma_plot_stats(data.frame(trait_BDT$TLP,trait_BDT$LMA),c("TLP","LMA"),nbtstrp,regression_type = regr_type)
  plot(trait_BDT$LMA[LMA_from_TLP$dataused],LMA_from_TLP$var_est,pch=16,xlab="LMA",ylab="LMA_est",main="LMA vs LMA_est")
  
  # LMA from LS
  LMA_from_LS <- sma_plot_stats(data.frame(trait_BDT$LS,trait_BDT$LMA),c("LS","LMA"),nbtstrp,regression_type = regr_type)
  
  # LMA from WD
  LMA_from_WD <- sma_plot_stats(data.frame(trait_BDT$WD,trait_BDT$LMA),c("WD","LMA"),nbtstrp,regression_type = regr_type)
  
  # Summarise statistics
  all_testnames_LMA <- c("LMA_from_TLP_LS_WD","LMA_from_TLP_LS","LMA_from_TLP_WD","LMA_from_LS_WD","LMA_from_TLP","LMA_from_LS","LMA_from_WD")
  all_R2_LMA <- c(LMA_from_TLP_LS_WD$R2,LMA_from_TLP_LS$R2,LMA_from_TLP_WD$R2,LMA_from_LS_WD$R2,LMA_from_TLP$R2,LMA_from_LS$R2,LMA_from_WD$R2)
  all_R2adj_LMA <- c(LMA_from_TLP_LS_WD$R2adj,LMA_from_TLP_LS$R2adj,LMA_from_TLP_WD$R2adj,LMA_from_LS_WD$R2adj,LMA_from_TLP$R2adj,LMA_from_LS$R2adj,LMA_from_WD$R2adj)
  all_rmse_LMA <- c(LMA_from_TLP_LS_WD$rmse,LMA_from_TLP_LS$rmse,LMA_from_TLP_WD$rmse,LMA_from_LS_WD$rmse,LMA_from_TLP$rmse,LMA_from_LS$rmse,LMA_from_WD$rmse)
  all_ndata_LMA <- c(LMA_from_TLP_LS_WD$ndata,LMA_from_TLP_LS$ndata,LMA_from_TLP_WD$ndata,LMA_from_LS_WD$ndata,LMA_from_TLP$ndata,LMA_from_LS$ndata,LMA_from_WD$ndata)
  
  all_LMA <- data.frame(all_testnames_LMA,all_R2_LMA,all_R2adj_LMA,all_rmse_LMA,all_ndata_LMA)
   if(view_stats==TRUE)View(all_LMA)
  
  # CHOICE: best model is LMA_from_TLP
  
  # Test MAT and PPT coverage of species for chosen model
  plot(trait_BDT$MAT[LMA_from_TLP$dataused],trait_BDT$MAP[LMA_from_TLP$dataused])
  # RESULT: WIDE CLIMATE COVERAGE
  
  # DECISION: LMA_from_TLP

  return_vals <- list("all_LMA"=all_LMA,
                      "LMA_from_TLP"=LMA_from_TLP)
  return(return_vals)
}

LMA_multivar_test_BE <- function(trait_BE,view_stats=FALSE, regr_type = 'lm' ) {
  # This version for BE only.
  
  # LMA from TLP, LS and WD
  LMA_from_TLP_LS_WD <- sma_plot_stats(vars=data.frame(trait_BE$TLP,trait_BE$LS,trait_BE$WD,trait_BE$LMA),labels=c("TLP","LS","WD","LMA"),nbtstrp, regression_type =  regr_type)
  plot(trait_BE$LMA[LMA_from_TLP_LS_WD$dataused],LMA_from_TLP_LS_WD$var_est,pch=16,xlab="LMA",ylab="LMA_est",main="LMA vs LMA_est")
  
  # LMA from TLP and LS
  LMA_from_TLP_LS <- sma_plot_stats(vars=data.frame(trait_BE$TLP,trait_BE$LS,trait_BE$LMA),labels=c("TLP","LS","LMA"),nbtstrp,T, regression_type =  regr_type)
  plot(trait_BE$LMA[LMA_from_TLP_LS$dataused],LMA_from_TLP_LS$var_est,pch=16,xlab="LMA",ylab="LMA_est",main="LMA vs LMA_est")
  
  # LMA from TLP and WD
  LMA_from_TLP_WD <- sma_plot_stats(vars=data.frame(trait_BE$TLP,trait_BE$WD,trait_BE$LMA),labels=c("TLP","WD","LMA"),nbtstrp,T, regression_type =  regr_type)
  plot(trait_BE$LMA[LMA_from_TLP_WD$dataused],LMA_from_TLP_WD$var_est,pch=16,xlab="LMA",ylab="LMA_est",main="LMA vs LMA_est")
  
  # LMA from LS and WD
  LMA_from_LS_WD <- sma_plot_stats(vars=data.frame(trait_BE$LS,trait_BE$WD,trait_BE$LMA),labels=c("LS","WD","LMA"),nbtstrp,T, regression_type =  regr_type)
  plot(trait_BE$LMA[LMA_from_LS_WD$dataused],LMA_from_LS_WD$var_est,pch=16,xlab="LMA",ylab="LMA_est",main="LMA vs LMA_est")
  
  # LMA from TLP
  LMA_from_TLP <- sma_plot_stats(vars=data.frame(trait_BE$TLP,trait_BE$LMA),labels=c("TLP","LMA"),nbtstrp , regression_type =  regr_type)
  plot(trait_BE$LMA[LMA_from_TLP$dataused],LMA_from_TLP$var_est,pch=16,xlab="LMA",ylab="LMA_est",main="LMA vs LMA_est")
  
  # LMA from LS
  LMA_from_LS <- sma_plot_stats(vars=data.frame(trait_BE$LS,trait_BE$LMA),labels=c("LS","LMA"),nbtstrp, regression_type =  regr_type)
  
  # LMA from WD
  LMA_from_WD <- sma_plot_stats(vars=data.frame(trait_BE$WD,trait_BE$LMA),labels=c("WD","LMA"),nbtstrp, regression_type =  regr_type)
  
  # Some of the above combinations have relatively high R2adj, but very reduced data points. Therefore test if the worse-performing combinations are only worse performing because they include a greater diversity of data.
  
  # LMA from LS and WD (same species as for LS, WD and TLP)
  LMA_from_LS_WD_limitTLPLSWD <- sma_plot_stats(data.frame(trait_BE$LS,trait_BE$WD,trait_BE$LMA),c("LS","WD","LMA"),nbtstrp,F,indin = LMA_from_TLP_LS_WD$dataused, regression_type =  regr_type)
  plot(trait_BE$LMA[LMA_from_LS_WD_limitTLPLSWD$dataused],LMA_from_LS_WD_limitTLPLSWD$var_est,pch=16,xlab="LMA",ylab="LMA_est",main="LMA vs LMA_est")
  
  # LMA from TLP and WD (same species as for LS, WD and TLP)
  LMA_from_TLP_WD_limitTLPLSWD <- sma_plot_stats(data.frame(trait_BE$TLP,trait_BE$WD,trait_BE$LMA),c("TLP","WD","LMA"),nbtstrp,F,indin = LMA_from_TLP_LS_WD$dataused, regression_type =  regr_type)
  plot(trait_BE$LMA[LMA_from_TLP_WD_limitTLPLSWD$dataused],LMA_from_TLP_WD_limitTLPLSWD$var_est,pch=16,xlab="LMA",ylab="LMA_est",main="LMA vs LMA_est")
  
  # Summarise statistics
  all_testnames_LMA <- c("LMA_from_TLP_LS_WD","LMA_from_TLP_LS","LMA_from_TLP_WD","LMA_from_LS_WD","LMA_from_TLP","LMA_from_LS","LMA_from_WD","LMA_from_LS_WD_limitTLPLSWD","LMA_from_TLP_WD_limitTLPLSWD")
  all_R2_LMA <- c(LMA_from_TLP_LS_WD$R2,LMA_from_TLP_LS$R2,LMA_from_TLP_WD$R2,LMA_from_LS_WD$R2,LMA_from_TLP$R2,LMA_from_LS$R2,LMA_from_WD$R2,LMA_from_LS_WD_limitTLPLSWD$R2,LMA_from_TLP_WD_limitTLPLSWD$R2)
  all_R2adj_LMA <- c(LMA_from_TLP_LS_WD$R2adj,LMA_from_TLP_LS$R2adj,LMA_from_TLP_WD$R2adj,LMA_from_LS_WD$R2adj,LMA_from_TLP$R2adj,LMA_from_LS$R2adj,LMA_from_WD$R2adj,LMA_from_LS_WD_limitTLPLSWD$R2adj,LMA_from_TLP_WD_limitTLPLSWD$R2adj)
  all_rmse_LMA <- c(LMA_from_TLP_LS_WD$rmse,LMA_from_TLP_LS$rmse,LMA_from_TLP_WD$rmse,LMA_from_LS_WD$rmse,LMA_from_TLP$rmse,LMA_from_LS$rmse,LMA_from_WD$rmse,LMA_from_LS_WD_limitTLPLSWD$rmse,LMA_from_TLP_WD_limitTLPLSWD$rmse)
  all_ndata_LMA <- c(LMA_from_TLP_LS_WD$ndata,LMA_from_TLP_LS$ndata,LMA_from_TLP_WD$ndata,LMA_from_LS_WD$ndata,LMA_from_TLP$ndata,LMA_from_LS$ndata,LMA_from_WD$ndata,LMA_from_LS_WD_limitTLPLSWD$ndata,LMA_from_TLP_WD_limitTLPLSWD$ndata)
  
  all_LMA <- data.frame(all_testnames_LMA,all_R2_LMA,all_R2adj_LMA,all_rmse_LMA,all_ndata_LMA)
   if(view_stats==TRUE)View(all_LMA)
  
  # CHOICE: best model in R2 terms is LMA_from_TLP_LS, but in RMSE is (marginally) LMA_from_TLP_LS_WD
  # Choose to go with LMA_from_TLP_LS on the basis that it is more parsimonious.
  
  # Test MAT and PPT coverage of species for chosen model
  plot(trait_BE$MAT[LMA_from_TLP_LS$dataused],trait_BE$MAP[LMA_from_TLP_LS$dataused])
  # RESULT: WIDE CLIMATE COVERAGE
  
  # DECISION: LMA_from_TLP_LS
  
  return_vals <- list("all_LMA"=all_LMA,
                      "LMA_from_TLP_LS_WD"=LMA_from_TLP_LS_WD)
  return(return_vals)
}

WD_multivar_test <- function(trait,view_stats=FALSE, leaf_type = NULL,regr_type = 'lm' )  {
  
  if(regr_type == 'lm' |regr_type == 'pcr'|regr_type == 'plsr'){# sma cannot handle more than 5 variables, so ecxlude from analysis here
  # WD from P50  slope Ks and LMA
  WD_from_P50_slope_Ks_LMA <- sma_plot_stats(data.frame(trait$P50,trait$slope,trait$Ks,trait$LMA,trait$WD),c("P50","slope","Ks","LMA","WD"),nbtstrp, regression_type = regr_type)
  plot(trait$WD[WD_from_P50_slope_Ks_LMA$dataused],WD_from_P50_slope_Ks_LMA$var_est,pch=16,xlab="WD",ylab="WD_est",main="WD vs WD_est")
  
  # WD from P50 Ks and LMA
  WD_from_P50_Ks_LMA_limitP50slopeKsLMA <- sma_plot_stats(data.frame(trait$P50,trait$Ks,trait$LMA,trait$WD),c("P50","Ks","LMA","WD"),nbtstrp,F, indin=WD_from_P50_slope_Ks_LMA$dataused, regression_type = regr_type)
  plot(trait$WD[WD_from_P50_Ks_LMA_limitP50slopeKsLMA$dataused],WD_from_P50_Ks_LMA_limitP50slopeKsLMA$var_est,pch=16,xlab="WD",ylab="WD_est",main="WD vs WD_est")
  
  }
  
  # WD from P50  slope and Ks
  WD_from_P50_slope_Ks <- sma_plot_stats(data.frame(trait$P50,trait$slope,trait$Ks,trait$WD),c("P50","slope","Ks","WD"),nbtstrp, regression_type = regr_type)
  plot(trait$WD[WD_from_P50_slope_Ks$dataused],WD_from_P50_slope_Ks$var_est,pch=16,xlab="WD",ylab="WD_est",main="WD vs WD_est")
  
  # WD from P50  LMA and Ks
  WD_from_P50_LMA_Ks <- sma_plot_stats(data.frame(trait$P50,trait$LMA,trait$Ks,trait$WD),c("P50","LMA","Ks","WD"),nbtstrp, regression_type = regr_type)
  plot(trait$WD[WD_from_P50_slope_Ks$dataused],WD_from_P50_slope_Ks$var_est,pch=16,xlab="WD",ylab="WD_est",main="WD vs WD_est")
  
  # WD from P50 and slope
  WD_from_P50_slope <- sma_plot_stats(data.frame(trait$P50,trait$slope,trait$WD),c("P50","slope","WD"),nbtstrp, regression_type = regr_type)
  plot(trait$WD[WD_from_P50_slope$dataused],WD_from_P50_slope$var_est,pch=16,xlab="WD",ylab="WD_est",main="WD vs WD_est")
  
  # WD from slope and P50*slope (interaction term) - NOTE: Testing the interactions because P50 and slope are so highly correlated with each other
  WD_from_slope_P50slope <- sma_plot_stats(data.frame(trait$slope,trait$P50*trait$slope,trait$WD),c("slope","P50*slope","WD"),nbtstrp,regression_type = regr_type)
  plot(trait$WD[WD_from_slope_P50slope$dataused],WD_from_slope_P50slope$var_est,pch=16,xlab="WD",ylab="WD_est",main="WD vs WD_est")
  
  # WD from slope Ks and P50*slope (interaction term) 
  WD_from_slope_Ks_P50slope <- sma_plot_stats(data.frame(trait$slope,trait$Ks,trait$P50*trait$slope,trait$WD),c("slope","Ks","P50*slope","WD"),nbtstrp,regression_type = regr_type)
  plot(trait$WD[WD_from_slope_Ks_P50slope$dataused],WD_from_slope_Ks_P50slope$var_est,pch=16,xlab="WD",ylab="WD_est",main="WD vs WD_est")
  
  # WD from P50 and P50*slope (interaction term)
  WD_from_P50_P50slope <- sma_plot_stats(data.frame(trait$P50,trait$P50*trait$slope,trait$WD),c("P50","P50*slope","WD"),nbtstrp,regression_type = regr_type)
  plot(trait$WD[WD_from_P50_P50slope$dataused],WD_from_P50_P50slope$var_est,pch=16,xlab="WD",ylab="WD_est",main="WD vs WD_est")
  
  if(regr_type=='lm'|regr_type == 'pcr'|regr_type == 'plsr'){
    #WD_from_slope_P50slope limit WD_from_slope_Ks_P50slope
    WD_from_slope_P50slope_limit_slope_Ks_P50slope <- sma_plot_stats(data.frame(trait$P50,trait$P50*trait$slope,trait$WD),c("P50","P50*slope","WD"),nbtstrp,T, WD_from_slope_Ks_P50slope$dataused,regression_type = regr_type)
  }
  # WD from P50 Ks
  WD_from_P50_Ks <- sma_plot_stats(data.frame(trait$P50,trait$Ks,trait$WD),c("P50","Ks","WD"),nbtstrp,regression_type = regr_type)
  plot(trait$WD[WD_from_P50_Ks$dataused],WD_from_P50_Ks$var_est,pch=16,xlab="WD",ylab="WD_est",main="WD vs WD_est")
  
  # WD from P50
  WD_from_P50 <- sma_plot_stats(data.frame(trait$P50,trait$WD),c("P50","WD"),nbtstrp,regression_type = regr_type)
  plot(trait$WD[WD_from_P50$dataused],WD_from_P50$var_est,pch=16,xlab="WD",ylab="WD_est",main="WD vs WD_est")
  
  # WD from slope
  WD_from_slope <- sma_plot_stats(data.frame(trait$slope,trait$WD),c("slope","WD"),nbtstrp,regression_type = regr_type)
  plot(trait$WD[WD_from_slope$dataused],WD_from_slope$var_est,pch=16,xlab="WD",ylab="WD_est",main="WD vs WD_est")
  
  # WD from Ks P50
  WD_from_Ks_P50 <- sma_plot_stats(data.frame(trait$Ks,trait$P50,trait$WD),c("Ks","P50","WD"),nbtstrp,regression_type = regr_type)
  plot(trait$WD[WD_from_Ks_P50$dataused],WD_from_Ks_P50$var_est,pch=16,xlab="WD",ylab="WD_est",main="WD vs WD_est")
  
  # WD_from P50 and Ks data limited to WD_from_P50slopeKs
  WD_from_P50_Ks_limitP50slopeKs  <- sma_plot_stats(data.frame(trait$P50,trait$Ks,trait$WD),c("P50","Ks","WD"),nbtstrp,F, indin= WD_from_P50_slope_Ks$dataused, regression_type = regr_type)
  plot(trait$WD[WD_from_P50_Ks_limitP50slopeKs$dataused],WD_from_P50_Ks_limitP50slopeKs$var_est,pch=16,xlab="WD",ylab="WD_est",main="WD vs WD_est")
  
  
  # Summarise statistics
  all_testnames_WD <- c("WD_from_P50_slope_Ks_LMA","WD_from_P50_Ks_LMA_limitP50slopeKsLMA","WD_from_P50_slope_Ks" ,  "WD_from_P50_LMA_Ks"          ,"WD_from_P50_slope",    "WD_from_slope_P50slope"," WD_from_slope_Ks_P50slope",     "WD_from_P50_P50slope"                ,"WD_from_P50_Ks"    ,"WD_from_P50"     ,"WD_from_slope",     "WD_from_Ks_P50","WD_from_P50_Ks_limitP50slopeKs")
  all_R2_WD <- c(WD_from_P50_slope_Ks_LMA$R2       ,WD_from_P50_Ks_LMA_limitP50slopeKsLMA$R2,WD_from_P50_slope_Ks$R2,  WD_from_P50_LMA_Ks$R2        ,WD_from_P50_slope$R2,         WD_from_slope_P50slope$R2, WD_from_slope_Ks_P50slope$R2,      WD_from_P50_P50slope$R2, WD_from_P50_Ks$R2   ,WD_from_P50$R2    ,WD_from_slope$R2,    WD_from_Ks_P50$R2,WD_from_P50_Ks_limitP50slopeKs$R2)
  all_R2adj_WD <- c(WD_from_P50_slope_Ks_LMA$R2adj,WD_from_P50_Ks_LMA_limitP50slopeKsLMA$R2adj,WD_from_P50_slope_Ks$R2adj,WD_from_P50_LMA_Ks$R2adj     ,WD_from_P50_slope$R2adj,WD_from_slope_P50slope$R2adj,  WD_from_slope_Ks_P50slope$R2adj, WD_from_P50_P50slope$R2adj    ,WD_from_P50_Ks$R2adj ,WD_from_P50$R2adj ,WD_from_slope$R2adj, WD_from_Ks_P50$R2adj,WD_from_P50_Ks_limitP50slopeKs$R2adj)
  all_rmse_WD <- c(WD_from_P50_slope_Ks_LMA$rmse,  WD_from_P50_Ks_LMA_limitP50slopeKsLMA$rmse,WD_from_P50_slope_Ks$rmse, WD_from_P50_LMA_Ks$rmse      ,WD_from_P50_slope$rmse,   WD_from_slope_P50slope$rmse, WD_from_slope_Ks_P50slope$rmse,   WD_from_P50_P50slope$rmse    , WD_from_P50_Ks$rmse, WD_from_P50$rmse   ,WD_from_slope$rmse,  WD_from_Ks_P50$rmse,WD_from_P50_Ks_limitP50slopeKs$rmse)
  all_ndata_WD <- c(WD_from_P50_slope_Ks_LMA$ndata,WD_from_P50_Ks_LMA_limitP50slopeKsLMA$ndata,WD_from_P50_slope_Ks$ndata,WD_from_P50_LMA_Ks$ndata     , WD_from_P50_slope$ndata,WD_from_slope_P50slope$ndata, WD_from_slope_Ks_P50slope$ndata,  WD_from_P50_P50slope$ndata , WD_from_P50_Ks$ndata  ,WD_from_P50$ndata  ,WD_from_slope$ndata, WD_from_Ks_P50$ndata,WD_from_P50_Ks_limitP50slopeKs$ndata)
  
  all_WD <- data.frame(all_testnames_WD,all_R2_WD,all_R2adj_WD,all_rmse_WD,all_ndata_WD)
  if(view_stats==TRUE)View(all_WD)
  
  # CHOICE: WD_from_slope_P50slope has the best combination of R2adj and RMSE
  
  # Test MAT and PPT coverage of species for chosen model
  plot(trait$MAT[WD_from_slope_P50slope$dataused],trait$MAP[WD_from_slope_P50slope$dataused])
  # WIDE CLIMATE COVERAGE
  
  # DECISION: WD_from_slope_P50slope
  
  return_vals <- list("all_WD"=all_WD,
                      "WD_from_P50_slope_Ks"=WD_from_P50_slope_Ks)
  #new:
  if(leaf_type== 'BDT'){
    # Test MAT and PPT coverage of species for chosen model
    plot(trait$MAT[WD_from_P50_slope_Ks$dataused],trait$MAP[WD_from_P50_slope_Ks$dataused],main = "Climate coverage")
    # WIDE CLIMATE COVERAGE
    
    # DECISION for BDT: :
    return_vals <- list("all_WD"=all_WD,
                        "WD_from_P50_Ks"= WD_from_P50_Ks)  # WD_from_P50_Ks 
  }
  
  if(leaf_type == 'BE'){
    # Test MAT and PPT coverage of species for chosen model
    plot(trait$MAT[WD_from_P50_slope_Ks_LMA$dataused],trait$MAP[WD_from_P50_slope_Ks_LMA$dataused],main = "Climate coverage")
    # WIDE CLIMATE COVERAGE
    
    # DECISION for BE: LS_from_LMA_TLP_Ks
    
    return_vals <- list("all_WD"= all_WD,
                        "WD_from_P50_LMA_Ks"= WD_from_P50_LMA_Ks) # WD_from_P50_LMA_Ks
  }
  
  return(return_vals)
}

slope_multivar_test <- function(trait,view_stats=FALSE, regr_type = 'lm') {
  
  # slope from P50, TLP, WD and Ks
  slope_from_P50_TLP_WD_Ks <- sma_plot_stats(data.frame(trait$P50,trait$TLP,trait$WD,trait$Ks,trait$slope),c("P50","TLP","WD","Ks","slope"),nbtstrp,regression_type = regr_type)
  plot(trait$slope[slope_from_P50_TLP_WD_Ks$dataused],slope_from_P50_TLP_WD_Ks$var_est,pch=16,xlab="slope",ylab="slope_est",main="slope vs slope_est")
  
  # slope from P50, TLP and WD
  slope_from_P50_TLP_WD <- sma_plot_stats(data.frame(trait$P50,trait$TLP,trait$WD,trait$slope),c("P50","TLP","WD","slope"),nbtstrp,regression_type = regr_type)
  plot(trait$slope[slope_from_P50_TLP_WD$dataused],slope_from_P50_TLP_WD$var_est,pch=16,xlab="slope",ylab="slope_est",main="slope vs slope_est")
  
  # slope from P50, TLP and Ks
  slope_from_P50_TLP_Ks <- sma_plot_stats(data.frame(trait$P50,trait$TLP,trait$Ks,trait$slope),c("P50","TLP","Ks","slope"),nbtstrp,regression_type = regr_type)
  plot(trait$slope[slope_from_P50_TLP_Ks$dataused],slope_from_P50_TLP_Ks$var_est,pch=16,xlab="slope",ylab="slope_est",main="slope vs slope_est")
  
  # slope from P50, WD and Ks
  slope_from_P50_WD_Ks <- sma_plot_stats(data.frame(trait$P50,trait$WD,trait$Ks,trait$slope),c("P50","WD","Ks","slope"),nbtstrp,regression_type = regr_type)
  plot(trait$slope[slope_from_P50_WD_Ks$dataused],slope_from_P50_WD_Ks$var_est,pch=16,xlab="slope",ylab="slope_est",main="slope vs slope_est")
  
  # slope from TLP, WD and Ks
  slope_from_TLP_WD_Ks <- sma_plot_stats(data.frame(trait$TLP,trait$WD,trait$Ks,trait$slope),c("TLP","WD","Ks","slope"),nbtstrp,regression_type = regr_type)
  plot(trait$slope[slope_from_TLP_WD_Ks$dataused],slope_from_TLP_WD_Ks$var_est,pch=16,xlab="slope",ylab="slope_est",main="slope vs slope_est")
  
  # slope from P50 and TLP
  slope_from_P50_TLP <- sma_plot_stats(data.frame(trait$P50,trait$TLP,trait$slope),c("P50","TLP","slope"),nbtstrp,regression_type = regr_type)
  plot(trait$slope[slope_from_P50_TLP$dataused],slope_from_P50_TLP$var_est,pch=16,xlab="slope",ylab="slope_est",main="slope vs slope_est")
  
  # slope from P50 and TLPP50 (interaction)
  slope_from_P50_TLPP50 <- sma_plot_stats(data.frame(trait$P50,trait$TLP*trait$P50,trait$slope),c("P50","TLP*P50","slope"),nbtstrp,regression_type = regr_type)
  plot(trait$slope[slope_from_P50_TLPP50$dataused],slope_from_P50_TLPP50$var_est,pch=16,xlab="slope",ylab="slope_est",main="slope vs slope_est")
  
  # slope from TLP and TLPP50 (interaction)
  slope_from_TLP_TLPP50 <- sma_plot_stats(data.frame(trait$TLP,trait$TLP*trait$P50,trait$slope),c("TLP","TLP*P50","slope"),nbtstrp,regression_type = regr_type)
  plot(trait$slope[slope_from_TLP_TLPP50$dataused],slope_from_TLP_TLPP50$var_est,pch=16,xlab="slope",ylab="slope_est",main="slope vs slope_est")
  
  # slope from TLP and Ks
  slope_from_TLP_Ks <- sma_plot_stats(data.frame(trait$TLP,trait$Ks,trait$slope),c("TLP","Ks","slope"),nbtstrp,regression_type = regr_type)
  plot(trait$slope[slope_from_TLP_Ks$dataused],slope_from_TLP_Ks$var_est,pch=16,xlab="slope",ylab="slope_est",main="slope vs slope_est")
  
  # slope from TLP and WD
  slope_from_TLP_WD <- sma_plot_stats(data.frame(trait$TLP,trait$WD,trait$slope),c("TLP","WD","slope"),nbtstrp,regression_type = regr_type)
  plot(trait$slope[slope_from_TLP_WD$dataused],slope_from_TLP_WD$var_est,pch=16,xlab="slope",ylab="slope_est",main="slope vs slope_est")
  
  # slope from TLP
  slope_from_TLP <- sma_plot_stats(data.frame(trait$TLP,trait$slope),c("TLP","slope"),nbtstrp,regression_type = regr_type)
  plot(trait$slope[slope_from_TLP$dataused],slope_from_TLP$var_est,pch=16,xlab="slope",ylab="slope_est",main="slope vs slope_est")
  
  # slope from P50
  slope_from_P50 <- sma_plot_stats(data.frame(trait$P50,trait$slope),c("P50","slope"),nbtstrp,regression_type = regr_type)
  plot(trait$slope[slope_from_P50$dataused],slope_from_P50$var_est,pch=16,xlab="slope",ylab="slope_est",main="slope vs slope_est")
  
  # slope from WD
  slope_from_WD <- sma_plot_stats(data.frame(trait$WD,trait$slope),c("WD","slope"),nbtstrp,regression_type = regr_type)
  plot(trait$slope[slope_from_WD$dataused],slope_from_WD$var_est,pch=16,xlab="slope",ylab="slope_est",main="slope vs slope_est")
  
  # slope from Ks
  slope_from_Ks <- sma_plot_stats(data.frame(trait$Ks,trait$slope),c("Ks","slope"),nbtstrp,regression_type = regr_type)
  plot(trait$slope[slope_from_Ks$dataused],slope_from_Ks$var_est,pch=16,xlab="slope",ylab="slope_est",main="slope vs slope_est")
  

  # Some of the above combinations have relatively high R2adj, but very reduced data points. Therefore test if the worse-performing combinations are only worse performing because they include a greater diversity of data.
  
  # slope from P50, WD and Ks (same species as slope_from_P50_TLP_WD_Ks) (NOTE: testing because the fit for slope_from_P50_WD_Ks is only marginally worse than for slope_from_P50_TLP_Ks, but it has many more species)
  slope_from_P50_WD_Ks_limitP50TLPWDKs <- sma_plot_stats(data.frame(trait$P50,trait$WD,trait$Ks,trait$slope),c("P50","WD","Ks","slope"),nbtstrp,F,indin = slope_from_P50_TLP_WD_Ks$dataused,regression_type = regr_type)
  plot(trait$slope[slope_from_P50_WD_Ks_limitP50TLPWDKs$dataused],slope_from_P50_WD_Ks_limitP50TLPWDKs$var_est,pch=16,xlab="slope",ylab="slope_est",main="slope vs slope_est")
  
  # slope from P50, TLP and Ks (same species as slope_from_P50_TLP_WD_Ks) (NOTE: testing because the fit for slope_from_P50_WD_Ks is only marginally worse than for slope_from_P50_TLP_Ks, but it has many more species)
  slope_from_P50_TLP_Ks_limitP50TLPWDKs <- sma_plot_stats(data.frame(trait$P50,trait$TLP,trait$Ks,trait$slope),c("P50","TLP","Ks","slope"),nbtstrp,F,indin = slope_from_P50_TLP_WD_Ks$dataused,regression_type = regr_type)
  plot(trait$slope[slope_from_P50_TLP_Ks_limitP50TLPWDKs$dataused],slope_from_P50_TLP_Ks_limitP50TLPWDKs$var_est,pch=16,xlab="slope",ylab="slope_est",main="slope vs slope_est")
  
  # slope from P50 (same species as slope_from_P50_TLP_WD_Ks)
  slope_from_P50_limitP50TLPWDKs <- sma_plot_stats(data.frame(trait$P50,trait$slope),c("P50","slope"),nbtstrp,F,indin = slope_from_P50_TLP_WD_Ks$dataused, regression_type = regr_type)
  plot(trait$slope[slope_from_P50_limitP50TLPWDKs$dataused],slope_from_P50_limitP50TLPWDKs$var_est,pch=16,xlab="slope",ylab="slope_est",main="slope vs slope_est")
  
  slope_from_P50_TLP_WD_limitP50TLPWDKs <- sma_plot_stats(data.frame(trait$P50,trait$TLP,trait$WD,trait$slope),c("P50","TLP","WD","slope"),nbtstrp,F, indin= slope_from_P50_TLP_WD_Ks$dataused, regression_type = regr_type)
  plot(trait$slope[slope_from_P50_TLP_WD_limitP50TLPWDKs$dataused],slope_from_P50_TLP_WD_limitP50TLPWDKs$var_est,pch=16,xlab="slope",ylab="slope_est",main="slope vs slope_est")
  
  slope_from_Ks_limitP50TLPWDKs <- sma_plot_stats(data.frame(trait$Ks,trait$slope),c("Ks","slope"),nbtstrp,F, indin= slope_from_P50_TLP_WD_Ks$dataused, regression_type = regr_type)
  plot(trait$slope[slope_from_Ks_limitP50TLPWDKs$dataused],slope_from_Ks_limitP50TLPWDKs$var_est,pch=16,xlab="slope",ylab="slope_est",main="slope vs slope_est")
  
  
  # Summarise statistics
  all_testnames_slope <- c("slope_from_P50_TLP_WD_Ks","slope_from_P50_TLP_WD","slope_from_P50_TLP_Ks","slope_from_P50_WD_Ks","slope_from_TLP_WD_Ks","slope_from_P50_TLP","slope_from_P50_TLPP50","slope_from_TLP_TLPP50","slope_from_TLP_Ks","slope_from_TLP_WD","slope_from_TLP","slope_from_P50","slope_from_WD","slope_from_Ks","slope_from_P50_WD_Ks_limitP50TLPWDKs","slope_from_P50_TLP_Ks_limitP50TLPWDKs","slope_from_P50_limitP50TLPWDKs","slope_from_P50_TLP_WD_limitP50TLPWDKs","slope_from_Ks_limitP50TLPWDKs")
  all_R2_slope <- c(slope_from_P50_TLP_WD_Ks$R2,slope_from_P50_TLP_WD$R2,slope_from_P50_TLP_Ks$R2,slope_from_P50_WD_Ks$R2,slope_from_TLP_WD_Ks$R2,slope_from_P50_TLP$R2,slope_from_P50_TLPP50$R2,slope_from_TLP_TLPP50$R2,slope_from_TLP_Ks$R2,slope_from_TLP_WD$R2,slope_from_TLP$R2,slope_from_P50$R2,slope_from_WD$R2,slope_from_Ks$R2,slope_from_P50_WD_Ks_limitP50TLPWDKs$R2,slope_from_P50_TLP_Ks_limitP50TLPWDKs$R2,slope_from_P50_limitP50TLPWDKs$R2,slope_from_P50_TLP_WD_limitP50TLPWDKs$R2,slope_from_Ks_limitP50TLPWDKs$R2)
  all_R2adj_slope <- c(slope_from_P50_TLP_WD_Ks$R2adj,slope_from_P50_TLP_WD$R2adj,slope_from_P50_TLP_Ks$R2adj,slope_from_P50_WD_Ks$R2adj,slope_from_TLP_WD_Ks$R2adj,slope_from_P50_TLP$R2adj,slope_from_P50_TLPP50$R2adj,slope_from_TLP_TLPP50$R2adj,slope_from_TLP_Ks$R2adj,slope_from_TLP_WD$R2adj,slope_from_TLP$R2adj,slope_from_P50$R2adj,slope_from_WD$R2adj,slope_from_Ks$R2adj,slope_from_P50_WD_Ks_limitP50TLPWDKs$R2adj,slope_from_P50_TLP_Ks_limitP50TLPWDKs$R2adj,slope_from_P50_limitP50TLPWDKs$R2adj,slope_from_P50_TLP_WD_limitP50TLPWDKs$R2adj,slope_from_Ks_limitP50TLPWDKs$R2adj)
  all_rmse_slope <- c(slope_from_P50_TLP_WD_Ks$rmse,slope_from_P50_TLP_WD$rmse,slope_from_P50_TLP_Ks$rmse,slope_from_P50_WD_Ks$rmse,slope_from_TLP_WD_Ks$rmse,slope_from_P50_TLP$rmse,slope_from_P50_TLPP50$rmse,slope_from_TLP_TLPP50$rmse,slope_from_TLP_Ks$rmse,slope_from_TLP_WD$rmse,slope_from_TLP$rmse,slope_from_P50$rmse,slope_from_WD$rmse,slope_from_Ks$rmse,slope_from_P50_WD_Ks_limitP50TLPWDKs$rmse,slope_from_P50_TLP_Ks_limitP50TLPWDKs$rmse,slope_from_P50_limitP50TLPWDKs$rmse,slope_from_P50_TLP_WD_limitP50TLPWDKs$rmse,slope_from_Ks_limitP50TLPWDKs$rmse)
  all_ndata_slope <- c(slope_from_P50_TLP_WD_Ks$ndata,slope_from_P50_TLP_WD$ndata,slope_from_P50_TLP_Ks$ndata,slope_from_P50_WD_Ks$ndata,slope_from_TLP_WD_Ks$ndata,slope_from_P50_TLP$ndata,slope_from_P50_TLPP50$ndata,slope_from_TLP_TLPP50$ndata,slope_from_TLP_Ks$ndata,slope_from_TLP_WD$ndata,slope_from_TLP$ndata,slope_from_P50$ndata,slope_from_WD$ndata,slope_from_Ks$ndata,slope_from_P50_WD_Ks_limitP50TLPWDKs$ndata,slope_from_P50_TLP_Ks_limitP50TLPWDKs$ndata,slope_from_P50_limitP50TLPWDKs$ndata,slope_from_P50_TLP_WD_limitP50TLPWDKs$ndata,slope_from_Ks_limitP50TLPWDKs$ndata)
  
  all_slope <- data.frame(all_testnames_slope,all_R2_slope,all_R2adj_slope,all_rmse_slope,all_ndata_slope)
   if(view_stats==TRUE)View(all_slope)
  
  # CHOICE: slope_from_P50_TLP_Ks
  
  # Test MAT and PPT coverage of species for chosen model
  plot(trait$MAT[slope_from_P50_TLP_WD$dataused],trait$MAP[slope_from_P50_TLP_WD$dataused])
  #plot(trait$MAT[slope_from_P50_TLP_WD_Ks$dataused],trait$MAP[slope_from_P50_TLP_WD_Ks$dataused])
  # WIDE CLIMATE COVERAGE
  
  # DECISION: slope_from_P50_TLP_WD
  
  return_vals <- list("all_slope"=all_slope,
                      "slope_from_P50_TLP_WD" = slope_from_P50_TLP_WD)
  
  return(return_vals)
}
  
Ks_multivar_test <- function(trait,view_stats=FALSE, regr_type = 'lm') {
  # not really a test, as the trait network is already established, but testing the thin lines in slide 10, too:
  
  # Ks from P50, TLP, LS
  Ks_from_P50_TLP_LS <- sma_plot_stats(data.frame(trait$P50,trait$TLP,trait$LS,trait$Ks),c("P50","TLP","LS","Ks"),nbtstrp,regression_type = regr_type)
  plot(trait$Ks[ Ks_from_P50_TLP_LS$dataused], Ks_from_P50_TLP_LS$var_est,pch=16,xlab="Ks",ylab="Ks_est",main="Ks vs Ks_est")
  
  # Ks from P50, LS, LMA
  Ks_from_P50_LS_LMA <- sma_plot_stats(data.frame(trait$P50,trait$LS,trait$LMA,trait$Ks),c("P50","LS","LMA","Ks"),nbtstrp,regression_type = regr_type)
  plot(trait$Ks[Ks_from_P50_LS_LMA$dataused],Ks_from_P50_LS_LMA$var_est,pch=16,xlab="Ks",ylab="Ks_est",main="Ks vs Ks_est")
  
  # Ks from P50, LS, slope
  Ks_from_P50_LS_slope <- sma_plot_stats(data.frame(trait$P50,trait$LS,trait$slope,trait$Ks),c("P50","LS","slope","Ks"),nbtstrp,regression_type = regr_type)
  plot(trait$Ks[Ks_from_P50_LS_slope$dataused],Ks_from_P50_LS_slope$var_est,pch=16,xlab="Ks",ylab="Ks_est",main="Ks vs Ks_est")
  
  # Ks from P50, LS, limited by data from Ks_from_P50_LS_slope only
  Ks_from_P50_LS_limitP50LSslope <- sma_plot_stats(data.frame(trait$P50,trait$LS,trait$Ks),c("P50","LS","Ks"),nbtstrp,makeplot=F,indin = Ks_from_P50_LS_slope$dataused, regression_type =  regr_type)
  plot(trait$Ks[Ks_from_P50_LS_limitP50LSslope$dataused],Ks_from_P50_LS_limitP50LSslope$var_est,pch=16,xlab="Ks",ylab="Ks_est",main="Ks vs Ks_est")
  
  # Ks from P50, LS, slope and WD
  Ks_from_P50_LS_slope_WD <- sma_plot_stats(data.frame(trait$P50,trait$LS,trait$slope,trait$WD,trait$Ks),c("P50","LS","slope","WD","Ks"),nbtstrp,regression_type = regr_type)
  plot(trait$Ks[Ks_from_P50_LS_slope_WD$dataused],Ks_from_P50_LS_slope_WD$var_est,pch=16,xlab="Ks",ylab="Ks_est",main="Ks vs Ks_est")
  
  # Ks from P50, LS and WD
  Ks_from_P50_LS_WD <- sma_plot_stats(data.frame(trait$P50,trait$LS,trait$WD,trait$Ks),c("P50","LS","WD","Ks"),nbtstrp,regression_type = regr_type)
  plot(trait$Ks[Ks_from_P50_LS_WD$dataused],Ks_from_P50_LS_WD$var_est,pch=16,xlab="Ks",ylab="Ks_est",main="Ks vs Ks_est")
  
  # Ks from P50, LS, slope limited by P50_LS_slope_WD
  Ks_from_P50_LS_slope_limitP50LSslopeWD <- sma_plot_stats(data.frame(trait$P50,trait$LS,trait$slope,trait$Ks),c("P50","LS","slope","Ks"),nbtstrp,indin=Ks_from_P50_LS_slope_WD$dataused,regression_type = regr_type)
  plot(trait$Ks[Ks_from_P50_LS_slope_limitP50LSslopeWD $dataused],Ks_from_P50_LS_slope_limitP50LSslopeWD $var_est,pch=16,xlab="Ks",ylab="Ks_est",main="Ks vs Ks_est")
  
  # Ks from P50, LS, WD limited by P50_LS_slope_WD
  Ks_from_P50_LS_WD_limitP50LSslopeWD <- sma_plot_stats(data.frame(trait$P50,trait$LS,trait$WD,trait$Ks),c("P50","LS","WD","Ks"),nbtstrp,indin=Ks_from_P50_LS_slope_WD$dataused,regression_type = regr_type)
  plot(trait$Ks[Ks_from_P50_LS_WD_limitP50LSslopeWD $dataused],Ks_from_P50_LS_WD_limitP50LSslopeWD $var_est,pch=16,xlab="Ks",ylab="Ks_est",main="Ks vs Ks_est")
  
  
  # Ks from P50, LS, LMA  limited by data from  Ks_from_P50_LS_WD only
  Ks_from_P50_LS_limitPSLSWD <- sma_plot_stats(data.frame(trait$P50,trait$LS,trait$Ks),c("P50","LS","Ks"),nbtstrp,makeplot=F,indin =  Ks_from_P50_LS_WD$dataused, regression_type =  regr_type)
  plot(trait$Ks[Ks_from_P50_LS_limitPSLSWD$dataused],Ks_from_P50_LS_limitPSLSWD$var_est,pch=16,xlab="Ks",ylab="Ks_est",main="Ks vs Ks_est")
  
  # Ks from P50, LS
  Ks_from_P50_LS <- sma_plot_stats(data.frame(trait$P50,trait$LS,trait$Ks),c("P50","LS","Ks"),nbtstrp,regression_type = regr_type)
  plot(trait$Ks[Ks_from_P50_LS$dataused],Ks_from_P50_LS$var_est,pch=16,xlab="Ks",ylab="Ks_est",main="Ks vs Ks_est")
  
  # Ks from LS
  Ks_from_LS <- sma_plot_stats(data.frame(trait$LS,trait$Ks),c("LS","Ks"),nbtstrp,regression_type = regr_type)
  plot(trait$Ks[ Ks_from_LS$dataused], Ks_from_LS$var_est,pch=16,xlab="KS",ylab="KS_est",main="KS vs Ks_est")
  
  # Ks from P50 
  Ks_from_P50 <- sma_plot_stats(data.frame(trait$P50,trait$Ks),c("P50","Ks"),nbtstrp,regression_type = regr_type)
  plot(trait$Ks[Ks_from_P50$dataused], Ks_from_P50$var_est,pch=16,xlab="Ks",ylab="Ks_est",main="Ks vs Ks_est")
  
  # Ks from LS
  Ks_from_LS_limitP50LS <- sma_plot_stats(data.frame(trait$LS,trait$Ks),c("LS","Ks"),nbtstrp,makeplot=F,indin = Ks_from_P50_LS$dataused,regression_type = regr_type)
  plot(trait$Ks[ Ks_from_LS$dataused], Ks_from_LS$var_est,pch=16,xlab="KS",ylab="KS_est",main="KS vs Ks_est")

  
  # Summarise statistics
  all_testnames_Ks <- c("Ks_from_P50_TLP_LS ","Ks_from_P50_LS_LMA",    "Ks_from_P50_LS_slope",   "Ks_from_P50_LS",    "Ks_from_P50_LS_WD","Ks_from_P50_LS_limitPSLSWD",          "Ks_from_P50_LS_slope_WD","Ks_from_P50_LS_limitP50LSslope",        "Ks_from_LS ",   "Ks_from_P50",      "Ks_from_LS_limitP50LS")
  all_R_Ks <- c(    Ks_from_P50_TLP_LS$R,     Ks_from_P50_LS_LMA$R,    Ks_from_P50_LS_slope$R,    Ks_from_P50_LS$R,   Ks_from_P50_LS_WD$R,Ks_from_P50_LS_limitPSLSWD$R,          Ks_from_P50_LS_slope_WD$R,  Ks_from_P50_LS_limitP50LSslope$R,         Ks_from_LS$R,     Ks_from_P50$R,     Ks_from_LS_limitP50LS$R)
  all_R2adj_Ks <- c(Ks_from_P50_TLP_LS$R2adj, Ks_from_P50_LS_LMA$R2adj,Ks_from_P50_LS_slope$R2,   Ks_from_P50_LS$R2adj, Ks_from_P50_LS_WD$R2adj,Ks_from_P50_LS_limitPSLSWD$R2adj,Ks_from_P50_LS_slope_WD$R2adj, Ks_from_P50_LS_limitP50LSslope$R2adj,         Ks_from_LS$R2adj, Ks_from_P50$R2adj, Ks_from_LS_limitP50LS$R2adj)
  all_rmse_Ks <- c( Ks_from_P50_TLP_LS$rmse,  Ks_from_P50_LS_LMA$rmse, Ks_from_P50_LS_slope$rmse, Ks_from_P50_LS$rmse,  Ks_from_P50_LS_WD$rmse,Ks_from_P50_LS_limitPSLSWD$rmse,  Ks_from_P50_LS_slope_WD$rmse, Ks_from_P50_LS_limitP50LSslope$rmse,         Ks_from_LS$rmse,  Ks_from_P50$rmse,  Ks_from_LS_limitP50LS$rmse)
  all_ndata_Ks <- c(Ks_from_P50_TLP_LS$ndata, Ks_from_P50_LS_LMA$ndata,Ks_from_P50_LS_slope$ndata,Ks_from_P50_LS$ndata,  Ks_from_P50_LS_WD$ndata,Ks_from_P50_LS_limitPSLSWD$ndata,Ks_from_P50_LS_slope_WD$ndata,Ks_from_P50_LS_limitP50LSslope$ndata,         Ks_from_LS$ndata, Ks_from_P50$ndata, Ks_from_LS_limitP50LS$ndata)
  
  all_Ks <- data.frame(all_testnames_Ks,all_R_Ks,all_R2adj_Ks,all_rmse_Ks,all_ndata_Ks)
   if(view_stats==TRUE)View(all_Ks)
  
  #CHOICE: Ks_from_LS has more data points, but marginally worse fit. Check whether both trait combination options
  # span the whole climate space. 
    
 # plot(trait$MAT[Ks_from_LS$dataused],trait$MAP[Ks_from_LS$dataused],main = "Climate coverage")  
  # Test MAT and PPT coverage of species for chosen model
  plot(trait$MAT[Ks_from_P50_LS_WD$dataused],trait$MAP[Ks_from_P50_LS_WD$dataused],main = "Climate coverage")
  # WIDE CLIMATE COVERAGE: The climate coverage is similar for both, despite less data for Ks_from_P50_LS.
  
  # DECISION: Ks_from_P50_LS_WD
  
  return_vals <- list("all_Ks"= all_Ks,
                      "Ks_from_P50_LS_WD"= Ks_from_P50_LS_WD)
  
  return(return_vals)
}

LS_multivar_test <- function(trait, view_stats = FALSE, leaf_type = NULL, regr_type = 'lm') {
# input:  leaf_type 'BDT' or 'BE', so that automatically the leaf type-specific model is selected. (BE: LS_from_TLP_Ks BDT: )
  
  if(is.null(leaf_type)){
    stop('must provide leaf_type e.g BDT or BE')
  }
  
  # not really a test, as the trait network is already established, but testing the thin lines in slide 10, too:
  
  if( regr_type=='lm'|regr_type == 'pcr'|regr_type == 'plsr'){  # smal regr cannot handle more 4 variables, so exclude there.
  # LS from LMA TLP, KS, P50
  LS_from_P50_TLP_Ks_LMA <- sma_plot_stats(data.frame(trait$P50,trait$TLP,trait$Ks,trait$LMA,trait$LS),c("P50","TLP","Ks","LMA","LS"),nbtstrp, regression_type = regr_type )
  plot(trait$LS[LS_from_P50_TLP_Ks_LMA$dataused], LS_from_P50_TLP_Ks_LMA$var_est,pch=16,xlab="LS",ylab="LS_est",main="LS vs LS_est")
  
  # LS from TLP, KS, LMA  ( samples from  LS_from_P50_TLP_Ks_LMA )
  LS_from_TLP_Ks_LMA_limitP50TLPKsMA <- sma_plot_stats(data.frame(trait$TLP,trait$Ks,trait$LMA,trait$LS),c("TLP","Ks","LMA","LS"),nbtstrp,F, LS_from_P50_TLP_Ks_LMA$dataused, regression_type =  regr_type)
  plot(trait$LS[LS_from_TLP_Ks_LMA_limitP50TLPKsMA$dataused], LS_from_TLP_Ks_LMA_limitP50TLPKsMA$var_est,pch=16,xlab="LS",ylab="LS_est",main="LS vs LS_est")
  }

  # LS from LMA TLP
  LS_from_LMA_TLP<- sma_plot_stats(data.frame(trait$LMA,trait$TLP,trait$LS),c("LMA","TLP","LS"),nbtstrp, regression_type =  regr_type)
  plot(trait$LS[LS_from_LMA_TLP$dataused], LS_from_LMA_TLP$var_est,pch=16,xlab="LS",ylab="LS_est",main="LS vs LS_est")
  
  # LS from P50, TLP and Ks
  LS_from_P50_TLP_Ks <- sma_plot_stats(data.frame(trait$P50,trait$TLP,trait$Ks,trait$LS),c("P50","TLP","Ks","LS"),nbtstrp, regression_type =  regr_type)
  plot(trait$LS[LS_from_P50_TLP_Ks$dataused],LS_from_P50_TLP_Ks$var_est,pch=16,xlab="LS",ylab="LS_est",main="LS vs LS_est")
  
  # LS from LMA, TLP and Ks
  LS_from_LMA_TLP_Ks <- sma_plot_stats(data.frame(trait$LMA,trait$TLP,trait$Ks,trait$LS),c("LMA","TLP","Ks","LS"),nbtstrp, regression_type =  regr_type)
  plot(trait$LS[LS_from_LMA_TLP_Ks$dataused],LS_from_LMA_TLP_Ks$var_est,pch=16,xlab="LS",ylab="LS_est",main="LS vs LS_est")
  
  if(regr_type=='lm'|regr_type == 'pcr'|regr_type == 'plsr'){
  # LS from TLP and Ks, limit to data from LS_from_LMA_TLP_Ks
  LS_from_TLP_Ks_limitLMA_TLP_Ks <- sma_plot_stats(data.frame(trait$TLP,trait$Ks,trait$LS),c("TLP","Ks","LS"),nbtstrp,  indin= LS_from_LMA_TLP_Ks$dataused, regression_type =  regr_type)
  plot(trait$LS[LS_from_TLP_Ks_limitLMA_TLP_Ks$dataused],LS_from_TLP_Ks_limitLMA_TLP_Ks$var_est,pch=16,xlab="LS",ylab="LS_est",main="LS vs LS_est")
  } # smal regr cannot handle more 4 variables

  # LS from TLP and Ks
  LS_from_TLP_Ks <- sma_plot_stats(data.frame(trait$TLP,trait$Ks,trait$LS),c("TLP","Ks","LS"),nbtstrp, regression_type =  regr_type)
  plot(trait$LS[LS_from_TLP_Ks$dataused],LS_from_TLP_Ks$var_est,pch=16,xlab="LS",ylab="LS_est",main="LS vs LS_est")
  
  # LS from TLP and LMA
  LS_from_TLP_LMA <- sma_plot_stats(data.frame(trait$TLP,trait$LMA,trait$LS),c("TLP","LMA","LS"),nbtstrp, regression_type =  regr_type)
  plot(trait$LS[LS_from_TLP_LMA$dataused],LS_from_TLP_LMA$var_est,pch=16,xlab="LS",ylab="LS_est",main="LS vs LS_est")
  
  # LS from P50, TLP
  LS_from_P50_TLP <- sma_plot_stats(data.frame(trait$P50,trait$TLP,trait$LS),c("P50","TLP","LS"),nbtstrp, regression_type =  regr_type)
  plot(trait$LS[LS_from_P50_TLP$dataused],LS_from_P50_TLP$var_est,pch=16,xlab="LS",ylab="LS_est",main="LS vs LS_est")
  
  # LS from P50, Ks
  LS_from_P50_Ks <- sma_plot_stats(data.frame(trait$P50,trait$Ks,trait$LS),c("P50","Ks","LS"),nbtstrp, regression_type =  regr_type)
  plot(trait$LS[LS_from_P50_Ks$dataused],LS_from_P50_Ks$var_est,pch=16,xlab="LS",ylab="LS_est",main="LS vs LS_est")
  
  # LS from P50
  LS_from_P50 <- sma_plot_stats(data.frame(trait$P50,trait$LS),c("P50","LS"),nbtstrp, regression_type =  regr_type)
  plot(trait$LS[LS_from_P50$dataused],LS_from_P50$var_est,pch=16,xlab="LS",ylab="LS_est",main="LS vs LS_est")
  
  # LS from TLP
  LS_from_TLP <- sma_plot_stats(data.frame(trait$TLP,trait$LS),c("TLP","LS"),nbtstrp, F, indin = LS_from_LMA_TLP_Ks$dataused, regression_type =  regr_type)
  plot(trait$LS[LS_from_TLP$dataused],LS_from_TLP$var_est,pch=16,xlab="LS",ylab="LS_est",main="LS vs LS_est")
  
  # LS from Ks
  LS_from_Ks <- sma_plot_stats(data.frame(trait$Ks,trait$LS),c("Ks","LS"),nbtstrp, regression_type =  regr_type)
  plot(trait$LS[LS_from_Ks$dataused],LS_from_Ks$var_est,pch=16,xlab="LS",ylab="LS_est",main="LS vs LS_est")
  
  # LS from P50
  LS_from_P50 <- sma_plot_stats(data.frame(trait$P50,trait$LS),c("P50","LS"),nbtstrp, regression_type =  regr_type)
  plot(trait$LS[LS_from_P50$dataused],LS_from_P50$var_est,pch=16,xlab="LS",ylab="LS_est",main="LS vs LS_est")
  
  # LS from LMA
  LS_from_LMA <- sma_plot_stats(data.frame(trait$LMA,trait$LS),c("LMA","LS"),nbtstrp, F, indin = LS_from_LMA_TLP_Ks$dataused, regression_type =  regr_type)
  plot(trait$LS[LS_from_LMA$dataused],LS_from_LMA$var_est,pch=16,xlab="LS",ylab="LS_est",main="LS vs LS_est")
  
  # LS from Ks (samples from TLP_Ks)
  LS_from_Ks_limitTLPKs <- sma_plot_stats(data.frame(trait$Ks,trait$LS),c("Ks","LS"),nbtstrp,F,indin = LS_from_TLP_Ks$dataused, regression_type =  regr_type)
  plot(trait$LS[LS_from_Ks$dataused],LS_from_Ks$var_est,pch=16,xlab="LS",ylab="LS_est",main="LS vs LS_est")
  
  
  # Summarise statistics
  if(regr_type == 'lm'|regr_type == 'pcr'|regr_type == 'plsr'){
    # Summarise statistics
    all_testnames_LS <-       c("LS_from_TLP_Ks_LMA_limitP50TLPKsMA","LS_from_P50_TLP_Ks_LMA"      ,"LS_from_P50_TLP_Ks","LS_from_LMA_TLP_Ks",       "LS_from_LMA_TLP"          ,"LS_from_TLP_Ks"    ,"LS_from_TLP_Ks_limitLMA_TLP_Ks" ,  "LS_from_P50_TLP","LS_from_P50_Ks","LS_from_P50","LS_from_TLP_limitLMA_TLP_Ks","LS_from_Ks",  "LS_from_LMA_limitLMA_TLP_Ks)", "LS_from_Ks_limitTLPKs")
    all_R2_LS       <- c(LS_from_TLP_Ks_LMA_limitP50TLPKsMA$R       , LS_from_P50_TLP_Ks_LMA$R     ,LS_from_P50_TLP_Ks$R,LS_from_LMA_TLP_Ks$R,         LS_from_LMA_TLP$R         ,LS_from_TLP_Ks$R   ,LS_from_TLP_Ks_limitLMA_TLP_Ks$R ,   LS_from_P50_TLP$R,LS_from_P50_Ks$R,LS_from_P50$R,LS_from_TLP$R,LS_from_Ks$R,  LS_from_LMA$R,LS_from_Ks_limitTLPKs$R)
    all_R2adj_LS <- c(LS_from_TLP_Ks_LMA_limitP50TLPKsMA$R2adj      , LS_from_P50_TLP_Ks_LMA$R2adj ,LS_from_P50_TLP_Ks$R2adj,LS_from_LMA_TLP_Ks$R2adj, LS_from_LMA_TLP$R2adj   , LS_from_TLP_Ks$R2adj,LS_from_TLP_Ks_limitLMA_TLP_Ks$R2,   LS_from_P50_TLP$R2adj,LS_from_P50_Ks$R2adj,LS_from_P50$R2adj,LS_from_TLP$R2adj,LS_from_Ks$R2adj,LS_from_LMA$R2adj, LS_from_Ks_limitTLPKs$R2adj)
    all_rmse_LS <- c(LS_from_TLP_Ks_LMA_limitP50TLPKsMA$rmse        ,LS_from_P50_TLP_Ks_LMA$rmse   ,LS_from_P50_TLP_Ks$rmse,LS_from_LMA_TLP_Ks$rmse   ,LS_from_LMA_TLP$rmse    ,LS_from_TLP_Ks$rmse  ,LS_from_TLP_Ks_limitLMA_TLP_Ks$rmse, LS_from_P50_TLP$rmse,LS_from_P50_Ks$rmse,LS_from_P50$rmse,LS_from_TLP$rmse,LS_from_Ks$rmse,LS_from_LMA$rmse,LS_from_Ks_limitTLPKs$rmse)
    all_ndata_LS <- c(LS_from_TLP_Ks_LMA_limitP50TLPKsMA$ndata,      LS_from_P50_TLP_Ks_LMA$ndata  ,LS_from_P50_TLP_Ks$ndata,LS_from_LMA_TLP_Ks$ndata ,LS_from_LMA_TLP$ndata   ,LS_from_TLP_Ks$ndata ,LS_from_TLP_Ks_limitLMA_TLP_Ks$ndata,LS_from_P50_TLP$ndata,LS_from_P50_Ks$ndata,LS_from_P50$ndata,LS_from_TLP$ndata,LS_from_Ks$ndata,LS_from_LMA$ndata,LS_from_Ks_limitTLPKs$ndata)
    
  }else{
    # Summarise statistics
    all_testnames_LS <-  c("LS_from_P50_TLP_Ks_LMA"      ,"LS_from_P50_TLP_Ks","LS_from_LMA_TLP_Ks",       "LS_from_LMA_TLP"          ,"LS_from_TLP_Ks"    ,"LS_from_TLP_Ks_limitLMA_TLP_Ks" ,  "LS_from_P50_TLP","LS_from_P50_Ks","LS_from_P50","LS_from_TLP_limitLMA_TLP_Ks","LS_from_Ks",  "LS_from_LMA_limitLMA_TLP_Ks)", "LS_from_Ks_limitTLPKs")
    all_R2_LS        <- c( LS_from_P50_TLP_Ks_LMA$R      ,LS_from_P50_TLP_Ks$R,LS_from_LMA_TLP_Ks$R,         LS_from_LMA_TLP$R         ,LS_from_TLP_Ks$R   ,LS_from_TLP_Ks_limitLMA_TLP_Ks$R ,   LS_from_P50_TLP$R,LS_from_P50_Ks$R,LS_from_P50$R,LS_from_TLP$R,LS_from_Ks$R,  LS_from_LMA$R,LS_from_Ks_limitTLPKs$R)
    all_R2adj_LS     <- c( LS_from_P50_TLP_Ks_LMA$R2adj  ,LS_from_P50_TLP_Ks$R2adj,LS_from_LMA_TLP_Ks$R2adj, LS_from_LMA_TLP$R2adj   , LS_from_TLP_Ks$R2adj,LS_from_TLP_Ks_limitLMA_TLP_Ks$R2,   LS_from_P50_TLP$R2adj,LS_from_P50_Ks$R2adj,LS_from_P50$R2adj,LS_from_TLP$R2adj,LS_from_Ks$R2adj,LS_from_LMA$R2adj, LS_from_Ks_limitTLPKs$R2adj)
    all_rmse_LS      <- c(LS_from_P50_TLP_Ks_LMA$rmse    ,LS_from_P50_TLP_Ks$rmse,LS_from_LMA_TLP_Ks$rmse   ,LS_from_LMA_TLP$rmse    ,LS_from_TLP_Ks$rmse  ,LS_from_TLP_Ks_limitLMA_TLP_Ks$rmse, LS_from_P50_TLP$rmse,LS_from_P50_Ks$rmse,LS_from_P50$rmse,LS_from_TLP$rmse,LS_from_Ks$rmse,LS_from_LMA$rmse,LS_from_Ks_limitTLPKs$rmse)
    all_ndata_LS     <- c( LS_from_P50_TLP_Ks_LMA$ndata  ,LS_from_P50_TLP_Ks$ndata,LS_from_LMA_TLP_Ks$ndata ,LS_from_LMA_TLP$ndata   ,LS_from_TLP_Ks$ndata ,LS_from_TLP_Ks_limitLMA_TLP_Ks$ndata,LS_from_P50_TLP$ndata,LS_from_P50_Ks$ndata,LS_from_P50$ndata,LS_from_TLP$ndata,LS_from_Ks$ndata,LS_from_LMA$ndata,LS_from_Ks_limitTLPKs$ndata)
    
  }
  
  all_LS <- data.frame(all_testnames_LS,all_R2_LS,all_R2adj_LS,all_rmse_LS,all_ndata_LS)
  if(view_stats==TRUE)View(all_LS)
  
  if(leaf_type== 'BDT'){
  # Test MAT and PPT coverage of species for chosen model
  plot(trait$MAT[LS_from_TLP_Ks$dataused],trait$MAP[LS_from_TLP_Ks$dataused],main = "Climate coverage")
  # WIDE CLIMATE COVERAGE
  
  # DECISION for BDT: :LS_from_TLP_Ks
    return_vals <- list("all_LS"=all_LS,
                        "LS_from_P50_TLP_Ks"= LS_from_P50_TLP_Ks)
  }
  
  if(leaf_type == 'BE'){
    # Test MAT and PPT coverage of species for chosen model
    plot(trait$MAT[LS_from_LMA_TLP_Ks$dataused],trait$MAP[LS_from_LMA_TLP_Ks$dataused],main = "Climate coverage")
    # WIDE CLIMATE COVERAGE
    
    # DECISION for BE: LS_from_LMA_TLP_Ks
    
    return_vals <- list("all_LS"=all_LS,
                        "LS_from_P50_TLP_Ks_LMA"= LS_from_P50_TLP_Ks_LMA)
  }
  
  return(return_vals)
}

# 
# ##LS_multivar_test_BDT <- function(trait_BDT,view_stats=FALSE, regr_type = 'lm') {
# ##  trait = trait_BDT
# # not really a test, as the trait network is already established, but testing the thin lines in slide 10, too:
# 
# 
# # LS from P50, TLP and Ks
# #  LS_from_P50_TLP_Ks <- sma_plot_stats(data.frame(trait$P50,trait$TLP,trait$Ks,trait$LS),c("P50","TLP","Ks","LS"),nbtstrp)
# #  plot(trait$LS[LS_from_P50_TLP_Ks$dataused],LS_from_P50_TLP_Ks$var_est,pch=16,xlab="LS",ylab="LS_est",main="LS vs LS_est")
# 
# # LS from TLP and Ks
# #  LS_from_TLP_Ks <- sma_plot_stats(data.frame(trait$TLP,trait$Ks,trait$LS),c("TLP","Ks","LS"),nbtstrp)
# plot(trait$LS[LS_from_TLP_Ks$dataused],LS_from_TLP_Ks$var_est,pch=16,xlab="LS",ylab="LS_est",main="LS vs LS_est")
# 
# # LS from P50, TLP
# LS_from_P50_TLP <- sma_plot_stats(data.frame(trait$P50,trait$TLP,trait$LS),c("P50","TLP","LS"),nbtstrp)
# plot(trait$LS[LS_from_P50_TLP$dataused],LS_from_P50_TLP$var_est,pch=16,xlab="LS",ylab="LS_est",main="LS vs LS_est")
# 
# # LS from P50, Ks
# LS_from_P50_Ks <- sma_plot_stats(data.frame(trait$P50,trait$Ks,trait$LS),c("P50","Ks","LS"),nbtstrp)
# plot(trait$LS[LS_from_P50_Ks$dataused],LS_from_P50_Ks$var_est,pch=16,xlab="LS",ylab="LS_est",main="LS vs LS_est")
# 
# # LS from P50
# LS_from_P50 <- sma_plot_stats(data.frame(trait$P50,trait$LS),c("P50","LS"),nbtstrp)
# plot(trait$LS[LS_from_P50$dataused],LS_from_P50$var_est,pch=16,xlab="LS",ylab="LS_est",main="LS vs LS_est")
# 
# # LS from TLP
# LS_from_TLP <- sma_plot_stats(data.frame(trait$TLP,trait$LS),c("TLP","LS"),nbtstrp)
# plot(trait$LS[LS_from_TLP$dataused],LS_from_TLP$var_est,pch=16,xlab="LS",ylab="LS_est",main="LS vs LS_est")
# 
# # LS from TLP and LMA
# LS_from_TLP_LMA <- sma_plot_stats(data.frame(trait$TLP,trait$LMA,trait$LS),c("TLP","LMA","LS"),nbtstrp)
# plot(trait$LS[LS_from_TLP_LMA$dataused],LS_from_TLP_LMA$var_est,pch=16,xlab="LS",ylab="LS_est",main="LS vs LS_est")
# 
# # LS from LMA
# LS_from_LMA<- sma_plot_stats(data.frame(trait$LMA,trait$LS),c("LMA","LS"),nbtstrp)
# plot(trait$LS[LS_from_LMA$dataused],LS_from_LMA$var_est,pch=16,xlab="LS",ylab="LS_est",main="LS vs LS_est")
# 
# # LS from Ks
# LS_from_Ks <- sma_plot_stats(data.frame(trait$Ks,trait$LS),c("Ks","LS"),nbtstrp)
# plot(trait$LS[LS_from_Ks$dataused],LS_from_Ks$var_est,pch=16,xlab="LS",ylab="LS_est",main="LS vs LS_est")
# 
# # LS from Ks (samples from TLP_Ks)
# LS_from_Ks_limitTLPKs <- sma_plot_stats(data.frame(trait$Ks,trait$LS),c("Ks","LS"),nbtstrp,F,LS_from_TLP_Ks$dataused)
# plot(trait$LS[LS_from_Ks$dataused],LS_from_Ks$var_est,pch=16,xlab="LS",ylab="LS_est",main="LS vs LS_est")
# 
# 
# 
# # Summarise statistics
# all_testnames_LS <- c("LS_from_TLP_Ks_LMA_limitP50TLPKsMA","LS_from_TLP_Ks_LMA","LS_from_TLP_LMA"     ,"LS_from_P50_TLP_Ks_LMA","LS_from_P50_TLP_Ks","LS_from_TLP_Ks","LS_from_P50_TLP","LS_from_P50_Ks","LS_from_P50","LS_from_TLP","LS_from_Ks","LS_from_Ks_limitTLPKs")
# all_R2_LS <- c( LS_from_TLP_Ks_LMA_limitP50TLPKsMA$R,LS_from_TLP_Ks_LMA$R,LS_from_TLP_LMA$R             ,LS_from_P50_TLP_Ks_LMA$R,LS_from_P50_TLP_Ks$R,LS_from_TLP_Ks$R,LS_from_P50_TLP$R,LS_from_P50_Ks$R,LS_from_P50$R,LS_from_TLP$R,LS_from_Ks$R,LS_from_Ks_limitTLPKs$R)
# all_R2adj_LS <- c( LS_from_TLP_Ks_LMA_limitP50TLPKsMA$R2adj,LS_from_TLP_Ks_LMA$R2adj,LS_from_TLP_LMA$R2adj  ,LS_from_P50_TLP_Ks_LMA$R2adj,LS_from_P50_TLP_Ks$R2adj,LS_from_TLP_Ks$R2adj,LS_from_P50_TLP$R2adj,LS_from_P50_Ks$R2adj,LS_from_P50$R2adj,LS_from_TLP$R2adj,LS_from_Ks$R2adj,LS_from_Ks_limitTLPKs$R2adj)
# all_rmse_LS <- c( LS_from_TLP_Ks_LMA_limitP50TLPKsMA$rmse, LS_from_TLP_Ks_LMA$rmse,LS_from_TLP_LMA$rmse    , LS_from_P50_TLP_Ks_LMA$rmse,LS_from_P50_TLP_Ks$rmse,LS_from_TLP_Ks$rmse,LS_from_P50_TLP$rmse,LS_from_P50_Ks$rmse,LS_from_P50$rmse,LS_from_TLP$rmse,LS_from_Ks$rmse,LS_from_Ks_limitTLPKs$rmse)
# all_ndata_LS <- c( LS_from_TLP_Ks_LMA_limitP50TLPKsMA$ndata,LS_from_TLP_Ks_LMA$ndata,LS_from_TLP_LMA$ndata  , LS_from_P50_TLP_Ks_LMA$ndata,LS_from_P50_TLP_Ks$ndata,LS_from_TLP_Ks$ndata,LS_from_P50_TLP$ndata,LS_from_P50_Ks$ndata,LS_from_P50$ndata,LS_from_TLP$ndata,LS_from_Ks$ndata,LS_from_Ks_limitTLPKs$ndata)
# 
# all_LS <- data.frame(all_testnames_LS,all_R2_LS,all_R2adj_LS,all_rmse_LS,all_ndata_LS)
# if(view_stats==TRUE)View(all_LS)
# 
# 
# # Test MAT and PPT coverage of species for chosen model
# plot(trait$MAT[LS_from_TLP_Ks_LMA$dataused],trait$MAP[LS_from_TLP_Ks_LMA$dataused],main = "Climate coverage")
# # WIDE CLIMATE COVERAGE
# 
# # DECISION: LS_from_TLP_Ks_LMA
# 
# return_vals <- list("all_LS"=all_LS,
#                     "LS_from_TLP_Ks_LMA"= LS_from_TLP_Ks_LMA)
# 
# return(return_vals)
# }