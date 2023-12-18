test_cor_signs <- function(trait_obs,trait_est){
  # test whether correlation between estimated traits and observed traits has the correct direction (sign)
  #excluding P88, as the result of slope and P50  
  # inputs:
  # trait = observed traits
  # trait_est = estimated traits (e.g. right after optimisation or before PCA)
  # outputs:
  # prints a table with TRUE FALSE to see whether observed and estimated sign of the relationship is correct(TRUE) or not(FALSE)
  # returns a table containing the correlation coefficients and their significance for both estimated and observed values
  # this is only a rough estimate, as lots of values are omitted, when one wants complete cases to do a correlation matrix for the optimisation.
  # in the cases against which I checked the complete cases against bivariate correlation with complete cases for two traits, (bivar$all_sma_bivar ['parson_vor']), 
  # the correlation signs seem to correspond
 
  
  #observed subset for which PCA will be done
  subs_obs <- na.omit(trait_obs[,c("WD" ,     "P50",             "LS",         "Ks",        "TLP",        "LMA", "slope")])
  cor_obs10 <- cor(subs_obs)
  
  # estimated (e.g. right after optimisation or before PCA) subset of traits
  cor_est10 <- subs_est <- cor(trait_est[,c("WD" ,     "P50",             "LS",         "Ks",        "TLP",        "LMA","slope")])
  
  # prepare for return object: correlations and significance for both observed and estimated traits ( below)
  # this bit is just used to get the right names. 
  # These correlations would be wrong, as too many incomplete.cases were removed for some bivariate relationships.
  z = cor_obs10
  z[lower.tri( z ,diag=TRUE)]=NA  #Prepare to drop duplicates and meaningless information
  z=as.data.frame(as.table( z ))  #Turn into a 3-column table
  z_obs =na.omit( z )  #Get rid of the junk we flagged above
  names(z_obs)[3] <- 'obs_cor'
  
  lgth <- dim(z_obs[,1:2])[1]
  # now populate with actual correlations:
  for(i in 1:lgth){
    z_obs$obs_cor[i] <- cor(na.omit(trait_obs[,c(as.matrix(z_obs[i,1:2])[1:2])]))[1,2]
}
  
  # 
  z = cor_est10
  z[lower.tri( z ,diag=TRUE)]=NA  #Prepare to drop duplicates and meaningless information
  z=as.data.frame(as.table( z ))  #Turn into a 3-column table
  z_est =na.omit( z )  #Get rid of the junk we flagged above
  names(z_est)[3] <- 'est_cor'
  
  lgth <- dim(z_est[,1:2])[1]
  # now populate with actual correlations:
  for(i in 1:lgth){
    z_est$est_cor[i] <- cor(na.omit(trait_est[,c(as.matrix(z_est[i,1:2])[1:2])]))[1,2]
  }
  
  cor_obsest <-merge(z_obs,z_est,by=c("Var1","Var2"))
  cor_obsest$traits <-  paste(cor_obsest$Var1,cor_obsest$Var2,sep='-')
  cor_obsest$Var1 <- NULL
  cor_obsest$Var2 <- NULL
  
  return_obj <- cor_obsest[order(cor_obsest$obs_cor,decreasing=TRUE),]
  
  # test for similar relationship ( as in: sign of the relationship ):
  return_obj$obs_sign <-NA
  return_obj$est_sign <-NA
  return_obj$relationship <- NA
  
  return_obj[which(return_obj$obs_cor < 0),]$obs_sign <- 0
  return_obj[which(return_obj$obs_cor> 0),]$obs_sign <- 1
  
  return_obj[which(return_obj$est_cor > 0),]$est_sign <- 1
  return_obj[which(return_obj$est_cor < 0),]$est_sign <- 0
  
  return_obj$relationship  <- return_obj$est_sign==return_obj$obs_sign
  
  return(return_obj)
}
