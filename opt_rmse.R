opt_rmse.R#-------------------------------------  
# Original code from T.Pugh
# Enabling of more flexible handling of inputs by Annemarie Eckes-Shephard May 2021
#trait: trait dataaframe from which to extract the observed trait values
#trait_names: list of traits (plain names, no _e suffix) for which to calculate summary statistics. Used to call the related objsects and perform stats on them.
#ind: indeces of the variables which were used to start the optimisation with. Used for subsetting trait dataframe
#spec_group_sel: ultimately subset fit calculations for the species groups for which the models and there fore predictions were made.
#output: dataframe with summary statistics on fit between predicted traits vs these traits' observations.
opt_rmse <- function(trait,
                     trait_names,
                     ind,
                     spec_group_sel) {

  nregress_opt=4 #Assume that there are 4 predictors for each of the 5 optimised traits (i.e. 5-1)
  nregress_other=5 #For the other traits (slope and WD) assume that all the optimised traits are predictors
  
  
  #prepare df to collect summary stats:
  df_sum_stats <- data.frame(all_names=trait_names,
                             all_e_rmse=rep(NA,length(trait_names)),
                             all_e_R =rep(NA,length(trait_names)),
                             all_e_R2 =rep(NA,length(trait_names)),
                             all_e_R2adj=rep(NA,length(trait_names)),
                             all_e_ndata=rep(NA,length(trait_names)))
  
  par(mfrow=c(2,3))
  par(mar=c(5,4,4,2))
  
  for (tn in trait_names){

    trait_sel       <- get(paste0(tn, '_e')) # call trait object with name tn. 
    # but for fit calculations we can only use the predicted traits, which can be paired with observations.
    # therefore, subset for complete cases of trait tn, from where optimisations with the complete starting pair (LS P50) have started (ind)
    trait_subs      <- trait[ind,] #  subset for complete predicted cases from where optimisation started from
    trait_subs$trait_sel <- trait_sel  # fill with predicted cases

    ind_trait_sel   =  which(!is.na(trait_subs[[tn]]) ) # subset for those cases which have observations only
    trait_subs      <-  trait_subs[ind_trait_sel,] 
    
    # calculate the residual and R of all available observed values of trait tn with the corresponding predicted trait value
    trait_subs$res <- trait_subs[[tn]] - trait_subs$trait_sel
    #trait_subs$R   <- cor(trait_subs[[tn]],trait_sel[ind_trait_sel,1])
    
    # but now make sure that we also take into account only the complete cases in the different species groups (BE and BDT), 
    # otherwise the comparisons against predictions would be incorrect ( e.g. (BE) model-predicted values would be compared against BDT observations ) We can only do this now, in order to avoid index confusions, when calculating the residuals above
    # so now filter out the residuals which actually matter for BDT or BE etc.:
    
    # Get index for selected species group
    if (spec_group_sel==1) {
      ind_spec_group = which( trait_subs$group=='BT' | trait_subs$group=='BD' )
      trait_sel_res  <- trait_subs[ind_spec_group,]$res
      trait_sel_R    <- cor(trait_subs[[tn]][ind_spec_group],trait_subs[ind_spec_group,]$trait_sel)
    } else if (spec_group_sel==2) {
      ind_spec_group = which(trait_subs$group=='BE')
      trait_sel_res  <- trait_subs[ind_spec_group,]$res
      trait_sel_R    <- cor(trait_subs[[tn]][ind_spec_group],trait_subs[ind_spec_group,]$trait_sel)
    } else if (spec_group_sel==3) {
      ind_spec_group = which(trait_subs$group=='BT')
      trait_sel_res  <- trait_subs[ind_spec_group,]$res
      trait_sel_R    <- cor(trait_subs[[tn]][ind_spec_group],trait_subs[ind_spec_group,]$trait_sel)
    } else if (spec_group_sel==4) {
      ind_spec_group = which(trait_subs$group=='BD')
      trait_sel_res  <- trait_subs[ind_spec_group,]$res
      trait_sel_R    <- cor(trait_subs[[tn]][ind_spec_group],trait_subs[ind_spec_group,]$trait_sel)
    }
    
    
    trait_sel_rmse  <- sqrt(mean( trait_sel_res^2, na.rm=T))  # standardised residuals
    trait_sel_ndata <- length(which(is.na(trait_sel_res)==F))
    trait_sel_R2    <- trait_sel_R^2
    # account for slope and WD having more predictors:
    if(tn == 'slope' | tn == 'WD'){
      trait_sel_R2adj <- 1 - ( ((1-trait_sel_R2)*(trait_sel_ndata-1))/(trait_sel_ndata-nregress_other-1) )
    }else{
    trait_sel_R2adj <- 1 - ( ((1-trait_sel_R2)*(trait_sel_ndata-1))/(trait_sel_ndata-nregress_opt-1) )
    }
    #plot(trait[[tn]][ind[ind_trait_sel]], trait_sel[ind_trait_sel,1],xlab=paste(tn),ylab=paste0(tn,'opt'))
    
     #collect output for each trait (tn) in dataframe
    df_sum_stats[which(df_sum_stats$all_names == tn),2:6] <- c(trait_sel_rmse,trait_sel_R,trait_sel_R2,trait_sel_R2adj, trait_sel_ndata)
         
  }
  
  # Set back to single plot
  par(mfrow=c(1,1))
  par(mar=c(5.1,4.1,4.1,2.1))
  
 
  return_vars <- df_sum_stats
  
  return(return_vars)
  
}
