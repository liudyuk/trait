#-------------------------------------  
# Original code from T.Pugh
# Enabling of more flexible handling of inputs by Annemarie Eckes-Shephard May 2021
#trait: trait dataaframe from which to extract the observed trait values
#trait_names: list of traits (plain names, no _e suffix) for which to calculate summary statistics. Used to call the related objsects and perform stats on them.
#ind: indeces of the variables which were used to start the optimisation with. Used for subsetting trait dataframe
#output: dataframe with summary statistics on fit between predicted traits vs these traits' observations.
opt_rmse <- function(trait,
                     trait_names,
                     ind) {

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

    trait_sel       <- get(paste0(tn, '_e')) # call trait object with name tn
    ind_trait_sel   = which(!is.na(trait[[tn]][ind]) & !is.na(trait_sel[,1]))
    trait_sel_res   <- trait[[tn]][ind[ind_trait_sel]]-trait_sel[ind_trait_sel,1]
    trait_sel_rmse  <- sqrt(mean(  trait_sel_res^2,na.rm=T))  # standardised residuals
    
    trait_sel_ndata <- length(which(is.na(trait_sel_res)==F))
    trait_sel_R     <- cor(trait[[tn]][ind[ind_trait_sel]],trait_sel[ind_trait_sel,1])
    trait_sel_R2    <- trait_sel_R^2
    # account for slope and WD having more predictors:
    if(tn == 'slope' | tn == 'WD'){
      trait_sel_R2adj <- 1 - ( ((1-trait_sel_R2)*(trait_sel_ndata-1))/(trait_sel_ndata-nregress_other-1) )
    }else{
    trait_sel_R2adj <- 1 - ( ((1-trait_sel_R2)*(trait_sel_ndata-1))/(trait_sel_ndata-nregress_opt-1) )
    }
    plot(trait[[tn]][ind[ind_trait_sel]], trait_sel[ind_trait_sel,1],xlab=paste(tn),ylab=paste0(tn,'opt'))
    
     #collect output for each trait (tn) in dataframe
    df_sum_stats[which(df_sum_stats$all_names == tn),2:6] <- c(trait_sel_rmse,trait_sel_R,trait_sel_R2,trait_sel_R2adj, trait_sel_ndata)
         
  }
  
  # Set back to single plot
  par(mfrow=c(1,1))
  par(mar=c(5.1,4.1,4.1,2.1))
  
 
  return_vars <- df_sum_stats
  
  return(return_vars)
  
}
