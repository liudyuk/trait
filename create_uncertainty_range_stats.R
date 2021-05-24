#-------------------------------------  
# Original code from T.Pugh main script
# Enabling of more flexible handling of input and automatised file naming by Annemarie Eckes-Shephard May 2021

# Stats defining the uncertainty range for each point
# input is the output from trait_optim. 
# The function will use the optimised (predicted) values and creates summary statistics (sum_stat) on their
# mean, median and upper and lower 5% of the uncertainty range.
# the output ( <variable>_e_<sum_stat>) is made a available to the global environment.
create_uncertainty_range_stats <- function(ins){
  
  predicted_names <- names(ins$predicted)
  extended_names <- c('_mean','_median','_5perc','_95perc')
  for (pn in predicted_names){
    collect <- list()
    collect[[1]] <- unname(apply(ins$predicted[[pn]], 1, mean,na.rm=T))
    collect[[2]] <- unname(apply(ins$predicted[[pn]], 1, median,na.rm=T))
    collect[[3]] <- unname(apply(ins$predicted[[pn]], 1, quantile,0.05,na.rm=T))
    collect[[4]] <- unname(apply(ins$predicted[[pn]], 1, quantile,0.95,na.rm=T))
    names(collect) <- paste0(pn,extended_names )
    list2env(collect , envir = .GlobalEnv)
  }
  
}
