unscale_traits <- function(varst, labels, nlabels, traits_mean, traits_sd){ 
  # for unscaling of single values during the network predictions, varst is a single value ( could be changed to a list, but done like this for now)
  #reverses the scaling of the traits in scale_traits.R after network prediction of the next value.
  
  if(labels == "TLP_e_last"  || labels == "TLP_e_start" || labels == "TLP_e" ){
    varst <- varst * traits_sd[,"TLP"]  
  }
  if(labels == "LS_e_last"   || labels == "LS_e_start" || labels == "LS_e" ){
    varst <- varst * traits_sd[,"LS"]  
  }
  if(labels == "LMA_e_last"   || labels == "LMA_e_start" || labels == "LMA_e" ){
    varst <- varst * traits_sd[,"LMA"] 
  }
  if(labels == "P50_e_last"   || labels == "P50_e_start" || labels == "P50_e" ){
    varst <- varst * traits_sd[,"P50"] 
  }
  if(labels == "Ks_e_last"   || labels == "Ks_e_start" || labels == "Ks_e" ){
    varst <- varst * traits_sd[,"Ks"]  
  }
  if(labels == "WD_e_last"   || labels == "WD_e_start" || labels == "WD_e" ){
    varst <- varst * traits_sd[,"WD"]
  }
  if(labels == "slope_e_last"   || labels == "slope_e_start" || labels == "slope_e" ){
    varst <- varst * traits_sd[,"slope"]  
  }
  
  return(varst)
}
