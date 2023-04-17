scale_traits <- function(varst, labels, nlabels, traits_mean, traits_sd){ 
  #to scale traits so that the correct value is used in the network predictions, when using plsr-model coefficients.
  # input:
  # varst dataframe of which the first nlables-1 variables will be scaled
  # labels =names of all traits in the dataframe
  # nlabels number of traits in the dataframe
  # traits_mean =  mean trait value from all traits in the dataset. will be subsetted for the correct one using 'labels'.
  # straits_sd = standard deviation from all traits in the dataset. will be subsetted for the correct one using 'labels'.
  # scaling is only performed on the explanatory variables X, because scale=TRUE only scales X variables. Therefore I do nlables-1
  # as index nlabels contains Y.
  
  # determine context in which the function is used ( e.g. in generating the model or predictions)
  # for scaling in generating the models, "varst" is multidimensional:
  if(!is.null(dim(varst))){ 
    ns <- names(varst)[1:(nlabels-1)]
    #print(paste('scaling',ns))
    for (n in 1:length(ns)){ # normalise predictor traits only
      varst[,ns[n]] <- (varst[,ns[n]] ) / traits_sd[,labels[n]] 
    }
  }else{  # for scaling of single values during the network predictions, varst is a single value ( could be changed to a list, but done like this for now)
    if(labels == "TLP_e_last"  || labels == "TLP_e_start" || labels == "TLP_e" ){
      varst <- (varst  ) / traits_sd[,"TLP"] 
    }
    if(labels == "LS_e_last"   || labels == "LS_e_start" || labels == "LS_e" ){
      varst <- (varst  ) / traits_sd[,"LS"] 
    }
    if(labels == "LMA_e_last"   || labels == "LMA_e_start" || labels == "LMA_e" ){
      varst <- (varst  ) / traits_sd[,"LMA"] 
    }
    if(labels == "P50_e_last"   || labels == "P50_e_start" || labels == "P50_e" ){
      varst <- (varst  ) / traits_sd[,"P50"] 
    }
    if(labels == "Ks_e_last"   || labels == "Ks_e_start" || labels == "Ks_e" ){
      varst <- (varst  ) / traits_sd[,"Ks"] 
    }
    if(labels == "WD_e_last"   || labels == "WD_e_start" || labels == "WD_e" ){
      varst <- (varst ) / traits_sd[,"WD"] 
    }
    if(labels == "slope_e_last"   || labels == "slope_e_start" || labels == "slope_e" ){
      varst <- (varst  ) / traits_sd[,"slope"] 
    }
  }
  return(varst)
}