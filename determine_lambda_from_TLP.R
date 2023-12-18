determine_lambda_from_TLP <- function(TLP_unlogged, plot=FALSE,single=FALSE){
  #input: TLP_unlogged  list; network-predicted TLP parameters
  #input single, default = FALSE, treating TLP_unlogged like a list, and determining min max lambda from there. However,
  # if TLP_unlogged not treated as list, lambda will be derived as was done before, using an empirical relationship between it and TLP_unlogged 
  #output: lambdas predicted lambdas, generated in this function
  
  ## create linear relationship between lambda and TLP with predicted values only, ranging between max and min value of lamba, 
  ##to make sure they spread across the whole range of (model-) possible lambda parameters
  # * knowing that there is a relationship between TLP and lambda
  # * using upper and lower ranges of isohydricity as permitted by the updated LPJ-GUESS-HYD model
  # see papastefanou et al 2020 for upper and lower ranges for Lambda:
  
  #upper and lower ranges for network-predicted parameters for TLP (y) and lambda (x)
  # to create slope using max and min values - upon which to project all TLPs, and then derive our lambdas
  
  # positive relationship between TLP and Lambda, therefore match min TLP with minlambda
  #(see Fu et al=:
  # large hydroscapes are associated with more anisohydric behaviour, and figure 3 a) shows that anisohydric behaviour is associated with low TLP)
  # https://academic-oup-com.ludwig.lub.lu.se/treephys/article/39/1/122/5107063
  y2 =  1.0  # extreme isohydric # 𝛌max
  y1 =  -0.3 # anisohydric   #𝛌min
  x1 =  min(TLP_unlogged)
  x2 =  max(TLP_unlogged)
  
  m = (y2-y1) / (x2-x1) 

  b = y2 - m * x2 
  
  lambdas = b + m * (TLP_unlogged) 

 if(single==TRUE){
   print("single set to true, using empirical relationship")
     lambdas =  -0.188+(-0.3*TLP_unlogged) # Based on https://docs.google.com/spreadsheets/d/1KvB97vt6OVdy3GPUdOdlxd8c9TdjEk9z4qIl2mcfPHk/edit#gid=51069844
  }
  return(lambdas)
}



