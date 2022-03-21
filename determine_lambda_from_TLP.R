determine_lambda_from_TLP <- function(TLP_unlogged, plot=FALSE){
  #input: TLP_unlogged  network-predicted TLP parameters
  #output: lambdas åredicted lambdas, generated in this function
  
  ## create linear relationship between lambda and TLP with predicted values only, 
  ##to make sure they spread across the whole range of (model-) possible lambda parameters
  # * knowing that there is a relationship between TLP and lambda
  # * using upper and lower ranges of isohydricity as permitted by the updated LPJ-GUESS-HYD model
  # see papastefanou et al 2020 for upper and lower ranges for Lambda:
  
  #upper and lower ranges for network-predicted parameters for TLP (y) and lambda (x)
  # to create slope using max and min values - upon which to project all TLPs, and then derive our lambdas
  
  # negative relationship between TLP and Lambda, therefore match max TLP with minlambda
  y1 =  1.0  # extreme isohydric
  y2 =  -0.3 # anisohydric
  x1 =  min(TLP_unlogged)
  x2 =  max(TLP_unlogged)
  
  m = (y2-y1) / (x2-x1) 

  b = y2 - m * x2
  
  if(plot==TRUE){
    points( -exp(TLP_e_mean),(m * -exp(TLP_e_mean) + b),col='purple')
  }
  
  lambdas = b + m * -exp(TLP_e_mean) 
  return(lambdas)
}



