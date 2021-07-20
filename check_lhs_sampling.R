# Function to check the distribution of sample-pairs from the latin hypercube sampling.
# Annemarie Eckes-Shephard, July 2021
# inputs:
# observed trait ( no NAs) 1
# observed trati ( no NAs) 2
# (latin hypercube) estimated trait 1
# (latin hypercube) estimated trait 2
# names of trait 1 and 2 as string for plotting
# returns:
# ggplot2 object that shows where the sampled values lie across the two-parameter space
check_lhs_sampling <- function(x,y,x_e,y_e,x_string,y_string){
  
  # most of this plootting solution is taken from:
  # https://stackoverflow.com/questions/8545035/scatterplot-with-marginal-histograms-in-ggplot2
  # Answer by user 'Ben'
  xy<-data.frame(x,y)
  # make the basic plot object
 p <- ggplot(xy, aes(x, y)) +        
    # set the locations of the x-axis labels as Tukey's five numbers   
    scale_x_continuous(limit=c(min(x), max(x)), 
                       breaks=round(fivenum(x),1)) +     
    # ditto for y-axis labels 
    scale_y_continuous(limit=c(min(y), max(y)),
                       breaks=round(fivenum(y),1)) +     
    # specify points
    geom_point() +
    # specify that we want the rug plot
    geom_rug(size=0.1) +   
    # improve the data/ink ratio
    theme_set(theme_minimal(base_size = 18)) +
    # add the hypervolume samples
    geom_point(data=data.frame(x_e,y_e),aes(x=x_e,y=y_e,colour="hq sample"), size=1) +
    
    theme_bw() +
    labs(x=paste0("Broadleaf ",x_string), y=paste0("Broadleaf ",y_string)) +
    scale_colour_manual(name='', values=c('data'='black', 'hq sample'='red'), guide='legend') 
  
  return(p)
}
