#Probably a useful function to make nicer later, for the publication. For now just a plot to visualise the distribution of the systematic sample across the hypervolume space
# hv hypervolume object
# var1 first subset of the systematic sample of the trait variable used in hypervolume e.g. Ks_e
# var2 second ubset of the systematic sample of the trait variable used in  hypervolume, e.g. LS_e
plot_Hypervolume <- function(hv,var1,var2){
  plot(hv@RandomPoints, col= 'light grey')
  points(hv@Data,pch=16, col= 'black',cex=0.5)
  points(var1,var2,col='red',pch=16)
}
