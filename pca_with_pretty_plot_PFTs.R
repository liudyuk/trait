pca_with_pretty_biplot_Pfts <- function(traits_PCA,center = T, scale = T){
  # performs pca and plots nice-lookin biplot with traits as loadings
  # trait scores are light grey dots
  # Annemarie Eckes-Shephard 12.10.2021
  # based on biplot function, but aesthetics are changed.
  # based on super helpful post by
  # https://stats.stackexchange.com/questions/276645/arrows-of-underlying-variables-in-pca-biplot-in-r
  library(RColorBrewer) # to plot the PFTs in different colours
  cols <- c(brewer.pal(9, "Set1") ,brewer.pal(8, "Set2"),brewer.pal(11, "Set3"),brewer.pal(3,'Paired'))
  CEN = scale(traits_PCA, center = center, scale = scale) # Centered and scaled
  PCA = prcomp(CEN)                                       # PCA analysis
  
  
  # Following getAnywhere(biplot.prcomp):
  
  choices = 1:2                                  # Selecting first two PC's
  scale = 1                                      # Default
  scores= PCA$x                                  # The scores
  lam = PCA$sdev[choices]                        # Sqrt e-vals (lambda) 2 PC's
  n = nrow(scores)                               # no. rows scores
  lam = lam * sqrt(n)                            # See below.
  s <- summary(PCA)                              # to plot PC axic labels with %ages (thanks https://www.benjaminbell.co.uk/2018/02/principal-components-analysis-pca-in-r.html )
  # deviation from above referenced post
  
  ############
  #
  # https://www.youtube.com/watch?v=EHb_kuw1GNU
  
  # make xtra sure that the placement of the arrows (eigenvectors) corresponds to the correlations between the traits
  # see 
  # https://stats.stackexchange.com/questions/104306/what-is-the-difference-between-loadings-and-correlation-loadings-in-pca-and
  #correlation loading plot
  #correlationloadings <- cor(traits_PCA, PCA$x)
  #plot(correlationloadings[,1], correlationloadings[,2],
  #     xlim=c(-1,1), ylim=c(-1,1),
  #     main='Correlation Loadings for PC1 vs. PC2')
  #############
  
  
  # at this point the following is called...
  # biplot.default(t(t(scores[,choices])      /  lam), 
  #                t(t(x$rotation[,choices]) *   lam))
  
  # Following from now on getAnywhere(biplot.default):
  
  x = t(t(scores[,choices])       / lam)         # scaled scores
  # "Scores that you get out of prcomp are scaled to have variance equal to      
  #  the eigenvalue. So dividing by the sq root of the eigenvalue (lam in 
  #  biplot) will scale them to unit variance. But if you want unit sum of 
  #  squares, instead of unit variance, you need to scale by sqrt(n)" (see comments).
  # > colSums(x^2)
  # PC1       PC2 
  # 0.9933333 0.9933333    # It turns out that the it's scaled to sqrt(n/(n-1)), 
  # ...rather than 1 (?) - 0.9933333=149/150
  
  y = t(t(PCA$rotation[,choices]) * lam)         # scaled eigenvecs (loadings)
  
  
  n = nrow(x)                                    # Same as dataset (150)
  p = nrow(y)                                    # Three var -> 3 rows
  
  # Names for the plotting:
  
  xlabs = 1L:n
  xlabs = as.character(xlabs)                    # no. from 1 to 150 
  dimnames(x) = list(xlabs, dimnames(x)[[2L]])   # no's and PC1 / PC2
  
  ylabs = dimnames(y)[[1L]]                      # trait names
  ylabs = as.character(ylabs)
  dimnames(y) <- list(ylabs, dimnames(y)[[2L]])  # Species and PC1/PC2
  
  # Function to get the range:
  unsigned.range = function(x) c(-abs(min(x, na.rm = TRUE)), 
                                 abs(max(x, na.rm = TRUE)))
  rangx1 = unsigned.range(x[, 1L])               # Range first col x
  # -0.1418269  0.1731236
  rangx2 = unsigned.range(x[, 2L])               # Range second col x
  # -0.2330564  0.2255037
  rangy1 = unsigned.range(y[, 1L])               # Range 1st scaled evec
  # -6.288626   11.986589
  rangy2 = unsigned.range(y[, 2L])               # Range 2nd scaled evec
  # -10.4776155   0.8761695
  
  (xlim = ylim = rangx1 = rangx2 = range(rangx1, rangx2))
  # range(rangx1, rangx2) = -0.2330564  0.2255037
  
  # And the critical value is the maximum of the ratios of ranges of 
  # scaled e-vectors / scaled scores:
  
  (ratio = max(rangy1/rangx1, rangy2/rangx2)) 
  # rangy1/rangx1   =   26.98328    53.15472
  # rangy2/rangx2   =   44.957418   3.885388
  # ratio           =   53.15472
  
  par(pty = "s")                                 # Calling a square plot
  
  # Plotting a box with x and y limits 
  # for the scaled scores:
  xlim =c(-0.3655831 , 0.3629699)
  ylim =c(-0.3655831 , 0.3629699)
  plot(x, type = "n", xlim = c(-0.3655831 , 0.3629699), ylim = c(-0.3655831 , 0.3629699),  # No points
       xlab=paste("PC 1 (", round(s$importance[2]*100, 1), "%)", sep = ""),
       ylab=paste("PC 2 (", round(s$importance[5]*100, 1), "%)", sep = ""))
  # Filling in the points as no's and the PC1 and PC2 labels:
  points(x,  pch=16 ,cex=2,col=cols) # exclude black
  text(x,labels=1:31)
  par(new = TRUE)                                # Avoids plotting what follows separately
  
  
  #par(new = TRUE)                                # Avoids plotting what follows separately
  
  # Setting now x and y limits for the arrows:
  
  (xlim = xlim * ratio)  # We multiply the original limits x ratio
  # -16.13617  15.61324
  (ylim = ylim * ratio)  # ... for both the x and y axis
  # -16.13617  15.61324
  
  # The following doesn't change the plot intially...
  plot(y, axes = FALSE, type = "n", 
       xlim = xlim, 
       ylim = ylim, xlab = "", ylab = "")
  abline(v=0, lty=2, col="grey50")
  abline(h=0, lty=2, col="grey50")
  # ... but it does now by plotting the ticks and new limits...
  # ... along the top margin (3) and the right margin (4)
  #axis(3); axis(4)
  text(y* 0.8, labels = ylabs, col = 2)  # This just prints the species
  
  arrow.len = 0.1                   # Length of the arrows about to plot.
  
  # The scaled e-vecs are further reduced to 60% of their value
  arrows(0, 0, y[, 1L] * 0.6, y[, 2L] * 0.6, 
         length = arrow.len, col = 2)
  
}