
## obtain approximate empirical distributions for Ks,LS,P50 and TLP for subsequent
# transforming of the hypercube stratified sample values. 
# This transformation is necessary to reflect the distributions of the observations.
# releases the objects fitLS, fitKs, fit 
estimate_trait_distr <- function(traits){


# not normally distributed
#shapiro.test(traits$LS[!is.na(traits$LS)]) 
fitLS <<- fit.st(traits$LS[!is.na(traits$LS)])
# Histogram, Kernel Density Plot and 'fitted' distribution.
hist(na.omit(traits$LS),prob = TRUE, main='Empirical and estimated density for broadleaf LS ',ylab='LS')
lines(density(na.omit(traits$LS)),main='empirical density for broadleaf LS observations')
x <- seq(from=-5, to=5, by=0.1)
lines(x,dst(x, mu=fitLS$par.ests[2]-1.2, sigma=fitLS$par.ests[3], nu=fitLS$par.ests[1], log=FALSE),col=3)
legend(legend=c('data','student t model'),fill=c(1,3),'topright')  

###not normally distributed
#shapiro.test(traits$Ks[!is.na(traits$Ks)])
fitKs <<- fit.st(na.omit(traits$Ks))
hist(na.omit(traits$Ks),prob = TRUE, main='Empirical and estimated density for broadleaf Ks ',ylab='Ks')
lines(density(na.omit(traits$Ks)),main='density distribution for broadleaf Ks observations')
x <- seq(from=-5, to=5, by=0.1)
lines(x,dst(x, mu=fitKs$par.ests[2], sigma=fitKs$par.ests[3], nu=fitKs$par.ests[1], log=FALSE),col=3)
legend(legend=c('data','student t model'),fill=c('grey',3),'topright')  


# not normally distributed: # If p> 0.05, normality can be assumed. Here p << 0.05
# shapiro.test(traits$P50[!is.na(traits$P50)]) 

fitP50 <<- fit.st(na.omit(traits$P50))
fitP50norm <- fitdistrplus::fitdist(traits$P50[!is.na(traits$P50)], "norm")
hist(na.omit(traits$P50),prob = TRUE, main='Empirical and estimated density of broadleaf P50 ',ylab='P50')
lines(density(na.omit(traits$P50)),main='density distribution for broadleaf Ks observations')
lines(x,dst(x,mu=fitP50$par.ests[2], sigma=fitP50$par.ests[3], nu=fitP50$par.ests[1], log=FALSE),col=3)
curve(dnorm(x, fitP50norm$estimate[1], fitP50norm$estimate[2]), col = 2, add = TRUE)
legend(legend=c('data','student t model','normal'),fill=c('grey',3,2),'topright')  


## _just_ unlikely to be normally distributed, according to this test..
#shapiro.test(traits$TLP[!is.na(traits$TLP)])
fitTLPnorm <- fitdistrplus::fitdist(traits$TLP[!is.na(traits$TLP)], "norm")
fitTLP <<- fit.st(na.omit(traits$TLP))
# check the fit is good:
# plot(fitTLP)
# good enough, proceeding with normal distribution of TLP
hist(na.omit(traits$TLP),prob = TRUE,main='density distribution for broadleaf Ks observations',ylab='TLP')
lines(density(na.omit(traits$TLP)))
lines(x,dst(x,mu=fitTLP$par.ests[2], sigma=fitTLP$par.ests[3], nu=fitTLP$par.ests[1], log=FALSE),col=3)
curve(dnorm(x, fitTLPnorm$estimate[1],fitTLPnorm$estimate[2]), col = 2, add = TRUE)
legend(legend=c('data','student t model','normal'),fill=c('grey',3,2),'topright')  

}
