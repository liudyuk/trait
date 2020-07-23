bivariate_sma_plot_stats <- function(var1,var2,label1,label2,nbtstrp) {

ind=which(!is.na(var1) & !is.na(var2))
ndata=length(ind)

mod <- sma_regress_multivar(data.frame(var1[ind],var2[ind]),nbtstrp,T)

var2_est <- mod$intercept_R + mod$slope_R.y1*var1[ind]

#Make plot
plot(var1[ind],var2[ind],pch=16,xlab=label1,ylab=label2,main=paste(label1," vs ",label2))
points(var1[ind],var2_est,col="red",pch=16)

#Calculate RMSE
res <- var2[ind]-var2_est
rmse <- sqrt(mean(res^2))

#Calculate R2
R <- cor(var2[ind],var2_est)
R2 <- R^2

return_vals <- list("mod"=mod,"rmse"=rmse,"R"=R,"R2"=R2)

return(return_vals)
}