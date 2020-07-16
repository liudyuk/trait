
#libarary
library(ggplot2)
library(cowplot)

#import the file
trait<-read.table('woody_trait.0625.txt',header=TRUE, stringsAsFactors=FALSE, sep="\t", dec=".",fileEncoding="latin1")
#2_0) all speccies for the bivariate regression
xy<-trait
names(xy)
summary(model<-lm(P50~TLP,data=xy))
P50_TLP<-ggplot(xy, aes(x=TLP, y=P50)) +geom_point(color="red") +labs(title = "y=0.89x+0.32, r2=0.20,p<0.001")+geom_smooth(method=lm, se=FALSE)+
  scale_x_continuous(name = "TLP",) +
  scale_y_continuous(name = "P50")
summary(model<-lm(slope~TLP,data=xy))
slope_TLP<-ggplot(xy, aes(x=TLP, y=slope)) +geom_point(color="red") +labs(title = "y=-0.91x+4.1, r2=0.11,p<0.001")+geom_smooth(method=lm, se=FALSE)+
  scale_x_continuous(name = "TLP",) +
  scale_y_continuous(name = "slope")
summary(model<-lm(slope~P50,data=xy))
slope_P50<-ggplot(xy, aes(x=P50, y=slope)) +geom_point(color="red") +labs(title = "y=-0.45x+3.8, r2=0.17,p<0.001")+geom_smooth(method=lm, se=FALSE)+
  scale_x_continuous(name = "P50",) +
  scale_y_continuous(name = "slope")
summary(model<-lm(WD~TLP,data=xy))
WD_TLP<-ggplot(xy, aes(x=WD, y=TLP)) +geom_point(color="red") +labs(title = "y=0.15x+0.49, r2=0.09,p<0.001")+geom_smooth(method=lm, se=FALSE)+
  scale_x_continuous(name = "TLP",) +
  scale_y_continuous(name = "WD")

summary(model<-lm(WD~P50,data=xy))
WD_P50<-ggplot(xy, aes(x=P50, y=WD)) +geom_point(color="red") +labs(title = "y=0.03x+0.54, r2=0.03,p<0.001")+geom_smooth(method=lm, se=FALSE)+
  scale_x_continuous(name = "P50",) +
  scale_y_continuous(name = "WD")
summary(model<-lm(LMA~TLP,data=xy))
LMA_TLP<-ggplot(xy, aes(x=TLP, y=LMA)) +geom_point(color="red") +labs(title = "y=0.86x+3.9, r2=0.19,p<0.001")+geom_smooth(method=lm, se=FALSE)+
  scale_x_continuous(name = "TLP",) +
  scale_y_continuous(name = "LMA")
summary(model<-lm(LMA~LS,data=xy))
LMA_LS<-ggplot(xy, aes(x=LMA, y=LS)) +geom_point(color="red") +labs(title = "y=-0.32x+4.4,r2=0.19,p<0.001")+geom_smooth(method=lm, se=FALSE)+
  scale_x_continuous(name = "LS",) +
  scale_y_continuous(name = "LMA")
summary(model<-lm(Ks~LS_Hmax,data=xy))
Ks_LS_Hmax<-ggplot(xy, aes(x=LS_Hmax, y=Ks)) +geom_point(color="red") +labs(title = "y=0.3x-0.19,r2=0.24,p<0.001")+geom_smooth(method=lm, se=FALSE)+
  scale_x_continuous(name = "LS_Hmax",) +
  scale_y_continuous(name = "Ks")

summary(model<-lm(P50~Ks,data=xy))
P50_Ks<-ggplot(xy, aes(x=Ks, y=P50)) +geom_point(color="red") +labs(title = "y=-0.29x+0.97,r2=0.2,p<0.001")+geom_smooth(method=lm, se=FALSE)+
  scale_x_continuous(name = "Ks",) +
  scale_y_continuous(name = "P50")
summary(model<-lm(TLP~LS,data=xy))
TLP_LS<-ggplot(xy, aes(x=LS, y=TLP)) +geom_point(color="red") +labs(title = "y=-0.16x+0.51,r2=0.16,p<0.001")+geom_smooth(method=lm, se=FALSE)+
  scale_x_continuous(name = "LS",) +
  scale_y_continuous(name = "TLP")
summary(model<-lm(LMA~WD,data=xy))
LMA_WD<-ggplot(xy, aes(x=WD, y=LMA)) +geom_point(color="red") +labs(title = "y=0.86x+4,r2=0.05,p<0.001")+geom_smooth(method=lm, se=FALSE)+
  scale_x_continuous(name = "WD",) +
  scale_y_continuous(name = "LMA")
summary(model<-lm(slope~Ks,data=xy))
slope_Ks<-ggplot(xy, aes(x=Ks, y=slope)) +geom_point(color="red") +labs(title = "y=-0.03x+3.2,r2=0.001,ns")+geom_smooth(method=lm, se=FALSE)+
  scale_x_continuous(name = "Ks",) +
  scale_y_continuous(name = "slope")

summary(model<-lm(Hmax~Ks,data=xy))
Hmax_Ks<-ggplot(xy, aes(x=Ks, y=Hmax)) +geom_point(color="red") +labs(title = "y=2.1x+18.9,r2=0.02,p<0.001")+geom_smooth(method=lm, se=FALSE)+
  scale_x_continuous(name = "Ks",) +
  scale_y_continuous(name = "Hmax")
summary(model<-lm(leafN~LMA,data=xy))
leafN_LMA<-ggplot(xy, aes(x=LMA, y=leafN)) +geom_point(color="red") +labs(title = "y=-0.39x+4.7,r2=0.31,p<0.001")+geom_smooth(method=lm, se=FALSE)+
  scale_x_continuous(name = "LMA",) +
  scale_y_continuous(name = "leafNleafN")
leafN_LMA

fig3<-plot_grid(P50_TLP,slope_TLP,slope_P50, WD_TLP,
                WD_P50, LMA_TLP,LMA_LS,Ks_LS_Hmax,
                P50_Ks,TLP_LS,LMA_WD,slope_Ks,
                Hmax_Ks,leafN_LMA,
                nrow=4)
fig3
ggsave(fig3,filename = '/Users/liudy/Downloads/ALL.tiff',width = 16,height = 16)

#2_1)broadleaf
trait_g<-trait[trait$group=="BT"|trait$group=="BD"|trait$group=="BE",]
names(trait_g)
xy<-trait_g
names(xy)

summary(model<-lm(P50~TLP,data=xy))
P50_TLP<-ggplot(xy, aes(x=TLP, y=P50)) +geom_point(color="red") +labs(title = "y=0.81x+0.33, r2=0.18,p<0.001")+geom_smooth(method=lm, se=FALSE)+
  scale_x_continuous(name = "TLP",) +
  scale_y_continuous(name = "P50")
summary(model<-lm(slope~TLP,data=xy))
slope_TLP<-ggplot(xy, aes(x=TLP, y=slope)) +geom_point(color="red") +labs(title = "y=-0.93x+4.1, r2=0.12,p<0.001")+geom_smooth(method=lm, se=FALSE)+
  scale_x_continuous(name = "TLP",) +
  scale_y_continuous(name = "slope")
summary(model<-lm(slope~P50,data=xy))
slope_P50<-ggplot(xy, aes(x=P50, y=slope)) +geom_point(color="red") +labs(title = "y=-0.49x+3.8, r2=0.16,p<0.001")+geom_smooth(method=lm, se=FALSE)+
  scale_x_continuous(name = "P50",) +
  scale_y_continuous(name = "slope")
summary(model<-lm(WD~TLP,data=xy))
WD_TLP<-ggplot(xy, aes(x=WD, y=TLP)) +geom_point(color="red") +labs(title = "y=0.18x+0.48, r2=0.14,p<0.001")+geom_smooth(method=lm, se=FALSE)+
  scale_x_continuous(name = "TLP",) +
  scale_y_continuous(name = "WD")

summary(model<-lm(WD~P50,data=xy))
WD_P50<-ggplot(xy, aes(x=P50, y=WD)) +geom_point(color="red") +labs(title = "y=0.06x+0.54, r2=0.07,p<0.001")+geom_smooth(method=lm, se=FALSE)+
  scale_x_continuous(name = "P50",) +
  scale_y_continuous(name = "WD")
summary(model<-lm(LMA~TLP,data=xy))
LMA_TLP<-ggplot(xy, aes(x=TLP, y=LMA)) +geom_point(color="red") +labs(title = "y=0.67x+3.9, r2=0.16,p<0.001")+geom_smooth(method=lm, se=FALSE)+
  scale_x_continuous(name = "TLP",) +
  scale_y_continuous(name = "LMA")
summary(model<-lm(LMA~LS,data=xy))
LMA_LS<-ggplot(xy, aes(x=LMA, y=LS)) +geom_point(color="red") +labs(title = "y=-0.22x+4.4,r2=0.11,p<0.001")+geom_smooth(method=lm, se=FALSE)+
  scale_x_continuous(name = "LS",) +
  scale_y_continuous(name = "LMA")
summary(model<-lm(Ks~LS_Hmax,data=xy))
Ks_LS_Hmax<-ggplot(xy, aes(x=LS_Hmax, y=Ks)) +geom_point(color="red") +labs(title = "y=0.32x-0.15,r2=0.27,p<0.001")+geom_smooth(method=lm, se=FALSE)+
  scale_x_continuous(name = "LS_Hmax",) +
  scale_y_continuous(name = "Ks")

summary(model<-lm(P50~Ks,data=xy))
P50_Ks<-ggplot(xy, aes(x=Ks, y=P50)) +geom_point(color="red") +labs(title = "y=-0.22x+0.84,r2=0.13,p<0.001")+geom_smooth(method=lm, se=FALSE)+
  scale_x_continuous(name = "Ks",) +
  scale_y_continuous(name = "P50")
summary(model<-lm(TLP~LS,data=xy))
TLP_LS<-ggplot(xy, aes(x=LS, y=TLP)) +geom_point(color="red") +labs(title = "y=-0.15x+0.51,r2=0.13,p<0.001")+geom_smooth(method=lm, se=FALSE)+
  scale_x_continuous(name = "LS",) +
  scale_y_continuous(name = "TLP")
summary(model<-lm(LMA~WD,data=xy))
LMA_WD<-ggplot(xy, aes(x=WD, y=LMA)) +geom_point(color="red") +labs(title = "y=1.06x+3.8,r2=0.1,p<0.001")+geom_smooth(method=lm, se=FALSE)+
  scale_x_continuous(name = "WD",) +
  scale_y_continuous(name = "LMA")
summary(model<-lm(slope~Ks,data=xy))
slope_Ks<-ggplot(xy, aes(x=Ks, y=slope)) +geom_point(color="red") +labs(title = "y=-0.05x+3.2,r2=0.005,ns")+geom_smooth(method=lm, se=FALSE)+
  scale_x_continuous(name = "Ks",) +
  scale_y_continuous(name = "slope")

summary(model<-lm(Hmax~Ks,data=xy))
Hmax_Ks<-ggplot(xy, aes(x=Ks, y=Hmax)) +geom_point(color="red") +labs(title = "y=3.7x+16.1,r2=0.08,p<0.001")+geom_smooth(method=lm, se=FALSE)+
  scale_x_continuous(name = "Ks",) +
  scale_y_continuous(name = "Hmax")
summary(model<-lm(leafN~LMA,data=xy))
leafN_LMA<-ggplot(xy, aes(x=LMA, y=leafN)) +geom_point(color="red") +labs(title = "y=-0.38x+4.7,r2=0.28,p<0.001")+geom_smooth(method=lm, se=FALSE)+
  scale_x_continuous(name = "LMA",) +
  scale_y_continuous(name = "leafNleafN")
leafN_LMA

fig3<-plot_grid(P50_TLP,slope_TLP,slope_P50, WD_TLP,
                WD_P50, LMA_TLP,LMA_LS,Ks_LS_Hmax,
                P50_Ks,TLP_LS,LMA_WD,slope_Ks,
                Hmax_Ks,leafN_LMA,
                nrow=4)
fig3

ggsave(fig3,filename = '/Users/liudy/Downloads/broadleaf.tiff',width = 16,height = 16)
#######################################################################################################################################
#2_3) Multi-bivariiate Broadleaf
#LMA~TLP #only one variable and do not test
###################################################################
#run the function from Tom's developed multivariable
#TLP ~other traits
trait_g<-trait[trait$group!="CC",]
xy<-trait_g
names(xy)
xy1<-xy[,c("LMA","LS","WD","P50","TLP")]
xy2<-na.omit(xy1)
model<-sma_regress_multivar(xy2)
pred<-model$intercept_R+model$slope_R.y1*xy2$LMA+
  model$slope_R.y2*xy2$LS+model$slope_R.y3*xy2$WD+model$slope_R.y4*(xy2$P50)
res<-xy2$TLP-pred
rmse<-sqrt(mean(res^2))
rmse
plot(pred,xy2$TLP)
abline(0, 1)

xy2$pred<-pred
model1<-sma(pred~TLP,data=xy2,method=c("SMA"))
model1
summary(model<-lm(pred~TLP,data=xy2))
tlp<-ggplot(xy2, aes(x=TLP, y=pred)) +geom_point(color="red") +labs(title = "y=0.5x+0.27, R2=0.3438,R2_adj=0.3408,p<0.001")+geom_smooth(method=lm, se=FALSE)+
  scale_x_continuous(name = "TLP_act",) +
  scale_y_continuous(name = "TLP_pred")+
  geom_abline(slope=1)


#3 out-smaple test (choose 150)
trait_g<-trait[trait$group!="CC",]
xy<-trait_g
names(xy)
xy1<-xy[,c("LMA","LS","WD","P50","TLP")]
xy2<-na.omit(xy1)
xy2_out<-xy2[c(1:150),]
model<-sma_regress_multivar(xy2_out)
model
pred<-model$intercept_R+model$slope_R.y1*xy2_out$LMA+
  model$slope_R.y2*xy2_out$LS+model$slope_R.y3*xy2_out$WD+model$slope_R.y4*(xy2_out$P50)
res<-xy2_out$TLP-pred
rmse<-sqrt(mean(res^2))
rmse

xy2_out$pred<-pred
summary(model<-lm(pred~TLP,data=xy2_out))
tlp<-ggplot(xy2_out, aes(x=TLP, y=pred)) +geom_point(color="red") +labs(title = "y=0.47x+0.27, R2=0.3074,R2_adj=0.3028,p<0.001")+geom_smooth(method=lm, se=FALSE)+
  scale_x_continuous(name = "TLP_act",) +
  scale_y_continuous(name = "TLP_pred")+
  geom_abline(slope=1)

#rest sp 73
xy2_out2<-xy2[c(151:223),]
pred<-model$intercept_R+model$slope_R.y1*xy2_out2$LMA+
  model$slope_R.y2*xy2_out2$LS+model$slope_R.y3*xy2_out2$WD+model$slope_R.y4*(xy2_out2$P50)
res<-xy2_out2$TLP-pred
rmse<-sqrt(mean(res^2))
rmse

xy2_out2$pred<-pred
summary(model<-lm(pred~TLP,data=xy2_out2))
tlp<-ggplot(xy2_out2, aes(x=TLP, y=pred)) +geom_point(color="red") +labs(title = "y=0.53x+0.26, R2=0.3658,R2_adj=0.3569,p<0.001")+geom_smooth(method=lm, se=FALSE)+
  scale_x_continuous(name = "TLP_act",) +
  scale_y_continuous(name = "TLP_pred")+
  geom_abline(slope=1)
tlp



#5test for bivariate regressioon wiith best correlatiom
trait_g<-trait[trait$group!="CC",]
xy<-trait_g
names(xy)
xy1<-na.omit(xy[,c("P50","TLP")])
model1<-sma(TLP~P50,data=xy1,method=c("SMA"))
model1
res<-residuals(model1)
rmse<-sqrt(mean(res^2))
rmse

xy1<-na.omit(xy[,c("LMA","TLP")])
model1<-sma(TLP~LMA,data=xy1,method=c("SMA"))
model1
res<-residuals(model1)
rmse<-sqrt(mean(res^2))
rmse

#for 223 sp get the TLP prediction
trait_g<-trait[trait$group!="CC",]
xy<-trait_g
names(xy)
xy1<-xy[,c("LMA","LS","WD","P50","TLP")]
xy2<-na.omit(xy1)
model1<-sma(TLP~LMA,data=xy2,method=c("SMA"))
model1
res<-residuals(model1)
rmse<-sqrt(mean(res^2))
rmse

result<-merge(xy1_pre,xy2_pred,by="P50")
names(result)
summary(model<-lm(pred~pre,data=result))
plot_pre<-ggplot(result, aes(x=pre, y=pred)) +geom_point(color="red") +labs(title = "y=0.5x+0.27, R2=0.28,R2_adj=0.27,p<0.001")+geom_smooth(method=lm, se=FALSE)+
  scale_x_continuous(name = "pre_bi",) +
  scale_y_continuous(name = "pre_mbi")+
  geom_abline(slope=1)
plot_pre
#####test the other traits that are working


