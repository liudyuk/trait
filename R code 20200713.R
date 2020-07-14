<<<<<<< HEAD
#libarary
library(ggplot2)
library(cowplot)

#import the file
trait<-read.table('woody_trait.0625.txt',header=TRUE, stringsAsFactors=FALSE, sep="\t", dec=".",fileEncoding="latin1")
#all speccies for the bivariate regression
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
=======
#import the file
trait<-read.table('woody_trait.0625.txt',header=TRUE, 
                  stringsAsFactors=FALSE, sep="\t", dec=".",fileEncoding="latin1")
>>>>>>> d9fe87cf1c0f27e896110ece9d41dcb117465faf
