# Make some plots that compare bivariate SMA regressions for broadleaf evergreen and broadleaf deciduous.
#T. Pugh
#27.08.20

source("trait_functions.R")

traits=read.csv("/Users/pughtam/Documents/TreeMort/Analyses/Hydraulic_modelling/Traits/woody_trait.0827.txt",sep="\t")

trait_allB<-subset(traits,group!="CC",drop = T)
trait_BE<-subset(traits,group=="BE",drop = T)
trait_BDT<-subset(traits,group=="BD" | group=="BT",drop = T)
trait_CC<-subset(traits,group=="CC",drop = T)

#source('trait_functions.R')
source('sma_multivar_regress.R')

nbtstrp=1000

#leafL_from_LMA <- sma_plot_stats_comp(data.frame(LMA,log(leafL)),c("LMA","leafL"),nbtstrp,T)

par(mfrow=c(4,4))
par(mar=c(2,2,2,2))

TLP_from_P50 <- sma_plot_stats_comp(data.frame(trait_BE$P50,trait_BE$TLP),data.frame(trait_BDT$P50,trait_BDT$TLP),c("P50","TLP"),nbtstrp,T)

TLP_from_slope <- sma_plot_stats_comp(data.frame(trait_BE$slope,trait_BE$TLP),data.frame(trait_BDT$slope,trait_BDT$TLP),c("slope","TLP"),nbtstrp,T)

P50_from_slope <- sma_plot_stats_comp(data.frame(trait_BE$slope,trait_BE$P50),data.frame(trait_BDT$slope,trait_BDT$P50),c("slope","P50"),nbtstrp,T)

TLP_from_WD <- sma_plot_stats_comp(data.frame(trait_BE$WD,trait_BE$TLP),data.frame(trait_BDT$WD,trait_BDT$TLP),c("WD","TLP"),nbtstrp,T)

P50_from_WD <- sma_plot_stats_comp(data.frame(trait_BE$WD,trait_BE$P50),data.frame(trait_BDT$WD,trait_BDT$P50),c("WD","P50"),nbtstrp,T)

TLP_from_LMA <- sma_plot_stats_comp(data.frame(trait_BE$LMA,trait_BE$TLP),data.frame(trait_BDT$LMA,trait_BDT$TLP),c("LMA","TLP"),nbtstrp,T)

LS_from_LMA <- sma_plot_stats_comp(data.frame(trait_BE$LMA,trait_BE$LS),data.frame(trait_BDT$LMA,trait_BDT$LS),c("LMA","LS"),nbtstrp,T)

Ks_from_LSHmax <- sma_plot_stats_comp(data.frame(trait_BE$LS_Hmax,trait_BE$Ks),data.frame(trait_BDT$LS_Hmax,trait_BDT$Ks),c("LS*Hmax","Ks"),nbtstrp,T)

Ks_from_P50 <- sma_plot_stats_comp(data.frame(trait_BE$P50,trait_BE$Ks),data.frame(trait_BDT$P50,trait_BDT$Ks),c("P50","Ks"),nbtstrp,T)

LS_from_TLP <- sma_plot_stats_comp(data.frame(trait_BE$TLP,trait_BE$LS),data.frame(trait_BDT$TLP,trait_BDT$LS),c("TLP","LS"),nbtstrp,T)

WD_from_LMA <- sma_plot_stats_comp(data.frame(trait_BE$LMA,trait_BE$WD),data.frame(trait_BDT$LMA,trait_BDT$WD),c("LMA","WD"),nbtstrp,T)

Ks_from_slope <- sma_plot_stats_comp(data.frame(trait_BE$slope,trait_BE$Ks),data.frame(trait_BDT$slope,trait_BDT$Ks),c("slope","Ks"),nbtstrp,T)

leafN_from_LMA <- sma_plot_stats_comp(data.frame(trait_BE$LMA,trait_BE$leafN),data.frame(trait_BDT$LMA,trait_BDT$leafN),c("LMA","leafN"),nbtstrp,T)

leafL_from_LMA <- sma_plot_stats_comp(data.frame(trait_BE$LMA,log(trait_BE$leafL)),data.frame(trait_BDT$LMA,log(trait_BDT$leafL)),c("LMA","leafL"),nbtstrp,T)

Ks_from_LMA <- sma_plot_stats_comp(data.frame(trait_BE$LMA,trait_BE$Ks),data.frame(trait_BDT$LMA,trait_BDT$Ks),c("LMA","Ks"),nbtstrp,T)
# Make some plots that compare bivariate SMA regressions for broadleaf evergreen and broadleaf deciduous.
#T. Pugh
#27.08.20

source("trait_functions.R")

traits=read.csv("/Users/pughtam/Documents/TreeMort/Analyses/Hydraulic_modelling/Traits/woody_trait.0827.txt",sep="\t")

trait_allB<-subset(traits,group!="CC",drop = T)
trait_BE<-subset(traits,group=="BE",drop = T)
trait_BDT<-subset(traits,group=="BD" | group=="BT",drop = T)
trait_CC<-subset(traits,group=="CC",drop = T)

#source('trait_functions.R')
source('sma_multivar_regress.R')

nbtstrp=1000

#leafL_from_LMA <- sma_plot_stats_comp(data.frame(LMA,log(leafL)),c("LMA","leafL"),nbtstrp,T)

par(mfrow=c(4,4))
par(mar=c(2,2,2,2))

TLP_from_P50 <- sma_plot_stats_comp(data.frame(trait_BE$P50,trait_BE$TLP),data.frame(trait_BDT$P50,trait_BDT$TLP),c("P50","TLP"),nbtstrp,T)

TLP_from_slope <- sma_plot_stats_comp(data.frame(trait_BE$slope,trait_BE$TLP),data.frame(trait_BDT$slope,trait_BDT$TLP),c("slope","TLP"),nbtstrp,T)

P50_from_slope <- sma_plot_stats_comp(data.frame(trait_BE$slope,trait_BE$P50),data.frame(trait_BDT$slope,trait_BDT$P50),c("slope","P50"),nbtstrp,T)

TLP_from_WD <- sma_plot_stats_comp(data.frame(trait_BE$WD,trait_BE$TLP),data.frame(trait_BDT$WD,trait_BDT$TLP),c("WD","TLP"),nbtstrp,T)

P50_from_WD <- sma_plot_stats_comp(data.frame(trait_BE$WD,trait_BE$P50),data.frame(trait_BDT$WD,trait_BDT$P50),c("WD","P50"),nbtstrp,T)

TLP_from_LMA <- sma_plot_stats_comp(data.frame(trait_BE$LMA,trait_BE$TLP),data.frame(trait_BDT$LMA,trait_BDT$TLP),c("LMA","TLP"),nbtstrp,T)

LS_from_LMA <- sma_plot_stats_comp(data.frame(trait_BE$LMA,trait_BE$LS),data.frame(trait_BDT$LMA,trait_BDT$LS),c("LMA","LS"),nbtstrp,T)

Ks_from_LSHmax <- sma_plot_stats_comp(data.frame(trait_BE$LS_Hmax,trait_BE$Ks),data.frame(trait_BDT$LS_Hmax,trait_BDT$Ks),c("LS*Hmax","Ks"),nbtstrp,T)

Ks_from_P50 <- sma_plot_stats_comp(data.frame(trait_BE$P50,trait_BE$Ks),data.frame(trait_BDT$P50,trait_BDT$Ks),c("P50","Ks"),nbtstrp,T)

LS_from_TLP <- sma_plot_stats_comp(data.frame(trait_BE$TLP,trait_BE$LS),data.frame(trait_BDT$TLP,trait_BDT$LS),c("TLP","LS"),nbtstrp,T)

WD_from_LMA <- sma_plot_stats_comp(data.frame(trait_BE$LMA,trait_BE$WD),data.frame(trait_BDT$LMA,trait_BDT$WD),c("LMA","WD"),nbtstrp,T)

Ks_from_slope <- sma_plot_stats_comp(data.frame(trait_BE$slope,trait_BE$Ks),data.frame(trait_BDT$slope,trait_BDT$Ks),c("slope","Ks"),nbtstrp,T)

leafN_from_LMA <- sma_plot_stats_comp(data.frame(trait_BE$LMA,trait_BE$leafN),data.frame(trait_BDT$LMA,trait_BDT$leafN),c("LMA","leafN"),nbtstrp,T)

leafL_from_LMA <- sma_plot_stats_comp(data.frame(trait_BE$LMA,log(trait_BE$leafL)),data.frame(trait_BDT$LMA,log(trait_BDT$leafL)),c("LMA","leafL"),nbtstrp,T)

Ks_from_LMA <- sma_plot_stats_comp(data.frame(trait_BE$LMA,trait_BE$Ks),data.frame(trait_BDT$LMA,trait_BDT$Ks),c("LMA","Ks"),nbtstrp,T)
