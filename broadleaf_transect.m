%Script to calculate community-weighted means of hydraulic trait data for inventory plots along an aridity gradient.
%
%T. Pugh
%03.08.20

%Plot summary data provided by Adriane
plotFUN=readtable('Stand-level-dynamics_v4.csv');
nplot=height(plotFUN);

%Find all species types in the dataset
[allspecies,ai]=unique(plotFUN.dom_50_agb_spp);
nspecies=length(allspecies);


%Check the availability of trait information for these species

traitfile='/Users/pughtam/Documents/TreeMort/Analyses/Hydraulic_modelling/Traits/mytrait-data/woody_trait.0803_NaN.txt';
traits=readtable(traitfile,'ReadVariableNames',true,'Delimiter','\t','Headerlines',0);

for nn=1:nplot
    aa=find(strcmp(traits.x,plotFUN.dom_50_agb_spp(nn)));
    if ~isempty(aa)
        plotFUN.group(nn)=traits.group(aa);
        plotFUN.P50(nn)=traits.P50(aa);
        plotFUN.TLP(nn)=traits.TLP(aa);
        plotFUN.LMA(nn)=traits.LMA(aa);
        plotFUN.Ks(nn)=traits.Ks(aa);
        plotFUN.WD(nn)=traits.WD(aa);
        plotFUN.slope(nn)=traits.slope(aa);
        plotFUN.LS(nn)=traits.LS(aa);
        plotFUN.Hmax(nn)=traits.Hmax(aa);
        plotFUN.Ks_Hmax(nn)=traits.Ks_Hmax(aa);
    end
end
clear nn aa

%Mark on a map the locations of the plots dominated by broadleaf species are found
figure
hold on
C = load('coast');
plot(C.long,C.lat,'k')
set(gca,'Xlim',[-20 20],'YLim',[25 65])

ind_toplot=find(strcmp(plotFUN.group,'BE') | strcmp(plotFUN.group,'BD') | strcmp(plotFUN.group,'BT'));
plot(plotFUN.longitude(ind_toplot),plotFUN.latitude(ind_toplot),'b.')

%Aggregate plot data to 0.25 degrees (bottom-left referencing)
%Calculate trait means per 0.25 degree grid cell for broadleaf species (currently based on plot dominance only)
lons_eur=-20:0.25:19.75;
lats_eur=25:0.25:64.75;
nlons_eur=length(lons_eur);
nlats_eur=length(lats_eur);

P50_mean=NaN(nlats_eur,nlats_eur);
TLP_mean=NaN(nlats_eur,nlats_eur);
LMA_mean=NaN(nlats_eur,nlats_eur);
Ks_mean=NaN(nlats_eur,nlats_eur);
WD_mean=NaN(nlats_eur,nlats_eur);
slope_mean=NaN(nlats_eur,nlats_eur);
LS_mean=NaN(nlats_eur,nlats_eur);
Hmax_mean=NaN(nlats_eur,nlats_eur);
Ks_Hmax_mean=NaN(nlats_eur,nlats_eur);
grid_lon=NaN;
grid_lat=NaN;
cc=0;
for xx=1:nlons_eur
    for yy=1:nlats_eur
        aa=find(plotFUN.longitude>=lons_eur(xx) & plotFUN.longitude<lons_eur(xx)+0.25 &...
            plotFUN.latitude>=lats_eur(yy) & plotFUN.latitude<lats_eur(yy)+0.25);

        if length(aa)>5 %Set a minimum of at least 5 plots to calculate stats over %NOTE: This assumption needs revisiting!
            P50_mean(yy,xx)=nanmean(plotFUN.P50(aa));
            TLP_mean(yy,xx)=nanmean(plotFUN.TLP(aa));
            LMA_mean(yy,xx)=nanmean(plotFUN.LMA(aa));
            Ks_mean(yy,xx)=nanmean(plotFUN.Ks(aa));
            WD_mean(yy,xx)=nanmean(plotFUN.WD(aa));
            slope_mean(yy,xx)=nanmean(plotFUN.slope(aa));
            LS_mean(yy,xx)=nanmean(plotFUN.LS(aa));
            Hmax_mean(yy,xx)=nanmean(plotFUN.Hmax(aa));
            Ks_Hmax_mean(yy,xx)=nanmean(plotFUN.Ks_Hmax(aa));
            
            cc=cc+1;
            grid_lon(cc)=lons_eur(xx);
            grid_lat(cc)=lats_eur(yy);
        end
    end
end
clear xx yy aa

%Read in the AI data at 30 second resolution
[AI,cmap,refmat,bbox]=geotiffread('/Users/pughtam/data/CGIARCSI_AI_annual/ai_et0/ai_et0.tif');
AI=flipud(AI);

%Aggregate to 0.25 degrees
AI_0p25=NaN(720,1440);
for xx=1:1440
    for yy=121:720
        yyy=yy-120;
        xx_s=(xx*30)-29;
        xx_e=xx*30;
        yy_s=(yyy*30)-29;
        yy_e=yyy*30;
        temp=AI(yy_s:yy_e,xx_s:xx_e);
        temp(temp<0)=NaN;
        AI_0p25(yy,xx)=nanmean(double(temp(:)));
    end
end
clear xx yy yyy xx_s xx_e yy_s yy_e temp
AI_0p25=AI_0p25/10000;

%Extract Europe only (20W-20E,25N-65N)
AI_eur=AI_0p25(460:619,640:799);
     
%Regress trait means against AI



%Output a gridlist for LPJ-GUESS using the plot locations
writetable(table(grid_lon'+0.125,grid_lat'+0.125),'gridlist_FUNDIV_0p25.txt','Delimiter',' ','WriteRowNames',false,'WriteVariableNames',false)


