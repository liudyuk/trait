%Script to calculate community-weighted means of hydraulic trait data for inventory plots along an aridity gradient.
%
%T. Pugh
%03.08.20

%Plot summary data provided by Adriane
plotFUN=readtable('/Users/pughtam/Documents/TreeMort/Analyses/Hydraulic_modelling/inv_Europe_spec/Stand-level-dynamics_v4.csv');
nplot=height(plotFUN);

%Find all species types in the dataset
[allspecies,ai]=unique(plotFUN.dom_50_agb_spp);
nspecies=length(allspecies);


%Check the availability of trait information for these species

traitfile='/Users/pughtam/Documents/TreeMort/Analyses/Hydraulic_modelling/Traits/mytrait-data/woody_trait.0803_NaN.txt';
traits=readtable(traitfile,'ReadVariableNames',true,'Delimiter','\t','Headerlines',0);

plotFUN.P50=NaN(height(plotFUN),1);
plotFUN.TLP=NaN(height(plotFUN),1);
plotFUN.LMA=NaN(height(plotFUN),1);
plotFUN.Ks=NaN(height(plotFUN),1);
plotFUN.WD=NaN(height(plotFUN),1);
plotFUN.slope=NaN(height(plotFUN),1);
plotFUN.LS=NaN(height(plotFUN),1);
plotFUN.Hmax=NaN(height(plotFUN),1);
plotFUN.Ks_Hmax=NaN(height(plotFUN),1);
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
    if (mod(nn,1000)==0)
        fprintf('nn is %d out of %d\n',nn,nplot)
    end
end
clear nn aa

%Convert traits to units used in LPJ-GUESS for easier comparison later
plotFUN.P50_LPJG=-exp(plotFUN.P50);
plotFUN.TLP_LPJG=-exp(plotFUN.TLP);
plotFUN.SLA_LPJG=(1./exp(plotFUN.LMA))*1000/2;
plotFUN.Ks_LPJG=exp(plotFUN.Ks);
plotFUN.WD_LPJG=(plotFUN.WD*1000)/2;
plotFUN.slope_LPJG=exp(plotFUN.slope);
plotFUN.LS_LPJG=exp(plotFUN.LS)*10000;

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
SLA_mean=NaN(nlats_eur,nlats_eur);
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
            P50_mean(yy,xx)=nanmean(plotFUN.P50_LPJG(aa));
            TLP_mean(yy,xx)=nanmean(plotFUN.TLP_LPJG(aa));
            LMA_mean(yy,xx)=nanmean(plotFUN.LMA(aa));
            SLA_mean(yy,xx)=nanmean(plotFUN.SLA_LPJG(aa));
            Ks_mean(yy,xx)=nanmean(plotFUN.Ks_LPJG(aa));
            WD_mean(yy,xx)=nanmean(plotFUN.WD_LPJG(aa));
            slope_mean(yy,xx)=nanmean(plotFUN.slope_LPJG(aa));
            LS_mean(yy,xx)=nanmean(plotFUN.LS_LPJG(aa));
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
[latsgrid,lonsgrid]=meshgrid(lats_eur,lons_eur);
indtoplot=find(latsgrid<55); %Exclude mid-northern Sweden

figure
subplot(2,4,1); hold on
plot(AI_eur(indtoplot),P50_mean(indtoplot),'r.')
title('P50')
subplot(2,4,2); hold on
plot(AI_eur(indtoplot),TLP_mean(indtoplot),'r.')
title('TLP')
subplot(2,4,3); hold on
plot(AI_eur(indtoplot),WD_mean(indtoplot),'r.')
title('WD')
subplot(2,4,4); hold on
plot(AI_eur(indtoplot),slope_mean(indtoplot),'r.')
title('slope')
subplot(2,4,5); hold on
plot(AI_eur(indtoplot),SLA_mean(indtoplot),'r.')
title('SLA')
subplot(2,4,6); hold on
plot(AI_eur(indtoplot),LS_mean(indtoplot),'r.')
title('LS')
subplot(2,4,7); hold on
plot(AI_eur(indtoplot),Ks_mean(indtoplot),'r.')
title('Ks')

%---
%Output a gridlist for LPJ-GUESS using the plot locations

writetable(table(grid_lon'+0.125,grid_lat'+0.125),'gridlist_FUNDIV_0p25.txt','Delimiter',' ','WriteRowNames',false,'WriteVariableNames',false)

%Coarsen to 0.5 degree resolution
%grid_lon_0p5=(round((grid_lon+0.125)*2)/2)+0.25;
%grid_lat_0p5=(round((grid_lat+0.125)*2)/2)+0.25;

%grids_0p5=table(grid_lon_0p5',grid_lat_0p5');
%grids_0p5=unique(grids_0p5);

%writetable(grids_0p5,'gridlist_FUNDIV_0p5.txt','Delimiter',' ','WriteRowNames',false,'WriteVariableNames',false)

aa=randsample(length(grid_lat),30);
writetable(table(grid_lon(aa)'+0.125,grid_lat(aa)'+0.125),'gridlist_FUNDIV_0p25_sample.txt','Delimiter',' ','WriteRowNames',false,'WriteVariableNames',false)

