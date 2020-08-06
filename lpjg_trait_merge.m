

fol='/Users/pughtam/Documents/TreeMort/Analyses/Hydraulic_modelling/LPJG/FUNDIV_TeBS_fullgridlist/';

lat=ncread([fol,'/AnnuallyOut.nc'],'Base/Latitude');
lon=ncread([fol,'/AnnuallyOut.nc'],'Base/Longitude');
cmass_in=squeeze(ncread([fol,'/AnnuallyOut.nc'],'Pft-Out/cmasstotal'));
pfts=h5read([fol,'/AnnuallyOut.nc'],'/Base/Pfts'); %Need to use h5read for strings in netcdf, see https://www.mathworks.com/matlabcentral/answers/65939-how-to-read-netcdf-files-containing-string-arrays
years=ncread([fol,'/AnnuallyOut.nc'],'Base/Time');
npft=length(pfts);
nyear=length(years);

%Format cmass into a grid

minlat=floor(min(lat));
maxlat=ceil(max(lat));
minlon=floor(min(lon));
maxlon=ceil(max(lon));
resolution=0.25;
lats=minlat:resolution:maxlat;
lons=minlon:resolution:maxlon;
nlons=length(lons);
nlats=length(lats);

cmass=NaN(nlats,nlons,npft,nyear);
for xx=1:nlons
    for yy=1:nlats
        aa=find(lat>=lats(yy) & lat<(lats(yy)+resolution) & lon>=lons(xx) & lon<(lons(xx)+resolution));
        if length(aa)==1
            cmass(yy,xx,:,:)=cmass_in(:,:,aa);
            aa
        elseif length(aa)>1
            fprintf('Warning: more than one candidate value for location %dE, %dN\n',xx,yy)
        end
    end
end
clear xx yy aa

%Extract tree data only (grass is first 2 PFTs, last is total)
cmass_tree=cmass(:,:,3:npft-1,:);

%Extract for last year
cmass_tree_lastyear=cmass_tree(:,:,:,nyear);

%Read in the trait data
traits=readtable('LPJG_PFT_summary_TeBS.csv');

%Calculated weighted trait means for each grid cell.
P50w=nansum(cmass_tree_lastyear.*permute(repmat(traits.P50,[1 nlats nlons]),[2 3 1]),3);
P50w=P50w./sum(cmass_tree_lastyear,3);
TLPw=nansum(cmass_tree_lastyear.*permute(repmat(traits.TLP,[1 nlats nlons]),[2 3 1]),3);
TLPw=TLPw./sum(cmass_tree_lastyear,3);
WDw=nansum(cmass_tree_lastyear.*permute(repmat(traits.WD,[1 nlats nlons]),[2 3 1]),3);
WDw=WDw./sum(cmass_tree_lastyear,3);
slopew=nansum(cmass_tree_lastyear.*permute(repmat(traits.slope,[1 nlats nlons]),[2 3 1]),3);
slopew=slopew./sum(cmass_tree_lastyear,3);
SLAw=nansum(cmass_tree_lastyear.*permute(repmat(traits.SLA,[1 nlats nlons]),[2 3 1]),3);
SLAw=SLAw./sum(cmass_tree_lastyear,3);
LSw=nansum(cmass_tree_lastyear.*permute(repmat(traits.LS,[1 nlats nlons]),[2 3 1]),3);
LSw=LSw./sum(cmass_tree_lastyear,3);
Ksw=nansum(cmass_tree_lastyear.*permute(repmat(traits.Ks,[1 nlats nlons]),[2 3 1]),3);
Ksw=Ksw./sum(cmass_tree_lastyear,3);

figure
subplot(2,4,1)
p1=pcolor(lons,lats,P50w); set(p1,'linestyle','none'); colorbar
title('P50')
subplot(2,4,2)
p2=pcolor(lons,lats,TLPw); set(p2,'linestyle','none'); colorbar
title('TLP')
subplot(2,4,3)
p3=pcolor(lons,lats,WDw); set(p3,'linestyle','none'); colorbar
title('WD')
subplot(2,4,4)
p4=pcolor(lons,lats,slopew); set(p4,'linestyle','none'); colorbar
title('slope')
subplot(2,4,5)
p5=pcolor(lons,lats,SLAw); set(p5,'linestyle','none'); colorbar
title('SLA')
subplot(2,4,6)
p6=pcolor(lons,lats,LSw); set(p6,'linestyle','none'); colorbar
title('LS')
subplot(2,4,7)
p7=pcolor(lons,lats,Ksw); set(p7,'linestyle','none'); colorbar
title('Ks')

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

%Extract Europe only
AI_lons=-179.875:0.25:179.875;
AI_lats=-89.875:0.25:89.875;
minlonind=find(AI_lons==minlon+(resolution/2));
maxlonind=find(AI_lons==maxlon+(resolution/2));
minlatind=find(AI_lats==minlat+(resolution/2));
maxlatind=find(AI_lats==maxlat+(resolution/2));
AI_eur=AI_0p25(minlatind:maxlatind,minlonind:maxlonind);

%Regress trait means against AI
[latsgrid,lonsgrid]=meshgrid(lats,lons);
indtoplot=find(latsgrid<55); %Exclude mid-northern Sweden, something weird is happening

figure
subplot(2,4,1); hold on
plot(AI_eur(indtoplot),P50w(indtoplot),'.')
title('P50')
subplot(2,4,2); hold on
plot(AI_eur(indtoplot),TLPw(indtoplot),'.')
title('TLP')
subplot(2,4,3); hold on
plot(AI_eur(indtoplot),WDw(indtoplot),'.')
title('WD')
subplot(2,4,4); hold on
plot(AI_eur(indtoplot),slopew(indtoplot),'.')
title('slope')
subplot(2,4,5); hold on
plot(AI_eur(indtoplot),SLAw(indtoplot),'.')
title('SLA')
subplot(2,4,6); hold on
plot(AI_eur(indtoplot),LSw(indtoplot),'.')
title('LS')
subplot(2,4,7); hold on
plot(AI_eur(indtoplot),Ksw(indtoplot),'.')
title('Ks')

