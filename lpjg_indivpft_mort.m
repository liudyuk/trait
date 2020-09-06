

%fol='/Users/pughtam/Documents/TreeMort/Analyses/Hydraulic_modelling/LPJG/FUNDIV_indivpft_v2_morepfts/';
%fol='/Users/pughtam/Documents/TreeMort/Analyses/Hydraulic_modelling/LPJG/FUNDIV_indivpft_v3_KsLS/';
%fol='/Users/pughtam/Documents/TreeMort/Analyses/Hydraulic_modelling/LPJG/FUNDIV_indivpft_v4/';
%fol='/Users/pughtam/Documents/TreeMort/Analyses/Hydraulic_modelling/LPJG/FUNDIV_indivpft_v5/';
%fol='/Users/pughtam/Documents/TreeMort/Analyses/Hydraulic_modelling/LPJG/FUNDIV_indivpft_v5_lowcton/';
%fol='/Users/pughtam/Documents/TreeMort/Analyses/Hydraulic_modelling/LPJG/FUNDIV_indivpft_v6/';
fol='/Users/pughtam/Documents/TreeMort/Analyses/Hydraulic_modelling/LPJG/FUNDIV_indivpft_v6_cton0p05/';
%fol='/Users/pughtam/Documents/TreeMort/Analyses/Hydraulic_modelling/LPJG/FUNDIV_indivpft_v6_ctonstan/';

pft='TeBS';
npft=28;

for nn=1:npft
    filein=['AnnuallyOut_',pft,'_pft',mat2str(nn),'.nc'];
    if nn==1
        lat=ncread([fol,'/',filein],'Base/Latitude');
        lon=ncread([fol,'/',filein],'Base/Longitude');
        pfts=h5read([fol,'/',filein],'/Base/Pfts'); %Need to use h5read for strings in netcdf, see https://www.mathworks.com/matlabcentral/answers/65939-how-to-read-netcdf-files-containing-string-arrays
        years=ncread([fol,'/',filein],'Base/Time');
        nyear=length(years);
        if length(pfts)~=4 %Two grass PFTs, 1 tree PFT and total.
            error('Number of tree pfts in file is not equal to 1')
        end
    end
    cmass_in=squeeze(ncread([fol,'/',filein],'Pft-Out/cmasstotal'));
    cmasslossbg_in=squeeze(ncread([fol,'/',filein],'Pft-Out/cmass_loss_bg'));
    cmasslossgreff_in=squeeze(ncread([fol,'/',filein],'Pft-Out/cmass_loss_greff'));
    cmasslosscav_in=squeeze(ncread([fol,'/',filein],'Pft-Out/cmass_loss_cav'));
    cavxylem_in=squeeze(ncread([fol,'/',filein],'Pft-Out/cav_xylem'));
    npp_in=squeeze(ncread([fol,'/',filein],'Pft-Out/npp'));
    height_in=squeeze(ncread([fol,'/',filein],'Pft-Out/height'));
    wnpp_in=squeeze(ncread([fol,'/',filein],'Pft-Out/wnpp'));
    cmassleaf_in=squeeze(ncread([fol,'/',filein],'Pft-Out/cmassleaf'));
    nmassleaf_in=squeeze(ncread([fol,'/',filein],'Pft-Out/nmassleaf'));
    lai_in=squeeze(ncread([fol,'/',filein],'Pft-Out/lai'));
    if nn==1
        dims=size(cmass_in);
        cmass_all=NaN(npft,nyear,dims(3));
        cmasslossbg_all=NaN(npft,nyear,dims(3));
        cmasslossgreff_all=NaN(npft,nyear,dims(3));
        cmasslosscav_all=NaN(npft,nyear,dims(3));
        cavxylem_all=NaN(npft,nyear,dims(3));
        npp_all=NaN(npft,nyear,dims(3));
        height_all=NaN(npft,nyear,dims(3));
        wnpp_all=NaN(npft,nyear,dims(3));
        cmassleaf_all=NaN(npft,nyear,dims(3));
        nmassleaf_all=NaN(npft,nyear,dims(3));
        lai_all=NaN(npft,nyear,dims(3));
    end
    cmass_all(nn,:,:)=cmass_in(3,:,:);
    cmasslossbg_all(nn,:,:)=cmasslossbg_in(3,:,:);
    cmasslossgreff_all(nn,:,:)=cmasslossgreff_in(3,:,:);
    cmasslosscav_all(nn,:,:)=cmasslosscav_in(3,:,:);
    cavxylem_all(nn,:,:)=cavxylem_in(3,:,:);
    npp_all(nn,:,:)=npp_in(3,:,:);
    height_all(nn,:,:)=height_in(3,:,:);
    wnpp_all(nn,:,:)=wnpp_in(3,:,:);
    cmassleaf_all(nn,:,:)=cmassleaf_in(3,:,:);
    nmassleaf_all(nn,:,:)=nmassleaf_in(3,:,:);
    lai_all(nn,:,:)=lai_in(3,:,:);
end
clear nn cmass_in dims filein cmasslossbg_in cmasslossgreff_in cmasslosscav_in cavxylem_in npp_in height_in wnpp_in

%Check for crazy lat or lon values
aa=find(lat>90 | lon>180 | lat<-90 | lon<-180);
if ~isempty(aa)
    lat(aa)=[];
    lon(aa)=[];
    cmass_all(:,:,aa)=[];
end
clear aa

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

cmass=NaN(npft,nyear,nlats,nlons);
cmasslossbg=NaN(npft,nyear,nlats,nlons);
cmasslossgreff=NaN(npft,nyear,nlats,nlons);
cmasslosscav=NaN(npft,nyear,nlats,nlons);
cavxylem=NaN(npft,nyear,nlats,nlons);
npp=NaN(npft,nyear,nlats,nlons);
height=NaN(npft,nyear,nlats,nlons);
wnpp=NaN(npft,nyear,nlats,nlons);
cmassleaf=NaN(npft,nyear,nlats,nlons);
nmassleaf=NaN(npft,nyear,nlats,nlons);
lai=NaN(npft,nyear,nlats,nlons);
for xx=1:nlons
    for yy=1:nlats
        aa=find(lat>=lats(yy) & lat<(lats(yy)+resolution) & lon>=lons(xx) & lon<(lons(xx)+resolution));
        if length(aa)==1
            cmass(:,:,yy,xx)=cmass_all(:,:,aa);
            cmasslossbg(:,:,yy,xx)=cmasslossbg_all(:,:,aa);
            cmasslossgreff(:,:,yy,xx)=cmasslossgreff_all(:,:,aa);
            cmasslosscav(:,:,yy,xx)=cmasslosscav_all(:,:,aa);
            cavxylem(:,:,yy,xx)=cavxylem_all(:,:,aa);
            npp(:,:,yy,xx)=npp_all(:,:,aa);
            height(:,:,yy,xx)=height_all(:,:,aa);
            wnpp(:,:,yy,xx)=wnpp_all(:,:,aa);
            cmassleaf(:,:,yy,xx)=cmassleaf_all(:,:,aa);
            nmassleaf(:,:,yy,xx)=nmassleaf_all(:,:,aa);
            lai(:,:,yy,xx)=lai_all(:,:,aa);
        elseif length(aa)>1
            fprintf('Warning: more than one candidate value for location %dE, %dN\n',xx,yy)
        end
    end
end
clear xx yy aa

%Extract for last year (or decade for fluxes)
cmass_lastyear=squeeze(cmass(:,nyear,:,:));
cmasslossbg_lastdec=squeeze(nanmean(cmasslossbg(:,nyear-9:nyear,:,:),2));
cmasslossgreff_lastdec=squeeze(nanmean(cmasslossgreff(:,nyear-9:nyear,:,:),2));
cmasslosscav_lastdec=squeeze(nanmean(cmasslosscav(:,nyear-9:nyear,:,:),2));
cavxylem_lastdec=squeeze(nanmean(cavxylem(:,nyear-9:nyear,:,:),2));
npp_lastdec=squeeze(nanmean(npp(:,nyear-9:nyear,:,:),2));
height_lastyear=squeeze(height(:,nyear,:,:));
wnpp_lastdec=squeeze(nanmean(wnpp(:,nyear-9:nyear,:,:),2));
cmassleaf_lastdec=squeeze(nanmean(cmassleaf(:,nyear-9:nyear,:,:),2));
nmassleaf_lastdec=squeeze(nanmean(nmassleaf(:,nyear-9:nyear,:,:),2));
lai_lastdec=squeeze(nanmean(lai(:,nyear-9:nyear,:,:),2));

%---
%Make maps for all pfts
nsub=ceil(sqrt(npft));
maxval=max((cmass_lastyear(:)));
maxval_npp=max((npp_lastdec(:)));

figure
for nn=1:npft
    subplot(nsub,nsub,nn)
    p1=pcolor(lons,lats,squeeze(cmass_lastyear(nn,:,:))); set(p1,'linestyle','none'); colorbar
    caxis([0 maxval])
    title(sprintf('%s PFT %d',pft,nn))
    
    hold on
    C = load('coast');
    plot(C.long,C.lat,'k')
    set(gca,'Xlim',[minlon maxlon],'YLim',[minlat maxlat])
end
clear nn

%---
%Get climate data

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

%Extract region of interest only
AI_lons=-179.875:0.25:179.875;
AI_lats=-89.875:0.25:89.875;
minlonind=find(AI_lons==minlon+(resolution/2));
maxlonind=find(AI_lons==maxlon+(resolution/2));
minlatind=find(AI_lats==minlat+(resolution/2));
maxlatind=find(AI_lats==maxlat+(resolution/2));
AI_eur=AI_0p25(minlatind:maxlatind,minlonind:maxlonind);

%Get MAT and MAP from WorldClim
addpath('/data/WorldClim/')
[MAT,MAP]=worldclim_MAT_MAP(true);
MAT_eur=MAT(minlatind:maxlatind,minlonind:maxlonind);
MAP_eur=MAP(minlatind:maxlatind,minlonind:maxlonind);

%Can also get the climate from LPJ-GUESS simulations
filein=['AnnuallyOut_',pft,'_pft1.nc'];
prec_lpjg_in=squeeze(ncread([fol,'/',filein],'Patch-Out/prec'));
tas_lpjg_in=squeeze(ncread([fol,'/',filein],'Patch-Out/temp_air'));
wcont_lpjg_in=squeeze(ncread([fol,'/',filein],'Patch-Out/wcont'));
mcwd_lpjg_in=squeeze(ncread([fol,'/',filein],'Patch-Out/mcwd_daily'));
prec_lpjg=NaN(nyear,nlats,nlons);
tas_lpjg=NaN(nyear,nlats,nlons);
wcont_lpjg=NaN(nyear,nlats,nlons);
mcwd_lpjg=NaN(nyear,nlats,nlons);
for xx=1:nlons
    for yy=1:nlats
        aa=find(lat>=lats(yy) & lat<(lats(yy)+resolution) & lon>=lons(xx) & lon<(lons(xx)+resolution));
        if length(aa)==1
            prec_lpjg(:,yy,xx)=prec_lpjg_in(:,aa);
            tas_lpjg(:,yy,xx)=tas_lpjg_in(:,aa);
            wcont_lpjg(:,yy,xx)=wcont_lpjg_in(:,aa);
            mcwd_lpjg(:,yy,xx)=mcwd_lpjg_in(:,aa);
        elseif length(aa)>1
            fprintf('Warning: more than one candidate value for location %dE, %dN\n',xx,yy)
        end
    end
end
clear xx yy aa
clear filein
%Extract for last decade
prec_lastdec=squeeze(nanmean(prec_lpjg(nyear-9:nyear,:,:),1));
tas_lastdec=squeeze(nanmean(tas_lpjg(nyear-9:nyear,:,:),1));
wcont_lastdec=squeeze(nanmean(wcont_lpjg(nyear-9:nyear,:,:),1));
mcwd_lastdec=squeeze(nanmean(mcwd_lpjg(nyear-9:nyear,:,:),1));

%---
%Regress against AI for all PFTs
figure
for nn=1:npft
    subplot(nsub,nsub,nn)
    cmass_ext=squeeze(cmass_lastyear(nn,:,:));
    plot(AI_eur(:),cmass_ext(:),'.')
    set(gca,'YLim',[0 maxval])
    title(sprintf('%s PFT %d',pft,nn))
end
clear nn cmass_ext

%Regress against MAT for all PFTs
figure
for nn=1:npft
    subplot(nsub,nsub,nn)
    cmass_ext=squeeze(cmass_lastyear(nn,:,:));
    plot(MAT_eur(:),cmass_ext(:),'.')
    set(gca,'YLim',[0 maxval])
    title(sprintf('%s PFT %d',pft,nn))
end
clear nn cmass_ext

%Regress against MAP for all PFTs
figure
for nn=1:npft
    subplot(nsub,nsub,nn)
    cmass_ext=squeeze(cmass_lastyear(nn,:,:));
    plot(MAP_eur(:),cmass_ext(:),'.')
    set(gca,'YLim',[0 maxval])
    title(sprintf('%s PFT %d',pft,nn))
end
clear nn cmass_ext

%Regress against CWD for all PFTs
figure
for nn=1:npft
    subplot(nsub,nsub,nn)
    cmass_ext=squeeze(cmass_lastyear(nn,:,:));
    plot(mcwd_lastdec(:),cmass_ext(:),'.')
    set(gca,'YLim',[0 maxval])
    title(sprintf('%s PFT %d',pft,nn))
end
clear nn cmass_ext

%Regress against wcont for all PFTs
figure
for nn=1:npft
    subplot(nsub,nsub,nn)
    cmass_ext=squeeze(cmass_lastyear(nn,:,:));
    plot(wcont_lastdec(:),cmass_ext(:),'.')
    set(gca,'YLim',[0 maxval])
    title(sprintf('%s PFT %d',pft,nn))
end
clear nn cmass_ext

%Regress NPP against MAP for all PFTs
figure
for nn=1:npft
    subplot(nsub,nsub,nn)
    npp_ext=squeeze(npp_lastdec(nn,:,:));
    plot(MAP_eur(:),npp_ext(:),'.')
    set(gca,'YLim',[0 maxval_npp])
    title(sprintf('%s PFT %d',pft,nn))
end
clear nn npp_ext

%---
%Check what is going on with all the locations with zero cmass. Suspect that this may be bioclimatic limits - check
%regression against temperature. Would also be good to plot against MAP.
% -> It is clearly bioclimatic limits - zero values fall to edge of temperature range.

%---
%Extract PFT numbers with biomass above a threshold for at least 5 sites
biomass_thres=0.5; %kg C m-2
pft_above_biomass_thres=false(npft,1);
for nn=1:npft
    cmass_ext=squeeze(cmass_lastyear(nn,:,:));
    cmass_ext=cmass_ext(isnan(cmass_ext(:))==0);
    n_above_thres=length(find(cmass_ext>biomass_thres));
    if n_above_thres>=5
        pft_above_biomass_thres(nn)=true;
    end
end
clear nn cmass_ext n_above_thres
%save(['pft_above_biomass_thres_',pft,'.mat'],'pft_above_biomass_thres')

%Prepare a legend for the next plots
pftnos=1:npft;
pftlabels=cell(npft,1);
for nn=1:npft
    pftlabels{nn}=['PFT',mat2str(pftnos(nn))];
end
clear nn

%For PFTs above the biomass threshold, plot their cmass to look for performance changes
figure
hold on
aa=find(pft_above_biomass_thres==true);
for nn=1:length(aa)
    cmass_ext=squeeze(cmass_lastyear(aa(nn),:,:));
    plot(MAP_eur(:),cmass_ext(:),'.','markersize',12)
end
clear nn cmass_ext
legend(pftlabels(aa))
clear aa

%Plot only the PFT with the highest cmass in each location
figure
hold on
aa=find(pft_above_biomass_thres==true);
[cmass_max,cmass_max_ind]=nanmax(cmass_lastyear,[],1);
cc=0; dd=0;
for nn=1:length(aa)
    pftind=pftnos(aa(nn));
    bb=find(cmass_max_ind==pftind);
    if ~isempty(bb)
        cc=cc+1;
        plot(MAP_eur(bb),cmass_max(bb),'.','markersize',12)
        dd(cc)=pftind;
    end
end
clear nn bb pftind
legend(pftlabels(dd))
clear aa cc dd

%Plot only the PFT with the highest NPP in each location
figure
hold on
aa=find(pft_above_biomass_thres==true);
[npp_max,npp_max_ind]=nanmax(npp_lastdec,[],1);
for nn=1:length(aa)
    pftind=pftnos(aa(nn));
    bb=find(npp_max_ind==pftind);
    plot(MAP_eur(bb),npp_max(bb),'.','markersize',12)
end
clear nn bb pftind
legend(pftlabels(aa))
clear aa

%Mean cmass, NPP and mortality values for all PFTs above a biomass threshold for at least 5 sites
cmass_pft=NaN(npft,1);
npp_pft=NaN(npft,1);
cmasslosscav_pft=NaN(npft,1);
cmasslossgreff_pft=NaN(npft,1);
cmasslossbg_pft=NaN(npft,1);
wnpp_pft=NaN(npft,1);
for nn=1:npft
    temp_cmass=cmass_lastyear(nn,:,:);
    temp_npp=npp_lastdec(nn,:,:);
    temp_cmasslosscav=cmasslosscav_lastdec(nn,:,:);
    temp_cmasslossgreff=cmasslossgreff_lastdec(nn,:,:);
    temp_cmasslossbg=cmasslossbg_lastdec(nn,:,:);
    temp_wnpp=wnpp_lastdec(nn,:,:);
    cmass_pft(nn)=nanmean(temp_cmass(:));
    npp_pft(nn)=nanmean(temp_npp(:));
    cmasslosscav_pft(nn)=nanmean(temp_cmasslosscav(:));
    cmasslossgreff_pft(nn)=nanmean(temp_cmasslossgreff(:));
    cmasslossbg_pft(nn)=nanmean(temp_cmasslossbg(:));
    wnpp_pft(nn)=nanmean(temp_wnpp(:));
end
clear nn temp_cmass temp_npp temp_cmasslosscav temp_cmasslossgreff temp_cmasslossbg
aa=find(pft_above_biomass_thres==true);
figure
hold on
subplot(7,1,1)
bar(cmass_pft(aa))
ylabel('Veg C')
set(gca,'XTickLabel',pftlabels(aa),'XTick',1:length(aa))
subplot(7,1,2)
bar(npp_pft(aa))
ylabel('NPP')
set(gca,'XTickLabel',pftlabels(aa),'XTick',1:length(aa))
subplot(7,1,3)
bar(wnpp_pft(aa)./npp_pft(aa))
ylabel('wNPP/NPP')
set(gca,'XTickLabel',pftlabels(aa),'XTick',1:length(aa))
subplot(7,1,4)
bar(cmasslosscav_pft(aa)./cmass_pft(aa))
ylabel('Mort rate cavitation')
set(gca,'XTickLabel',pftlabels(aa),'XTick',1:length(aa))
subplot(7,1,5)
bar(cmasslosscav_pft(aa)./cmass_pft(aa))
ylabel('Mort rate cavitation')
set(gca,'XTickLabel',pftlabels(aa),'XTick',1:length(aa))
set(gca,'YLim',[0 0.01])
subplot(7,1,6)
bar(cmasslossgreff_pft(aa)./cmass_pft(aa))
ylabel('Mort rate greff')
set(gca,'XTickLabel',pftlabels(aa),'XTick',1:length(aa))
set(gca,'YLim',[0 0.01])
subplot(7,1,7)
bar(cmasslossbg_pft(aa)./cmass_pft(aa))
ylabel('Mort rate background')
set(gca,'XTickLabel',pftlabels(aa),'XTick',1:length(aa))
set(gca,'YLim',[0 0.01])
clear aa

%Calculate diversity - number of species above a biomass threshold for each site
%For biomass threshold use 10% of maximum simulated for any PFT at that site.
biomass_thres=0.5; %kg C m-2 (set a minimum threshold for consideration at all)
cmass_div_count=NaN(nlats,nlons);
for xx=1:nlons
    for yy=1:nlats
        maxbiom=max(cmass_lastyear(:,yy,xx),[],1);
        if maxbiom>biomass_thres
            biomass_thres_div=max(maxbiom*0.1,biomass_thres);
            cmass_div_count(yy,xx)=length(find(cmass_lastyear(:,yy,xx)>biomass_thres_div));
        end
    end
end
clear xx yy aa biomass_thres_div maxbiom

figure
hold on
%plot(MAP_eur(:),cmass_div_count(:),'.','markersize',12)
plot(mcwd_lastdec(:),cmass_div_count(:),'.','markersize',12)
xlabel('MCWD')
ylabel('No. viable PFTs')

figure
hold on
plot(MAT_eur(:),cmass_div_count(:),'.','markersize',12)
xlabel('MAT')
ylabel('No. viable PFTs')

%Calculate diversity - number of species above an NPP threshold for each site
%For NPP threshold use 50% of maximum simulated for any PFT at that site.
biomass_thres=0.5; %kg C m-2 (set a minimum threshold for consideration at all)
npp_div_count=NaN(nlats,nlons);
for xx=1:nlons
    for yy=1:nlats
        maxbiom=max(cmass_lastyear(:,yy,xx),[],1);
        maxnpp=max(npp_lastdec(:,yy,xx),[],1);
        if maxbiom>biomass_thres
            npp_thres_div=maxnpp*0.5;
            npp_div_count(yy,xx)=length(find(npp_lastdec(:,yy,xx)>npp_thres_div & cmass_lastyear(:,yy,xx)>biomass_thres));
        end
    end
end
clear xx yy aa npp_thres_div maxbiom maxnpp

figure
hold on
plot(MAT_eur(:),npp_div_count(:),'.','markersize',12)

%Plot mortality rates by climate gradient for selected PFTs
pftind=2;
cmass_sel=squeeze(cmass_lastyear(pftind,:,:));
cmasslossbg_sel=squeeze(cmasslossbg_lastdec(pftind,:,:));
cmasslossgreff_sel=squeeze(cmasslossgreff_lastdec(pftind,:,:));
cmasslosscav_sel=squeeze(cmasslosscav_lastdec(pftind,:,:));

figure
hold on
plot(MAP_eur(:),cmasslossbg_sel(:)./cmass_sel(:),'.','markersize',12)
plot(MAP_eur(:),cmasslossgreff_sel(:)./cmass_sel(:),'.','markersize',12)
plot(MAP_eur(:),cmasslosscav_sel(:)./cmass_sel(:),'.','markersize',12)
xlabel('MAP')
title(['PFT',mat2str(pftind)])

figure
hold on
plot(AI_eur(:),cmasslossbg_sel(:)./cmass_sel(:),'.','markersize',12)
plot(AI_eur(:),cmasslossgreff_sel(:)./cmass_sel(:),'.','markersize',12)
plot(AI_eur(:),cmasslosscav_sel(:)./cmass_sel(:),'.','markersize',12)
xlabel('AI')

figure
hold on
plot(MAT_eur(:),cmasslossbg_sel(:)./cmass_sel(:),'.','markersize',12)
plot(MAT_eur(:),cmasslossgreff_sel(:)./cmass_sel(:),'.','markersize',12)
plot(MAT_eur(:),cmasslosscav_sel(:)./cmass_sel(:),'.','markersize',12)
xlabel('MAT')

%Plot xylem cavitation by climate gradient for selected PFTs
pftind=3;
cavxylem_sel=squeeze(cavxylem_lastdec(pftind,:,:));

figure
subplot(1,3,1)
hold on
plot(MAP_eur(:),cavxylem_sel(:),'.','markersize',12)
xlabel('MAP')
title(['PFT',mat2str(pftind)])

subplot(1,3,2)
hold on
plot(AI_eur(:),cavxylem_sel(:),'.','markersize',12)
xlabel('AI')

subplot(1,3,3)
hold on
plot(MAT_eur(:),cavxylem_sel(:),'.','markersize',12)
xlabel('MAT')

%Plot height by climate gradient for selected PFTs
pftind=23;
height_sel=squeeze(height_lastyear(pftind,:,:));

figure
subplot(1,3,1)
hold on
plot(MAP_eur(:),height_sel(:),'.','markersize',12)
xlabel('MAP')
title(['PFT',mat2str(pftind)])

subplot(1,3,2)
hold on
plot(AI_eur(:),height_sel(:),'.','markersize',12)
xlabel('AI')

subplot(1,3,3)
hold on
plot(MAT_eur(:),height_sel(:),'.','markersize',12)
xlabel('MAT')

%---
%Merge with traits

%Read in the trait data
%traits=readtable('LPJG_PFT_summary_TrBE.csv');
traits=readtable('LPJG_PFT_summary_TeBS.csv');
%traits=readtable('morePFTs/LPJG_PFT_summary_TeBS.csv');

%Calculated weighted trait means for each grid cell.
P50w=squeeze(nansum(cmass_lastyear.*repmat(traits.P50,[1 nlats nlons]),1));
P50w=P50w./squeeze(sum(cmass_lastyear,1));
TLPw=squeeze(nansum(cmass_lastyear.*repmat(traits.TLP,[1 nlats nlons]),1));
TLPw=TLPw./squeeze(sum(cmass_lastyear,1));
WDw=squeeze(nansum(cmass_lastyear.*repmat(traits.WD,[1 nlats nlons]),1));
WDw=WDw./squeeze(sum(cmass_lastyear,1));
slopew=squeeze(nansum(cmass_lastyear.*repmat(traits.slope,[1 nlats nlons]),1));
slopew=slopew./squeeze(sum(cmass_lastyear,1));
SLAw=squeeze(nansum(cmass_lastyear.*repmat(traits.SLA,[1 nlats nlons]),1));
SLAw=SLAw./squeeze(sum(cmass_lastyear,1));
LSw=squeeze(nansum(cmass_lastyear.*repmat(traits.LS,[1 nlats nlons]),1));
LSw=LSw./squeeze(sum(cmass_lastyear,1));
Ksw=squeeze(nansum(cmass_lastyear.*repmat(traits.Ks,[1 nlats nlons]),1));
Ksw=Ksw./squeeze(sum(cmass_lastyear,1));

%Regress trait means against AI
[latsgrid,lonsgrid]=meshgrid(lats,lons);
indtoplot=find(latsgrid<80); 

figure
subplot(2,4,1); hold on
plot(AI_eur(indtoplot),P50w(indtoplot),'r.')
title('P50')
subplot(2,4,2); hold on
plot(AI_eur(indtoplot),TLPw(indtoplot),'r.')
title('TLP')
subplot(2,4,3); hold on
plot(AI_eur(indtoplot),WDw(indtoplot),'r.')
title('WD')
subplot(2,4,4); hold on
plot(AI_eur(indtoplot),slopew(indtoplot),'r.')
title('slope')
subplot(2,4,5); hold on
plot(AI_eur(indtoplot),SLAw(indtoplot),'r.')
title('SLA')
subplot(2,4,6); hold on
plot(AI_eur(indtoplot),LSw(indtoplot),'r.')
title('LS')
subplot(2,4,7); hold on
plot(AI_eur(indtoplot),Ksw(indtoplot),'r.')
title('Ks')

figure
subplot(2,4,1); hold on
plot(prec_lastdec(indtoplot),P50w(indtoplot),'b.')
title('P50')
subplot(2,4,2); hold on
plot(prec_lastdec(indtoplot),TLPw(indtoplot),'b.')
title('TLP')
subplot(2,4,3); hold on
plot(prec_lastdec(indtoplot),WDw(indtoplot),'b.')
title('WD')
subplot(2,4,4); hold on
plot(prec_lastdec(indtoplot),slopew(indtoplot),'b.')
title('slope')
subplot(2,4,5); hold on
plot(prec_lastdec(indtoplot),SLAw(indtoplot),'b.')
title('SLA')
subplot(2,4,6); hold on
plot(prec_lastdec(indtoplot),LSw(indtoplot),'b.')
title('LS')
subplot(2,4,7); hold on
plot(prec_lastdec(indtoplot),Ksw(indtoplot),'b.')
title('Ks')

%---
%Plot trait values as ranges
figure
for xx=1:nlons
    for yy=1:nlats
        temp=cmass_lastyear(:,yy,xx);
        indall=find(temp>2);
        indmax=find(temp==nanmax(temp));
        if ~isempty(indall)
            subplot(3,3,1); 
            hold on
            plot([mcwd_lastdec(yy,xx),mcwd_lastdec(yy,xx)],[min(traits.TLP(indall)),max(traits.TLP(indall))],'color',[0.5 0.5 0.5])
            plot(mcwd_lastdec(yy,xx),traits.TLP(indmax),'k.','markersize',14)
            ylabel('TLP (MPa)')
            
            subplot(3,3,2); 
            hold on
            plot([mcwd_lastdec(yy,xx),mcwd_lastdec(yy,xx)],[min(traits.P50(indall)),max(traits.P50(indall))],'color',[0.5 0.5 0.5])
            plot(mcwd_lastdec(yy,xx),traits.P50(indmax),'k.','markersize',14)
            ylabel('P50 (MPa)')
            
            subplot(3,3,3); 
            hold on
            plot([mcwd_lastdec(yy,xx),mcwd_lastdec(yy,xx)],[min(traits.P88(indall)),max(traits.P88(indall))],'color',[0.5 0.5 0.5])
            plot(mcwd_lastdec(yy,xx),traits.P88(indmax),'k.','markersize',14)
            ylabel('P88 (MPa)')
            
            subplot(3,3,4); 
            hold on
            plot([mcwd_lastdec(yy,xx),mcwd_lastdec(yy,xx)],[min(traits.slope(indall)),max(traits.slope(indall))],'color',[0.5 0.5 0.5])
            plot(mcwd_lastdec(yy,xx),traits.slope(indmax),'k.','markersize',14)
            ylabel('Cavitation slope (% MPa^{-1})')
            
            subplot(3,3,5); 
            hold on
            plot([mcwd_lastdec(yy,xx),mcwd_lastdec(yy,xx)],[min(traits.Ks(indall)),max(traits.Ks(indall))],'color',[0.5 0.5 0.5])
            plot(mcwd_lastdec(yy,xx),traits.Ks(indmax),'k.','markersize',14)
            ylabel('Ks (kg m^{-1} s^{-1} MPa^{-1})')
            
            subplot(3,3,6); 
            hold on
            plot([mcwd_lastdec(yy,xx),mcwd_lastdec(yy,xx)],[min(traits.LS(indall)),max(traits.LS(indall))],'color',[0.5 0.5 0.5])
            plot(mcwd_lastdec(yy,xx),traits.LS(indmax),'k.','markersize',14)
            ylabel('LS (m^{2} cm^{-2})')
            
            subplot(3,3,7); 
            hold on
            plot([mcwd_lastdec(yy,xx),mcwd_lastdec(yy,xx)],[min(traits.WD(indall)),max(traits.WD(indall))],'color',[0.5 0.5 0.5])
            plot(mcwd_lastdec(yy,xx),traits.WD(indmax),'k.','markersize',14)
            ylabel('WD (kg C m^{-2})')
            
            subplot(3,3,8); 
            hold on
            plot([mcwd_lastdec(yy,xx),mcwd_lastdec(yy,xx)],[min(traits.SLA(indall)),max(traits.SLA(indall))],'color',[0.5 0.5 0.5])
            plot(mcwd_lastdec(yy,xx),traits.SLA(indmax),'k.','markersize',14)
            ylabel('SLA (m^2 kg^{-1} C)')
            
            subplot(3,3,9); 
            hold on
            plot([mcwd_lastdec(yy,xx),mcwd_lastdec(yy,xx)],[min(traits.leafN_LPJG(indall)),max(traits.leafN_LPJG(indall))],'color',[0.5 0.5 0.5])
            plot(mcwd_lastdec(yy,xx),traits.leafN_LPJG(indmax),'k.','markersize',14)
            ylabel('leafN (mg N g^{-1})')
        end
    end
end

%---
% %Make a plot of model C:N ratio against SLA
% %CN_lastdec=cmassleaf_lastdec./nmassleaf_lastdec;
% CN_pft=NaN(npft,1);
% for nn=1:npft
%     CN_temp=(lai_lastdec(nn,:,:)*traits.SLA(nn))./nmassleaf_lastdec(nn,:,:);
%     CN_pft(nn)=nanmean(CN_temp(cmass_lastyear(nn,:,:)>2));
% end
% figure
% plot(traits.SLA,CN_pft,'.')

