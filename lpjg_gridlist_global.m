
grid_file='/Users/pughtam/LPJG/trunk_r8278/benchmarks/global/config/gridlist.txt';

grid_in=dlmread(grid_file,'',0,0);
grid_lon=grid_in(:,1);
grid_lat=grid_in(:,2);

%Read year 2000 canopy cover map (Hansen et al., 2013, Science) and make a mask for gridcells with > 5% canopy cover
forest_file='/Users/pughtam/Documents/TreeMort/Analyses/Temperate_dist/TempBoreal/hansen_forested_canopy_frac_0p5deg.nc4';

can_frac=ncread(forest_file,'canopy_cover_frac')';
can_lats=ncread(forest_file,'Latitude');
can_lons=ncread(forest_file,'Longitude');
can_mask=false(size(can_frac));
can_mask(can_frac>5)=true;

%Remove all cells in gridlist that lie outside the canopy cover mask
mask_out=false(size(grid_lat));
for nn=1:length(grid_lat)
    lat_ind=find(can_lats==grid_lat(nn));
    lon_ind=find(can_lons==grid_lon(nn));
    if can_mask(lat_ind,lon_ind)==false
        mask_out(nn)=true;
    end
end
clear nn lat_ind lon_ind

grid_lat_masked=grid_lat;
grid_lon_masked=grid_lon;
grid_lat_masked(mask_out)=[];
grid_lon_masked(mask_out)=[];

figure
plot(grid_lon_masked,grid_lat_masked,'.')

%Randomly make a selection of gridcells
nsel=2000;
aa=randsample(length(grid_lat_masked),nsel);

figure
hold on
plot(grid_lon_masked,grid_lat_masked,'.')
plot(grid_lon_masked(aa),grid_lat_masked(aa),'r.')

%Split into tropical and extratropical sections
grid_lon_masked_trop=grid_lon_masked(aa);
grid_lat_masked_trop=grid_lat_masked(aa);
grid_lon_masked_extratrop=grid_lon_masked(aa);
grid_lat_masked_extratrop=grid_lat_masked(aa);

grid_lon_masked_trop(abs(grid_lat_masked(aa))>23)=[];
grid_lat_masked_trop(abs(grid_lat_masked(aa))>23)=[];
grid_lon_masked_extratrop(abs(grid_lat_masked(aa))<=23)=[];
grid_lat_masked_extratrop(abs(grid_lat_masked(aa))<=23)=[];

%Write out new gridlists
writetable(table(grid_lon_masked_trop+0.125,grid_lat_masked_trop+0.125),'gridlist_global_0p25_sample2000_trop.txt','Delimiter',' ','WriteRowNames',false,'WriteVariableNames',false)
writetable(table(grid_lon_masked_extratrop+0.125,grid_lat_masked_extratrop+0.125),'gridlist_global_0p25_sample2000_extratrop.txt','Delimiter',' ','WriteRowNames',false,'WriteVariableNames',false)

