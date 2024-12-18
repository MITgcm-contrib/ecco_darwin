clear
close all;

dataDir = '/Users/carrolld/Documents/research/bathy/raw_data/Schaffer_2016/';
saveDir = '/Users/carrolld/Documents/research/bathy/mat/bathy/Schaffer_2016_raw/';

%% 

fileName = 'RTopo-2.0.1_30sec_ice_base_topography.nc';

ncdisp([dataDir fileName]);

%% 

lon = ncread([dataDir fileName],'lon');
lat = ncread([dataDir fileName],'lat');

%[xx yy] = meshgrid(lon,lat);

%% 

base = ncread([dataDir fileName],'ice_base_topography');

temp = base';

%% 

imagesc(lon,lat,temp);

ax = gca;
ax.YDir= 'normal';
 
caxis([-5000 0]);

colormap jet

colorbar

%% 

save([saveDir 'Schaffer_2016_ice_lon_lat.mat'],'lon','lat','base','-v7.3');

%%
