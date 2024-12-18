clear
close all;

dataDir = '/Users/carrolld/Documents/research/bathy/raw_data/GEBCO_2024/';
saveDir = '/Users/carrolld/Documents/research/bathy/mat/bathy/GEBCO_2024_raw/';

%% 

fileName = 'GEBCO_2024_sub_ice_topo.nc';

ncdisp([dataDir fileName]);

%% 

lon = ncread([dataDir fileName],'lon');
lat = ncread([dataDir fileName],'lat');

%[xx yy] = meshgrid(lon,lat);

%% 

elevation = ncread([dataDir fileName],'elevation');

temp = elevation';

%% 

imagesc(lon,lat,temp);

ax = gca;
ax.YDir= 'normal';
 
caxis([-5000 0]);

colormap jet

colorbar

%% 

save([saveDir 'GEBCO_2024_lon_lat.mat'],'lon','lat','elevation','-v7.3');

%%
