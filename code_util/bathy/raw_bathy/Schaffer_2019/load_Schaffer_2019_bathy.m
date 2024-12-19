clear
close all;

dataDir = '/Users/carrolld/Documents/research/bathy/raw_data/Schaffer_2019/';
saveDir = '/Users/carrolld/Documents/research/bathy/mat/bathy/Schaffer_2019_raw/';

%% 

fileName = 'RTopo-2.0.4_30sec_bedrock_topography.nc';

ncdisp([dataDir fileName]);

%% 

lon = double(ncread([dataDir fileName],'lon'));
lat = double(ncread([dataDir fileName],'lat'));

%[xx yy] = meshgrid(lon,lat);

%% 

elevation = ncread([dataDir fileName],'bedrock_topography');

temp = elevation';

%% 

imagesc(lon,lat,temp);

ax = gca;
ax.YDir= 'normal';
 
caxis([-5000 0]);

colormap jet

colorbar

%% 

save([saveDir 'Schaffer_2019_bathy_lon_lat.mat'],'lon','lat','elevation','-v7.3');

%%
