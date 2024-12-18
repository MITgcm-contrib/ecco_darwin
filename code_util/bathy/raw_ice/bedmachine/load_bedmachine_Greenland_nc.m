clear
close all;

dataDir = '/Users/carrolld/Documents/research/bathy/raw_data/bedmachine/';
saveDir = '/Users/carrolld/Documents/research/bathy/mat/bedmachine/';

%% 

fileName = 'BedMachineGreenland-2021-04-20.nc';

ncdisp([dataDir fileName]);

%% 

x = double(ncread([dataDir fileName],'x'));
y = double(ncread([dataDir fileName],'y'));

%mask (0 = ocean, 1 = ice-free land, 2 = grounded ice, 3 = floating ice, 4 = non-Greenland land)
mask = ncread([dataDir fileName],'mask'); 

%bed = ncread([dataDir fileName],'bed');
%thickness = ncread([dataDir fileName],'thickness'); 
%surface = ncread([dataDir fileName],'surface'); 
%firn = ncread([dataDir fileName],'firn'); 
error = ncread([dataDir fileName],'errbed');

%% 

% geoid = 'eigen-6c4'
% grid_mapping_name = 'polar_stereographic'
% latitude_of_projection_origin = 90
% standard_parallel = 70
% straight_vertical_longitude_from_pole = -45
% semi_major_axis = 6378273
% inverse_flattening = 298.2794
% false_easting = 0
% false_northing = 0

earthRadius = 6378273;
eccentricity = 0.08181919;

lat_true = 70;
lon_posy = -45;

[x y] = meshgrid(x,y);

[lat lon]= polarstereo_inv(x,y,earthRadius,eccentricity,lat_true,lon_posy);

lon = lon';
lat = lat';

%%

figure

pcolorcen(lon,lat,mask);

colorbar

%%
