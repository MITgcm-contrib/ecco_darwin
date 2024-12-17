clear
close all;

saveDir = '/Users/carrolld/Documents/research/bathy/mat/bedmachine/';

tic

%%

scale = 8;

dx = 1/(48*scale);
dy = 1/(48*scale);

lon = -180:dy:180;
lat = 60:dx:84.5;

lat = fliplr(lat);

%%

[xx yy] = meshgrid(lon,lat);

%%

%base, meters base of the ice sheet (bottom of ice shelves, but same as bed over grounded ice)
base = bedmachine_interp('base',yy,xx,'greenland','method','nearest');
mask = bedmachine_interp('mask',yy,xx,'greenland','method','nearest');

%%

%0 = ocean
%1 = ice-free land
%2 = grounded ice
%3 = floating ice
%4 = non-Greenland land or Vostok

%%

base(base >= 0) = 0;

%%

figure

subplot(121);

mypcolor(lon,lat,base);

colorbar

subplot(122);

mypcolor(lon,lat,mask);

colorbar

%%

save([saveDir 'bedmachine_greenland_lon_lat_interp.mat'],'lon','lat','base','mask','-v7.3');

%%

toc

%%
