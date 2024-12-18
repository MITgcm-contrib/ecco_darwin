clear
close all;

dataDir = '/Users/carrolld/Documents/research/bathy/raw_data/GEBCO_2024/';
saveDir = '/Users/carrolld/Documents/research/bathy/mat/bathy/GEBCO_2024_raw/';

%% 

fileName = 'GEBCO_2024_sub_ice_topo.nc';

ncdisp([dataDir fileName]);

%% 

lon = single(ncread([dataDir fileName],'lon'));
lat = single(ncread([dataDir fileName],'lat'));

%% 

elevation = single(ncread([dataDir fileName],'elevation'));

temp = elevation';

close all

hFig1 = figure(1);
set(hFig1,'units','normalized','outerposition',[0 0 1 1]);
set(gcf,'color',[1 1 1]);

fs = 24;
lw = 2;

colors = cmocean('deep',1000);

imagesc(lon,lat,temp);

ax = gca;
ax.YDir= 'normal';
 
caxis([-5000 0]);

colormap(colors);

hcb = colorbar;

set(get(hcb,'ylabel'),'String',{'Depth (m)'});

xlabel('Longitude (deg)');
ylabel('Latitude (deg)');

box on
grid off
set(gca,'GridLineStyle','-.');
set(gca,'LineWidth',lw);
set(gca,'FontSize',fs);

title('GEBCO 2024');

%%
