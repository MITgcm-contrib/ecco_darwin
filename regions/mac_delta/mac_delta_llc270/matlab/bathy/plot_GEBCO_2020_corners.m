clear
close all;

plotGrid = 1;
savePlot = 0;

gridDir = '/Users/carrolld/Documents/research/mackenzie/grid/LLC_270/';

dataDir1 = '/Users/carrolld/Documents/research/LLC_540/raw_data/gebco_2020/';
dataDir2 = '/Users/carrolld/Documents/research/mackenzie/mat/corners/';
dataDir3 = '/Users/carrolld/Documents/research/mackenzie/mat/bathy/';

figureDir = '/Users/carrolld/Documents/research/mackenzie/figures/bathy/';  

bathyFileName = {'Depth.data','GEBCO_mackenzie_bathy_dustin_median.bin'};

bathyTitle = {'Hong','Dustin Median'};

%%

nx = 46;
ny = 68;

kx = 1;
prec = 'real*4';

XG = readbin([gridDir 'XC.data'],[nx ny]);
YG = readbin([gridDir 'YC.data'],[nx ny]);

fileName = 'GEBCO_2020.nc';

lon = ncread([dataDir1 fileName],'lon');
lat = ncread([dataDir1 fileName],'lat');

elevation = ncread([dataDir1 fileName],'elevation');

load([dataDir2 'cell_corners.mat']);

xPoly = [XGsw(:) XGse(:) XGne(:) XGnw(:) XGsw(:)];
yPoly = [YGsw(:) YGse(:) YGne(:) YGnw(:) YGsw(:)];

%%

fs = 14;
lw = 2;
bgColor = [0 0 0];
gridColor = [1 1 1];

if plotGrid
    
    suffix = 'grid';
    
else
    
    suffix = 'no_grid';
    
end

%%

close all

lonBounds = [min(min(xPoly)) max(max(xPoly))];
latBounds = [min(min(yPoly)) max(max(yPoly))];

maxDepth = 1;

%% 

cc = -maxDepth:0;
colors = flipud(cbrewer('div','Spectral',(maxDepth + 1)));

dx = 0;
dy = 0;

ix1 = find(lon > lonBounds(1)- dx & lon < lonBounds(2) + dx);
iy1 = find(lat > latBounds(1) - dy & lat < latBounds(2) + dy);

ic = find(min(xPoly') >= lonBounds(1) - dx & max(xPoly') <= lonBounds(2)+ dx & ...
    min(yPoly') >= latBounds(1) - dx & max(yPoly') <= latBounds(2) + dy);

GEBCOBathy = elevation(ix1,iy1);

imAlpha = ones(size(GEBCOBathy'));
imAlpha(GEBCOBathy' >= 0) = 0;

bathy1 = -readbin([dataDir3 bathyFileName{1}],[nx ny]);
bathy2 = -readbin([dataDir3 bathyFileName{2}],[nx ny]);

%bathy1(bathy1 == 0) = nan;
%bathy2(bathy2 == 0) = nan;

%%

hFig1 = figure(1);
set(hFig1,'units','normalized','outerposition',[0 0 1 1]);
set(gcf,'color',[1 1 1]);

cc1 = subplot(1,3,1);

hold on

imagesc(lon(ix1),lat(iy1),GEBCOBathy','AlphaData',imAlpha);

set(gca,'color',bgColor);

if plotGrid
    
    for i =1:length(ic)
        
        line(xPoly(ic(i),:),yPoly(ic(i),:),'Color',gridColor);
        
    end
    
end

caxis([-maxDepth 0]);

colormap(cc1,colors);

%hcb = colorbar('horizontal');
%ylabel(hcb,'Depth (m)');

xlim([lonBounds(1) lonBounds(2)]);
ylim([latBounds(1) latBounds(2)]);

xlabel('Longitude');
ylabel('Latitude');

set(gca,'yDir','normal');
set(gca,'FontSize',fs);
set(gca,'LineWidth',lw);

title('GEBCO');

cc2 = subplot(132);

hold on

set(gca,'color',bgColor);

pcolorcen(XG,YG,bathy1);
shading flat

if plotGrid
    
    for i =1:length(ic)
        
        line(xPoly(ic(i),:),yPoly(ic(i),:),'Color',gridColor);
        
    end
    
end

caxis([-maxDepth 0]);

colormap(cc2,colors);

%hcb = colorbar('horizontal');
%ylabel(hcb,'Depth (m)');

xlim([lonBounds(1) lonBounds(2)]);
ylim([latBounds(1) latBounds(2)]);

xlabel('Longitude');
ylabel('Latitude');

set(gca,'yDir','normal');
set(gca,'FontSize',fs);
set(gca,'LineWidth',lw);

title(bathyTitle{1});

cc3 = subplot(133);
    
hold on

set(gca,'color',bgColor);

pcolorcen(XG,YG,bathy2);
shading flat

if plotGrid
    
    for i =1:length(ic)
        
        line(xPoly(ic(i),:),yPoly(ic(i),:),'Color',gridColor);
        
    end
    
end

caxis([-maxDepth 0]);

colormap(cc3,colors);

%hcb = colorbar('horizontal');
%ylabel(hcb,'Depth (m)');

xlim([lonBounds(1) lonBounds(2)]);
ylim([latBounds(1) latBounds(2)]);

xlabel('Longitude');
ylabel('Latitude');

set(gca,'yDir','normal');
set(gca,'FontSize',fs);
set(gca,'LineWidth',lw);

title(bathyTitle{2});

%% 
