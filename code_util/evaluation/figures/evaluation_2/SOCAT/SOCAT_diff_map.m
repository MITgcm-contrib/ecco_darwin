function [] = SOCAT_diff_map(gridType,startTime,endTime,caption,savePPT)

dataDir1 = '../../data/model_data_pairs/';
dataDir2 = '../../data/indices/';
dataDir3 = '../../data/coastline/';

%%

if savePPT
    
        slideNum = exportToPPTX('addslide');
        exportToPPTX('addtext',caption,'HorizontalAlignment','center','VerticalAlignment','middle');

end

%%

load([dataDir1 'SOCAT_model_data_pairs.mat']);
load([dataDir2 'ECCO_darwin_combined_biome_indices_LLC_' num2str(gridType) '.mat']);
load([dataDir3 'coastline_map.mat']);

%%

fs = 26;
lw = 2;
mw = 15;

colors = cmocean('balance',1000);

%%

obsAll = vertcat(observations.fCO2{:});
modelAll = vertcat(model.fCO2{:});

%%

obs.time = vertcat(obsAll(:).time);

obs.lon = vertcat(obsAll(:).lon);
obs.lat = vertcat(obsAll(:).lat);
obs.gridIndex = vertcat(obsAll(:).gridIndex);

obs.data = vertcat(obsAll(:).data2);
model.data = vertcat(modelAll(:).data);

[years, ~, ~, ~, ~, ~] = datevec(obs.time);
[~, months, ~, ~, ~, ~] = datevec(obs.time);

obs.years = years;
obs.months = months;

%%

xxi = find(obs.time >= startTime & obs.time <= endTime);

lon = obs.lon(xxi);

%Wunsch centering
lon(lon < 0) = lon(lon < 0) + 360;
lon(lon>113) = lon(lon>113)-360;

lat = obs.lat(xxi);
gridIndex = obs.gridIndex(xxi);

obsData = obs.data(xxi);
modelData = model.data(xxi);

obsDif = modelData - obsData;

clear xxi

%%

uGridIndex = unique(gridIndex);

for j = 1:length(uGridIndex)
    
    gi = find(gridIndex == uGridIndex(j));
    
    mLon(j) = nanmean(lon(gi));
    mLat(j) = nanmean(lat(gi));
    mGridIndex(j) = uGridIndex(j);
    
    meanObs(j) = nanmean(obsData(gi));
    meanModel(j) = nanmean(modelData(gi));
    
    meanDif(j) = meanModel(j) - meanObs(j);
    
    clear gi
    
end

%%

sigma = 2;

temp = abs(obsDif);

cMin = -(nanmean(temp) + (sigma .* nanstd(temp)));
cMax = nanmean(temp) + (sigma .* nanstd(temp));

%%

hFig1 = figure(1);
set(hFig1,'units','normalized','outerposition',[0 0 1 1]);
set(gcf,'color',[1 1 1]);

hold on

contour(coastline.lon,coastline.lat,coastline.map,'Color','k','LineWidth',1.5);
    
scatter(mLon,mLat,mw,meanDif,'filled','MarkerEdgeColor','none');

colormap(colors);

caxis([cMin cMax]);

hcb = colorbar;
set(get(hcb,'ylabel'),'String','Model-data Difference (\muAtm)');

xlim([-247 113]);
ylim([-90 90]);

set(gca,'xtick',[-240:60:60]);
set(gca,'ytick',[-80:20:80]);

set(gca,'xticklabel',{'120\circE','180\circ','120\circW','60\circW','0\circ','60\circE'});
set(gca,'yticklabel',{'80\circS','60\circS','40\circS','20\circS','0\circ', ...
    '20\circN','40\circN','60\circN','80\circN'});

box on
grid off
set(gca,'GridLineStyle','--');
set(gca,'LineWidth',lw);
set(gca,'FontSize',fs);

drawnow

title('Surface-ocean fCO_2 Model-data Difference');

if savePPT
    
    slideNum = exportToPPTX('addslide');
    exportToPPTX('addpicture',hFig1,'maxscale');
   
end

clear obsAll obs.time obs.lon obs.lat obs.depth obs.gridIndex obs.depthIndex obs.data modelAll model.data ...
    years months obs.years obs.months meanObs stdObs lon lat depth gridIndex depthIndex ...
    obsData modelData obsDif uGridIndex mLon mLat mGridIndex meanObs meanModel meanDif

close all;

end

%%
