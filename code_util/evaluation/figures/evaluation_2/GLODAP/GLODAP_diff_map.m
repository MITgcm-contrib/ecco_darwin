function [] = GLODAP_diff_map(gridType,startTime,endTime,minDepth,maxDepth,caption,savePPT)

gridDir = '../../data/grid/';

dataDir1 = '../../data/model_data_pairs/';
dataDir2 = '../../data/indices/';
dataDir3 = '../../data/coastline/';

%%

if savePPT
    
    slideNum = exportToPPTX('addslide');
    exportToPPTX('addtext',caption,'HorizontalAlignment','center','VerticalAlignment','middle');
    
end

%%

nz = 50;

DRF = readbin([gridDir 'DRF.data'],[nz],1,'real*4');

%%

load([dataDir1 'GLODAP_model_data_pairs.mat']);
load([dataDir2 'ECCO_darwin_combined_biome_indices_LLC_' num2str(gridType) '.mat']);
load([dataDir3 'coastline_map.mat']);

%%

bgcVar{1} = 'THETA';
bgcVar{2} = 'SAL';
bgcVar{3} = 'DIC';
bgcVar{4} = 'ALK';
bgcVar{5} = 'PH';
bgcVar{6} = 'NO3';
bgcVar{7} = 'PO4';
bgcVar{8} = 'SiO2';
bgcVar{9} = 'O2';

bgcName{1} = 'Pot. Temp.';
bgcName{2} = 'Salinity';
bgcName{3} = 'DIC';
bgcName{4} = 'Alkalinity';
bgcName{5}  = 'pH';
bgcName{6} = 'NO_3';
bgcName{7} = 'PO_4';
bgcName{8} = 'SiO_2';
bgcName{9} = 'O_2';

bgcUnits{1} = '(\circC)';
bgcUnits{2} = '(psu)';
bgcUnits{3} = '(mmol m^-^3)';
bgcUnits{4} = '(meq. m^-^3)';
bgcUnits{5}  = '';
bgcUnits{6} = '(mmol m^-^3)';
bgcUnits{7} = '(mmol m^-^3)';
bgcUnits{8} = '(mmol m^-^3)';
bgcUnits{9} = '(mmol m^-^3)';

%%

fs = 26;
lw = 2;
mw = 30;

colors = cmocean('balance',1000);

for i = 1:length(bgcVar)
    
    bVar = bgcVar{i};
    bName = bgcName{i};
    bUnits = bgcUnits{i};
    
    %%
    
    eval(['obsAll = vertcat(observations.' bVar '{:});']);
    eval(['modelAll = vertcat(model.' bVar '{:});']);
    
    obs.time = vertcat(obsAll(:).time);
    obs.lon = vertcat(obsAll(:).lon);
    obs.lat = vertcat(obsAll(:).lat);
    obs.depth = vertcat(obsAll(:).depth);
    obs.gridIndex = vertcat(obsAll(:).gridIndex);
    obs.depthIndex = vertcat(obsAll(:).depthIndex);
    
    obs.data = vertcat(obsAll(:).data);
    model.data = vertcat(modelAll(:).data);
    
    [years, ~, ~, ~, ~, ~] = datevec(obs.time);
    
    %%
    
    xxi = find(obs.time >= startTime & obs.time <= endTime & ...
        obs.depth >= minDepth-eps & obs.depth <= maxDepth+eps);
    
    time = obs.time(xxi);
    lon = obs.lon(xxi);
    
    %Wunsch centering
    lon(lon < 0) = lon(lon < 0) + 360;
    lon(lon>113) = lon(lon>113)-360;
    
    lat = obs.lat(xxi);
    depth = obs.depth(xxi);
    gridIndex = obs.gridIndex(xxi);
    depthIndex = obs.depthIndex(xxi);
    obsThickness = DRF(depthIndex);
    
    obsData = obs.data(xxi);
    modelData = model.data(xxi);
    
    obsDif = modelData - obsData;
    
    clear xxi
    
    %%
    %depth average
    
    uGridIndex = unique(gridIndex);
    
    for j = 1:length(uGridIndex)
        
        gi = find(gridIndex == uGridIndex(j));
        
        mLon(j) = nanmean(lon(gi));
        mLat(j) = nanmean(lat(gi));
        mGridIndex(j) = uGridIndex(j);
        
        meanObs(j) = nansum(obsData(gi) .* obsThickness(gi)) ./ nansum(obsThickness(gi));
        meanModel(j) = nansum(modelData(gi) .* obsThickness(gi)) ./ nansum(obsThickness(gi));
        
        meanDif(j) = meanModel(j) - meanObs(j);
        
        clear gi
        
    end
    
    %%
    
    sigma = 2;
    
    temp = abs(obsDif);
    
    cMin = -(nanmean(temp) + (sigma .* nanstd(temp)));
    cMax = nanmean(temp) + (sigma .* nanstd(temp));
    
    hFig1 = figure(1);
    set(hFig1,'units','normalized','outerposition',[0 0 1 1]);
    set(gcf,'color',[1 1 1]);
    
    hold on
    
    contour(coastline.lon,coastline.lat,coastline.map,'Color','k','LineWidth',1.5);
    
    scatter(mLon,mLat,mw,meanDif,'filled','MarkerEdgeColor','none');
    
    colormap(colors);
    
    caxis([cMin cMax]);
    
    hcb = colorbar;
    set(get(hcb,'ylabel'),'String',['Model-data Difference ' bgcUnits{i}]);
    
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
    
    title([bgcName{i} ' Model-data Difference, ' num2str(minDepth) ' to ' num2str(maxDepth) '-m Depth']);
    
    if savePPT
        
        slideNum = exportToPPTX('addslide');
        exportToPPTX('addpicture',hFig1,'maxscale');
        
    end
    
    close all
    
    clear obsAll obs.time obs.lon obs.lat obs.depth obs.gridIndex obs.depthIndex obs.data modelAll model.data ...
        years months obs.years obs.months meanObs stdObs lon lat depth gridIndex depthIndex obsThickness ...
        obsData modelData obsDif uGridIndex mLon mLat mGridIndex meanObs meanModel meanDif
    
end

close all;

end

%%
