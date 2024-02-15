function [] = GLODAP_seasonal_climatology(gridType,startTime,endTime,minDepth,maxDepth,caption,savePPT)

gridDir = '../../data/grid/';

dataDir1 = '../../data/model_data_pairs/';
dataDir2 = '../../data/indices/';

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

biomeString = {'NH';'NS';'EQ';'SS';'SH'; ...
    'NP1'; ...
    'NP2'; ...
    'NP3'; ...
    'NP4'; ...
    'WPE'; ...
    'EPE'; ...
    'SP'; ...
    'NA1'; ...
    'NA2'; ...
    'NA3'; ...
    'NA4'; ...
    'AE'; ...
    'SA'; ...
    'IO'; ...
    'SO1'; ...
    'SO2'; ...
    'SO3'; ...
    'All Biomes';
    'Global';
    };

biomeNames{1} = 'fay18';
biomeNames{2} = 'fay19';
biomeNames{3} = 'fay20';
biomeNames{4} = 'fay21';
biomeNames{5} = 'fay22';

c = 1;

for i = 6:22
    
    biomeNames{i} = ['fay' num2str(c)];
    
    c = c + 1;
    
end

biomeNames{i+1} = 'allBiomes';
biomeNames{i+2} = 'global';

xi = [];

%%

fs = 16;
lw = 2;
mw = 5;

colors = cbrewer('qual','Paired',6);
bgColor = [1  1 1];

marg = 0.05;

%%

for i = 1:length(bgcVar)
    
    bVar = bgcVar{i};
    bName = bgcName{i};
    bUnits = bgcUnits{i};
    
    hFig1 = figure(1);
    set(hFig1,'units','normalized','outerposition',[0 0 1 1]);
    set(gcf,'color',[1 1 1]);
    
    for i = 1:length(biomeNames)
        
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
        
        %%
        
        eval(['xi = [index.' biomeNames{i} '];']);
        
        xi = unique(xi);
        xi = ismember(obs.gridIndex,xi);
        
        obs.time = obs.time(xi);
        
        [~, months, ~, ~, ~, ~] = datevec(obs.time);
        
        obs.months = months;
        obs.gridIndex = obs.gridIndex(xi);
        obs.depthIndex = obs.depthIndex(xi);
        obs.depth = obs.depth(xi);
        
        obs.data = obs.data(xi);
        model.data = model.data(xi);
        
        clear xi
        
        %%
        
        xxi = find(obs.time >= startTime & obs.time <= endTime & ...
            obs.depth >= minDepth-eps & obs.depth <= maxDepth+eps);
        
        time = obs.time(xxi);
        lon = obs.lon(xxi);
        lat = obs.lat(xxi);
        depth = obs.depth(xxi);
        gridIndex = obs.gridIndex(xxi);
        depthIndex = obs.depthIndex(xxi);
        obsThickness = DRF(depthIndex);
        obsMonths = obs.months(xxi);
        
        obsData = obs.data(xxi);
        modelData = model.data(xxi);
        
        clear xxi
        
        %%
        
        for j = 1:12
            
            im = find(obsMonths == j);
            
            seasTime(j) = j;
            
            if ~isempty(im)
                
                seasObs(j) = nansum(obsData(im) .* obsThickness(im)) ./ nansum(obsThickness(im));
                seasModel(j) = nansum(modelData(im) .* obsThickness(im)) ./ nansum(obsThickness(im));
                
            else
                
                seasObs(j) = nan;
                seasModel(j) = nan;
                
            end
            
            clear im
            
        end
        
        %%
        
        eval(['subplot_tight(6,5,' num2str(i) ',marg);']);
        
        hold on
        
        set(gca,'color',bgColor);
        
        p1 = scatter(obsMonths,obsData,mw,'MarkerFaceColor',colors(2,:),'MarkerEdgeColor',colors(2,:));
        p2 = scatter(obsMonths,modelData,mw,'MarkerFaceColor',colors(6,:)','MarkerEdgeColor',colors(6,:));
        
        plot(seasTime,seasObs,'Color',colors(2,:),'LineWidth',lw);
        plot(seasTime,seasModel,'Color',colors(6,:),'LineWidth',lw);
        
        if i == 1
            
            ylabel([bName ' ' bUnits]);
            
            legend([p1 p2],{'Obs.','Model'},'Location','Best');
            
        end
        
        if i == 21
            
            xlabel(['Month']);
            
        end
        
        set(gca,'xtick',[1:12]);
        
        xlim([1 12]);
        
        axis tight
        
        box on
        grid on
        set(gca,'GridLineStyle','--');
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs-2);
        
        title(['\bf' biomeString{i} '\rm'],'FontWeight','Normal','FontSize',fs-2);
        
        drawnow
        
        clear obsAll obs.time obs.lon obs.lat obs.depth obs.gridIndex obs.depthIndex obs.data modelAll model.data ...
            years months obs.years obs.months meanObs stdObs lon lat depth gridIndex depthIndex obsThickness ...
            obsData modelData obsDif uGridIndex mLon mLat mGridIndex meanObs meanModel meanDif ...
            seasTime seasObs seasModel
        
    end
    
    if savePPT
        
        slideNum = exportToPPTX('addslide');
        exportToPPTX('addpicture',hFig1,'maxscale');
        
    end
    
    close all;
    
end

end

%%
