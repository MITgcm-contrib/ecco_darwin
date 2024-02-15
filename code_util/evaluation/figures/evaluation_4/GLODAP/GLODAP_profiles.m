function [] = GLODAP_profiles(gridType,startTime,endTime,caption,savePPT);

dataDir1 = '../../data/model_data_pairs/';
dataDir2 = '../../data/indices/';

%%

if savePPT
    
    slideNum = exportToPPTX('addslide');
    exportToPPTX('addtext',caption,'HorizontalAlignment','center','VerticalAlignment','middle');
    
end

%%

load([dataDir1 'GLODAP_model_data_pairs.mat']);
load([dataDir2 'ECCO_darwin_combined_biome_indices_LLC_' num2str(gridType) '.mat']);

%%

bgcVar{1} = 'DIC';
bgcVar{2} = 'ALK';
bgcVar{3} = 'PH';
bgcVar{4} = 'NO3';
bgcVar{5} = 'PO4';
bgcVar{6} = 'SiO2';
bgcVar{7} = 'O2';

bgcName{1} = 'DIC';
bgcName{2} = 'Alkalinity';
bgcName{3}  = 'pH';
bgcName{4} = 'NO_3';
bgcName{5} = 'PO_4';
bgcName{6} = 'SiO_2';
bgcName{7} = 'O_2';

bgcUnits{1} = '(mmol m^-^3)';
bgcUnits{2} = '(meq. m^-^3)';
bgcUnits{3}  = '';
bgcUnits{4} = '(mmol m^-^3)';
bgcUnits{5} = '(mmol m^-^3)';
bgcUnits{6} = '(mmol m^-^3)';
bgcUnits{7} = '(mmol m^-^3)';

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

marg = 0.05;

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
        obs.gridIndex = vertcat(obsAll(:).gridIndex);
        obs.depth = vertcat(obsAll(:).depth);
        
        obs.data = vertcat(obsAll(:).data);
        model.data = vertcat(modelAll(:).data);
        
        %%
        
        eval(['xi = [index.' biomeNames{i} '];']);
        
        xi = unique(xi);
        xi = ismember(obs.gridIndex,xi);
        
        obs.time = obs.time(xi);
        obs.depth = obs.depth(xi);
        
        obs.data = obs.data(xi);
        model.data = model.data(xi);
        
        clear xi
        
        %%
        
        xxi = find(obs.time >= startTime & obs.time <= endTime);
        
        lon = obs.lon(xxi);
        lat = obs.lat(xxi);
        gridIndex = obs.gridIndex(xxi);
        depth = obs.depth(xxi);
        
        obsData = obs.data(xxi);
        modelData = model.data(xxi);
        
        clear xxi
        
        obsDif = modelData - obsData;
        
        %%
        
        eval(['subplot_tight(6,5,' num2str(i) ',marg);']);
        
        hold on
        
        p1 = scatter(obsData,-depth,mw,'MarkerFaceColor',colors(2,:),'MarkerEdgeColor',colors(2,:));
        p2 = scatter(modelData,-depth,mw,'MarkerFaceColor',colors(6,:),'MarkerEdgeColor',colors(6,:));
        
        if i == 1
            
            ylabel('Depth (m)');
            legend([p1 p2],{'Obs.','Model'},'Location','SouthWest');
            
        end
        
        if i == 21
            
            xlabel([bName ' ' bUnits]);
            
        end
        
        axis tight
        
        ylim([-6000 500]);
        
        box on
        grid on
        set(gca,'GridLineStyle','--');
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs-2);
        
        title(['\bf' biomeString{i} '\rm'],'FontWeight','Normal','FontSize',fs-2);
        
        drawnow
        
        clear obsAll obs.time obs.lon obs.lat obs.depth obs.gridIndex obs.depthIndex obs.data modelAll model.data ...
            years months obs.years obs.months meanObs stdObs lon lat depth gridIndex depthIndex obsThickness ...
            obsData modelData obsDif uGridIndex mLon mLat mGridIndex meanObs meanModel meanDif
        
        
    end
    
    if savePPT
        
        slideNum = exportToPPTX('addslide');
        exportToPPTX('addpicture',hFig1,'maxscale');
        
    end
    
    close all;
    
end

end

%%
