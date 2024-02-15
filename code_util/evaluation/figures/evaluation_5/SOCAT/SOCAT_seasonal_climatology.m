function [] = SOCAT_seasonal_climatology(gridType,startTime,endTime,caption,savePPT)

dataDir1 = '../../data/model_data_pairs/';
dataDir2 = '../../data/indices/';

%%

load([dataDir1 'SOCAT_model_data_pairs.mat']);
load([dataDir2 'ECCO_darwin_combined_biome_indices_LLC_' num2str(gridType) '.mat']);

%%

if savePPT
    
    slideNum = exportToPPTX('addslide');
    exportToPPTX('addtext',caption,'HorizontalAlignment','center','VerticalAlignment','middle');
    
end

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
    'All Biomes'; ...
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

%%

fs = 16;
lw = 2;
mw = 10;

colors = cbrewer('qual','Paired',6);
bgColor = [1  1 1];

marg = 0.05;

%%

close all

hFig1 = figure(1);
set(hFig1,'units','normalized','outerposition',[0 0 1 1]);
set(gcf,'color',[1 1 1]);

for i = 1:length(biomeNames)
    
    obsAll = vertcat(observations.fCO2{:});
    modelAll = vertcat(model.fCO2{:});
    
    obs.time = vertcat(obsAll(:).time);
    obs.lon = vertcat(obsAll(:).lon);
    obs.lat = vertcat(obsAll(:).lat);
    obs.gridIndex = vertcat(obsAll(:).gridIndex);
    
    obs.data = vertcat(obsAll(:).data2);
    model.data = vertcat(modelAll(:).data);
    
    %%
    
    eval(['xi = [index.' biomeNames{i} '];']);
    
    xi = unique(xi);
    xi = ismember(obs.gridIndex,xi);
    
    obs.time = obs.time(xi);
    obs.lon = obs.lon(xi);
    obs.lat = obs.lat(xi);
    obs.gridIndex = obs.gridIndex(xi);
    
    obs.data = obs.data(xi);
    model.data = model.data(xi);
    
    clear xi
    
    %%
    
    [~, months, ~, ~, ~, ~] = datevec(obs.time);
    
    obs.months = months;
    
    xxi = find(obs.time >= startTime & obs.time <= endTime);
    
    obsMonths = obs.months(xxi);
    obsData = obs.data(xxi);
    modelData = model.data(xxi);
    
    clear xxi
    
    %%
    
    for j = 1:12
        
        im = find(obsMonths == j);
        
        seasTime(j) = j;
        
        if ~isempty(im)
            
            seasObs(j) = nanmean(obsData(im));
            seasModel(j) = nanmean(modelData(im));
            
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
    p2 = scatter(obsMonths,modelData,mw,'MarkerFaceColor',colors(6,:),'MarkerEdgeColor',colors(6,:));
    
    plot(seasTime,seasObs,'Color',colors(2,:),'LineWidth',lw);
    plot(seasTime,seasModel,'Color',colors(6,:),'LineWidth',lw);
    
    if i == 1
        
        ylabel('fCO_2 (\muAtm)');
        
        legend([p1 p2],{'Obs.','Model'},'Location','Best');
        
    end
    
    if i == 21
        
        xlabel('Month');
        
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
        years months obs.years obs.months meanObs stdObs lon lat depth gridIndex depthIndex ...
        obsData modelData obsDif uGridIndex mLon mLat mGridIndex meanObs meanModel meanDif ...
        seasTime seasObs seasModel
    
end

if savePPT
    
    slideNum = exportToPPTX('addslide');
    exportToPPTX('addpicture',hFig1,'maxscale');
    
end

close all;

end

%%
