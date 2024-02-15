function [] = BGC_Argo_time_series(gridType,startTime,endTime,minDepth,maxDepth,caption,savePPT)

gridDir = '../../data/grid/';

dataDir1 = '../../data/model_data_pairs/';
dataDir2 = '../../data/indices/';

%%

if savePPT
    
    slideNum = exportToPPTX('addslide');
    exportToPPTX('addtext',caption,'HorizontalAlignment','center','VerticalAlignment','middle');
    
end

%%

load([dataDir1 'BGC-Argo_model_data_pairs.mat']);
load([dataDir2 'ECCO_darwin_combined_biome_indices_LLC_' num2str(gridType) '.mat']);

%%

nz = 50;

DRF = readbin([gridDir 'DRF.data'],[nz],1,'real*4');

%%

bgcVar{1} = 'PH';
bgcVar{2} = 'NO3';
bgcVar{3} = 'O2';

bgcName{1}  = 'pH';
bgcName{2} = 'NO_3';
bgcName{3} = 'O_2';

bgcUnits{1}  = '';
bgcUnits{2} = '(mmol m^-^3)';
bgcUnits{3} = '(mmol m^-^3)';

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

xtick = [datenum(1995,1,1,0,0,0) datenum(2000,1,1,0,0,0) ...
    datenum(2005,1,1,0,0,0) datenum(2010,1,1,0,0,0) datenum(2015,1,1,0,0,0) ...
    datenum(2020,1,1,0,0,0)];

%%

for i = 1:length(bgcVar)
    
    bVar = bgcVar{i};
    bName = bgcName{i};
    bUnits = bgcUnits{i};
    
    hFig1 = figure(1);
    set(hFig1,'units','normalized','outerposition',[0 0 1 1]);
    set(gcf,'color',[1 1 1]);
    
    for j = 1:length(biomeNames)
        
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
        
        eval(['xi = [index.' biomeNames{j} '];']);
        
        xi = unique(xi);
        xi = ismember(obs.gridIndex,xi);
        
        obs.time = obs.time(xi);
        
        [years, ~, ~, ~, ~, ~] = datevec(obs.time);
        
        obs.years = years;
        obs.gridIndex = obs.gridIndex(xi);
        obs.depthIndex = obs.depthIndex(xi);
        obs.depth = obs.depth(xi);
        
        obs.data = obs.data(xi);
        model.data = model.data(xi);
        
        clear xi
        
        %%
        
        xxi = find(obs.time >= startTime & obs.time <= endTime & ...
            obs.depth >= minDepth-eps & obs.depth <= maxDepth+eps);
        
        obsTime = obs.time(xxi);
        obsYears = obs.years(xxi);
        depthIndex = obs.depthIndex(xxi);
        obsThickness = DRF(depthIndex);
        
        uniqueYears = unique(obsYears);
        
        obsData = obs.data(xxi);
        modelData = model.data(xxi);
        
        %%
        
        clear xxi
        
        [obsTime ti] = sort(obsTime);
        
        obsYears = obsYears(ti);
        obsThickness = obsThickness(ti);
        obsData = obsData(ti);
        modelData = modelData(ti);
        
        %%
        
        for k = 1:length(uniqueYears)
            
            im = find(obsYears == uniqueYears(k));
            
            annualTime(k) = datenum(uniqueYears(k),6,15,0,0,0);
            
            annualObs(k) = nansum(obsData(im) .* obsThickness(im)) ./ nansum(obsThickness(im));
            annualModel(k) = nansum(modelData(im) .* obsThickness(im)) ./ nansum(obsThickness(im));
            
        end
        
        %%
        
        eval(['subplot_tight(6,5,' num2str(j) ',marg);']);
        
        hold on
        
        set(gca,'color',bgColor);
        
        if ~isempty(k)
            
            p1 = scatter(obsTime,obsData,mw,'MarkerFaceColor',colors(1,:),'MarkerEdgeColor',colors(1,:));
            p2 = scatter(obsTime,modelData,mw,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5]);
            
            plot(annualTime,annualObs,'Color',colors(2,:),'LineWidth',lw+2);
            plot(annualTime,annualModel,'Color',colors(6,:),'LineWidth',lw+2);
            
        end
        
        if j == 1
            
            ylabel([bName ' ' bUnits]);
            
            legend([p1 p2],{'Obs','Model'},'Location','Best');
            
        end
        
        if j == 21
            
            xlabel(['Year']);
            
        end
        
        set(gca,'xtick',xtick);
        datetick('x',11,'keepticks');
        
        axis tight
        
        xlim([startTime endTime]);
        
        box on
        grid on
        set(gca,'GridLineStyle','--');
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs-2);
        
        title(['\bf' biomeString{j} '\rm'],'FontWeight','Normal','FontSize',fs-2);
        
        drawnow
        
        clear obsAll obs.time obs.lon obs.lat obs.depth obs.gridIndex obs.depthIndex obs.data modelAll model.data ...
            years months obs.years obs.months meanObs stdObs lon lat depth gridIndex depthIndex obsThickness ...
            obsTime obsYears obsData modelData obsDif uGridIndex mLon mLat mGridIndex meanObs meanModel meanDif ...
            annualTime annualObs annualModel
        
    end
    
    if savePPT
        
        slideNum = exportToPPTX('addslide');
        exportToPPTX('addpicture',hFig1,'maxscale');
        
    end
    
    close all;
    
end

end

%%
