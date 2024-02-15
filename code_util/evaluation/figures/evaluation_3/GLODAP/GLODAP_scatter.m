function [] = GLODAP_scatter(gridType,startTime,endTime,numBins,caption,savePPT)

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

colors = flipud(cbrewer('div','RdYlBu',1000));

bgColor = [0.65 0.65 0.65];

marg = 0.05;

c = 1;

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
        
        time = obs.time(xxi);
        lon = obs.lon(xxi);
        lat = obs.lat(xxi);
        gridIndex = obs.gridIndex(xxi);
        depth = obs.depth(xxi);
        
        obsData = obs.data(xxi);
        modelData = model.data(xxi);
        
        obsDif = modelData - obsData;
        
        %%
        
        if ~isempty(obsData)
            
            mdl = fitlm(obsData,modelData)
            
            R2 = mdl.Rsquared.Ordinary;
            
            xMin = min(obsData);
            xMax = max(obsData);
            
            yMin = min(modelData);
            yMax = max(modelData);
            
            bin = max(xMax-xMin,yMax-yMin) ./ numBins;
            
            vec = min([xMin yMin]):bin:max([xMax yMax]);
            
            [xx yy] = meshgrid(vec);
            
            for j = 1:length(vec)
                
                for k = 1:length(vec)
                    
                    id = find(obs.data >= (vec(j) - bin/2) & obs.data < (vec(j) + bin/2) ...
                        & model.data >= (vec(k) - bin/2) & model.data < (vec(k) + bin/2));
                    
                    if ~isempty(id)
                        
                        dens(j,k) = length(id);
                        
                    else
                        
                        dens(j,k) = nan;
                        
                    end
                    
                end
                
            end
            
        else
            
            R2 = 0;
            
            xx = ones(2,2) .* nan;
            yy = xx .* nan;
            
            dens = xx .* nan;
            
        end
        
        dens = log10(dens);
        dens = dens ./ nanmax(dens(:));
        
        %%
        
        eval(['subplot_tight(6,5,' num2str(i) ',marg);']);
        
        hold on
        
        set(gca,'color',bgColor);
        
        pcolorcen(xx,yy,dens');
        
        colormap(colors);
        
        caxis([0 1]);
        
        line([0 10^8],[0 10^8],'Color','k','LineWidth',1.5,'LineStyle','--');
        
        xlim([min([xMin yMin]) max([xMax yMax])]);
        ylim([min([xMin yMin]) max([xMax yMax])]);
        
        if i == 1
            
            ylabel(['\bfModel\rm ' bName ' ' bUnits]);
            
        end
        
        if i == 21
            
            xlabel(['\bfData\rm ' bName ' ' bUnits]);
            
        end
        
        if i == 24
            
            hcb = colorbar;
            
            set(get(hcb,'ylabel'),'String',{'log_1_0(Density)'});
            pos = get(hcb,'Position');
            
            pos = pos + [0.0325 0 0 0];
            set(hcb,'Position',pos);
            
        end
        
        box on
        grid on
        set(gca,'GridLineStyle','--');
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs-2);
        
        title(['\bf' biomeString{i} '\rm, n = ' num2str(length(obsData)) ', r^2 = ' num2str(R2)],'FontWeight','Normal','FontSize',fs-2);
        
        drawnow
        
        clear obsAll obs.time obs.lon obs.lat obs.depth obs.gridIndex obs.depthIndex obs.data modelAll model.data ...
            years months obs.years obs.months meanObs stdObs lon lat depth gridIndex depthIndex obsThickness ...
            obsData modelData obsDif uGridIndex mLon mLat mGridIndex meanObs meanModel meanDif ...
            xx yy dens
        
        
    end
    
    if savePPT
        
        slideNum = exportToPPTX('addslide');
        exportToPPTX('addpicture',hFig1,'maxscale');
        
    end
    
    close all
    
end

end

%%