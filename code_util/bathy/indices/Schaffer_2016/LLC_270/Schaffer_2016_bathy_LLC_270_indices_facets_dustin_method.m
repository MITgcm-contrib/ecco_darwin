clear
close all;

tic

plotPoly = 0;
maskDryCells = 0;

dataDir1 = '/Users/carrolld/Documents/research/bathy/mat/cell_corners/LLC_270/';
dataDir2 = '/Users/carrolld/Documents/research/bathy/raw_data/Schaffer_2016/';
saveDir = '/Users/carrolld/Documents/research/bathy/mat/indices/Schaffer_2016/LLC_270/';

%%

load([dataDir1 'cell_corners_facets.mat']);

numFacets = length(facet);

dx = 1; %bounding box for polygons
dy = 1;

%%

fileName = 'RTopo-2.0.1_30sec_bedrock_topography.nc';

Schaffer_2016.lon = ncread([dataDir2 fileName],'lon');
Schaffer_2016.lat = ncread([dataDir2 fileName],'lat');

Schaffer_2016.lon2 = Schaffer_2016.lon;
Schaffer_2016.lon2(Schaffer_2016.lon2 < 0) = Schaffer_2016.lon2(Schaffer_2016.lon2 < 0) + 360;

elevation = -ncread([dataDir2 fileName],'bedrock_topography');
elevation(elevation <= 0) = 0;

%%

for i = 1:numFacets
    
    [m n] = size(facet{i}.XGne);
    
    XGne = facet{i}.XGne;
    XGnw = facet{i}.XGnw;
    XGse = facet{i}.XGse;
    XGsw = facet{i}.XGsw;
    
    YGne = facet{i}.YGne;
    YGnw = facet{i}.YGnw;
    YGse = facet{i}.YGse;
    YGsw = facet{i}.YGsw;
    
    %wrapped set
    XG2ne = XGne;
    XG2nw = XGnw;
    XG2se = XGse;
    XG2sw = XGsw;
    
    XG2ne(find(XGne < 0)) = XGne(find(XGne < 0)) + 360;
    XG2nw(find(XGnw < 0)) = XGnw(find(XGnw < 0)) + 360;
    XG2se(find(XGse < 0)) = XGse(find(XGse < 0)) + 360;
    XG2sw(find(XGsw< 0)) = XGsw(find(XGsw < 0)) + 360;
    
    bathy.minDepth = ones(m,n) * -9999;
    bathy.maxDepth = ones(m,n) * -9999;
    
    bathy.meanDepth = ones(m,n) * -9999;
    bathy.medianDepth = ones(m,n) * -9999;
    bathy.numWetCells = ones(m,n) * -9999;
    
    c = 1;
    
    for j = 1:m
        
        for k = 1:n
            
            %compute y points from corners
            xPoly = [XGsw(j,k) XGse(j,k) XGne(j,k) XGnw(j,k)];
            yPoly = [YGsw(j,k) YGse(j,k) YGne(j,k) YGnw(j,k)];
            
            minY = min(yPoly);
            maxY = max(yPoly);
            
            iy = find(Schaffer_2016.lat >= (minY - dy) & Schaffer_2016.lat <= (maxY + dy));
            
            %process 90W-90E using XG
            if min(xPoly) >= -90 & max(xPoly) <= 90
                
                wrap = 0;
                
                minX = min(xPoly);
                maxX = max(xPoly);
                
                ix = find(Schaffer_2016.lon >= (minX - dx) & Schaffer_2016.lon <= (maxX + dx));
                
            else %process 90E-270E using XG2 (wrap)
                
                wrap = 1;
                
                xPoly = [XG2sw(j,k) XG2se(j,k) XG2ne(j,k) XG2nw(j,k)];
                
                minX = min(xPoly);
                maxX = max(xPoly);
                
                ix = find(Schaffer_2016.lon2 >= (minX - dx) & Schaffer_2016.lon2 <= (maxX + dx));
                
            end
            
            %create sub grid
            
            bathyPoly = double(elevation(ix,iy));
            
            [indX indY] = meshgrid(ix,iy); %Schaffer_2016 grid indices
            
            if wrap
                
                [x y] = meshgrid(Schaffer_2016.lon2(ix),Schaffer_2016.lat(iy));
                
            else
                
                [x y] = meshgrid(Schaffer_2016.lon(ix),Schaffer_2016.lat(iy));
                
            end
            
            x = x';
            y = y';
            
            in = inpolygon(x,y,xPoly,yPoly);
            
            indXX = indX(in == 1);
            indYY = indY(in == 1);
            
            bathyPoly(in == 0) = nan; %mask out region outside grid-cell polygon
            
            numTotalCells = length(find(~isnan(bathyPoly)));
            
            numDryCells = length(find(bathyPoly == 0));
            
            %if Schaffer_2016 search region contains >= 90% land, set model grid cell to land.
            %otherwise, exclude Schaffer_2016 land cells and compute median.
            
            if (numDryCells >= (numTotalCells .* 0.9))
                
                bathy.numWetCells(j,k) = nan;
                bathy.medianDepth(j,k) = nan;
                
            else
                
                if plotPoly
                    
                    subplot(121);
                    
                    hold on
                    
                    scatter(x(:),y(:),'r');
                    scatter(xPoly(:),yPoly(:),'b','filled');
                    line([xPoly xPoly(1)],[yPoly yPoly(1)]);
                    
                    scatter(x(in == 1),y(in == 1),'m');
                    
                    subplot(122);
                    
                    hold on
                    
                    pcolorcen(x,y,bathyPoly);
                    scatter(x(in == 1),y(in == 1),'m');
                    
                    drawnow
                    
                    pause
                    
                end
                
                if(maskDryCells)
                    
                    bathyPoly(bathyPoly == 0) = nan;
                    
                end
                
                bathyPoly(isnan(bathyPoly)) = [];
                
                bathy.numWetCells(j,k) = length(find(bathyPoly ~= 0));
                bathy.medianDepth(j,k) = median(bathyPoly);
                
                if(isnan(bathy.medianDepth(j,k)))
                    
                    disp('NaN found');
                    pause
                    
                end
                
                c = c + 1;
                
            end
            
            clear xPoly yPoly iy ix bathyPoly indX indY x y in indXX indYY numTotalCells numDryCells
            
        end
        
        disp(num2str(j));
        
    end
    
    %%
    
    if maskDryCells
        
        suffix = 'wet_dustin';
        
    else
        
        suffix = 'all_dustin';
        
    end
    
    save([saveDir 'Schaffer_2016_bathy_LLC_270_indices_facet_' num2str(i) '_' suffix  '.mat'],'bathy','-v7.3');
    
    clear bathy
    
end

toc

%%
