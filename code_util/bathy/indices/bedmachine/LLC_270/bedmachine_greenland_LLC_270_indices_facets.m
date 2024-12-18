clear
close all;

tic

maskDryCells = 0;

dataDir1 = '/Users/carrolld/Documents/research/bathy/mat/cell_corners/LLC_270/';
dataDir2 = '/Users/carrolld/Documents/research/bathy/mat/bedmachine/';
saveDir = '/Users/carrolld/Documents/research/bathy/mat/ice/bedmachine/greenland/LLC_270/';

%%

load([dataDir1 'cell_corners_facets.mat']);

numFacets = length(facet);

dx = 0.2; %bounding box for polygons
dy = 0.2;

%%

bedmachine = load([dataDir2 'bedmachine_greenland_lon_lat_interp.mat']);

bedmachine.lon2 = bedmachine.lon;
bedmachine.lon2(bedmachine.lon2 < 0) = bedmachine.lon2(bedmachine.lon2 < 0) + 360;

elevation = bedmachine.base';
mask = bedmachine.mask' + 1;

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
    
    bathy.medianDepth = ones(m,n);
    bathy.medianMask = ones(m,n);
    
    %%
    
    c = 1;
    
    for j = 1:m
        
        for k = 1:n
            
            %compute y points from corners
            xPoly = [XGsw(j,k) XGse(j,k) XGne(j,k) XGnw(j,k)];
            yPoly = [YGsw(j,k) YGse(j,k) YGne(j,k) YGnw(j,k)];
            
            minY = min(yPoly);
            maxY = max(yPoly);
            
            iy = find(bedmachine.lat >= (minY - dy) & bedmachine.lat <= (maxY + dy));
            
            %process 90W-90E using XG
            if min(xPoly) >= -90 & max(xPoly) <= 90
                
                wrap = 0;
                
                minX = min(xPoly);
                maxX = max(xPoly);
                
                ix = find(bedmachine.lon >= (minX - dx) & bedmachine.lon <= (maxX + dx));
                
            else %process 90E-270E using XG2 (wrap)
                
                wrap = 1;
                
                xPoly = [XG2sw(j,k) XG2se(j,k) XG2ne(j,k) XG2nw(j,k)];
                
                minX = min(xPoly);
                maxX = max(xPoly);
                
                ix = find(bedmachine.lon2 >= (minX - dx) & bedmachine.lon2 <= (maxX + dx));
                
            end
            
            %%
            %create sub grid
            
            bathyPoly = double(elevation(ix,iy));
            maskPoly = double(mask(ix,iy));
            
            [indX indY] = meshgrid(ix,iy); %bedmachine grid indices
            
            if wrap
                
                [x y] = meshgrid(bedmachine.lon2(ix),bedmachine.lat(iy));
                
            else
                
                [x y] = meshgrid(bedmachine.lon(ix),bedmachine.lat(iy));
                
            end
            
            x = x';
            y = y';
            
            in = inpolygon(x,y,xPoly,yPoly);
            
            indXX = indX(in == 1);
            indYY = indY(in == 1);
            
            bathyPoly(in == 0) = nan; %mask out region outside grid-cell polygon
            maskPoly(in == 0) = nan; %mask out region outside grid-cell polygon
            
            numTotalCells1 = length(find(~isnan(bathyPoly)));
            numTotalCells2 = length(find(~isnan(maskPoly)));
            
            numDryCells1 = length(find(bathyPoly == 0));
            numDryCells2 = length(find(maskPoly == 0));
            
            %if bedmachine search region contains >= 90% land, set model grid cell to land.
            %otherwise, exclude bedmachine land cells and compute median.
            
            if (numDryCells1 >= (numTotalCells1 * 0.9))
                
                bathy.medianDepth(j,k) = 0;
                
            else
                
                bathyPoly(isnan(bathyPoly)) = [];
                bathy.medianDepth(j,k) = nanmedian(bathyPoly(:));
                
            end
            
            if (numDryCells2 >= (numTotalCells2 * 0.9))
                
                bathy.medianMask(j,k) = 0;
                
            else
                
                maskPoly(isnan(maskPoly)) = [];
                
                bathy.medianMask(j,k) = nanmedian(maskPoly(:));
                
            end

            clear xPoly yPoly iy ix bathyPoly maskPoly indX indY x y in indXX indYY numTotalCells1 numTotalCells2 numDryCells1 numDryCells2

        end
        
    disp(num2str(j));
    
    end
    
    %%
    
    figure
    
    subplot(121);
    
    hold on
    
    set(gca,'color',[0.5 0.5 0.5]);
    
    mypcolor(bathy.medianDepth');
    
    caxis([-2000 0]);
    
    colorbar
    
    title('Median Ice');
    
    subplot(122);
    
    hold on
    
    set(gca,'color',[0.5 0.5 0.5]);
    
    mypcolor(bathy.medianMask');
    
    colorbar
    
    title('Median Mask');
    
    drawnow
    
    %pause
    
    %%
    
    if maskDryCells
        
        suffix = 'wet_dustin';
        
    else
        
        suffix = 'all_dustin';
        
    end
    
    save([saveDir 'bedmachine_greenland_LLC_270_indices_facet_' num2str(i) '_' suffix  '.mat'],'bathy','-v7.3');
    
    clear bathy
    
end

toc

%%
