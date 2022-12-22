clear
close all;

%addpath(genpath('/nobackup/dcarrol2/MATLAB'));
addpath ~dmenemen/matlab

computeGrid = 0;

gridDir = '/nobackup/dcarrol2/LOAC/grid/LLC_270_raw/';
dataDir = '/nobackup/dcarrol2/LOAC/mat/jra55_do/coast_mask/';

%%

load([dataDir 'LLC_270_coastMask_orig.mat']);

binDir1 = '/nobackup/dcarrol2/LOAC/bin/jra55_do/v1.4.0/';
binDir2 = '/nobackup/dmenemen/forcing/jra55_do/GlobalNEWS/GlobalNEWS2_on_jra55v1.4.0/';
%b2 = binDir;

saveDir = '/nobackup/rsavelli/LOAC/mat/jra55_do/';
writeDir = '/nobackup/rsavelli/LOAC/write_bin/jra55_do/v1.4.0/LLC_270/';

%%

numFaces = 13;

nx = 270;
ny = numFaces .* nx;

xc = readbin([gridDir 'XC.data'],[nx ny],1,'real*4');
yc = readbin([gridDir 'YC.data'],[nx ny],1,'real*4');
RAC = readbin([gridDir 'RAC.data'],[nx ny],1,'real*4');
hFacC = readbin([gridDir 'hFacC.data'],[nx ny],1,'real*4');

xc(xc < 0) = xc(xc < 0) + 360;

%%

nLon = 1440;
nLat = 720;

lonInc = ones(1,nLon-1) .* 0.25;
latInc = ones(1,nLat-1) .* 0.25;

lon0 = 0.125;
lat0 = -89.875;

lon = cumsum([lon0 lonInc]);
lat = cumsum([lat0 latInc]);

[xx yy] = meshgrid(lon,lat);

%%

latLim = [lat(1) lat(end)];
lonLim = [lon(1) lon(end)];
rasterSize = [nLat nLon];

%R = georefcells(latLim,lonLim,rasterSize,'ColumnsStartFrom','north');
%sphericalEarth = referenceSphere('earth','m');

%[a aVec] = areamat(xx.*0 + 1,R,sphericalEarth);

%da = repmat(aVec,[1 nLon]);

%%

if computeGrid
    
    files = dir([binDir1 'jra55_do_runoff_*']); 
    
    cellArea = readbin([binDir1 'cellarea.bin'], [nLon, nLat], 1,'real*4',0)';
    
    for i = 1:length(files)
        
        years(i) = str2num(files(i).name(end-3:end));
                        
        fileName = files(i).name;
        
        numDays = sum(eomday(years(i),1:12),'omitnan');

        runoff(:,:,i) = sum(readbin([binDir1 fileName], [nLon, nLat numDays], 1,'real*4'),3,'omitnan');
        
    end
   
     runoff = sum(runoff,3,'omitnan')';
    
    %% 
    
    wetXc = xc;
    wetYc = yc;
    
    wetXc(~isnan(coastMask)) = nan;
    wetYc(~isnan(coastMask)) = nan;
    
    ir = find(runoff ~= 0);
    
    runoffLon = xx(ir);
    runoffLat = yy(ir);
    
    for i = 1:length(runoffLon)
        
        distWet = ((runoffLon(i) - wetXc).^2 .* (cosd(runoffLat(i)))) + ...
            ((runoffLat(i) - wetYc).^2);
        
        ix = find(distWet == min(distWet(:)),1); %index of closest wet grid cell to release location
        
        jra55.index(i) = ir(i);
        jra55.area(i) = cellArea(ir(i));
        
        jra55.llc270Index(i) = ix;
        jra55.llc270Area(i) = RAC(ix);
        
        jra55.llc270Weight(i) = RAC(ix) ./ cellArea(ir(i));
        
        clear ix
        
    end
    
    jra55.years = years;
    
    save([saveDir 'jra55_do_LLC_270_grid_orig.mat']);
    
else
    load([saveDir 'jra55_do_LLC_270_grid_orig.mat']);
    
end

%%

disp(binDir2);

files = dir([binDir2 '*_*']); %go through all nutrient runoff files

for i = 1:length(files)
    years(i) = str2num(files(i).name(end-3:end));
end

for i = 1:length(years)
    
    fileName = files(i).name;
    
    numDays = sum(eomday(years(i),1:12),'omitnan');

    runoffAll = readbin([binDir2 files(i).name], [nLon, nLat numDays],1,'real*4');
    
    for j = 1:numDays
        
        nutrient_runoff = runoffAll(:,:,j)';
        
        jra55.nutrient_runoff = nutrient_runoff(jra55.index);
        
        uniqueIndex = unique(jra55.llc270Index);
        
        for k = 1:length(uniqueIndex)
            
            ix = find(jra55.llc270Index == uniqueIndex(k));
            
            jra55.sumRunoff(k) = sum(jra55.nutrient_runoff(ix),'omitnan');
            jra55.sumLlc270Weight(k) = mean(jra55.llc270Weight(ix),'omitnan');
            
            clear ix
            
        end
        
        llc270Runoff = (hFacC .* 0);
        
        llc270Runoff(uniqueIndex) = (jra55.sumRunoff ./ jra55.sumLlc270Weight);
    
    
    if length(fileName) == 8
        writebin([writeDir fileName(1:3) '_LLC_270_' fileName(end-3:end)],llc270Runoff,1,'real*4',j-1);
    else
        writebin([writeDir fileName(1:2) '_LLC_270_' fileName(end-3:end)],llc270Runoff,1,'real*4',j-1);
    end
        disp(num2str(j));
        
        clear llc270Runoff
        
    end
    
    disp(years(i));
    
end

%%
