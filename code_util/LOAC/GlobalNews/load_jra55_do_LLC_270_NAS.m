clear
close all;

addpath(genpath('/nobackup/dcarrol2/MATLAB'));

computeGrid = 1;

gridDir = '/nobackup/dcarrol2/LOAC/grid/LLC_270_raw/';
dataDir = '/nobackup/dcarrol2/LOAC/mat/jra55_do/coast_mask/';

%%

load([dataDir 'LLC_270_coastMask_orig.mat']);

binDir = '/nobackup/dcarrol2/LOAC/bin/jra55_do/v1.4.0/';
b2 = binDir;

saveDir = '/nobackup/dcarrol2/LOAC/mat/jra55_do/';
writeDir = '/nobackup/dcarrol2/LOAC/write_bin/jra55_do/v1.4.0/LLC_270/';

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
    
    files = dir([binDir 'jra55_do_runoff_*']);
    cellArea = readbin([binDir 'cellarea.bin'], [nLon, nLat], 1,'real*4',0)';
    
    for i = 1:length(files)
        
        years(i) = str2num(files(i).name(end-3:end));
        
        fileName = files(i).name;
        
        numDays = sum(eomday(years(i),1:12),'omitnan');

        runoff(:,:,i) = sum(readbin([binDir fileName], [nLon, nLat numDays], 1,'real*4'),3,'omitnan');
        
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

disp(binDir);

for i = 1:length(years)
    
    fileName = files(i).name;
    
    numDays = sum(eomday(years(i),1:12),'omitnan');

    runoffAll = readbin([b2 files(i).name], [nLon, nLat numDays],1,'real*4');
    
    for j = 1:numDays
        
        runoff = runoffAll(:,:,j)';
        
        jra55.runoff = runoff(jra55.index);
        
        uniqueIndex = unique(jra55.llc270Index);
        
        for k = 1:length(uniqueIndex)
            
            ix = find(jra55.llc270Index == uniqueIndex(k));
            
            jra55.sumRunoff(k) = sum(jra55.runoff(ix),'omitnan');
            jra55.sumLlc270Weight(k) = mean(jra55.llc270Weight(ix),'omitnan');
            
            clear ix
            
        end
        
        llc270Runoff = (hFacC .* 0);
        
        llc270Runoff(uniqueIndex) = (jra55.sumRunoff ./ jra55.sumLlc270Weight);
    
    	%spread amazon runoff 
    	%ix = 44;
    	%iy = 3413;

    	%runoff1 = llc270Runoff(ix,iy);
    
    	%runoffArea1 = RAC(ix,iy);
    	%runoffArea2 = RAC(ix,iy+1);
    	%runoffArea3 = RAC(ix-1,iy);
    	%runoffArea4 = RAC(ix-1,iy+1);
    
    	%y1 = runoff1 .* runoffArea1;
    
    	%x1 = 0:10^-10:10^-2;
    
    	%y2 = x1.*runoffArea1 + x1.*runoffArea2 + x1.*runoffArea3 + x1.*runoffArea4;

    	%ir = find(y2 >= y1,1);
    
    	%llc270Runoff(ix,iy) = x1(ir);
    	%llc270Runoff(ix,iy+1) = x1(ir);
    	%llc270Runoff(ix-1,iy) = x1(ir);
	%llc270Runoff(ix-1,iy+1) = x1(ir);
    
        writebin([writeDir fileName(1:end-5) '_LLC_270_' fileName(end-3:end)],llc270Runoff,1,'real*4',j-1);
       
        disp(num2str(j));
        
        clear llc270Runoff
        
    end
    
    disp(years(i));
    
end

%%
