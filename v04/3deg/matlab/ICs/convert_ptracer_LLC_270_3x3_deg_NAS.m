clear
close all;

addpath(genpath('/nobackup/dcarrol2/MATLAB'));

writePickup = 1;

gridDir = '/nobackup/dcarrol2/grid/LLC_270/';
interpWeightsDir = '/nobackup/dcarrol2/v4_3deg/mat/interp/'

pickupDir = '/nobackup/ojahn/ecco_darwin/v4_llc270/darwin3/';
pickupSaveDir = '/nobackup/dcarrol2/v4_3deg/pickup/';

cd(pickupSaveDir)
delete *

load([interpWeightsDir 'interp_weights_llc270_3x3_deg.mat']);

numTracers = 31;

%% 

numFaces = 13;

nx = 270;
ny = nx .* numFaces;
nz = 50;

input.XC = readbin([gridDir 'XC.data'],[nx ny],1,'real*4');
input.YC = readbin([gridDir 'YC.data'],[nx ny],1,'real*4');

input.hFacC = readbin([gridDir 'hFacC.data'],[nx ny nz],1,'real*4');
input.z = abs(readbin([gridDir 'RF.data'],[nz],1,'real*4'));

for i = 1:nz
    
    tmpin = input.hFacC(:,:,i);
    z = tmpin(:).';
    zi = sum(z(tri) .* w,2);
    
    tmpout = reshape(zi,siz);
    
    output.hFacC(:,:,i) = tmpout;
    
    clear tmpout
    
    disp(num2str(i));
    
end

%% 

lat0 = -90 + 2.8125/2;
lon0 = 2.8125/2;

lon_inc = ones(1,128) * 2.8125;
lat_inc = ones(1,64) * 2.8125;

lonGrid = [lon0 (lon0 + cumsum(lon_inc(1:end-1)))];
latGrid = [lat0 (lat0 + cumsum(lat_inc(1:end-1)))];

[output.XC output.YC] = meshgrid(lonGrid,latGrid);

output.delZ =  [50 70 100 140 190 240 290 340 390 440 490 540 590 640 690];
output.z = cumsum(output.delZ);

%% 

for i = 1:length(output.z)
    
    zIndices(i) = find(input.z >= output.z(i),1);
    
end

%%

for i = 1:numTracers
    
    for j = 1:length(zIndices)
        
	pickupName =  ['ptracers_optimized_' num2str(i,'%02.f') '.0000000001'];
        
        %read input pickup
        tracer = readbin([pickupDir 'pickup_ptracers_optimized.0000000001.data'],[nx ny],1,'real*8',(i-1)*50+zIndices(j)-1);
     
        tmpin = tracer;
        z = tmpin(:).';
        zi = sum(z(tri) .* w,2);
    
        interpTracer = reshape(zi,siz);
    
        if writePickup
            
            writebin([pickupSaveDir 'v4_3deg_pickup_ptracers_optimized_darwin3.0000000001.data'],interpTracer,1,'real*8',(i-1)*50+j-1);
            
	    writebin([pickupSaveDir pickupName],interpTracer,1,'real*4',j-1);
            
        end
   
     
    end
    
    disp(['Tracer Number: ' num2str(i)]);
    
end

