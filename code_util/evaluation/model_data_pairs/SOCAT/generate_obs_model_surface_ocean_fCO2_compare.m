clear
close all;

addpath(genpath('/nobackup/dcarrol2/MATLAB'));

gridDir = '/nobackup/dcarrol2/grid/LLC_270/';

dataDir = '/nobackup/dcarrol2/evaluation/mat/observations/SOCAT/gridded/';
modelDir = '/nobackup/dcarrol2/v05_latest/darwin3/run/diags/monthly/';
saveDir = '/nobackup/dcarrol2/evaluation/mat/model_obs_compare/SOCAT/gridded/';

%% 

numFaces = 13;

nx = 270;
ny = 270 .* numFaces;
nz = 50; 

XC = readbin([gridDir 'XC.data'],[nx ny],1,'real*4');
YC = readbin([gridDir 'YC.data'],[nx ny],1,'real*4');
hFacC = readbin([gridDir 'hFacC.data'],[nx ny],1,'real*4');

%% 

startYear = 1992;

SOCAT = load([dataDir 'SOCAT_2023_gridded.mat']);

startDate = datenum(startYear,1,1,0,0,0);
endDate = datenum(2022,12,15,0,0,0);

SOCAT.time = SOCAT.tmnth + datenum(1970,1,1,0,0,0);

SOCAT.timeLowerBound = squeeze(SOCAT.tmnth_bnds(1,:)) + datenum(1970,1,1,0,0,0);
SOCAT.timeUpperBound = squeeze(SOCAT.tmnth_bnds(2,:)) + datenum(1970,1,1,0,0,0);

si = find(SOCAT.time >= startDate,1);
ei = find(SOCAT.time >= endDate,1);

atm_to_uAtm = 10^6;

%% 

obs.time = SOCAT.time(si:ei);

obs.timeLowerBound = SOCAT.timeLowerBound(si:ei);
obs.timeUpperBound = SOCAT.timeUpperBound(si:ei);

obs.data_unweighted = SOCAT.fco2_ave_unwtd(:,:,si:ei);
obs.data_weighted = SOCAT.fco2_ave_weighted(:,:,si:ei);

%obs.lon = SOCAT.xlon;
%obs.lat = SOCAT.ylat;

%use the following for SOCATv2023

[xx yy] = meshgrid(SOCAT.xlon,SOCAT.ylat);

xx = xx';
yy = yy';

obs.lon = xx;
obs.lat = yy;

%%

fCO2.files = dir([modelDir 'fCO2*.*data']);
SIarea.files = dir([modelDir 'SIarea*.*data']);
mldDepth.files = dir([modelDir 'mldDepth*.*data']);
wspeed.files = dir([modelDir 'wspeed*.*data']);
apCO2.files = dir([modelDir 'apCO2.*data']);

%%

for i = 1:length(fCO2.files)
    
    fileName = fCO2.files(i).name;
    
    temp = datenum(ts2dte(str2num(fileName(end-14:end-5)),1200,startYear));

    temp = addtodate(temp, -1, 'month');
    monthlyModelTime(i) = addtodate(temp, 14, 'day');
    
end

%% 

cFCO2 = 1;

for i = 1:length(obs.time)

    obsField1 = obs.data_unweighted(:,:,i);
    obsField2 = obs.data_weighted(:,:,i);
    
    %find model output within monthly SOCAT time
    it = find(monthlyModelTime >= obs.timeLowerBound(i) & monthlyModelTime <= obs.timeUpperBound(i),1);
    
    modelField1 = readbin([modelDir fCO2.files(it).name],[nx ny],1,'real*4')  .* atm_to_uAtm;
    modelField1(hFacC == 0) = nan;
  
    modelField2 = readbin([modelDir SIarea.files(it).name],[nx ny],1,'real*4');
    modelField2(hFacC == 0) = nan;

    modelField3 = readbin([modelDir mldDepth.files(it).name],[nx ny],1,'real*4');
    modelField3(hFacC == 0) = nan;

    modelField4 = readbin([modelDir wspeed.files(it).name],[nx ny],1,'real*4');
    modelField4(hFacC == 0) = nan;

    modelField5 = readbin([modelDir apCO2.files(it).name],[nx ny],1,'real*4');
    modelField5(hFacC == 0) = nan;

    %% 
    
    iObs = find(~isnan(obsField1)); %index of available observations

    if ~isempty(iObs)
    
    	for j = 1:length(iObs)
        
            ds = (((XC(:) - obs.lon(iObs(j))) .^2) .* cosd(obs.lat(iObs(j)))) ...
                + ((YC(:) - obs.lat(iObs(j))) .^2);
        
            iModel(j) = find(ds == min(ds),1); %indices of closest grid cells to observations
        
    	end
    
    	%remove model-obs pairs on dry cells
    	il = find(hFacC(iModel) == 0); 
    
    	iModel(il) = [];
    	iObs(il) = [];
     
    	if(~isempty(iObs) && ~isempty(iModel))
        
            observations.fCO2{cFCO2}.time = (ones(1,length(obs.lat(iObs))) .* obs.time(i))';
            observations.fCO2{cFCO2}.lon = obs.lon(iObs);
            observations.fCO2{cFCO2}.lat = obs.lat(iObs);
            observations.fCO2{cFCO2}.gridLon = XC(iModel)';
            observations.fCO2{cFCO2}.gridLat = YC(iModel)';   
            observations.fCO2{cFCO2}.depth = (ones(1,length(obs.lat(iObs))) .* 0)';
            observations.fCO2{cFCO2}.gridIndex = iModel';
            observations.fCO2{cFCO2}.type = (ones(1,length(obs.lat(iObs))) .* 1)'; %1 = SOCAT

            observations.fCO2{cFCO2}.data1 = obsField1(iObs); %unweighted
            observations.fCO2{cFCO2}.data2 = obsField2(iObs); %weighted
        
            model.fCO2{cFCO2}.data = modelField1(iModel)';
            model.fCO2{cFCO2}.SIarea = modelField2(iModel)';
            model.fCO2{cFCO2}.mldDepth = modelField3(iModel)';
            model.fCO2{cFCO2}.wspeed = modelField4(iModel)';
            model.fCO2{cFCO2}.apCO2 = modelField5(iModel)';

    	end
    
    	cFCO2 = cFCO2 + 1;
 
    end
    
    disp(datestr(obs.time(i)));

    clear iModel iObs it il modelField1 modelField2 modelField3 modelField4 modelField5
    
end

observations.metadata.type = '1 = SOCAT, 2 = GLODAP, 3 = BGC-Argo';
observations.metadata.fCO2 = 'Surface-ocean fCO2 (uAtm), data1 = unweighted, data2 = weighted';

save([saveDir 'SOCAT_model_data_pairs.mat'],'observations','model','-v7.3');

%% 
