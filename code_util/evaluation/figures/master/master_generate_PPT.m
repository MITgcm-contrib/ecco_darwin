clear
close all;

rootDir = '/Users/carrolld/Documents/research/evaluation/m_files/manuscript/';

cd(rootDir);

%%

%step 1: biome/superbiome maps
%step 2: misfit maps
%step 3: scatter plots
%step 4: vertical profiles
%step 5: seasonal climatology
%step 6: interannual time series

doStep = [1 1 1 1 1 1]; %set flag to 1 to generate desired analysis step

gridType = 270; %LLC grid 

startTime1 = datenum(1992,1,1,0,0,0);
startTime2 = datenum(1995,1,1,0,0,0);

endTime = datenum(2022,12,31,23,59,59);

numBins = 100; %number of bins (nx by ny) for scatter plots

%depth levels for vertical averaging
depthBin1 = [0 100];
depthBin2 = [100 500];
depthBin3 = [500 6000];
depthBin4 = [500 2000];

%%

fileName = 'v05_evaluation.ppt';

eval(['delete ' fileName]);

sOpen  = exportToPPTX();

exportToPPTX('close');

exportToPPTX('new','Dimensions',[12 6], ...
    'Title','ECCO-Darwin Evaluation', ...
    'Author','Dustin Carroll', ...
    'Subject','', ...
    'Comments','');

step = {'evaluation_1/';'evaluation_2/';'evaluation_3/';'evaluation_4/';'evaluation_5/';'evaluation_6/'};

obsType = {'SOCAT';'GLODAP';'BGC-Argo'};

%%

if doStep(1)
    
    %biome map
    caption = ['Superbiome and biome regions'];
    
    cd([rootDir step{1}]);
    
    biome_map(gridType,caption,doStep(1));
    
end

%%

if doStep(2)
    
    %SOCAT misfit map
    cd([rootDir step{2} obsType{1}]);
    
    caption = ['ECCO-Darwin vs. SOCAT model-data difference map: surface ocean'];
    SOCAT_diff_map(gridType,startTime2,endTime,caption,doStep(2));
    
    %%
    
    %GLODAP misfit map
    cd([rootDir step{2} obsType{2}]);
    
    caption = ['ECCO-Darwin vs. GLODAP model-data difference map: ' num2str(depthBin1(1)) ' to ' num2str(depthBin1(2)) '-m depth'];
    GLODAP_diff_map(gridType,startTime2,endTime,depthBin1(1),depthBin1(2),caption,doStep(2));
    
    caption = ['ECCO-Darwin vs. GLODAP model-data difference map: ' num2str(depthBin2(1)) ' to ' num2str(depthBin2(2)) '-m depth'];
    GLODAP_diff_map(gridType,startTime2,endTime,depthBin2(1),depthBin2(2),caption,doStep(2));
    
    caption = ['ECCO-Darwin vs. GLODAP model-data difference map: ' num2str(depthBin3(1)) ' to ' num2str(depthBin3(2)) '-m depth'];
    GLODAP_diff_map(gridType,startTime2,endTime,depthBin3(1),depthBin3(2),caption,doStep(2));
    
    %%
    
    %BGC-Argo misfit map
    cd([rootDir step{2} obsType{3}]);
    
    caption = ['ECCO-Darwin vs. BGC-Argo model-data difference map: ' num2str(depthBin1(1)) ' to ' num2str(depthBin1(2)) '-m depth'];
    BGC_Argo_diff_map(gridType,startTime2,endTime,depthBin1(1),depthBin1(2),caption,doStep(2));
    
    caption = ['ECCO-Darwin vs. BGC-Argo model-data difference map: ' num2str(depthBin2(1)) ' to ' num2str(depthBin2(2)) '-m depth'];
    BGC_Argo_diff_map(gridType,startTime2,endTime,depthBin2(1),depthBin2(2),caption,doStep(2));
    
    caption = ['ECCO-Darwin vs. BGC-Argo model-data difference map: ' num2str(depthBin4(1)) ' to ' num2str(depthBin4(2)) '-m depth'];
    BGC_Argo_diff_map(gridType,startTime2,endTime,depthBin4(1),depthBin4(2),caption,doStep(2));
    
end

%%

if doStep(3)
    
    %SOCAT scatter
    cd([rootDir step{3} obsType{1}]);
    
    caption = 'ECCO-Darwin vs. SOCAT scatter: surface ocean';
    SOCAT_scatter(gridType,startTime2,endTime,numBins,caption,doStep(3));
    
    %%
    
    %GLODAP scatter
    cd([rootDir step{3} obsType{2}]);
    
    caption = 'ECCO-Darwin vs. GLODAP scatter: all depths';
    GLODAP_scatter(gridType,startTime2,endTime,numBins,caption,doStep(3));
    
    %%
    
    %BGC-Argo scatter
    cd([rootDir step{3} obsType{3}]);
    
    caption = 'ECCO-Darwin vs. BGC-Argo scatter: all depths';
    BGC_Argo_scatter(gridType,startTime2,endTime,numBins,caption,doStep(3));
    
end

%%

if doStep(4)
    
    %GLODAP profiles
    cd([rootDir step{4} obsType{2}]);
    
    caption = 'ECCO-Darwin vs. GLODAP vertical profiles: all depths';
    GLODAP_profiles(gridType,startTime2,endTime,caption,doStep(4));
    
    %%
    
    %BGC-Argo profiles
    cd([rootDir step{4} obsType{3}]);
    
    caption = 'ECCO-Darwin vs. BGC-Argo vertical profiles: all depths';
    BGC_Argo_profiles(gridType,startTime2,endTime,caption,doStep(4));
    
end

%%

if doStep(5)
    
    %SOCAT seasonal climatology
    cd([rootDir step{5} obsType{1}]);
    
    caption = 'ECCO-Darwin vs. SOCAT seasonal climatology: surface ocean';
    SOCAT_seasonal_climatology(gridType,startTime2,endTime,caption,doStep(5));
    
    %%
    
    %GLODAP seasonal climatology
    cd([rootDir step{5} obsType{2}]);
    
    caption = ['ECCO-Darwin vs. GLODAP seasonal climatology: ' num2str(depthBin1(1)) ' to ' num2str(depthBin1(2)) '-m depth'];
    GLODAP_seasonal_climatology(gridType,startTime2,endTime,depthBin1(1),depthBin1(2),caption,doStep(5));
    
    caption = ['ECCO-Darwin vs. GLODAP seasonal climatology: ' num2str(depthBin2(1)) ' to ' num2str(depthBin2(2)) '-m depth'];
    GLODAP_seasonal_climatology(gridType,startTime2,endTime,depthBin2(1),depthBin2(2),caption,doStep(5));
    
    caption = ['ECCO-Darwin vs. GLODAP seasonal climatology: ' num2str(depthBin3(1)) ' to ' num2str(depthBin3(2)) '-m depth'];
    GLODAP_seasonal_climatology(gridType,startTime2,endTime,depthBin3(1),depthBin3(2),caption,doStep(5));
    
    %%
    
    %BGC-Argo seasonal climatology
    cd([rootDir step{5} obsType{3}]);
    
    caption = ['ECCO-Darwin vs. BGC-Argo seasonal climatology: ' num2str(depthBin1(1)) ' to ' num2str(depthBin1(2)) '-m depth'];
    BGC_Argo_seasonal_climatology(gridType,startTime2,endTime,depthBin1(1),depthBin1(2),caption,doStep(5));
    
    caption = ['ECCO-Darwin vs. BGC-Argo seasonal climatology: ' num2str(depthBin2(1)) ' to ' num2str(depthBin2(2)) '-m depth'];
    BGC_Argo_seasonal_climatology(gridType,startTime2,endTime,depthBin2(1),depthBin2(2),caption,doStep(5));
    
    caption = ['ECCO-Darwin vs. BGC-Argo seasonal climatology: ' num2str(depthBin4(1)) ' to ' num2str(depthBin4(2)) '-m depth'];
    BGC_Argo_seasonal_climatology(gridType,startTime2,endTime,depthBin4(1),depthBin4(2),caption,doStep(5));
    
end

%%

if doStep(6)
    
    %SOCAT interannual time series
    caption = 'SOCAT time series: surface ocean';
    
    cd([rootDir step{6} obsType{1}]);
    
    SOCAT_time_series(gridType,startTime1,endTime,caption,doStep(6));
    
    %%
    %GLODAP interannual time series
    cd([rootDir step{6} obsType{2}]);
    
    caption = ['GLODAP time series: ' num2str(depthBin1(1)) ' to ' num2str(depthBin1(2)) '-m depth'];
    GLODAP_time_series(gridType,startTime1,endTime,depthBin1(1),depthBin1(2),caption,doStep(6));
    
    caption = ['GLODAP time series: ' num2str(depthBin2(1)) ' to ' num2str(depthBin2(2)) '-m depth'];
    GLODAP_time_series(gridType,startTime1,endTime,depthBin2(1),depthBin2(2),caption,doStep(6));
    
    caption = ['GLODAP time series: ' num2str(depthBin3(1)) ' to ' num2str(depthBin3(2)) '-m depth'];
    GLODAP_time_series(gridType,startTime1,endTime,depthBin3(1),depthBin3(2),caption,doStep(6));
    
    %%
    %BGC-Argo interannual time series
    cd([rootDir step{6} obsType{3}]);
    
    caption = ['BGC-Argo time series: ' num2str(depthBin1(1)) ' to ' num2str(depthBin1(2)) '-m depth'];
    BGC_Argo_time_series(gridType,startTime1,endTime,depthBin1(1),depthBin1(2),caption,doStep(6));
    
    caption = ['BGC-Argo time series: ' num2str(depthBin2(1)) ' to ' num2str(depthBin2(2)) '-m depth'];
    BGC_Argo_time_series(gridType,startTime1,endTime,depthBin2(1),depthBin2(2),caption,doStep(6));
    
    caption = ['BGC-Argo time series: ' num2str(depthBin4(1)) ' to ' num2str(depthBin4(2)) '-m depth'];
    BGC_Argo_time_series(gridType,startTime1,endTime,depthBin4(1),depthBin4(2),caption,doStep(6));
    
end

%%

cd(rootDir);

newFile = exportToPPTX('saveandclose',fileName(1:end-4));

%%
