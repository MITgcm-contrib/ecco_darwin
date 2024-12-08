clear 
close all

addpath(genpath('/nobackup/dcarrol2/MATLAB'));

gridDir = '/nobackup/dcarrol2/grid/LLC_270/';

dataDir1 = '/nobackup/dcarrol2/bin_average/mat/';
saveDir = '/home6/dmenemen/public/llc_270/ecco_darwin_v5/output/bin_average/';

%% 

load([dataDir1 'bin_average_LLC_270_to_1deg.mat']);

nx = 270;
ny = nx .* 13;

input.hFacC = readbin([gridDir 'hFacC.data'],[nx ny]);
output.hFacC = reshape(bin_average_cons*input.hFacC(:),[360 180]);

%% 

dataDir2 = '/home6/dmenemen/public/llc_270/ecco_darwin_v5/output/monthly/CO2_flux/';
dataDir3 = '/home6/dmenemen/public/llc_270/ecco_darwin_v5/output/monthly/pCO2/';
dataDir4 = '/home6/dmenemen/public/llc_270/ecco_darwin_v5/output/monthly/apCO2/';
dataDir5 = '/home6/dmenemen/public/llc_270/ecco_darwin_v5/output/monthly/mldDepth/';
dataDir6 = '/home6/dmenemen/public/llc_270/ecco_darwin_v5/output/monthly/SST/';
dataDir7 = '/home6/dmenemen/public/llc_270/ecco_darwin_v5/output/monthly/SSSanom/';
dataDir8 = '/home6/dmenemen/public/llc_270/ecco_darwin_v5/output/monthly/Chl1/';
dataDir9 = '/home6/dmenemen/public/llc_270/ecco_darwin_v5/output/monthly/Chl2/';
dataDir10 = '/home6/dmenemen/public/llc_270/ecco_darwin_v5/output/monthly/Chl3/';
dataDir11 = '/home6/dmenemen/public/llc_270/ecco_darwin_v5/output/monthly/Chl4/';
dataDir12 = '/home6/dmenemen/public/llc_270/ecco_darwin_v5/output/monthly/Chl5/';
dataDir13 = '/home6/dmenemen/public/llc_270/ecco_darwin_v5/output/monthly/SIarea/';
dataDir14 = '/home6/dmenemen/public/llc_270/ecco_darwin_v5/output/monthly/wspeed/';

files1 = dir([dataDir2 'CO2_flux.*data']);
files2 = dir([dataDir3 'pCO2.*data']);
files3 = dir([dataDir4 'apCO2.*data']);
files4 = dir([dataDir5 'mldDepth.*data']);
files5 = dir([dataDir6 'SST.*data']);
files6 = dir([dataDir7 'SSSanom.*data']);
files7 = dir([dataDir8 'Chl1.*data']);
files8 = dir([dataDir9 'Chl2.*data']);
files9 = dir([dataDir10 'Chl3.*data']);
files10 = dir([dataDir11 'Chl4.*data']);
files11 = dir([dataDir12 'Chl5.*data']);
files12 = dir([dataDir13 'SIarea.*data']);
files13 = dir([dataDir14 'wspeed.*data']);

%% 

c = 1;

%36:311, Jan 1995 to Dec 2017
for i = 36:311

    fileName1 = files1(i).name;
    
    time(c) = datenum(ts2dte(str2num(fileName1(10:end-5)),1200,1992));
    
    field1 = readbin([dataDir2 files1(i).name],[nx ny],1,'real*4');
    field2 = readbin([dataDir3 files2(i).name],[nx ny],1,'real*4');
    field3 = readbin([dataDir4 files3(i).name],[nx ny],1,'real*4');
    field4 = readbin([dataDir5 files4(i).name],[nx ny],1,'real*4');
    field5 = readbin([dataDir6 files5(i).name],[nx ny],1,'real*4');
    field6 = readbin([dataDir7 files6(i).name],[nx ny],1,'real*4') + 35;
    field7 = readbin([dataDir8 files7(i).name],[nx ny],1,'real*4');
    field8 = readbin([dataDir9 files8(i).name],[nx ny],1,'real*4');
    field9 = readbin([dataDir10 files9(i).name],[nx ny],1,'real*4');
    field10 = readbin([dataDir11 files10(i).name],[nx ny],1,'real*4');
    field11 = readbin([dataDir12 files11(i).name],[nx ny],1,'real*4');
    field12 = readbin([dataDir13 files12(i).name],[nx ny],1,'real*4');
    field13 = readbin([dataDir14 files13(i).name],[nx ny],1,'real*4');
    
    field1(input.hFacC == 0) = 0;
    field2(input.hFacC == 0) = 0;
    field3(input.hFacC == 0) = 0;
    field4(input.hFacC == 0) = 0;
    field5(input.hFacC == 0) = 0;
    field6(input.hFacC == 0) = 0;
    field7(input.hFacC == 0) = 0;
    field8(input.hFacC == 0) = 0;
    field9(input.hFacC == 0) = 0;
    field10(input.hFacC == 0) = 0;
    field11(input.hFacC == 0) = 0;
    field12(input.hFacC == 0) = 0;
    field13(input.hFacC == 0) = 0;
    
    binField1 = -reshape(bin_average_cons*field1(:),[360 180]); %CO2 flux
    binField2 = reshape(bin_average_cons*field2(:),[360 180]) ./ output.hFacC;
    binField3 = reshape(bin_average_cons*field3(:),[360 180]) ./ output.hFacC;
    binField4 = reshape(bin_average_cons*field4(:),[360 180]) ./ output.hFacC;
    binField5 = reshape(bin_average_cons*field5(:),[360 180]) ./ output.hFacC;
    binField6 = reshape(bin_average_cons*field6(:),[360 180]) ./ output.hFacC;
    binField7 = reshape(bin_average_cons*field7(:),[360 180]) ./ output.hFacC;
    binField8 = reshape(bin_average_cons*field8(:),[360 180]) ./ output.hFacC;
    binField9 = reshape(bin_average_cons*field9(:),[360 180]) ./ output.hFacC;
    binField10 = reshape(bin_average_cons*field10(:),[360 180]) ./ output.hFacC;
    binField11 = reshape(bin_average_cons*field11(:),[360 180]) ./ output.hFacC;
    binField12 = reshape(bin_average_cons*field12(:),[360 180]) ./ output.hFacC;
    binField13 = reshape(bin_average_cons*field13(:),[360 180]) ./ output.hFacC;
      
    %phytoplankton check for negative values
    binField7(binField7 < 0) = 0;
    binField8(binField8 < 0) = 0;
    binField9(binField9 < 0) = 0;
    binField10(binField10 < 0) = 0;
    binField11(binField11 < 0) = 0;
    
    CO2_flux(:,:,c) = binField1;
    pCO2(:,:,c) = binField2;
    apCO2(:,:,c) = binField3;
    mldDepth(:,:,c) = binField4;
    SST(:,:,c) = binField5;
    SSS(:,:,c) = binField6;
    
    Chl1(:,:,c) = binField7;
    Chl2(:,:,c) = binField8;
    Chl3(:,:,c) = binField9;
    Chl4(:,:,c) = binField10;
    Chl5(:,:,c) = binField11;
    
    seaIceArea(:,:,c) = binField12;
    windSpeed(:,:,c) = binField13;
    
    disp(num2str(c));
    
    c = c + 1;
    
end

%% 

lon = output.XG;
lat = output.YG;
area = output.RAC;

%% 

[m n t] = size(CO2Flux);

ncFileName = ['v05_ECCO-Darwin_bin_average_1x1_deg.nc'];

cd(saveDir)

delete 'ECCO-Darwin_V4.nc'

%time
nccreate([saveDir ncFileName],'time','Dimensions', {'time',inf});
ncwrite([saveDir ncFileName],'time',time);
ncwriteatt([saveDir ncFileName],'time','long_name','Time');
ncwriteatt([saveDir ncFileName],'time','units','MATLAB datenum');
ncwriteatt([saveDir ncFileName],'time','description','number of days from January 0, 0000');

%longitude
nccreate([saveDir ncFileName],'lon','Dimensions',{'x',m,'y',n});
ncwrite([saveDir ncFileName],'lon',lon);
ncwriteatt([saveDir ncFileName],'lon', 'long_name', 'Longitude');
ncwriteatt([saveDir ncFileName],'lon', 'units', 'degrees');

%latitude
nccreate([saveDir ncFileName],'lat','Dimensions',{'x',m,'y',n});
ncwrite([saveDir ncFileName],'lat',lat);
ncwriteatt([saveDir ncFileName],'lat','long_name','Latitude');
ncwriteatt([saveDir ncFileName],'lat','units','degrees');

%area
nccreate([saveDir ncFileName],'area','Dimensions',{'x',m,'y',n});
ncwrite([saveDir ncFileName],'area',area);
ncwriteatt([saveDir ncFileName],'area','long_name','area');
ncwriteatt([saveDir ncFileName],'area','units','m^-2');

%CO2Flux
nccreate([saveDir ncFileName],'CO2_flux','Dimensions', ...
    {'x',m,'y',n,'time',length(time)});

ncwrite([saveDir ncFileName],'CO2_flux',CO2_flux);
ncwriteatt([saveDir ncFileName],'CO2_flux','long_name','CO2 Flux');
ncwriteatt([saveDir ncFileName],'CO2_flux','units','mol C m^-2 s^-1');
ncwriteatt([saveDir ncFileName],'CO2_flux','description','Monthly time-averaged air-sea CO2 flux (negative values are flux into the ocean)');

%sw_pCO2
nccreate([saveDir ncFileName],'pCO2','Dimensions', ...
    {'x',m,'y',n,'time',length(time)});

ncwrite([saveDir ncFileName],'pCO2',pCO2);
ncwriteatt([saveDir ncFileName],'pCO2','long_name','Seawater pCO2');
ncwriteatt([saveDir ncFileName],'pCO2','units','atm');
ncwriteatt([saveDir ncFileName],'pCO2','description','Monthly time-averaged seawater pCO2');

%atm_pCO2
nccreate([saveDir ncFileName],'apCO2','Dimensions', ...
    {'x',m,'y',n,'time',length(time)});

ncwrite([saveDir ncFileName],'apCO2',apCO2);
ncwriteatt([saveDir ncFileName],'apCO2','long_name','Atmospheric pCO2');
ncwriteatt([saveDir ncFileName],'apCO2','units','atm');
ncwriteatt([saveDir ncFileName],'apCO2','description','Monthly time-averaged atmospheric pCO2');

%MLD
nccreate([saveDir ncFileName],'mldDepth','Dimensions', ...
    {'x',m,'y',n,'time',length(time)});

ncwrite([saveDir ncFileName],'mldDepth',mldDepth);
ncwriteatt([saveDir ncFileName],'mldDepth','long_name','Mixed layer depth');
ncwriteatt([saveDir ncFileName],'mldDepth','units','m');
ncwriteatt([saveDir ncFileName],'mldDepth','description','Monthly time-averaged mixed layer depth (Boyer formulation)');

%SST
nccreate([saveDir ncFileName],'SST','Dimensions', ...
    {'x',m,'y',n,'time',length(time)});

ncwrite([saveDir ncFileName],'SST',SST);
ncwriteatt([saveDir ncFileName],'SST','long_name','Sea-surface temperature');
ncwriteatt([saveDir ncFileName],'SST','units','deg C');
ncwriteatt([saveDir ncFileName],'SST','description','Monthly time-averaged sea-surface temperature');

%SSS
nccreate([saveDir ncFileName],'SSS','Dimensions', ...
    {'x',m,'y',n,'time',length(time)});

ncwrite([saveDir ncFileName],'SSS',SSS);
ncwriteatt([saveDir ncFileName],'SSS','long_name','Sea-surface salinity');
ncwriteatt([saveDir ncFileName],'SSS','units','psu');
ncwriteatt([saveDir ncFileName],'SSS','description','Monthly time-averaged sea-surface salinity');

%CHL1
nccreate([saveDir ncFileName],'Chl1','Dimensions', ...
    {'x',m,'y',n,'time',length(time)});

ncwrite([saveDir ncFileName],'Chl1',Chl1);
ncwriteatt([saveDir ncFileName],'Chl1','long_name','Chl a 1');
ncwriteatt([saveDir ncFileName],'Chl1','units','mg Chl a m^-3');
ncwriteatt([saveDir ncFileName],'Chl1','description','Monthly time-averaged Chl a (from phytoplankton type 1)');

%CHL2
nccreate([saveDir ncFileName],'Chl2','Dimensions', ...
    {'x',m,'y',n,'time',length(time)});

ncwrite([saveDir ncFileName],'Chl2',Chl2);
ncwriteatt([saveDir ncFileName],'Chl2','long_name','Chl a 2');
ncwriteatt([saveDir ncFileName],'Chl2','units','mg Chl a m^-3');
ncwriteatt([saveDir ncFileName],'Chl2','description','Monthly time-averaged Chl a (from phytoplankton type 2)');

%CHL3
nccreate([saveDir ncFileName],'Chl3','Dimensions', ...
    {'x',m,'y',n,'time',length(time)});

ncwrite([saveDir ncFileName],'Chl3',Chl3);
ncwriteatt([saveDir ncFileName],'Chl3','long_name','Chl a 3');
ncwriteatt([saveDir ncFileName],'Chl3','units','mg Chl a m^-3');
ncwriteatt([saveDir ncFileName],'Chl3','description','Monthly time-averaged Chl a (from phytoplankton type 3)');

%CHL4
nccreate([saveDir ncFileName],'Chl4','Dimensions', ...
    {'x',m,'y',n,'time',length(time)});

ncwrite([saveDir ncFileName],'Chl4',Chl4);
ncwriteatt([saveDir ncFileName],'Chl4','long_name','Chl a 4');
ncwriteatt([saveDir ncFileName],'Chl4','units','mg Chl a m^-3');
ncwriteatt([saveDir ncFileName],'Chl4','description','Monthly time-averaged Chl a (from phytoplankton type 4)');

%CHL5
nccreate([saveDir ncFileName],'Chl5','Dimensions', ...
    {'x',m,'y',n,'time',length(time)});

ncwrite([saveDir ncFileName],'Chl5',Chl5);
ncwriteatt([saveDir ncFileName],'Chl5','long_name','Chl a 5');
ncwriteatt([saveDir ncFileName],'Chl5','units','mg Chl a m^-3');
ncwriteatt([saveDir ncFileName],'Chl5','description','Monthly time-averaged Chl a (from phytoplankton type 5)');

%seaIceArea
nccreate([saveDir ncFileName],'seaIceArea','Dimensions', ...
    {'x',m,'y',n,'time',length(time)});

ncwrite([saveDir ncFileName],'seaIceArea',seaIceArea);
ncwriteatt([saveDir ncFileName],'seaIceArea','long_name','Sea-ice area');
ncwriteatt([saveDir ncFileName],'seaIceArea','units','none');
ncwriteatt([saveDir ncFileName],'seaIceArea','description','Monthly time-averaged sea-ice area (from 0 to 1)');

%wind speed
nccreate([saveDir ncFileName],'windSpeed','Dimensions', ...
    {'x',m,'y',n,'time',length(time)});

ncwrite([saveDir ncFileName],'windSpeed',windSpeed);
ncwriteatt([saveDir ncFileName],'windSpeed','long_name','Wind speed');
ncwriteatt([saveDir ncFileName],'windSpeed','units','m s^-1');
ncwriteatt([saveDir ncFileName],'windSpeed','description','Monthly time-averaged wind speed');

%global attributes
ncwriteatt([saveDir ncFileName],'/','creation date',datestr(now));
ncwriteatt([saveDir ncFileName],'/','file description','ECCO-Darwin v05, bin averaged onto 1 x 1 degree grid');

%% 
