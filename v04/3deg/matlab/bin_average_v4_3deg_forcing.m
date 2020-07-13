clear
close all

codeDir = '/Users/carrolld/Documents/research/carbon/m_files/bin_average/standard/';
dataDir = '/Users/carrolld/Documents/research/carbon/MITgcm_setup/forcing/exf/raw/';
gridDir = '/Users/carrolld/Documents/research/carbon/simulations/grid/LLC_270/';
saveDir1 = '/Users/carrolld/Documents/research/carbon/mat/bin_average/standard/';
saveDir2 =  '/Users/carrolld/Documents/research/carbon/MITgcm_setup/forcing/exf/bin_average/';
saveDir3 = '/Users/carrolld/Documents/research/carbon/mat/bin_average/exf/';

cd(saveDir2)

delete *

%% 

binLonInc = 3;
binLatInc = 3;

%% 

numFaces = 13;

nx = 270;
ny = nx .* numFaces;
nz = 50;

LLC270.XC = readbin([gridDir 'XC.data'],[nx ny],1,'real*4');
LLC270.YC = readbin([gridDir 'YC.data'],[nx ny],1,'real*4');

LLC270.dXG = readbin([gridDir 'dXG.data'],[nx ny]);
LLC270.dYG = readbin([gridDir 'dYG.data'],[nx ny]);

LLC270.hFacC = readbin([gridDir 'hFacC.data'],[nx ny],1,'real*4');
LLC270.RAC = readbin([gridDir 'RAC.data'],[nx ny],1,'real*4');

rSphere = mmax(LLC270.dYG).*4.*nx/2/pi;
rDeg = 2.*pi.*rSphere/360;

%% 

%  atemp_lon0        =   0.0000000D0,
%  atemp_lon_inc     =   0.7031250D0,
%  atemp_lat0        = -89.4628220D0,
%  atemp_lat_inc     = 0.6958694,0.6999817,0.7009048,0.7012634,0.7014313,
%                      245*0.7017418,
%                      0.7014313,0.7012634,0.7009048,0.6999817,0.6958694
%  atemp_nlon        = 512,
%  atemp_nlat        = 256,
%

%  atempstartdate1   = 19920101,
%  atempstartdate2   = 000000,
%  atempperiod       = 21600.0

atemp.filename = 'EIG_tmp2m_degC_plus_ECCO_v4r1_ctrl_2000';

atemp.lon0 = 0;
atemp.lon_inc =  0.7031250;
atemp.lat0 = -89.4628220;
atemp.lat_inc = [0.6958694,0.6999817,0.7009048,0.7012634,0.7014313, ...
    (ones(1,245) * 0.7017418), 0.7014313,0.7012634,0.7009048,0.6999817,0.6958694];
atemp.nlon = 512;
atemp.nlat = 256;
atemp.lonVec = [atemp.lon0 (atemp.lon0 + cumsum(ones(1,atemp.nlon-1)*atemp.lon_inc))];
atemp.latVec = [atemp.lat0 (atemp.lat0 + cumsum(atemp.lat_inc))];

[atemp.lon atemp.lat] = meshgrid(atemp.lonVec,atemp.latVec);

atemp.lon = atemp.lon';
atemp.lat = atemp.lat';

atemp.area = rDeg^2*cosd(atemp.lat);

atemp.time = datenum(1992,1,1,0,0,0):1/4:datenum(1992,12,31,23,0,0);
[atemp.years atemp.months atemp.days atemp.hours atemp.minutes atemp.seconds] = datevec(atemp.time);

atemp.field = readbin([dataDir atemp.filename],[atemp.nlon atemp.nlat length(atemp.time)],1,'real*4');

saveFilename = 'bin_average_atemp.mat'

cd(codeDir);
[blon blat bin_average] = compute_bin_average(atemp.lon,atemp.lat,atemp.area,binLonInc,binLatInc,saveDir1,saveFilename);

atemp.uniqueMonths = unique(atemp.months);

for i = 1:length(atemp.uniqueMonths)
    
    im = find(atemp.months == atemp.uniqueMonths(i));
    
    meanField = squeeze(nanmean(atemp.field(:,:,im),3));
    
    atemp.meanField(:,:,i) = meanField;
    atemp.binMeanField(:,:,i) = reshape(bin_average*double(meanField(:)),length(blon),length(blat));
    
    clear im meanField
    
end

writebin([saveDir2 atemp.filename],atemp.binMeanField,1,'real*4');
save([saveDir3 'atemp.mat'],'atemp','blon','blat','-v7.3');

%% 

%  aqh_lon0        =   0.0000000D0,
%  aqh_lon_inc     =   0.7031250D0,
%  aqh_lat0        = -89.4628220D0,
%  aqh_lat_inc     = 0.6958694,0.6999817,0.7009048,0.7012634,0.7014313,
%                      245*0.7017418,
%                      0.7014313,0.7012634,0.7009048,0.6999817,0.6958694
%  aqh_nlon        = 512,
%  aqh_nlat        = 256,

%  aqhstartdate1     = 19920101,
%  aqhstartdate2     = 000000,
%  aqhperiod         = 21600.0,
 
aqh.filename = 'EIG_spfh2m_plus_ECCO_v4r1_ctrl_2000';

aqh.lon0 = 0;
aqh.lon_inc =  0.7031250;
aqh.lat0 = -89.4628220;
aqh.lat_inc = [0.6958694,0.6999817,0.7009048,0.7012634,0.7014313, ...
    (ones(1,245) * 0.7017418), 0.7014313,0.7012634,0.7009048,0.6999817,0.6958694];
aqh.nlon = 512;
aqh.nlat = 256;
aqh.lonVec = [aqh.lon0 (aqh.lon0 + cumsum(ones(1,aqh.nlon-1)*aqh.lon_inc))];
aqh.latVec = [aqh.lat0 (aqh.lat0 + cumsum(aqh.lat_inc))];

[aqh.lon aqh.lat] = meshgrid(aqh.lonVec,aqh.latVec);

aqh.lon = aqh.lon';
aqh.lat = aqh.lat';

aqh.area = rDeg^2*cosd(aqh.lat);

aqh.time = datenum(1992,1,1,0,0,0):1/4:datenum(1992,12,31,23,0,0);
[aqh.years aqh.months aqh.days aqh.hours aqh.minutes aqh.seconds] = datevec(aqh.time);

aqh.field = readbin([dataDir aqh.filename],[aqh.nlon aqh.nlat length(aqh.time)],1,'real*4');

saveFilename = 'bin_average_aqh.mat'

cd(codeDir);
[blon blat bin_average] = compute_bin_average(aqh.lon,aqh.lat,aqh.area,binLonInc,binLatInc,saveDir1,saveFilename);

aqh.uniqueMonths = unique(aqh.months);

for i = 1:length(aqh.uniqueMonths)
    
    im = find(aqh.months == aqh.uniqueMonths(i));
    
    meanField = squeeze(nanmean(aqh.field(:,:,im),3));
    
    aqh.meanField(:,:,i) = meanField;
    aqh.binMeanField(:,:,i) = reshape(bin_average*double(meanField(:)),length(blon),length(blat));
    
    clear im meanField
    
end

writebin([saveDir2 aqh.filename],aqh.binMeanField,1,'real*4');
save([saveDir3 'aqh.mat'],'aqh','blon','blat','-v7.3');

%% 

%  precip_lon0        =   0.0000000D0,
%  precip_lon_inc     =   0.7031250D0,
%  precip_lat0        = -89.4628220D0,
%  precip_lat_inc     = 0.6958694,0.6999817,0.7009048,0.7012634,0.7014313,
%                      245*0.7017418,
%                      0.7014313,0.7012634,0.7009048,0.6999817,0.6958694
%  precip_nlon        = 512,
%  precip_nlat        = 256,

%  precipstartdate1  = 19920101,
%  precipstartdate2  = 030000,
%  precipperiod      = 21600.0,

precip.filename = 'EIG_rain_plus_ECCO_v4r1_ctrl_2000';

precip.lon0 = 0;
precip.lon_inc =  0.7031250;
precip.lat0 = -89.4628220;
precip.lat_inc = [0.6958694,0.6999817,0.7009048,0.7012634,0.7014313, ...
    (ones(1,245) * 0.7017418), 0.7014313,0.7012634,0.7009048,0.6999817,0.6958694];
precip.nlon = 512;
precip.nlat = 256;
precip.lonVec = [precip.lon0 (precip.lon0 + cumsum(ones(1,precip.nlon-1)*precip.lon_inc))];
precip.latVec = [precip.lat0 (precip.lat0 + cumsum(precip.lat_inc))];

[precip.lon precip.lat] = meshgrid(precip.lonVec,precip.latVec);

precip.lon = precip.lon';
precip.lat = precip.lat';

precip.area = rDeg^2*cosd(precip.lat);

precip.time = datenum(1992,1,1,0,0,0):1/4:datenum(1992,12,31,23,0,0);
[precip.years precip.months precip.days precip.hours precip.minutes precip.seconds] = datevec(precip.time);

precip.field = readbin([dataDir precip.filename],[precip.nlon precip.nlat length(precip.time)],1,'real*4');

saveFilename = 'bin_average_precip.mat'

cd(codeDir);
[blon blat bin_average] = compute_bin_average(precip.lon,precip.lat,precip.area,binLonInc,binLatInc,saveDir1,saveFilename);

precip.uniqueMonths = unique(precip.months);

for i = 1:length(precip.uniqueMonths)
    
    im = find(precip.months == precip.uniqueMonths(i));
    
    meanField = squeeze(nanmean(precip.field(:,:,im),3));
    
    precip.meanField(:,:,i) = meanField;
    precip.binMeanField(:,:,i) = reshape(bin_average*double(meanField(:)),length(blon),length(blat));
    
    clear im meanField
    
end

writebin([saveDir2 precip.filename],precip.binMeanField,1,'real*4');
save([saveDir3 'precip.mat'],'precip','blon','blat','-v7.3');

%% 

%  uwind_lon0        =   0.0000000D0,
%  uwind_lon_inc     =   0.7031250D0,
%  uwind_lat0        = -89.4628220D0,
%  uwind_lat_inc     = 0.6958694,0.6999817,0.7009048,0.7012634,0.7014313,
%                      245*0.7017418,
%                      0.7014313,0.7012634,0.7009048,0.6999817,0.6958694
%  uwind_nlon        = 512,
%  uwind_nlat        = 256,

%  uwindstartdate1   = 19920101,
%  uwindstartdate2   = 000000,
%  uwindperiod       = 21600.0,

uwind.filename = 'EIG_u10m_2000';

uwind.lon0 = 0;
uwind.lon_inc =  0.7031250;
uwind.lat0 = -89.4628220;
uwind.lat_inc = [0.6958694,0.6999817,0.7009048,0.7012634,0.7014313, ...
    (ones(1,245) * 0.7017418), 0.7014313,0.7012634,0.7009048,0.6999817,0.6958694];
uwind.nlon = 512;
uwind.nlat = 256;
uwind.lonVec = [uwind.lon0 (uwind.lon0 + cumsum(ones(1,uwind.nlon-1)*uwind.lon_inc))];
uwind.latVec = [uwind.lat0 (uwind.lat0 + cumsum(uwind.lat_inc))];

[uwind.lon uwind.lat] = meshgrid(uwind.lonVec,uwind.latVec);

uwind.lon = uwind.lon';
uwind.lat = uwind.lat';

uwind.area = rDeg^2*cosd(uwind.lat);

uwind.time = datenum(1992,1,1,0,0,0):1/4:datenum(1992,12,31,23,0,0);
[uwind.years uwind.months uwind.days uwind.hours uwind.minutes uwind.seconds] = datevec(uwind.time);

uwind.field = readbin([dataDir uwind.filename],[uwind.nlon uwind.nlat length(uwind.time)],1,'real*4');

saveFilename = 'bin_average_uwind.mat'

cd(codeDir);
[blon blat bin_average] = compute_bin_average(uwind.lon,uwind.lat,uwind.area,binLonInc,binLatInc,saveDir1,saveFilename);

uwind.uniqueMonths = unique(uwind.months);

for i = 1:length(uwind.uniqueMonths)
    
    im = find(uwind.months == uwind.uniqueMonths(i));
    
    meanField = squeeze(nanmean(uwind.field(:,:,im),3));
    
    uwind.meanField(:,:,i) = meanField;
    uwind.binMeanField(:,:,i) = reshape(bin_average*double(meanField(:)),length(blon),length(blat));
    
    clear im meanField
    
end

writebin([saveDir2 uwind.filename],uwind.binMeanField,1,'real*4');
save([saveDir3 'uwind.mat'],'uwind','blon','blat','-v7.3');

%% 

% vwind_lon0        =   0.0000000D0,
% vwind_lon_inc     =   0.7031250D0,
% vwind_lat0        = -89.4628220D0,
% vwind_lat_inc     = 0.6958694,0.6999817,0.7009048,0.7012634,0.7014313,
%                     245*0.7017418,
%                     0.7014313,0.7012634,0.7009048,0.6999817,0.6958694
% vwind_nlon        = 512,
% vwind_nlat        = 256,

%  vwindstartdate1   = 19920101,
%  vwindstartdate2   = 000000,
%  vwindperiod       = 21600.0,

vwind.filename = 'EIG_v10m_2000';

vwind.lon0 = 0;
vwind.lon_inc =  0.7031250;
vwind.lat0 = -89.4628220;
vwind.lat_inc = [0.6958694,0.6999817,0.7009048,0.7012634,0.7014313, ...
    (ones(1,245) * 0.7017418), 0.7014313,0.7012634,0.7009048,0.6999817,0.6958694];
vwind.nlon = 512;
vwind.nlat = 256;
vwind.lonVec = [vwind.lon0 (vwind.lon0 + cumsum(ones(1,vwind.nlon-1)*vwind.lon_inc))];
vwind.latVec = [vwind.lat0 (vwind.lat0 + cumsum(vwind.lat_inc))];

[vwind.lon vwind.lat] = meshgrid(vwind.lonVec,vwind.latVec);

vwind.lon = vwind.lon';
vwind.lat = vwind.lat';

vwind.area = rDeg^2*cosd(vwind.lat);

vwind.time = datenum(1992,1,1,0,0,0):1/4:datenum(1992,12,31,23,0,0);
[vwind.years vwind.months vwind.days vwind.hours vwind.minutes vwind.seconds] = datevec(vwind.time);

vwind.field = readbin([dataDir vwind.filename],[vwind.nlon vwind.nlat length(vwind.time)],1,'real*4');

saveFilename = 'bin_average_vwind.mat'

cd(codeDir);
[blon blat bin_average] = compute_bin_average(vwind.lon,vwind.lat,vwind.area,binLonInc,binLatInc,saveDir1,saveFilename);

vwind.uniqueMonths = unique(vwind.months);

for i = 1:length(vwind.uniqueMonths)
    
    im = find(vwind.months == vwind.uniqueMonths(i));
    
    meanField = squeeze(nanmean(vwind.field(:,:,im),3));
    
    vwind.meanField(:,:,i) = meanField;
    vwind.binMeanField(:,:,i) = reshape(bin_average*double(meanField(:)),length(blon),length(blat));
    
    clear im meanField
    
end

writebin([saveDir2 vwind.filename],vwind.binMeanField,1,'real*4');
save([saveDir3 'vwind.mat'],'vwind','blon','blat','-v7.3');

%% 

%  swdown_lon0        =   0.0000000D0,
%  swdown_lon_inc     =   0.7031250D0,
%  swdown_lat0        = -89.4628220D0,
%  swdown_lat_inc     = 0.6958694,0.6999817,0.7009048,0.7012634,0.7014313,
%                      245*0.7017418,
%                      0.7014313,0.7012634,0.7009048,0.6999817,0.6958694
%  swdown_nlon        = 512,
%  swdown_nlat        = 256,

%  swdownstartdate1  = 19920101,
%  swdownstartdate2  = 030000,
%  swdownperiod      = 21600.0,

swdown.filename = 'EIG_dsw_plus_ECCO_v4r1_ctrl_2000';

swdown.lon0 = 0;
swdown.lon_inc =  0.7031250;
swdown.lat0 = -89.4628220;
swdown.lat_inc = [0.6958694,0.6999817,0.7009048,0.7012634,0.7014313, ...
    (ones(1,245) * 0.7017418), 0.7014313,0.7012634,0.7009048,0.6999817,0.6958694];
swdown.nlon = 512;
swdown.nlat = 256;
swdown.lonVec = [swdown.lon0 (swdown.lon0 + cumsum(ones(1,swdown.nlon-1)*swdown.lon_inc))];
swdown.latVec = [swdown.lat0 (swdown.lat0 + cumsum(swdown.lat_inc))];

[swdown.lon swdown.lat] = meshgrid(swdown.lonVec,swdown.latVec);

swdown.lon = swdown.lon';
swdown.lat = swdown.lat';

swdown.area = rDeg^2*cosd(swdown.lat);

swdown.time = datenum(1992,1,1,0,0,0):1/4:datenum(1992,12,31,23,0,0);
[swdown.years swdown.months swdown.days swdown.hours swdown.minutes swdown.seconds] = datevec(swdown.time);

swdown.field = readbin([dataDir swdown.filename],[swdown.nlon swdown.nlat length(swdown.time)],1,'real*4');

saveFilename = 'bin_average_swdown.mat'

cd(codeDir);
[blon blat bin_average] = compute_bin_average(swdown.lon,swdown.lat,swdown.area,binLonInc,binLatInc,saveDir1,saveFilename);

swdown.uniqueMonths = unique(swdown.months);

for i = 1:length(swdown.uniqueMonths)
    
    im = find(swdown.months == swdown.uniqueMonths(i));
    
    meanField = squeeze(nanmean(swdown.field(:,:,im),3));
    
    swdown.meanField(:,:,i) = meanField;
    swdown.binMeanField(:,:,i) = reshape(bin_average*double(meanField(:)),length(blon),length(blat));
    
    clear im meanField
    
end

writebin([saveDir2 swdown.filename],swdown.binMeanField,1,'real*4');
save([saveDir3 'swdown.mat'],'swdown','blon','blat','-v7.3');

%% 

%  lwdown_lon0        =   0.0000000D0,
%  lwdown_lon_inc     =   0.7031250D0,
%  lwdown_lat0        = -89.4628220D0,
%  lwdown_lat_inc     = 0.6958694,0.6999817,0.7009048,0.7012634,0.7014313,
%                      245*0.7017418,
%                      0.7014313,0.7012634,0.7009048,0.6999817,0.6958694
%  lwdown_nlon        = 512,
%  lwdown_nlat        = 256,

%  lwdownstartdate1  = 19920101,
%  lwdownstartdate2  = 030000,
%  lwdownperiod      = 21600.0,

lwdown.filename = 'EIG_dlw_plus_ECCO_v4r1_ctrl_2000';

lwdown.lon0 = 0;
lwdown.lon_inc =  0.7031250;
lwdown.lat0 = -89.4628220;
lwdown.lat_inc = [0.6958694,0.6999817,0.7009048,0.7012634,0.7014313, ...
    (ones(1,245) * 0.7017418), 0.7014313,0.7012634,0.7009048,0.6999817,0.6958694];
lwdown.nlon = 512;
lwdown.nlat = 256;
lwdown.lonVec = [lwdown.lon0 (lwdown.lon0 + cumsum(ones(1,lwdown.nlon-1)*lwdown.lon_inc))];
lwdown.latVec = [lwdown.lat0 (lwdown.lat0 + cumsum(lwdown.lat_inc))];

[lwdown.lon lwdown.lat] = meshgrid(lwdown.lonVec,lwdown.latVec);

lwdown.lon = lwdown.lon';
lwdown.lat = lwdown.lat';

lwdown.area = rDeg^2*cosd(lwdown.lat);

lwdown.time = datenum(1992,1,1,0,0,0):1/4:datenum(1992,12,31,23,0,0);
[lwdown.years lwdown.months lwdown.days lwdown.hours lwdown.minutes lwdown.seconds] = datevec(lwdown.time);

lwdown.field = readbin([dataDir lwdown.filename],[lwdown.nlon lwdown.nlat length(lwdown.time)],1,'real*4');

saveFilename = 'bin_average_lwdown.mat'

cd(codeDir);
[blon blat bin_average] = compute_bin_average(lwdown.lon,lwdown.lat,lwdown.area,binLonInc,binLatInc,saveDir1,saveFilename);

lwdown.uniqueMonths = unique(lwdown.months);

for i = 1:length(lwdown.uniqueMonths)
    
    im = find(lwdown.months == lwdown.uniqueMonths(i));
    
    meanField = squeeze(nanmean(lwdown.field(:,:,im),3));
    
    lwdown.meanField(:,:,i) = meanField;
    lwdown.binMeanField(:,:,i) = reshape(bin_average*double(meanField(:)),length(blon),length(blat));
    
    clear im meanField
    
end

writebin([saveDir2 lwdown.filename],lwdown.binMeanField,1,'real*4');
save([saveDir3 'lwdown.mat'],'lwdown','blon','blat','-v7.3');

%% 

%  apco2_lon0         = 0.0D0,
%  apco2_lon_inc      = 360.0D0,
%  apco2_lat0        = -89.4628220D0,
%  apco2_lat_inc     = 0.6958694,0.6999817,0.7009048,0.7012634,0.7014313,
%                      245*0.7017418,
%                      0.7014313,0.7012634,0.7009048,0.6999817,0.6958694
%  apco2_nlon         = 2,
%  apco2_nlat         = 256,

%  apco2startdate1  = 19920101,
%  apco2startdate2  = 030000,
%  apco2period      = 86400.0,

apco2.filename = 'apCO2_2000';

apco2.lon0 = 0;
apco2.lon_inc = 0.7031250;
apco2.lat0 = -89.4628220;
apco2.lat_inc = [0.6958694,0.6999817,0.7009048,0.7012634,0.7014313, ...
    (ones(1,245) * 0.7017418), 0.7014313,0.7012634,0.7009048,0.6999817,0.6958694];
apco2.nlon = 512;
apco2.nlat = 256;

apco2.lonVec = [apco2.lon0 (apco2.lon0 + cumsum(ones(1,apco2.nlon-1)*apco2.lon_inc))];
apco2.latVec = [apco2.lat0 (apco2.lat0 + cumsum(apco2.lat_inc))];

[apco2.lon apco2.lat] = meshgrid(apco2.lonVec,apco2.latVec);

apco2.lon = apco2.lon';
apco2.lat = apco2.lat';

apco2.area = rDeg^2*cosd(apco2.lat);

apco2.time = datenum(1992,1,1,0,0,0):1:datenum(1992,12,31,23,0,0);
[apco2.years apco2.months apco2.days apco2.hours apco2.minutes apco2.seconds] = datevec(apco2.time);

apco2.field = readbin([dataDir apco2.filename],[1 apco2.nlat length(apco2.time)],1,'real*4');

saveFilename = 'bin_average_apco2.mat'

cd(codeDir);
[blon blat bin_average] = compute_bin_average(apco2.lon,apco2.lat,apco2.area,binLonInc,binLatInc,saveDir1,saveFilename);

apco2.uniqueMonths = unique(apco2.months);

for i = 1:length(apco2.uniqueMonths)
    
    im = find(apco2.months == apco2.uniqueMonths(i));
    
    meanField = squeeze(nanmean(apco2.field(:,:,im),3));
    meanField = repmat(meanField,[apco2.nlon 1]);
    
    apco2.meanField(:,:,i) = meanField;
    apco2.binMeanField(:,:,i) = reshape(bin_average*double(meanField(:)),length(blon),length(blat));
    
    clear im meanField
    
end

writebin([saveDir2 apco2.filename],apco2.binMeanField,1,'real*4');
save([saveDir3 'apco2.mat'],'apco2','blon','blat','-v7.3');

%% 

%  runofffile        = 'runoff-2d-Fekete-1deg-mon-V4-SMOOTH.bin',
%  runoffperiod      = -12,

runoff.filename = 'runoff-2d-Fekete-1deg-mon-V4-SMOOTH.bin';

runoff.lon = LLC270.XC;
runoff.lat = LLC270.YC;
runoff.area = LLC270.RAC;

runoff.field = readbin([dataDir runoff.filename],[nx ny 12],1,'real*4');

saveFilename = 'bin_average_runoff.mat'

cd(codeDir);
[blon blat bin_average] = compute_bin_average(runoff.lon,runoff.lat,runoff.area,binLonInc,binLatInc,saveDir1,saveFilename);

for i = 1:12
    
    meanField = squeeze(runoff.field(:,:,i));
    
    runoff.meanField(:,:,i) = meanField;
    runoff.binMeanField(:,:,i) = reshape(bin_average*double(meanField(:)),length(blon),length(blat));
    
    clear im meanField
    
end

writebin([saveDir2 runoff.filename],runoff.binMeanField,1,'real*4');
save([saveDir3 'runoff.mat'],'runoff','blon','blat','-v7.3');

%% 

% llc270_Mahowald_2009_soluble_iron_dust.bin

irondust.filename = 'llc270_Mahowald_2009_soluble_iron_dust.bin';

irondust.lon = LLC270.XC;
irondust.lat = LLC270.YC;
irondust.area = LLC270.RAC;

irondust.field = readbin([dataDir irondust.filename],[nx ny 12],1,'real*4');

saveFilename = 'bin_average_irondust.mat'

cd(codeDir);
[blon blat bin_average] = compute_bin_average(irondust.lon,irondust.lat,irondust.area,binLonInc,binLatInc,saveDir1,saveFilename);

for i = 1:12
    
    meanField = squeeze(irondust.field(:,:,i));
    
    irondust.meanField(:,:,i) = meanField;
    irondust.binMeanField(:,:,i) = reshape(bin_average*double(meanField(:)),length(blon),length(blat));
    
    clear im meanField
    
end

writebin([saveDir2 irondust.filename],irondust.binMeanField,1,'real*4');
save([saveDir3 'irondust.mat'],'irondust','blon','blat','-v7.3');

%% 
