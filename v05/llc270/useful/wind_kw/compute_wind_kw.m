close all
clear all

gridDir = '/Users/carrolld/Documents/research/carbon/simulations/grid/LLC_270/';

dataDir = '/Users/carrolld/Documents/research/darwin3/mat/wind_kw/';

%% 

numFaces = 13;
nx = 270;
ny = nx .* numFaces;
nz = 50;

dt = 1200;

XC = readbin([gridDir 'XC.data'],[nx ny],1,'real*4');
YC = readbin([gridDir 'YC.data'],[nx ny],1,'real*4');
hFacC = readbin([gridDir 'hFacC.data'],[nx ny nz],1,'real*4');
DRF = readbin([gridDir 'DRF.data'],[nz],1,'real*4');
RAC = readbin([gridDir 'RAC.data'],[nx ny],1,'real*4');

%% 

%load one year of year 2000 3-hourly mean wind speed, wind speed std, mean SST, and
%mean sea-ice fraction
load('windSpeedStats.mat'); 

%scaling = the global estimate of 14C bomb flux. The default (16 cm/hr) is from
%Wanninkhof (2013). But Sweeney et al. (2007) suggest 14.6 cm/hr
        
scaling = 16;

Uavg = meanField1;
Ustd = stdField1;

T = meanField2;
ice = meanField3;

Uavg2 = Uavg .^ 2;
U2 = (Uavg2 + Ustd) .^ 2;

Sc = schmidt_number(T);

Sc660 = (Sc ./ 660) .^ -0.5;

a = calculate_scaled_alpha(U2,Sc660,ice,scaling);

kw = a .* Uavg2 .* Sc660;

%% 

meanKw = nanmean(nanmean(kw .* RAC)) ./ nanmean(nanmean(RAC));

hFig1 = figure(1);
set(hFig1,'units','normalized','outerposition',[0 0 1 1]);
set(gcf,'color',[1 1 1]);
set(gca,'color',[0.5 0.5 0.5]);

colors = cmocean('speed',1000);

lw = 2;
fs = 24;

hold on

set(gca,'Color',[0.5 0.5 0.5]);

quikplot_llc(kw);

colormap(colors);

hcb = colorbar

set(get(hcb,'ylabel'),'String',{'Kw'});

axis tight

box on
set(gca,'LineWidth',lw);
set(gca,'FontSize',fs);

title(['Global Area-weighted Mean Kw = ' num2str(meanKw)]);