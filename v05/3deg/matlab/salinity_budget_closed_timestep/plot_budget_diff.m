clear
close all

modelDir = '../../../darwin3/run/';
dataDir = '/Users/carrolld/Documents/research/v05_budget/mat/salinity_budget/';

B1 = load([dataDir 'budget_closed.mat']);
B2 = load([dataDir 'budget_equation_12.mat']);

%%
%plotting params.

fs = 12;
lw = 2;

colors = flipud(cbrewer('div','RdBu',500));

%%
%constants

secPerDay = 86400;
secPerHour = 3600;
hoursPerDay = 24;
deltaT = 3600;

rhoConst = 1029;
mmol_to_mol = 1 ./ 1000;

%%
%load grid

nx = 128;
ny = 64;
nz = 15;

XC = readbin([modelDir 'XC.data'],[nx ny],1,'real*8');
YC = readbin([modelDir 'YC.data'],[nx ny],1,'real*8');
hFacC = readbin([modelDir 'hFacC.data'],[nx ny nz],1,'real*8');
RAC = readbin([modelDir 'RAC.data'],[nx ny],1,'real*8');
DXG = readbin([modelDir 'DXG.data'],[nx ny],1,'real*8');
DYG = readbin([modelDir 'DYG.data'],[nx ny],1,'real*8');
DRF = readbin([modelDir 'dRF.data'],[nz],1,'real*8');
RC = readbin([modelDir 'RC.data'],[nz],1,'real*8');

numLevels = numel(RC);

depth = readbin([modelDir 'Depth.data'],[nx ny],1,'real*8');

dzMatF = repmat(DRF,1,nx,ny);
dzMatF = permute(dzMatF,[2 3 1]);
dzMat = dzMatF .* hFacC;

RACMat = repmat(RAC,1,1,nz);

mskC = hFacC .* 1;
mskC(hFacC == 0) = nan;
mskC(~isnan(mskC)) = 1;

VVV = mskC .* RACMat .* dzMatF .* hFacC;

%% 

intTendSal = B1.budget.intTendSal - B2.budget.intTendSal;

intHAdvSal = B1.budget.intHAdvSal - B2.budget.intHAdvSal;
intVAdvSal = B1.budget.intVAdvSal - B2.budget.intVAdvSal;

intHDifSal = B1.budget.intHDifSal -B2.budget.intHDifSal;
intVDifSal = B1.budget.intVDifSal - B2.budget.intVDifSal;
intForcSal = B1.budget.intForcSal - B2.budget.intForcSal;

%% 

hFig1 = figure(1);
set(hFig1,'units','normalized','outerposition',[0 0 1 0.5]);
set(gcf,'color',[1 1 1]);

cMin = -10^-6;
cMax = 10^-6;

timeStep = 2;

subplot(181);

hold on

set(gca,'color',[0.5 0.5 0.5]);

pcolor(XC,YC,intTendSal);
shading flat

caxis([cMin cMax]);

colormap(colors);

colorbar

axis tight

xlim([0 360]);
ylim([-90 90]);

set(gca,'xtick',[0:180:360]);
set(gca,'ytick',[-90:30:90]);

xlabel('Longitude (deg)');
ylabel('Latitude (deg)');

box on
set(gca,'LineWidth',lw);
set(gca,'FontSize',fs);

title({'Salinity Tendency (psu m s^-^1)'},'FontWeight','Bold','FontSize',fs);

subplot(1,8,2);

hold on

set(gca,'color',[0.5 0.5 0.5]);

pcolor(XC,YC,intHAdvSal);
shading flat

caxis([cMin cMax]);

colormap(colors);

colorbar

axis tight

xlim([0 360]);
ylim([-90 90]);

set(gca,'xtick',[0:180:360]);
set(gca,'ytick',[-90:30:90]);

xlabel('Longitude (deg)');

box on
set(gca,'LineWidth',lw);
set(gca,'FontSize',fs);

title({'Horizontal';'Advection';'(psu m s^-^1)'},'FontWeight','Bold','FontSize',fs);

subplot(1,8,3);

hold on

set(gca,'color',[0.5 0.5 0.5]);

pcolor(XC,YC,intVAdvSal);
shading flat

caxis([cMin cMax]);

colormap(colors);

colorbar

axis tight

xlim([0 360]);
ylim([-90 90]);

set(gca,'xtick',[0:180:360]);
set(gca,'ytick',[-90:30:90]);

xlabel('Longitude (deg)');

box on
set(gca,'LineWidth',lw);
set(gca,'FontSize',fs);

title({'Vertical';'Advection';'(psu m s^-^1)'},'FontWeight','Bold','FontSize',fs);

subplot(1,8,4);

hold on

set(gca,'color',[0.5 0.5 0.5]);

pcolor(XC,YC,intHDifSal);
shading flat

caxis([cMin cMax]);

colormap(colors);

colorbar

axis tight

xlim([0 360]);
ylim([-90 90]);

set(gca,'xtick',[0:180:360]);
set(gca,'ytick',[-90:30:90]);

xlabel('Longitude (deg)');

box on
set(gca,'LineWidth',lw);
set(gca,'FontSize',fs);

title({'Horizontal';'Diffusion';'(psu m s^-^1)'},'FontWeight','Bold','FontSize',fs);

subplot(1,8,5);

hold on

set(gca,'color',[0.5 0.5 0.5]);

pcolor(XC,YC,intVDifSal);
shading flat

caxis([cMin cMax]);

colormap(colors);

colorbar

axis tight

xlim([0 360]);
ylim([-90 90]);

set(gca,'xtick',[0:180:360]);
set(gca,'ytick',[-90:30:90]);

xlabel('Longitude (deg)');

box on
set(gca,'LineWidth',lw);
set(gca,'FontSize',fs);

title({'Vertical';'Diffusion';'(psu m s^-^1)'},'FontWeight','Bold','FontSize',fs);

subplot(1,8,6);

hold on

set(gca,'color',[0.5 0.5 0.5]);

pcolor(XC,YC,intForcSal);
shading flat

caxis([cMin cMax]);

colormap(colors);

colorbar

axis tight

xlim([0 360]);
ylim([-90 90]);

set(gca,'xtick',[0:180:360]);
set(gca,'ytick',[-90:30:90]);

xlabel('Longitude (deg)');

box on
set(gca,'LineWidth',lw);
set(gca,'FontSize',fs);

title({'Surface';'Volume';'Forcing';'(psu m s^-^1)'},'FontWeight','Bold','FontSize',fs);

subplot(1,8,7);

hold on

set(gca,'color',[0.5 0.5 0.5]);

pcolor(XC,YC,(intHAdvSal + intVAdvSal + intHDifSal + intVDifSal + intForcSal));
shading flat

caxis([cMin cMax]);

colormap(colors);

colorbar

axis tight

xlim([0 360]);
ylim([-90 90]);

set(gca,'xtick',[0:180:360]);
set(gca,'ytick',[-90:30:90]);

xlabel('Longitude (deg)');

box on
set(gca,'LineWidth',lw);
set(gca,'FontSize',fs);

title({'Budget';'Total';'(psu m s^-^1)'},'FontWeight','Bold','FontSize',fs);

subplot(1,8,8);

hold on

set(gca,'color',[0.5 0.5 0.5]);

pcolor(XC,YC,intTendSal - (intHAdvSal + intVAdvSal + intHDifSal + intVDifSal + intForcSal));
shading flat

caxis([cMin cMax]);

colormap(colors);

colorbar

axis tight

xlim([0 360]);
ylim([-90 90]);

set(gca,'xtick',[0:180:360]);
set(gca,'ytick',[-90:30:90]);

xlabel('Longitude (deg)');

box on
set(gca,'LineWidth',lw);
set(gca,'FontSize',fs);

title({'Budget Total - Tend';'(psu m s^-^1)'},'FontWeight','Bold','FontSize',fs);

drawnow

%%

hFig2 = figure(2);
set(hFig2,'units','normalized','outerposition',[0 0 1 1]);
set(gcf,'color',[1 1 1]);

cMin = -10^-4;
cMax = 10^-4;

term1 = B1.budget.test;
term2 = B2.budget.test;

term3 = term1 - term2;

timeStep = 2;

subplot(131);

hold on

set(gca,'color',[0.5 0.5 0.5]);

pcolor(XC,YC,term1);
shading flat

caxis([cMin cMax]);

colormap(colors);

colorbar

axis tight

xlim([0 360]);
ylim([-90 90]);

set(gca,'xtick',[0:180:360]);
set(gca,'ytick',[-90:30:90]);

xlabel('Longitude (deg)');

box on
set(gca,'LineWidth',lw);
set(gca,'FontSize',fs);

subplot(132);

hold on

set(gca,'color',[0.5 0.5 0.5]);

pcolor(XC,YC,term2);
shading flat

caxis([cMin cMax]);

colormap(colors);

colorbar

axis tight

xlim([0 360]);
ylim([-90 90]);

set(gca,'xtick',[0:180:360]);
set(gca,'ytick',[-90:30:90]);

xlabel('Longitude (deg)');

box on
set(gca,'LineWidth',lw);
set(gca,'FontSize',fs);

subplot(133);

hold on

set(gca,'color',[0.5 0.5 0.5]);

pcolor(XC,YC,term3);
shading flat

caxis([cMin cMax]);

colormap(colors);

colorbar

axis tight

xlim([0 360]);
ylim([-90 90]);

set(gca,'xtick',[0:180:360]);
set(gca,'ytick',[-90:30:90]);

xlabel('Longitude (deg)');

box on
set(gca,'LineWidth',lw);
set(gca,'FontSize',fs);
