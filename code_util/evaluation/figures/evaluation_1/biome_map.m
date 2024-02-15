function [] = biome_map(gridType,caption,savePPT)

gridDir = '../data/grid/';

dataDir = '../data/indices/';

%%

fileName = ['biome_map.ppt'];

eval(['delete ' fileName]);

if savePPT
    
    slideNum = exportToPPTX('addslide');
    exportToPPTX('addtext',caption,'HorizontalAlignment','center','VerticalAlignment','middle');
    
end

%%

load([dataDir 'ECCO_darwin_combined_biome_indices_LLC_' num2str(gridType) '.mat']);

%%

numFaces = 13;
nx = gridType;
ny = nx .* numFaces;

input.XC = readbin([gridDir 'XC.data'],[nx ny],1,'real*4');
input.YC = readbin([gridDir 'YC.data'],[nx ny],1,'real*4');
input.RAC = readbin([gridDir 'RAC.data'],[nx ny],1,'real*4');
input.hFacC = readbin([gridDir 'hFacC.data'],[nx ny],1,'real*4');
input.hFacC(input.hFacC == 0) = nan;

lon = -180:0.25:180;
lat = -90:0.25:90;

[xx yy] = meshgrid(lon,lat);
xx = xx';
yy = yy';

output.XC = xx;
output.YC = yy;
output.hFacC = (output.XC .* 0) + 1;

%%

tempColors1 = cbrewer('qual','Accent',20);

colors1 = tempColors1([6 6 7 8 9 10 11],:);

colors1(1,:) = [0.45 0.45 0.45];
colors1(2,:) = [1 1 1];

colors2 = tempColors1([1 1:18],:);
colors2(1,:) = [0.45 0.45 0.45];
colors2(2,:) = [1 1 1];

numColors1 = length(colors1)+1;
numColors2 = length(colors2)+1;

%%

biomeMap1 = (input.hFacC(:,:,1).*0) + 1;
biomeMap2 = (input.hFacC(:,:,1).*0) + 1;

%%

% %superbiomes
% biomeMap1(index.fay18) = 2; %North Hemisphere High Latitudes (NH HL)
% biomeMap1(index.fay19) = 3; %North Hemisphere Subtropics (NH SUB)
% biomeMap1(index.fay20) = 4; %Equatorial (EQU)
% biomeMap1(index.fay21) = 5; %Southern Hemisphere Subtropics (SH SUB)
% biomeMap1(index.fay22) = 6; %Southern Hemisphere High Latitudes (SH HL)
%
% %Fay et al. biomes
% biomeMap2(index.fay1) = 2; %NP ICE, North Pacific Ice
% biomeMap2(index.fay2) = 3; %NP SPSS, North Pacific Subpolar Seasonally Stratified
% biomeMap2(index.fay3) = 4; %NP STSS, North Pacific Subtropical Seasonally Stratified
% biomeMap2(index.fay4) = 5; %NP STPS, North Pacific Subtropical Permanently Stratified
% biomeMap2(index.fay5) = 6; %PEQU-W, West Pacific Equatorial
% biomeMap2(index.fay6) = 7; %PEQU-E, East Pacific Equatorial
% biomeMap2(index.fay7) = 8; %SP STPS, South Pacific Subtropical Permanently Stratified
% biomeMap2(index.fay8) = 9; %NA ICE, North Atlantic Ice
% biomeMap2(index.fay9) = 10; %NA SPSS, North Atlantic Subpolar Seasonally Stratified
% biomeMap2(index.fay10) = 11; %NA STSS, North Atlantic Subtropical Seasonally Stratified
% biomeMap2(index.fay11) = 12; %NA STPS, North Atlantic Subtropical Permanently Stratified
% biomeMap2(index.fay12) = 13; %AEQU, Atlantic Equatorial
% biomeMap2(index.fay13) = 14; %SA STPS, South Atlantic Subtropical Permanently Stratified
% biomeMap2(index.fay14) = 15; %IND STPS, Indian Ocean Subtropical Permanently Stratified
% biomeMap2(index.fay15) = 16; %SO STSS, Southern Ocean Subtropical Seasonally Stratified
% biomeMap2(index.fay16) = 17; %SO SPSS, Southern Ocean Subpolar Seasonally Stratified
% biomeMap2(index.fay17) = 18; %SO ICE, Southern Ocean Ice

%%

%superbiomes
biomeMap1(index.fay18) = 2; %sNH, North Hemisphere High Latitudes
biomeMap1(index.fay19) = 3; %sNS, North Hemisphere Subtropics
biomeMap1(index.fay20) = 4; %sEQ, Equatorial
biomeMap1(index.fay21) = 5; %sSS, Southern Hemisphere Subtropics
biomeMap1(index.fay22) = 6; %sSH, Southern Hemisphere High Latitudes

%biomes
biomeMap2(index.fay1) = 2; %NP1, North Pacific Ice
biomeMap2(index.fay2) = 3; %NP2, North Pacific Subpolar Seasonally Stratified
biomeMap2(index.fay3) = 4; %NP3, North Pacific Subtropical Seasonally Stratified
biomeMap2(index.fay4) = 5; %NP4, North Pacific Subtropical Permanently Stratified
biomeMap2(index.fay5) = 6; %WPE, West Pacific Equatorial
biomeMap2(index.fay6) = 7; %EPE, East Pacific Equatorial
biomeMap2(index.fay7) = 8; %SP, South Pacific Subtropical Permanently Stratified

biomeMap2(index.fay8) = 9; %NA1, North Atlantic Ice
biomeMap2(index.fay9) = 10; %NA2, North Atlantic Subpolar Seasonally Stratified
biomeMap2(index.fay10) = 11; %NA3, North Atlantic Subtropical Seasonally Stratified
biomeMap2(index.fay11) = 12; %NA4, North Atlantic Subtropical Permanently Stratified
biomeMap2(index.fay12) = 13; %AE, Atlantic Equatorial
biomeMap2(index.fay13) = 14; %SA, South Atlantic Subtropical Permanently Stratified

biomeMap2(index.fay14) = 15; %IO, Indian Ocean Subtropical Permanently Stratified

biomeMap2(index.fay15) = 16; %SO1, Southern Ocean Subtropical Seasonally Stratified
biomeMap2(index.fay16) = 17; %SO2, Southern Ocean Subpolar Seasonally Stratified
biomeMap2(index.fay17) = 18; %SO3, Southern Ocean Ice

%%

fay1String = '\bfsNH\rm';
fay2String = '\bfsNS\rm';
fay3String = '\bfsEQ\rm';
fay4String = '\bfsSS\rm';
fay5String = '\bfsSH\rm';

fay6String = '\bfNP1\rm';
fay7String = '\bfNP2\rm';
fay8String = '\bfNP3\rm';
fay9String = '\bfNP4\rm';
fay10String = '\bfWPE\rm';
fay11String = '\bfEPE\rm';
fay12String = '\bfSP\rm';

fay13String = '\bfNA1\rm';
fay14String = '\bfNA2\rm';
fay15String = '\bfNA3\rm';
fay16String = '\bfNA4\rm';
fay17String = '\bfAE\rm';
fay18String = '\bfSA\rm';

fay19String = '\bfIO\rm';

fay20String = '\bfSO1\rm';
fay21String = '\bfSO2\rm';
fay22String = '\bfSO3\rm';

%%

locString1 = {'No Superbiome', ...
    fay1String, ...
    fay2String, ...
    fay3String, ...
    fay4String, ...
    fay5String};

locString2 = {'No Biome', ...
    fay6String, ...
    fay7String, ...
    fay8String, ...
    fay9String, ...
    fay10String, ...
    fay11String, ...
    fay12String, ...
    fay13String, ...
    fay14String, ...
    fay15String, ...
    fay16String, ...
    fay17String, ...
    fay18String, ...
    fay19String, ...
    fay20String, ...
    fay21String, ...
    fay22String};

%%

[ngi ogi] = nninterp(input,output,1,1);

%%

tempField = biomeMap1;

interpField = output.XC * 0;

for i = 1:length(ngi)
    
    if ~isnan(ngi(i)) && ~isnan(ogi(i))
        
        interpField(ngi(i)) = tempField(ogi(i));
        
    else
        
        interpField(ngi(i)) = nan;
        
    end
    
end

%%

hFig1 = figure(1);
set(hFig1,'units','normalized','outerposition',[0 0 1 1]);
set(gcf,'color',[1 1 1]);

fs = 22;
gFs = 14;
cFs = 14;
tFs = 30;
conLw = 1;

lon(lon < 0) = lon(lon < 0) + 360;
[sortLon ix] = sort(lon);
interpField = interpField(ix,:)';

[xx yy] = meshgrid(sortLon,lat);
[LG,LT] = meshgrid(sortLon,lat);

ind=[find(sortLon>113) find(sortLon<=113)];
interpField2 = interpField(:,ind);
LT=LT(:,ind);
LG=LG(:,ind);
LG(LG>113) = LG(LG>113)-360;

interpField2(isnan(interpField2)) = 0;

ccc1 = subplot(121);

hold on

set(gca,'color',[0.65 0.65 0.65]);

m_proj('robinson','lon',[-247 113]);

m_pcolor(LG,LT,interpField2);

shading flat

hold on

[CS CH] = m_contour(LG,LT,interpField2,[2:6],'Color','k','LineWidth',conLw);

set(gca,'color',[0.65 0.65 0.65]);

m_grid('tickdir','out','linewi',2,'FontSize',gFs,'ticklen',0.005,'xlabeldir','middle');

set(gcf,'color','w');

colormap(colors1);

caxis([0 numColors1-2]);

colorBarTick = 1.5:1:7.5;

h = m_contfbar(1.05,[0.3 0.7],CS,1:numColors1-1,'endpiece','none');

set(h,'Ylim',[1 7],'YTick',colorBarTick','YTickLabel',locString1, ...
    'FontSize',cFs+10,'YDir','Reverse');

grid off
box on
set(gca,'LineWidth',2);
set(gca,'FontSize',fs);

title('Superbiomes','FontSize',tFs);

if savePPT
    
    slideNum = exportToPPTX('addslide');
    exportToPPTX('addpicture',hFig1,'Scale','maxfixed');
    
end

%%

tempField = biomeMap2;

interpField = output.XC * 0;

for i = 1:length(ngi)
    
    if ~isnan(ngi(i)) && ~isnan(ogi(i))
        
        interpField(ngi(i)) = tempField(ogi(i));
        
    else
        
        interpField(ngi(i)) = nan;
        
    end
    
end

hFig2 = figure(2);
set(hFig2,'units','normalized','outerposition',[0 0 1 1]);
set(gcf,'color',[1 1 1]);

lon(lon < 0) = lon(lon < 0) + 360;
[sortLon ix] = sort(lon);
interpField = interpField(ix,:)';

[xx yy] = meshgrid(sortLon,lat);
[LG,LT] = meshgrid(sortLon,lat);

ind=[find(sortLon>113) find(sortLon<=113)];
interpField2 = interpField(:,ind);
LT=LT(:,ind);
LG=LG(:,ind);
LG(LG>113) = LG(LG>113)-360;

interpField2(isnan(interpField2)) = 0;

ccc1 = subplot(121);

hold on

set(gca,'color',[0.65 0.65 0.65]);

m_proj('robinson','lon',[-247 113]);

m_pcolor(LG,LT,interpField2);

shading flat

hold on

[CS CH] = m_contour(LG,LT,interpField2,[2:18],'Color','k','LineWidth',conLw);

set(gca,'color',[0.65 0.65 0.65]);

m_grid('tickdir','out','linewi',2,'FontSize',gFs,'ticklen',0.005,'xlabeldir','middle');

set(gcf,'color','w');

colormap(colors2);

caxis([0 numColors2-2]);

colorBarTick = 1.5:1:18.5;

h = m_contfbar(1.05,[0.3 0.7],CS,1:numColors2-1,'endpiece','none');

set(h,'Ylim',[1 19],'YTick',colorBarTick','YTickLabel',locString2, ...
    'FontSize',cFs+4,'YDir','Reverse');

grid off
box on
set(gca,'LineWidth',2);
set(gca,'FontSize',fs);

title('Biomes','FontSize',tFs);

if savePPT
    
    slideNum = exportToPPTX('addslide');
    exportToPPTX('addpicture',hFig2,'Scale','maxfixed');
    
end

close all;

end

%%
