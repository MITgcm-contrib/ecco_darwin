%% =========================================================================
%  OAEMIP Gaussian Alkalinity Injection Forcing Generator for ECCO-Darwin
%
%  This script generates daily Ocean Alkalinity Enhancement (OAE) forcing
%  files for the ECCO v4r5 1° grid.
%
%  The injection follows the OAEMIP protocol:
%    - Gaussian distribution centered on a specified location
%    - Standard deviation = sigma_km
%    - Cutoff radius = 3 × sigma
%    - Land cells receive no alkalinity
%    - Remaining ocean cells are renormalized so that the total injected
%      alkalinity equals the prescribed flux (mol s-1)
%
%  Output:
%    Binary forcing files (lon × lat × day)
%    Units: mmol m-2 s-1
%
% =========================================================================

clear all
close all

%% -------------------------------------------------------------------------
%  Add MATLAB toolboxes
% -------------------------------------------------------------------------

addpath /Users/rsavelli/Documents/MATLAB/dmenem/
addpath /Users/rsavelli/Documents/MATLAB/

p = genpath('/Users/rsavelli/Documents/MATLAB/gcmfaces/'); addpath(p);
p = genpath('/Users/rsavelli/Documents/MATLAB/m_map/'); addpath(p);
grid_load;
gcmfaces_global;
%% -------------------------------------------------------------------------
%  User-defined injection parameters
% -------------------------------------------------------------------------
plot_control = 1; % quality control plots

% Total alkalinity addition (mol s-1)
totalFlux = 3.17e4;

% Injection period
start_date = '2003-01-01';
stop_date  = '2003-12-31';

% Injection location
lat = 26.5;
lon = 134.5;

% Gaussian width (standard deviation)
sigma_km = 250;

% Injection radius (OAEMIP = 3 sigma)
cutoff = 3*sigma_km;

%% -------------------------------------------------------------------------
%  Directories
% -------------------------------------------------------------------------

gridDir  = '/Users/rsavelli/Documents/Models/grid/ECCO_V4r5/';
writeDir = '/Users/rsavelli/Documents/OAEMIP/forcings/test/ECCO_V4r5/';

%% -------------------------------------------------------------------------
%  Load ECCO grid
% -------------------------------------------------------------------------

numFaces = 13;
nx = 90;
ny = nx*numFaces;
nz = 50;

% Horizontal cell area (m2)
RAC = readbin([gridDir 'RAC.data'],[nx ny],1,'real*4');

% Ocean depth (m)
Depth = readbin([gridDir 'Depth.data'],[nx ny],1,'real*4');

% Longitude and latitude
XC = readbin([gridDir 'XC.data'],[nx ny],1,'real*4');
YC = readbin([gridDir 'YC.data'],[nx ny],1,'real*4');

% Ocean fraction (surface layer)
hFacC = readbin([gridDir 'hFacC.data'],[nx ny],1,'real*4');

% Ocean mask
ocean = hFacC > 0;

%% -------------------------------------------------------------------------
%  Compute great-circle distance from every grid cell to injection center
% -------------------------------------------------------------------------

R = 6371;                 % Earth radius (km)

% Wrapped longitude difference
dlon = deg2rad(XC-lon);
dlon = atan2(sin(dlon),cos(dlon));

dlat = deg2rad(YC-lat);

lat1 = deg2rad(lat);
lat2 = deg2rad(YC);

% Haversine formula
a = sin(dlat/2).^2 + ...
    cos(lat1).*cos(lat2).*sin(dlon/2).^2;

dist = 2*R*asin(sqrt(a));

%% -------------------------------------------------------------------------
%  Build Gaussian footprint
% -------------------------------------------------------------------------

% Classical Gaussian
raw = exp(-0.5*(dist/sigma_km).^2);

% Remove land cells
raw(~ocean) = 0;

% Remove cells outside 3 sigma
raw(dist > cutoff) = 0;

%% -------------------------------------------------------------------------
%  Normalize Gaussian
%
%  The normalization guarantees:
%
%      sum(Flux * Area) = totalFlux
%
%  even if land cells were removed.
% -------------------------------------------------------------------------

weightedArea = sum(raw(:).*RAC(:));

% Flux in mol m-2 s-1
flux_mol = raw * totalFlux / weightedArea;

% Convert to mmol m-2 s-1 (required by ECCO-Darwin)
flux_mmol = 1000*flux_mol;

if plot_control
%quick plot for quality control
figure
test = convert2gcmfaces(flux_mmol);
[X,Y,FLD]=convert2pcol(mygrid.XC, mygrid.YC, test);
m_proj('robinson','lon',[-179.9 179.9]);
m_pcolor(X,Y,FLD); hold on
m_coast('patch',[.7 .7 .7],'edgecolor','none');
m_grid('tickdir','out','linewi',2,'linestyle','none');

cmap = brewermap([],'Blues');                                                                                                                                                                
colormap(cmap);
h=colorbar('southoutside');
title('ALK Injection');
set(get(h,'ylabel'),'String','Injection (mmol m^{-2} s^{-1})');
set(gcf,'color','w');   % Need to do this otherwise 'print' turns the lakes black
set(gca,'FontSize',20) % Creates an axes and sets its FontSize to 18

end

%% -------------------------------------------------------------------------
%  Generate forcing files
%
%  One binary file is produced for each year.
%  The forcing is constant throughout each day.
% -------------------------------------------------------------------------

% Years containing OAE forcing
for year = 1992:2025

    disp(year)

    ndays = double(days(datetime(year+1,1,1)-datetime(year,1,1)));

    ALKforcing = repmat(flux_mmol,[1 1 ndays]);

    % Year injection (Tmol yr-1)
    sum(ALKforcing.*RAC,'all')*86400*1e-15

    writebin([writeDir ...
        'nonzero/Ainjection_ECCO_V4r5_' num2str(year)], ...
        ALKforcing,1,'real*4');

end

%% -------------------------------------------------------------------------
%  Generate zero-forcing files
%
%  Useful before/after the injection period.
% -------------------------------------------------------------------------

for year = 1992:2025

    disp(year)

    ndays = double(days(datetime(year+1,1,1)-datetime(year,1,1)));

    ALKforcing = zeros(nx,ny,ndays);

    writebin([writeDir ...
        'zero/Ainjection0_ECCO_V4r5_' num2str(year)], ...
        ALKforcing,1,'real*4');

end

%% -------------------------------------------------------------------------
%  Plot Gaussian transect for quality control
%
%  Extract the smallest rectangle containing the Gaussian footprint and
%  plot a zonal transect through its centre.
% -------------------------------------------------------------------------

[row,col] = find(raw>0);

imin = min(row);
imax = max(row);

jmin = min(col);
jmax = max(col);

dist_samp = dist(imin:imax,jmin:jmax);

ALKforcing_samp = flux_mmol(imin:imax,jmin:jmax);

% Row passing through the centre
row = 8;

signedDist = zeros(size(dist_samp(row,:)));

centerCol = find(dist_samp(row,:)==min(dist_samp(row,:)),1);

% Create signed distance (negative west, positive east)
signedDist(1:centerCol) = -dist_samp(row,1:centerCol);
signedDist(centerCol:end) = dist_samp(row,centerCol:end);

if plot_control
figure
plot(signedDist,ALKforcing_samp(row,:),'LineWidth',2)

xlabel('Distance from injection centre (km)')
ylabel('Alkalinity flux (mmol m^{-2} s^{-1})')

title('Gaussian OAE injection profile')

grid on
end
