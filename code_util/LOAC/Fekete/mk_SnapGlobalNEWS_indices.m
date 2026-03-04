% Snap GlobalNEWS2 2000 time-mean Qact runoff to Fekete locations
clear, close all

addpath(genpath('/nobackup/dcarrol2/MATLAB'));

% Load LLC90 grid files
gridDir = '/nobackup/dcarrol2/LOAC/grid/ECCO_V4r5_raw/';
nx = 90;
ny = 1170;
nz = 50;
RAC = readbin([gridDir 'RAC.data'],[nx ny],1,'real*4');
XC = readbin([gridDir 'XC.data'],[nx ny],1,'real*4');
YC = readbin([gridDir 'YC.data'],[nx ny],1,'real*4');
hFacC = readbin([gridDir 'hFacC.data'],[nx ny 1],1,'real*4');
wet_mask = hFacC(:,:) > 0;   % surface ocean mask
wet_vec = wet_mask(:);
XC_vec = XC(:);
YC_vec = YC(:);
wet_XC = XC_vec(wet_vec);
wet_YC = YC_vec(wet_vec);
wet_indices = find(wet_vec);

% Load GlobalNEWS2 mouth_lon, mouth_lat, Qact (km3/yr) and end of the basin
% (land or ocean)
gns=xlsread('globalnews');
glon=gns(:,1);  % GlobalNEWS2 longitude (deg);
glat=gns(:,2);  % GlobalNEWS2 latitude (deg)
gQact=gns(:,3); % GlobalNEWS2 actual discharge (km^3/yr)
gDIN=gns(:,4);  % GlobalNEWS2 load DIN (Mg/yr)
gDIP=gns(:,5);  % GlobalNEWS2 load DIP (Mg/yr)
gDON=gns(:,6);  % GlobalNEWS2 load DON (Mg/yr)
gDOP=gns(:,7);  % GlobalNEWS2 load DON (Mg/yr)
gDOC=gns(:,8);  % GlobalNEWS2 load DOC (Mg/yr)
gDSi=gns(:,9);  % GlobalNEWS2 load DSi (Mg/yr)
gPN=gns(:,10);  % GlobalNEWS2 load PN (Mg/yr)
gPP=gns(:,11);  % GlobalNEWS2 load PP (Mg/yr)
gPOC=gns(:,12); % GlobalNEWS2 load POC (Mg/yr)
gTSS=gns(:,13); % GlobalNEWS2 load TSS (Mg/yr)
clear gns ix

% File DIC_final_globalnews.xlsx is derived from 
% computation of DIC fluxes from BDIC_GlobalNEWS.xlsx
gns=xlsread('DIC_final_globalnews.xlsx');
gDIC=gns(:,1); % GlobalNEWS2 load DIC (Mg/yr)
clear gns ix

% Correct DIC inputs from the Amazon
% relation from Li et al., 2017 overestimates DIC export so here average of
% literature (da Cuha et al., 2013, Probst et al, 1994, Li et al 2017)
% 2.54 Tmol yr-1 so 30480000 Mg yr-1
gDIC(1) = 30480000;

% load basins id
load GlobalNEWS_basinsID.mat
gbasins = GlobalNEWS_basinsID;

% Remove rivers ending on land
load GlobalNEWS_basins_end.mat
gbasins_end = GlobalNEWSbasinsend;
ix = strcmp(cellstr(gbasins_end),"Land");
for f={'lat','lon','Qact','DIN','DIP','DON','DOP','DOC','DSi','PN','PP','POC','TSS','DIC','basins'}
    eval(['g' f{1} '(ix) = [];'])
end

% remove rivers ending in Aral or Black Sea
load GlobalNEWS_basins_end2.mat
gbasins_end2 = GlobalNEWSbasinsend2;
ix = find(strcmp(cellstr(gbasins_end2),"Aral Sea") | strcmp(cellstr(gbasins_end2),"Black Sea"));
for f={'lat','lon','Qact','DIN','DIP','DON','DOP','DOC','DSi','PN','PP','POC','TSS','DIC','basins'}
    eval(['g' f{1} '(ix) = [];'])
end

nR = length(gQact);

% sort Global NEWS freshwater discharge from lowest to highest
[gQact,sort_idx]= sort(gQact,'descend');
for f={'lat','lon','DIN','DIP','DON','DOP','DOC','DSi','PN','PP','POC','TSS','DIC','basins'}
    eval(['g' f{1} ' = g' f{1} '(sort_idx);'])
end

% load Fekete monthly climatology and convert to km3/yr
fin = '/nobackup/hzhang1/pub/Release5/input_bin/runoff-2d-Fekete-1deg-mon-V4-SMOOTH_S60scalving_v3.bin';
fekete = readbin(fin,[90 1170 12],1,'real*4');
[nx,ny,nm] = size(fekete);

% Annual Fekete discharge
Fekete_annual = sum(fekete,3).*RAC*30.5*24*60*60/1e9; %km3/yr

fprintf('Total Fekete discharge: %.2f km3/yr\n', sum(Fekete_annual,'all'));
fprintf('Total NEWS discharge:   %.2f km3/yr\n', sum(gQact));


% Structure to store mapping
GN2fekete(nR).LLC90     = [];
GN2fekete(nR).Q_GN      = [];
GN2fekete(nR).Q_FekSum  = [];
GN2fekete(nR).riverID   = [];
GN2fekete(nR).lat   = [];
GN2fekete(nR).lon   = [];
GN2fekete(nR).DIN   = [];
GN2fekete(nR).DIP   = [];
GN2fekete(nR).DON   = [];
GN2fekete(nR).DOP   = [];
GN2fekete(nR).DOC   = [];
GN2fekete(nR).DSi   = [];
GN2fekete(nR).PN   = [];
GN2fekete(nR).PP   = [];
GN2fekete(nR).POC   = [];
GN2fekete(nR).TSS   = [];
GN2fekete(nR).DIC   = [];

% set distance limit for exploration
Qmax = max(gQact);
Rmin = 50e3;      % 50 km
RmaxAmazon = 1000e3;  % 2000 km

%% Loop through rivers
    for r = 1:nR

    river_lon = glon(r);
    river_lat = glat(r);
    river_Q   = gQact(r);
    river_ID = gbasins(r);
    river_DIN   = gDIN(r);
    river_DIP   = gDIP(r);
    river_DON   = gDON(r);
    river_DOP   = gDOP(r);
    river_DOC   = gDOC(r);
    river_DSi   = gDSi(r);
    river_PN   = gPN(r);
    river_PP   = gPP(r);
    river_POC   = gPOC(r);
    river_TSS   = gTSS(r);
    river_DIC   = gDIC(r);

    if river_Q <= 0
        GN2fekete(r).Q_GN      = river_Q;
        GN2fekete(r).riverID   = river_ID;
        GN2fekete(r).lon   = river_lon;
        GN2fekete(r).lat   = river_lat;
        continue
    end

    % Compute radius scaling
    Rmax = Rmin + (RmaxAmazon - Rmin) * (river_Q / Qmax).^0.25;

    % Snap river mouth to closest wet grid cells
    d_snap = haversine_distance(river_lat, river_lon, wet_YC, wet_XC);

    [~, snap_id] = min(d_snap);

    snap_index = wet_indices(snap_id);

    % Replace river location with snapped grid location
    river_lat = YC_vec(snap_index);
    river_lon = XC_vec(snap_index);
    
    % Compute distance from all cells
    d = haversine_distance(river_lat, river_lon, YC, XC);

    d_vec = d(:);
    F_vec = Fekete_annual(:);

    % Ocean mask
    ocean_mask = ~isnan(F_vec);

    valid = wet_vec & (d_vec <= Rmax);

    if ~any(valid)
        fprintf('River %d skipped (no cells in Rmax)\n', river_ID);
        GN2fekete(r).Q_GN      = river_Q;
        GN2fekete(r).riverID   = river_ID;
        GN2fekete(r).lon   = river_lon;
        GN2fekete(r).lat   = river_lat
        continue
    end

    % ---------------------------------------------------------
    % Sort by distance (radial expansion)
    % ---------------------------------------------------------
    [dist_sorted, idx_sorted] = sort(d_vec(valid),'ascend');
    valid_idx = find(valid);
    sorted_idx = valid_idx(idx_sorted);

    cumF = cumsum(F_vec(sorted_idx));

    match_idx = find(cumF >= river_Q,1,'first');

    if isempty(match_idx)
        fprintf('River %d not satisfied within Rmax\n', river_ID);
        match_idx = length(sorted_idx);
    end

    selected_idx = sorted_idx(1:match_idx);
    
    Fsum = sum(F_vec(selected_idx));

    % ---------------------------------------------------------
    % Store
    % ---------------------------------------------------------
    GN2fekete(r).LLC90     = selected_idx;
    GN2fekete(r).Q_GN      = river_Q;
    GN2fekete(r).Q_FekSum  = Fsum;
    GN2fekete(r).riverID   = river_ID;
    GN2fekete(r).lon   = river_lon;
    GN2fekete(r).lat   = river_lat;
    
    GN2fekete(r).DIN   = river_DIN;
    GN2fekete(r).DIP   = river_DIP;
    GN2fekete(r).DON   = river_DON;
    GN2fekete(r).DOP   = river_DOP;
    GN2fekete(r).DOC   = river_DOC;
    GN2fekete(r).DSi   = river_DSi;
    GN2fekete(r).PN   = river_PN;
    GN2fekete(r).PP   = river_PP;
    GN2fekete(r).POC   = river_POC;
    GN2fekete(r).TSS   = river_TSS;
    GN2fekete(r).DIC   = river_DIC;
end


fprintf('Assigned Fekete discharge: %.2f km3/yr\n', sum([GN2fekete.Q_FekSum]));
fprintf('Assigned NEWS discharge:   %.2f km3/yr\n', sum([GN2fekete.Q_GN]));

save('River_Fekete_Mapping_Continuous_RadiusLimited_x2.mat',...
     'GN2fekete','-v7.3');

