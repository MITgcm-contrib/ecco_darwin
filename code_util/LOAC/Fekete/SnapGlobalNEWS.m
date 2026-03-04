% ============================================================
% Compute nutrient & carbon fluxes on LLC grid
% Units of output: mmol m-2 s-1
% ============================================================

clear; close all;

addpath(genpath('/nobackup/dcarrol2/MATLAB'));

% output directory
pout='/nobackup/rsavelli/LOAC/Fekete/bgc_runoff/';

% Load mapping
load River_Fekete_Mapping_Continuous_RadiusLimited_x2.mat

% Load grid & runoff
gridDir = '/nobackup/dcarrol2/LOAC/grid/ECCO_V4r5_raw/';
nx = 90; ny = 1170; nm = 12;

fin = '/nobackup/hzhang1/pub/Release5/input_bin/runoff-2d-Fekete-1deg-mon-V4-SMOOTH_S60scalving_v3.bin';
fekete = readbin(fin,[nx ny nm],1,'real*4');

nR = length(GN2fekete);

% Conversion factors (g → mol)
gP_to_molP  = 0.03228539149637;
gN_to_molN  = 0.071394404106606;
gC_to_molC  = 0.083259093974539;
gSi_to_molSi = 0.03560556158872;

% LOOP OVER ELEMENTS

for f={'DIN','DIP','DON','DOP','DOC','DSi','PN','PP','POC','TSS','DIC'}

    fout = [pout f{1} '_Fekete_ECCO_V4r5.bin'];

    fprintf('\nProcessing %s\n',f{1});

    FLUX = zeros(nx,ny,nm,'double');

    % Loop over rivers
    for r = 1:nR
        
	fprintf('\nProcessing river %d\n',GN2fekete(r).riverID);

        if isempty(GN2fekete(r).LLC90)
            continue
        end

        river_Q = GN2fekete(r).Q_GN;     % km3/yr

        eval(['load_Mg=GN2fekete(r).' f{1} ';'])
    
    	if river_Q <= 0
            continue
        end

        % Convert load to mmol/yr
        Load_g = load_Mg * 1e6;

        switch f{1}
            case {'DIN','DON','PN'}
                mol = Load_g * gN_to_molN;
            case {'DIP','DOP','PP'}
                mol = Load_g * gP_to_molP;
            case {'DOC','POC','DIC'}
                mol = Load_g * gC_to_molC;
            case 'DSi'
                mol = Load_g * gSi_to_molSi;
            case 'TSS'
                continue % skip or define separately
        end

        Load_mmol = mol * 1000;

        % Convert discharge to m3/yr
        Q_m3yr = river_Q * 1e9;

	% Total Fekete discharge for this river
        QF_sum = (GN2fekete(r).Q_FekSum)*1e9;

        % Correction factor
        scale = Q_m3yr / QF_sum;
	
	% River concentration (mmol m-3)
        C_r = Load_mmol / Q_m3yr;

        % Apply to each associated grid cell
        idx = GN2fekete(r).LLC90;

        for m = 1:nm
            runoff_m = fekete(:,:,m);
            runoff_vec = runoff_m(:);
            
	    % Scale runoff so total matches GN discharge
            runoff_scaled = runoff_vec(idx) * scale;

            % mmol m-2 s-1
            flux_local = C_r .* runoff_scaled;

            tmp = FLUX(:,:,m);
            tmp(idx) = tmp(idx) + flux_local;
            FLUX(:,:,m) = tmp;
        end

    end
    writebin(fout, FLUX, 1, 'real*4', 0);
    fprintf('%s written\n',fout);

end
