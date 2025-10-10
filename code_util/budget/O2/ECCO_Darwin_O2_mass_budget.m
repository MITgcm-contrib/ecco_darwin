clear
close all;

% Based on ECCO-Darwin DIC budget code from 
% https://github.com/MITgcm-contrib/ecco_darwin/blob/master/code_util/budget/DIC_only/fields/mass/v05_biogeochem_DIC_mass_budget_fields.m
% modified by J Koelling, 2025/08/25

%%
%settings, modify as needed

saveMat = 0; %save budget .mat files

useVol = 1; %integrate w/ volume
useLLC270 = 1; %use LLC270 grid

startIntLevel = 1; %vertical integration start k level
endIntLevel = [37 50]; %1900m and botton

nanString = 'omitnan';

%%
gridDir = '~/Documents/MATLAB/ECCO-Darwin/gridDir/';


figureDir = 'figures/';
saveDir = 'mat/';

mkdir figures
mkdir mat

%%
%constants

mmol_to_mol = 10^-3;
deltaT = 1200; % time step in seconds


%%
%load grid

global mygrid
mygrid = [];

grid_load(gridDir,5,'compact');

numLevels = numel(mygrid.RC);

dzMatF = mk3D(mygrid.DRF, mygrid.hFacC);
dzMat = dzMatF .* mygrid.hFacC;
VVV = mygrid.mskC .* mygrid.hFacC .* mk3D(mygrid.RAC,mygrid.mskC) .* mk3D(mygrid.DRF,mygrid.mskC);

vol = convert2gcmfaces(VVV);

%%
%filenames

diagDir = '~/Documents/MATLAB/data/ECCO2/LLC270/ECCO-Darwin_extension/budget/';

filename_2d = 'average_2d';
filename_O2 = 'average_O2_3d';
filename_snap2 = 'snap_2d';
filename_snap3 = 'snap_3d';
filename_O2flx = 'surfO2_tend';

%%

fn = dir([diagDir filename_2d '/' filename_2d '*.*.data']);
tt0 = zeros(length(fn),1);

for i = 1:length(tt0)
    nme = fn(i).name;
    tt0(i) = str2num(nme(end-14:end-5));
end

% Different steps for budget terms and snapshots
tt = tt0(2:end);
ttS = tt0(1:end);

numFiles = numel(tt);
dt = diff(ttS'); 
dt = dt .* deltaT; %seconds
%dt = dt ./ (secPerDay ./ deltaT) .* hoursPerDay; %hours

numTimeSteps = numel(tt);

disp(['Number of timesteps: ' num2str(numTimeSteps)]);

%%

for timeStep = 1:numFiles
    
    disp(num2str(timeStep));
    
    %set file numbers to load in
    ttAverage = tt(timeStep);
    ttSnap = [ttS(timeStep) ttS(timeStep+1)];
        
    
    disp(sprintf('ttA: %i, ttS: %i, %i', ttAverage, ttSnap))
    %two-dimensional time-averaged fields
    %
    ETAN = convert2gcmfaces(rdmds([diagDir filename_2d '/' filename_2d],ttAverage,'rec',1));
    
    %O2 content
    O2 = convert2gcmfaces(rdmds([diagDir filename_O2 '/' filename_O2],ttAverage,'rec',1)) .* mmol_to_mol; %mol m^-3
    ADVx_O2 = convert2gcmfaces(rdmds([diagDir filename_O2 '/' filename_O2],ttAverage,'rec',2)) .* mmol_to_mol; %mol s^-1
    ADVy_O2 = convert2gcmfaces(rdmds([diagDir filename_O2 '/' filename_O2],ttAverage,'rec',3)) .* mmol_to_mol; %mol s^-1
    ADVr_O2 = convert2gcmfaces(rdmds([diagDir filename_O2 '/' filename_O2],ttAverage,'rec',4)) .* mmol_to_mol; %mol s^-1
    DFxE_O2 = convert2gcmfaces(rdmds([diagDir filename_O2 '/' filename_O2],ttAverage,'rec',5)) .* mmol_to_mol; %mol s^-1
    DFyE_O2 = convert2gcmfaces(rdmds([diagDir filename_O2 '/' filename_O2],ttAverage,'rec',6)) .* mmol_to_mol; %mol s^-1
    DFrE_O2 = convert2gcmfaces(rdmds([diagDir filename_O2 '/' filename_O2],ttAverage,'rec',7)) .* mmol_to_mol; %mol s^-1
    DFrI_O2 = convert2gcmfaces(rdmds([diagDir filename_O2 '/' filename_O2],ttAverage,'rec',8)) .* mmol_to_mol; %mol s^-1
    gDAR_O2 = convert2gcmfaces(rdmds([diagDir filename_O2 '/' filename_O2],ttAverage,'rec',9)) .* mmol_to_mol; %mol s^-1
    bioC_O2 = convert2gcmfaces(rdmds([diagDir filename_O2 '/' filename_O2],ttAverage,'rec',10)) .* mmol_to_mol; %mol m^-3 s^-1
    bioS_O2 = convert2gcmfaces(rdmds([diagDir filename_O2 '/' filename_O2],ttAverage,'rec',11)) .* mmol_to_mol; %mol m^-3 s^-1
    fluxO2 = convert2gcmfaces(rdmds([diagDir '../monthly/' filename_O2flx '/' filename_O2flx],ttAverage,'rec',1)) .* mmol_to_mol; %mol m^-3 s^-1
    %%
    %load snapshots
    
    ETAN_SNAP = rdmds([diagDir filename_snap2 '/' filename_snap2],ttSnap,'rec',1);
    O2_SNAP = rdmds([diagDir filename_snap3 '/' filename_snap3],ttSnap,'rec',11) .* mmol_to_mol; %mol m^-3 s^-1
           
    
    ETAN_SNAP = convert2gcmfaces(ETAN_SNAP);
    O2_SNAP = convert2gcmfaces(O2_SNAP);
   
    %%
    %volume budget, s^-1
    
    sStarMean = (1+mk3D(ETAN ./ mygrid.Depth,dzMat));
    sStarSnap = 0 .* O2_SNAP;
    
    for nt = 1:2
        
        sStarSnap(:,:,:,nt) = (1+mk3D(ETAN_SNAP(:,:,nt) ./ mygrid.Depth,dzMat));
        
    end

    %O2 tendency due to air-sea flux, mol m^-3 s^-1
    %Same dimensions as O2 array
    forcO2 = 0 .* O2;
    forcO2(:,:,1) = fluxO2;
    

    gDAR_O2 = (gDAR_O2 - (forcO2)) ./ mygrid.hFacC; %remove air-sea O2 flux from gDAR, so it is just biology
    bioS_O2 = (bioS_O2 - (forcO2)) ./ mygrid.hFacC;

    bioC_O2 = (bioC_O2) ./ mygrid.hFacC;

    forcO2 = mygrid.mskC .* (forcO2); 
    
    %O2 budget
    
    O_snap = 0 .* O2_SNAP;
    
    for nt = 1:2
        
        O_snap(:,:,:,nt) = sStarSnap(:,:,:,nt) .* O2_SNAP(:,:,:,nt);
        
    end
    
    %tendency, mol m^-3 s^-1
    tendO2 = (O_snap(:,:,:,2) - O_snap(:,:,:,1)) ./ (dt(timeStep));
    
    %horizontal divergences, mol m^-3 s^-1
    adv_hConvO2 = calc_UV_conv(ADVx_O2,ADVy_O2) ./ VVV;
    dif_hConvO2 = calc_UV_conv(DFxE_O2,DFyE_O2) ./ VVV;
    
    %vertical divergences, mol m^-3 s^-1
    adv_vConvO2 = 0 .* ADVx_O2;
    dif_vConvO2 = 0 .* ADVx_O2;
    
    for nz = 1:numLevels
        
        nzp1 = min([nz+1,numLevels]);
        
        adv_vConvO2(:,:,nz) = squeeze(ADVr_O2(:,:,nzp1) .* double(nz~=numLevels) - ADVr_O2(:,:,nz));
        dif_vConvO2(:,:,nz) = squeeze(DFrI_O2(:,:,nzp1) .* double(nz~=numLevels) - DFrI_O2(:,:,nz) ...
            + DFrE_O2(:,:,nzp1) .* double(nz~=numLevels) - DFrE_O2(:,:,nz));
        
    end
    
    adv_vConvO2 = adv_vConvO2 ./ VVV;
    dif_vConvO2 = dif_vConvO2 ./ VVV;

    %%
    %convert gcmfaces objects to matrices
    
    sStarMean(mygrid.hFacC == 0) = nan;
    sStar = repmat(sStarMean,[1 1 nz]);
    
    sStarMean = convert2gcmfaces(sStarMean);
    sStarSnap = convert2gcmfaces(sStarSnap);
    
    O2 = convert2gcmfaces(O2);
    tendO2 = convert2gcmfaces(tendO2);
    adv_hConvO2 = convert2gcmfaces(adv_hConvO2);
    adv_vConvO2 = convert2gcmfaces(adv_vConvO2);
    dif_hConvO2 = convert2gcmfaces(dif_hConvO2);
    dif_vConvO2 = convert2gcmfaces(dif_vConvO2);

    gDARO2 = convert2gcmfaces(gDAR_O2);
    forcO2 = convert2gcmfaces(forcO2);
    bioC_O2 = convert2gcmfaces(bioC_O2);
    bioS_O2 = convert2gcmfaces(bioS_O2);
    for ii = 1:length(endIntLevel) % Loop over end int levels
        
    %O2; multiplied by dt to get monthly integral
    intO2 = dt(timeStep).*sum((O2(:,:,startIntLevel:endIntLevel(ii)) .* vol(:,:,startIntLevel:endIntLevel(ii))),3,nanString);
    intTendO2 = dt(timeStep).*sum((tendO2(:,:,startIntLevel:endIntLevel(ii)) .* vol(:,:,startIntLevel:endIntLevel(ii))),3,nanString);
    intHAdvO2 = dt(timeStep).*sum((adv_hConvO2(:,:,startIntLevel:endIntLevel(ii)) .* vol(:,:,startIntLevel:endIntLevel(ii))),3,nanString);
    intVAdvO2 = dt(timeStep).*sum((adv_vConvO2(:,:,startIntLevel:endIntLevel(ii)) .* vol(:,:,startIntLevel:endIntLevel(ii))),3,nanString);
    intHDifO2 = dt(timeStep).*sum((dif_hConvO2(:,:,startIntLevel:endIntLevel(ii)) .* vol(:,:,startIntLevel:endIntLevel(ii))),3,nanString);
    intVDifO2 = dt(timeStep).*sum((dif_vConvO2(:,:,startIntLevel:endIntLevel(ii)) .* vol(:,:,startIntLevel:endIntLevel(ii))),3,nanString);
    intForcO2 = dt(timeStep).*sum((forcO2(:,:,startIntLevel:endIntLevel(ii)) .* vol(:,:,startIntLevel:endIntLevel(ii))),3,nanString);%fluxO2.*mygrid.RAC*10;
    
    intBioC_O2 = dt(timeStep).*sum((bioC_O2(:,:,startIntLevel:endIntLevel(ii)) .* vol(:,:,startIntLevel:endIntLevel(ii))),3,nanString);
    intBioS_O2 = dt(timeStep).*sum((bioS_O2(:,:,startIntLevel:endIntLevel(ii)) .* vol(:,:,startIntLevel:endIntLevel(ii))),3,nanString);
    
    intTotalO2 = intHAdvO2 + intVAdvO2 ...
        + intHDifO2 + intVDifO2 + intBioC_O2 + intBioS_O2;
    
    intResidualO2 = intTendO2 - intTotalO2;

    %%
    %save budget
    
    if saveMat
        
        save([saveDir num2str(startIntLevel) '_' num2str(endIntLevel(ii)) '/'...
            'O2_' num2str(startIntLevel) '_' num2str(endIntLevel(ii)) '_' num2str(tt(timeStep)) '_budget.mat'], ...
            'intTendO2','intHAdvO2','intVAdvO2','intHDifO2','intVDifO2','intTotalO2','intResidualO2', ...
            'intForcO2', 'intBioC_O2', 'intBioS_O2', '-v7.3');
        
    end
    end


    %%
    %clear gcmfaces objects to avoid memory leaks
    
    clear ETAN  ADVx_O2 ADVy_O2 ADVr_O2 DFxE_O2 DFyE_O2 DFrE_O2 DFrI_O2 DICTFLX ...
        ETAN_SNAP O2_SNAP
     
    clear sStarMean sStarSnap
    clear O2 tendO2 adv_hConvO2 adv_vConvO2 dif_hConvO2 dif_vConvO2 virtualFluxO2 gDAR_O2 bioC_O2 bioS_O2
    
    
end
