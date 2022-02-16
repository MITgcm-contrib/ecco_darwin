clear
close all;

%%
%settings, modify as needed

runOnline = 1;
savePlot = 1;
saveMat = 1; %save budget .mat files

useVol = 1; %integrate w/ volume
useLLC270 = 1; %use LLC270 grid

startIntLevel = 1; %vertical integration start k level
endIntLevel = 10;

nanString = 'omitnan';

%%

if runOnline
    
    %addpath(genpath('/nobackup/dcarrol2/MATLAB'));
    
    gridDir = '/nobackup/dcarrol2/grid/LLC_270/';
    modelDir = '/nobackup/dcarrol2/v05_latest/darwin3/run/'; 

else
    
    gridDir = '../../../../../darwin3/run/';
    modelDir = '../../../../../darwin3/run/';
    
end

figureDir = 'figures/';
saveDir = 'mat/';

mkdir figures
mkdir mat

%%
%plotting parameters

fs = 12;
lw = 2;

colors = cmocean('balance',1000);

%%
%constants

secPerDay = 86400;
secPerHour = 3600;
hoursPerDay = 24;

rhoConst = 1029;
mmol_to_mol = 10^-3;
umol_to_mol = 10^-6;

if useLLC270
    
    deltaT = 1200;
    
else
    
    deltaT = 3600;
    
end

%%
%load grid

global mygrid
mygrid = [];

if useLLC270
    
    grid_load(gridDir,5,'compact');
    plotString = 'quikplot_llc';
    flipString = [];
else
    
    plotString = 'pcolor';
    
    nF = 1;
    fileFormat = 'straight';
    
    grid_load(gridDir,nF,fileFormat);
    
    mygrid.domainPeriodicity = [1 0];
    
    flipString = '''';
    
end

nx = mygrid.ioSize(1);
ny = mygrid.ioSize(2);
numLevels = numel(mygrid.RC);

dzMatF = mk3D(mygrid.DRF, mygrid.hFacC);
dzMat = dzMatF .* mygrid.hFacC;
RACMat = mk3D(mygrid.RAC, mygrid.hFacC);
VVV = mygrid.mskC .* mygrid.hFacC .* mk3D(mygrid.RAC,mygrid.mskC) .* mk3D(mygrid.DRF,mygrid.mskC);
VVV2 = mygrid.mskC .* mk3D(mygrid.RAC,mygrid.mskC) .* mk3D(mygrid.DRF,mygrid.mskC);

vol = convert2gcmfaces(VVV);

if ~useVol
    
    vol(:) = 1;
    
end

%%
%filenames

diagDir = [modelDir 'diags/budget/'];

filename1 = 'average_2d';
filename2 = 'average_velmass_3d';
filename3 = 'average_salt_3d';
filename4 = 'average_DIC_3d';
filename5 = 'average_ALK_3d';
filename6 = 'average_NO3_3d';
filename7 = 'average_NO2_3d';
filename8 = 'average_NH4_3d';
filename9 = 'average_PO4_3d';
filename10 = 'average_Fe_3d';
filename11 = 'average_Fe_darwin_2d';
filename12 = 'average_SiO2_3d';
filename13 = 'snap_2d';
filename14 = 'snap_3d';

%%

fn = dir([diagDir filename1 '.*.data']);
tt = zeros(length(fn),1);

for i = 1:length(tt)
    
    nme = fn(i).name;
    tt(i) = str2num(nme(end-14:end-5));
    
end

numFiles = numel(tt);
dt = diff([0 tt']);
dt = dt ./ (secPerDay ./ deltaT) .* hoursPerDay; %hours

numTimeSteps = numel(tt);

disp(['Number of timesteps: ' num2str(numTimeSteps)]);

%%

for timeStep = 1:numFiles
    
    disp(num2str(timeStep));
    
    %set file numbers to load in
    if timeStep ==1
        
        ttAverage = tt(timeStep);
        ttSnap = [1 tt(timeStep)];
        
    else %timeStep~=1
        
        ttAverage = tt(timeStep);
        ttSnap = [tt(timeStep-1) tt(timeStep)];
        
    end
    
    %two-dimensional time-averaged fields
    %
    ETAN = convert2gcmfaces(rdmds([diagDir filename1],ttAverage,'rec',1));
    oceFWflx = convert2gcmfaces(rdmds([diagDir filename1],ttAverage,'rec',2));
    SFLUX = convert2gcmfaces(rdmds([diagDir filename1],ttAverage,'rec',3));
    oceSPflx = convert2gcmfaces(rdmds([diagDir filename1],ttAverage,'rec',4));
    DICTFLX = convert2gcmfaces(rdmds([diagDir filename1],ttAverage,'rec',5)) .* mmol_to_mol; %mol m^-3 s^-1
    fluxCO2 = convert2gcmfaces(rdmds([diagDir filename1],ttAverage,'rec',6)) .* mmol_to_mol; %mol m^-2 s^-1
    
    %two-dimensional surface forcing
    DIC_Epr = convert2gcmfaces(rdmds([diagDir filename1],ttAverage,'rec',7)) .* mmol_to_mol; %mol m^-3 s^-1
   
    %three-dimensional time-averaged fields
    %
    UVELMASS = convert2gcmfaces(rdmds([diagDir filename2],ttAverage,'rec',1));
    VVELMASS = convert2gcmfaces(rdmds([diagDir filename2],ttAverage,'rec',2));
    WVELMASS = convert2gcmfaces(rdmds([diagDir filename2],ttAverage,'rec',3));
    
    %salt content
    SALT = convert2gcmfaces(rdmds([diagDir filename3],ttAverage,'rec',1));
    ADVx_SLT = convert2gcmfaces(rdmds([diagDir filename3],ttAverage,'rec',2));
    ADVy_SLT = convert2gcmfaces(rdmds([diagDir filename3],ttAverage,'rec',3));
    ADVr_SLT = convert2gcmfaces(rdmds([diagDir filename3],ttAverage,'rec',4));
    DFxE_SLT = convert2gcmfaces(rdmds([diagDir filename3],ttAverage,'rec',5));
    DFyE_SLT = convert2gcmfaces(rdmds([diagDir filename3],ttAverage,'rec',6));
    DFrE_SLT = convert2gcmfaces(rdmds([diagDir filename3],ttAverage,'rec',7));
    DFrI_SLT = convert2gcmfaces(rdmds([diagDir filename3],ttAverage,'rec',8));
    oceSPtnd = convert2gcmfaces(rdmds([diagDir filename3],ttAverage,'rec',9));
    
    %DIC content
    DIC = convert2gcmfaces(rdmds([diagDir filename4],ttAverage,'rec',1)) .* mmol_to_mol; %mol m^-3
    ADVx_DIC = convert2gcmfaces(rdmds([diagDir filename4],ttAverage,'rec',2)) .* mmol_to_mol; %mol s^-1
    ADVy_DIC = convert2gcmfaces(rdmds([diagDir filename4],ttAverage,'rec',3)) .* mmol_to_mol; %mol s^-1
    ADVr_DIC = convert2gcmfaces(rdmds([diagDir filename4],ttAverage,'rec',4)) .* mmol_to_mol; %mol s^-1
    DFxE_DIC = convert2gcmfaces(rdmds([diagDir filename4],ttAverage,'rec',5)) .* mmol_to_mol; %mol s^-1
    DFyE_DIC = convert2gcmfaces(rdmds([diagDir filename4],ttAverage,'rec',6)) .* mmol_to_mol; %mol s^-1
    DFrE_DIC = convert2gcmfaces(rdmds([diagDir filename4],ttAverage,'rec',7)) .* mmol_to_mol; %mol s^-1
    DFrI_DIC = convert2gcmfaces(rdmds([diagDir filename4],ttAverage,'rec',8)) .* mmol_to_mol; %mol s^-1
    gDAR_DIC = convert2gcmfaces(rdmds([diagDir filename4],ttAverage,'rec',9)) .* mmol_to_mol; %mol m^-3 s^-1
    
    %bio decomposition
    cDIC = convert2gcmfaces(rdmds([diagDir filename4],ttAverage,'rec',10)) .* mmol_to_mol; %mol m^-3 s^-1
    cDIC_PIC = convert2gcmfaces(rdmds([diagDir filename4],ttAverage,'rec',11)) .* mmol_to_mol; %mol m^-3 s^-1
    respDIC = convert2gcmfaces(rdmds([diagDir filename4],ttAverage,'rec',12)) .* mmol_to_mol; %mol m^-3 s^-1
    rDIC_DOC = convert2gcmfaces(rdmds([diagDir filename4],ttAverage,'rec',13)) .* mmol_to_mol; %mol m^-3 s^-1
    rDIC_POC = convert2gcmfaces(rdmds([diagDir filename4],ttAverage,'rec',14)) .* mmol_to_mol; %mol m^-3 s^-1
    dDIC_PIC = convert2gcmfaces(rdmds([diagDir filename4],ttAverage,'rec',15)) .* mmol_to_mol; %mol m^-3 s^-1
     
    %%
    %load snapshots
    
    ETAN_SNAP = nan*ones(nx,ny,2);
    SALT_SNAP = nan .* ones(nx,ny,numLevels,2);
    DIC_SNAP = nan .* ones(nx,ny,numLevels,2);
    
    if timeStep == 1
        
        %no pre-initial snapshot
        ETAN_SNAP(:,:,2) = rdmds([diagDir filename13],tt(1),'rec',1);
        
        SALT_SNAP(:,:,:,2) = rdmds([diagDir filename14],tt(1),'rec',2);
        DIC_SNAP(:,:,:,2) = rdmds([diagDir filename14],tt(1),'rec',3) .* mmol_to_mol; %mol m^-3 s^-1
        
    elseif timeStep == numFiles %no final snapshot
        
        ETAN_SNAP(:,:,1) = rdmds([diagDir filename13],tt(timeStep-1),'rec',1);
        
        SALT_SNAP(:,:,:,1) = rdmds([diagDir filename14],tt(timeStep-1),'rec',2);
        DIC_SNAP(:,:,:,1) = rdmds([diagDir filename14],tt(timeStep-1),'rec',3) .* mmol_to_mol; %mol m^-3 s^-1
         
    else %timeStep~=1 & timeStep~=numFiles
        
        ETAN_SNAP = rdmds([diagDir filename13],ttSnap,'rec',1);
        
        SALT_SNAP = rdmds([diagDir filename14],ttSnap,'rec',2);
        DIC_SNAP = rdmds([diagDir filename14],ttSnap,'rec',3) .* mmol_to_mol; %mol m^-3 s^-1
           
    end
    
    ETAN_SNAP = convert2gcmfaces(ETAN_SNAP);
    SALT_SNAP = convert2gcmfaces(SALT_SNAP);
    DIC_SNAP = convert2gcmfaces(DIC_SNAP);
   
    %%
    %volume budget, s^-1
    
    sStarMean = (1+mk3D(ETAN ./ mygrid.Depth,dzMat));
    sStarSnap = 0 .* SALT_SNAP;
    
    for nt = 1:2
        
        sStarSnap(:,:,:,nt) = (1+mk3D(ETAN_SNAP(:,:,nt) ./ mygrid.Depth,dzMat));
        
    end

    %total tendency
    tendV = (1 ./ mk3D(mygrid.Depth,mygrid.mskC)) .* mk3D((ETAN_SNAP(:,:,2) - ETAN_SNAP(:,:,1)) ...
        ./ (secPerHour .* dt(timeStep)),mygrid.mskC);
    
    %horizontal divergence
    adv_hConvV = mygrid.mskC .* calc_UV_conv(UVELMASS,VVELMASS,{'dh'}) ./ (RACMat.*mygrid.hFacC);
    
    %vertical divergence
    adv_vConvV = 0*adv_hConvV;
    
    for nz = 1:numLevels
        
        nzp1 = min([nz+1,numLevels]);
        
        adv_vConvV(:,:,nz) = squeeze(WVELMASS(:,:,nzp1) .* double(nz~=numLevels) - WVELMASS(:,:,nz) .* double(nz~=1)) ...
            ./ (dzMat(:,:,nz));
        
    end
    
    %surface forcing
    forcV = mygrid.mskC .* mk3D(oceFWflx,mygrid.mskC) ./ (dzMat .* rhoConst);
    forcV(:,:,2:numLevels) = 0 .* mygrid.mskC(:,:,2:numLevels);
    
    %%
    %salt budget, psu s^-1
    
    S_snap = 0 .* SALT_SNAP;
    
    for nt = 1:2
        
        S_snap(:,:,:,nt) = sStarSnap(:,:,:,nt) .* SALT_SNAP(:,:,:,nt);
        
    end
    
    tendS = (S_snap(:,:,:,2) - S_snap(:,:,:,1)) ./ (secPerHour .* dt(timeStep));
    
    %horizontal divergences
    adv_hConvS = calc_UV_conv(ADVx_SLT,ADVy_SLT) ./ VVV;
    dif_hConvS = calc_UV_conv(DFxE_SLT,DFyE_SLT) ./ VVV;
    
    %vertical divergences
    adv_vConvS = 0 .* ADVx_SLT;
    dif_vConvS = 0 .* ADVx_SLT;
    
    for nz = 1:numLevels
        
        nzp1 = min([nz+1,numLevels]);
        
        adv_vConvS(:,:,nz) = squeeze(ADVr_SLT(:,:,nzp1) .* double(nz~=numLevels) - ADVr_SLT(:,:,nz));
        dif_vConvS(:,:,nz) = squeeze(DFrI_SLT(:,:,nzp1) .* double(nz~=numLevels) - DFrI_SLT(:,:,nz) ...
            + DFrE_SLT(:,:,nzp1) .* double(nz~=numLevels) - DFrE_SLT(:,:,nz));
        
    end
    
    adv_vConvS = adv_vConvS ./ VVV;
    dif_vConvS = dif_vConvS ./ VVV;
    
    %surface salt flux
    forcS = 0 .* oceSPtnd;
    
    for nz = 1:numLevels
        
        if nz == 1
            
            forcS(:,:,nz) = SFLUX ./ rhoConst;
            
        end
        
        forcS(:,:,nz) = forcS(:,:,nz) + (oceSPtnd(:,:,nz) ./ rhoConst);
        
    end
    
    forcS = forcS ./ dzMat;
    
    %%
    %salinity budget
    
    %tendency
    tendSal = mygrid.mskC .* (SALT_SNAP(:,:,:,2) - SALT_SNAP(:,:,:,1)) ...
        ./ (secPerHour .* dt(timeStep));
    
    %advection
    adv_hConvSal = ((-SALT .* adv_hConvV) + adv_hConvS) ./ sStarMean;
    adv_vConvSal = ((-SALT .* adv_vConvV) + adv_vConvS) ./ sStarMean;
    
    %diffusion
    dif_hConvSal = dif_hConvS ./ sStarMean;
    dif_vConvSal = dif_vConvS ./ sStarMean;
    
    %forcing
    forcSal = ((-SALT .* forcV) + forcS) ./ sStarMean;
    
    %%
    %DIC budget
    
    S_snap = 0 .* SALT_SNAP;
    
    for nt = 1:2
        
        S_snap(:,:,:,nt) = sStarSnap(:,:,:,nt) .* DIC_SNAP(:,:,:,nt);
        
    end
    
    %tendency, mol m^-3 s^-1
    tendDIC = (S_snap(:,:,:,2) - S_snap(:,:,:,1)) ./ (secPerHour .* dt(timeStep));
   
    DICSnap1 = S_snap(:,:,:,1);
    DICSnap2 = S_snap(:,:,:,2);

    DICdt = (secPerHour .* dt(timeStep));

    adv_hConvDIC = calc_UV_conv(ADVx_DIC,ADVy_DIC) ./ VVV;
    dif_hConvDIC = calc_UV_conv(DFxE_DIC,DFyE_DIC) ./ VVV;
    
    %vertical divergences, mol m^-3 s^-1
    adv_vConvDIC = 0 .* ADVx_DIC;
    dif_vConvDIC = 0 .* ADVx_DIC;
    
    for nz = 1:numLevels
        
        nzp1 = min([nz+1,numLevels]);
        
        adv_vConvDIC(:,:,nz) = squeeze(ADVr_DIC(:,:,nzp1) .* double(nz~=numLevels) - ADVr_DIC(:,:,nz));
        dif_vConvDIC(:,:,nz) = squeeze(DFrI_DIC(:,:,nzp1) .* double(nz~=numLevels) - DFrI_DIC(:,:,nz) ...
            + DFrE_DIC(:,:,nzp1) .* double(nz~=numLevels) - DFrE_DIC(:,:,nz));
        
    end
    
    adv_vConvDIC = adv_vConvDIC ./ VVV;
    dif_vConvDIC = dif_vConvDIC ./ VVV;
    
    %DIC tendency due to air-sea CO2 flux, mol m^-3 s^-1
    forcDIC = 0 .* oceSPtnd;
    
    %DIC tendency due to E/P/runoff, mol m^-3 s^-1
    virtualFluxDIC = 0 .* oceSPtnd;
    
    for nz = 1:numLevels
        
        if nz == 1
            
            forcDIC(:,:,1) = fluxCO2;
            virtualFluxDIC(:,:,1) = DIC_Epr;
            
        else
            
            forcDIC(:,:,nz) = 0;
            virtualFluxDIC(:,:,nz) = 0;
            
        end
        
    end
    
    gDAR_DIC = (gDAR_DIC - (forcDIC ./ dzMat)) ./ mygrid.hFacC; %remove air-sea CO2 flux from gDAR, so it is just biology
    
    forcDIC = mygrid.mskC .* (forcDIC ./ dzMat);
    
    virtualFluxDIC = mygrid.mskC .* virtualFluxDIC .* 0;
    
    %individual biology terms
    bioCons_DIC = mygrid.mskC .* (-cDIC ./ mygrid.hFacC);
    bioCons_DIC_PIC = mygrid.mskC .* (-cDIC_PIC ./ mygrid.hFacC);
    bioResp_DIC = mygrid.mskC .* (respDIC ./ mygrid.hFacC);
    bioRemin_DIC_DOC = mygrid.mskC .* (rDIC_DOC ./ mygrid.hFacC);
    bioRemin_DIC_POC = mygrid.mskC .* (rDIC_POC ./ mygrid.hFacC);
    bioDissc_DIC_PIC = mygrid.mskC .* (dDIC_PIC ./ mygrid.hFacC);
    DIC = mygrid.mskC .* DIC;
 
    %%
    %convert gcmfaces objects to matrices
    
    sStarMean(mygrid.hFacC == 0) = nan;
    sStar = repmat(sStarMean,[1 1 nz]);
    
    sStarMean = convert2gcmfaces(sStarMean);
    sStarSnap = convert2gcmfaces(sStarSnap);
    
    tendV = convert2gcmfaces(tendV);
    adv_hConvV = convert2gcmfaces(adv_hConvV);
    adv_vConvV = convert2gcmfaces(adv_vConvV);
    forcV = convert2gcmfaces(forcV);
    
    tendS = convert2gcmfaces(tendS);
    adv_hConvS = convert2gcmfaces(adv_hConvS);
    adv_vConvS = convert2gcmfaces(adv_vConvS);
    dif_hConvS = convert2gcmfaces(dif_hConvS);
    dif_vConvS = convert2gcmfaces(dif_vConvS);
    forcS = convert2gcmfaces(forcS);
    
    tendSal = convert2gcmfaces(tendSal);
    adv_hConvSal = convert2gcmfaces(adv_hConvSal);
    adv_vConvSal = convert2gcmfaces(adv_vConvSal);
    dif_hConvSal = convert2gcmfaces(dif_hConvSal);
    dif_vConvSal = convert2gcmfaces(dif_vConvSal);
    forcSal = convert2gcmfaces(forcSal);
    
    DIC = convert2gcmfaces(DIC);
    DICSnap1 = convert2gcmfaces(DICSnap1);
    DICSnap2 = convert2gcmfaces(DICSnap2);
    tendDIC = convert2gcmfaces(tendDIC);
    adv_hConvDIC = convert2gcmfaces(adv_hConvDIC);
    adv_vConvDIC = convert2gcmfaces(adv_vConvDIC);
    dif_hConvDIC = convert2gcmfaces(dif_hConvDIC);
    dif_vConvDIC = convert2gcmfaces(dif_vConvDIC);
    forcDIC = convert2gcmfaces(forcDIC);
    virtualFluxDIC = convert2gcmfaces(virtualFluxDIC);
    gDARDIC = convert2gcmfaces(gDAR_DIC);
    
    bioConsDIC = convert2gcmfaces(bioCons_DIC);
    bioConsDIC_PIC = convert2gcmfaces(bioCons_DIC_PIC);
    bioRespDIC = convert2gcmfaces(bioResp_DIC);
    bioReminDIC_DOC = convert2gcmfaces(bioRemin_DIC_DOC);
    bioReminDIC_POC = convert2gcmfaces(bioRemin_DIC_POC);
    bioDisscDIC_PIC = convert2gcmfaces(bioDissc_DIC_PIC);
    
    %%
    %vertical integration
    
    %volume
    intTendV = sum((tendV(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intHConvV = sum((adv_hConvV(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intVConvV = sum((adv_vConvV(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intForcV = sum((forcV(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    
    intTotalV = intHConvV + intVConvV + intForcV;
    intResidualV = intTendV - intTotalV;
    
    %salt
    intTendS = sum((tendS(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intHAdvS = sum((adv_hConvS(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intVAdvS = sum((adv_vConvS(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intHDifS = sum((dif_hConvS(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intVDifS = sum((dif_vConvS(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intForcS = sum((forcS(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    
    intTotalS = intHAdvS + intVAdvS + intHDifS + intVDifS + intForcS;
    intResidualS = intTendS - intTotalS;
    
    %salinity
    intTendSal = sum((tendSal(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intHAdvSal = sum((adv_hConvSal(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intVAdvSal = sum((adv_vConvSal(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intHDifSal = sum((dif_hConvSal(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intVDifSal = sum((dif_vConvSal(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intForcSal = sum((forcSal(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    
    intTotalSal = intHAdvSal + intVAdvSal + intHDifSal + intVDifSal + intForcSal;
    intResidualSal = intTendSal - intTotalSal;
    
    %DIC
    intDIC = sum((DIC(:,:,startIntLevel:endIntLevel) .* sStarMean(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intDICSnap1 = sum(DICSnap1(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel),3,nanString);
    intDICSnap2 = sum(DICSnap2(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel),3,nanString);
    intTendDIC = sum((tendDIC(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intHAdvDIC = sum((adv_hConvDIC(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intVAdvDIC = sum((adv_vConvDIC(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intHDifDIC = sum((dif_hConvDIC(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intVDifDIC = sum((dif_vConvDIC(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intForcDIC = sum((forcDIC(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intVirtualFluxDIC = sum((virtualFluxDIC(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intGDARDIC = sum((gDARDIC(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    
    %DIC biology decomposition
    intBioConsDIC = sum((bioConsDIC(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intBioConsDIC_PIC = sum((bioConsDIC_PIC(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intBioRespDIC = sum((bioRespDIC(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intBioReminDIC_DOC = sum((bioReminDIC_DOC(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intBioReminDIC_POC = sum((bioReminDIC_POC(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intBioDisscDIC_PIC = sum((bioDisscDIC_PIC(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    
    intBioDIC = intBioConsDIC + intBioConsDIC_PIC + intBioRespDIC + ...
        intBioReminDIC_DOC + intBioReminDIC_POC + intBioDisscDIC_PIC;
    
    intTotalDIC = intHAdvDIC + intVAdvDIC ...
        + intHDifDIC + intVDifDIC + intForcDIC + intVirtualFluxDIC.*0 + intBioDIC;
    
    intResidualDIC = intTendDIC - intTotalDIC;
    
    %%
    %save budget
    
    if saveMat

        save([saveDir 'sStar_' num2str(startIntLevel) '_' num2str(endIntLevel) '_' num2str(tt(timeStep)) '_budget.mat'], ...
            'sStarMean','sStarSnap','-v7.3');

        save([saveDir 'volume_' num2str(startIntLevel) '_' num2str(endIntLevel) '_' num2str(tt(timeStep)) '_budget.mat'], ...
            'intTendV','intHConvV','intVConvV','intForcV','intTotalV','intResidualV','-v7.3');
        
        save([saveDir 'salt_' num2str(startIntLevel) '_' num2str(endIntLevel) '_' num2str(tt(timeStep)) '_budget.mat'], ...
            'intTendS','intHAdvS','intVAdvS','intHDifS','intVDifS','intForcS','intTotalS','intResidualS','-v7.3');
        
        save([saveDir 'salinity_' num2str(startIntLevel) '_' num2str(endIntLevel) '_' num2str(tt(timeStep)) '_budget.mat'], ...
            'intTendSal','intHAdvSal','intVAdvSal','intHDifSal','intVDifSal','intForcSal','intTotalSal','intResidualSal','-v7.3');
        
        save([saveDir 'DIC_' num2str(startIntLevel) '_' num2str(endIntLevel) '_' num2str(tt(timeStep)) '_budget.mat'], ...
            'DICdt','intDIC','intDICSnap1','intDICSnap2','intTendDIC','intHAdvDIC','intVAdvDIC','intHDifDIC','intVDifDIC','intForcDIC', ...
            'intBioConsDIC','intBioConsDIC_PIC','intBioRespDIC','intBioReminDIC_DOC','intBioReminDIC_POC','intBioDisscDIC_PIC','intTotalDIC','intResidualDIC','-v7.3');

    end
    
    %%
    %clear gcmfaces objects to avoid memory leaks
    
    clear ETAN oceFWflx SFLUX oceSPflx UVELMASS VVELMASS WVELMASS ...
        SALT ADVr_SLT ADVx_SLT ADVy_SLT DFrI_SLT DFrE_SLT DFxE_SLT DFyE_SLT oceSPtnd ...
        DIC ADVx_DIC ADVy_DIC ADVr_DIC DFxE_DIC DFyE_DIC DFrE_DIC DFrI_DIC DICTFLX ...
        CONSUMP_DIC CONSUMP_DIC_PIC REMIN_DOC REMIN_POC DISSC_PIC ...
        ETAN_SNAP SALT_SNAP DIC_SNAP ALK_SNAP NO3_SNAP NO2_SNAP NH4_SNAP PO4_SNAP Fe_SNAP SiO2_SNAP ...
     
    clear sStarMean sStarSnap
    clear tendV adv_hConvV adv_vConvV forcV
    clear tendS adv_hConvS adv_vConvS dif_hConvS dif_vConvS forcS
    clear tendSal adv_hConvSal adv_vConvSal dif_hConvSal dif_vConvSal forcSal
    clear DIC tendDIC adv_hConvDIC adv_vConvDIC dif_hConvDIC dif_vConvDIC forcDIC virtualFluxDIC gDAR_DIC bioDIC
    clear bioCons_DIC bioCons_DIC_PIC bioResp_DIC bioRemin_DIC_DOC bioRemin_DIC_POC bioDissc_DIC_PIC
    clear ALK tendALK adv_hConvALK adv_vConvALK dif_hConvALK dif_vConvALK forcALK virtualFluxALK gDAR_ALK bioC_ALK bioS_ALK
    clear NO3 tendNO3 adv_hConvNO3 adv_vConvNO3 dif_hConvNO3 dif_vConvNO3 virtualFluxNO3 gDAR_NO3 bioC_NO3 bioS_NO3
    clear NO2 tendNO2 adv_hConvNO2 adv_vConvNO2 dif_hConvNO2 dif_vConvNO2 virtualFluxNO2 gDAR_NO2 bioC_NO2 bioS_NO2
    clear NH4 tendNH4 adv_hConvNH4 adv_vConvNH4 dif_hConvNH4 dif_vConvNH4 virtualFluxNH4 gDAR_NH4 bioC_NH4 bioS_NH4
    clear PO4 tendPO4 adv_hConvPO4 adv_vConvPO4 dif_hConvPO4 dif_vConvPO4 virtualFluxPO4 gDAR_PO4 bioC_PO4 bioS_PO4
    clear Fe tendFe adv_hConvFe adv_vConvFe dif_hConvFe dif_vConvFe forcFe virtualFluxFe gDAR_Fe bioFe sedFe freeFe
    clear SiO2 tendSiO2 dv_hConvSiO2 adv_vConvSiO2 dif_hConvSiO2 dif_vConvSiO2 virtualFluxSiO2 gDAR_SiO2 bioC_SiO2 bioS_SiO2
    
    if(savePlot)
        
        close all
        
    end
    
end
