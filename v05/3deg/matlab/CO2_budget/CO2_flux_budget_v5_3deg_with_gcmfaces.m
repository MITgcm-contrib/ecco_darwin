clear
close all

%%
%settings, modify as needed

intLevel = 1; %integration k level

%set to 1 to plot budget terms
plotCO2Budget = 1;

gridDir = '../../../../../darwin3/run/';
modelDir = '../../../../../darwin3/run/';

%%
%plotting params.

fs = 12;
lw = 2;

%colors1 = parula(500);
colors1 = flipud(cbrewer('div','RdBu',500));

%%
%constants

secPerDay = 86400;
secPerHour = 3600;
hoursPerDay = 24;
deltaT = 3600;

rhoConst = 1029;
mmol_to_mol = 1 ./ 1000;
atm_2_uatm = 10^6;
molPerM3_2_umolPerKg = 10^6 ./ rhoConst;

%%
%load grid

nF = 1;
fileFormat='straight';

global mygrid
grid_load(gridDir,nF,fileFormat);

nx = mygrid.ioSize(1);
ny = mygrid.ioSize(2);
numLevels = numel(mygrid.RC);

dzMatF = mk3D(mygrid.DRF, mygrid.hFacC);
dzMat = dzMatF .* mygrid.hFacC;
RACMat = mk3D(mygrid.RAC, mygrid.hFacC);
VVV = mygrid.mskC .* mygrid.hFacC .* mk3D(mygrid.RAC,mygrid.mskC) .* mk3D(mygrid.DRF,mygrid.mskC);

%%

diagDir = [modelDir 'diags/budget/'];

filename1 = 'average_2d';
filename2 = 'average_velmass_3d';
filename3 = 'average_salt_3d';
filename4 = 'average_dic_3d';
filename5 = 'co2_flux_budget_2d';
filename6 = 'co2_flux_budget_3d';
filename7 = 'snap_2d';
filename8 = 'snap_3d';

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
    
        %load two-dimensional time-averaged fields
    ETAN = convert2gcmfaces(rdmds([diagDir filename1],ttAverage,'rec',1));    
    oceFWflx = convert2gcmfaces(rdmds([diagDir filename1],ttAverage,'rec',2));
    SFLUX = convert2gcmfaces(rdmds([diagDir filename1],ttAverage,'rec',3));
    oceSPflx = convert2gcmfaces(rdmds([diagDir filename1],ttAverage,'rec',4));
    DICTFLX = convert2gcmfaces(rdmds([diagDir filename1],ttAverage,'rec',5)) .* mmol_to_mol; %mol m^-3 s^-1

    %load three-dimensional time-averaged fields
    UVELMASS = convert2gcmfaces(rdmds([diagDir filename2],ttAverage,'rec',1));
    VVELMASS = convert2gcmfaces(rdmds([diagDir filename2],ttAverage,'rec',2));
    WVELMASS = convert2gcmfaces(rdmds([diagDir filename2],ttAverage,'rec',3));
    
    SALT = convert2gcmfaces(rdmds([diagDir filename3],ttAverage,'rec',1));
    ADVx_SLT = convert2gcmfaces(rdmds([diagDir filename3],ttAverage,'rec',2));
    ADVy_SLT = convert2gcmfaces(rdmds([diagDir filename3],ttAverage,'rec',3));
    ADVr_SLT = convert2gcmfaces(rdmds([diagDir filename3],ttAverage,'rec',4));
    DFxE_SLT = convert2gcmfaces(rdmds([diagDir filename3],ttAverage,'rec',5));
    DFyE_SLT = convert2gcmfaces(rdmds([diagDir filename3],ttAverage,'rec',6));
    DFrE_SLT = convert2gcmfaces(rdmds([diagDir filename3],ttAverage,'rec',7));
    DFrI_SLT = convert2gcmfaces(rdmds([diagDir filename3],ttAverage,'rec',8));
    oceSPtnd = convert2gcmfaces(rdmds([diagDir filename3],ttAverage,'rec',9));
    
    DIC = convert2gcmfaces(rdmds([diagDir filename4],ttAverage,'rec',1)) .* mmol_to_mol; %mol m^-3
    ADVx_DIC = convert2gcmfaces(rdmds([diagDir filename4],ttAverage,'rec',2)) .* mmol_to_mol; %mol s^-1
    ADVy_DIC = convert2gcmfaces(rdmds([diagDir filename4],ttAverage,'rec',3)) .* mmol_to_mol; %mol s^-1
    ADVr_DIC = convert2gcmfaces(rdmds([diagDir filename4],ttAverage,'rec',4)) .* mmol_to_mol; %mol s^-1
    DFxE_DIC = convert2gcmfaces(rdmds([diagDir filename4],ttAverage,'rec',5)) .* mmol_to_mol; %mol s^-1
    DFyE_DIC = convert2gcmfaces(rdmds([diagDir filename4],ttAverage,'rec',6)) .* mmol_to_mol; %mol s^-1
    DFrE_DIC = convert2gcmfaces(rdmds([diagDir filename4],ttAverage,'rec',7)) .* mmol_to_mol; %mol s^-1
    DFrI_DIC = convert2gcmfaces(rdmds([diagDir filename4],ttAverage,'rec',8)) .* mmol_to_mol; %mol s^-1
    BIO_DIC = convert2gcmfaces(rdmds([diagDir filename4],ttAverage,'rec',9)) .* mmol_to_mol;  %mol m^-3 s^-1

    APCO2 = mygrid.mskC(:,:,1) .* convert2gcmfaces(rdmds([diagDir filename5],tt(timeStep),'rec',3)) .* atm_2_uatm; %uatm
    THETA = mygrid.mskC .* convert2gcmfaces(rdmds([diagDir filename6],tt(timeStep),'rec',2));
    SALT = mygrid.mskC .* convert2gcmfaces(rdmds([diagDir filename6],tt(timeStep),'rec',3));
    ALK = mygrid.mskC .* convert2gcmfaces(rdmds([diagDir filename6],tt(timeStep),'rec',4)) .* mmol_to_mol;  %mol m^-3 s^-1
    PO4 = mygrid.mskC .* convert2gcmfaces(rdmds([diagDir filename6],tt(timeStep),'rec',5)) .* mmol_to_mol;  %mol m^-3 s^-1
    SIO2 = mygrid.mskC .* convert2gcmfaces(rdmds([diagDir filename6],tt(timeStep),'rec',6)) .* mmol_to_mol;  %mol m^-3 s^-1
    
    %%
    %snapshots
    
    ETAN_SNAP = nan*ones(nx,ny,2);
    SALT_SNAP = nan .* ones(nx,ny,numLevels,2);
    DIC_SNAP = nan .* ones(nx,ny,numLevels,2);
    
    deltaTHETA = convert2gcmfaces(nan .* ones(nx,ny,numLevels,1));
    deltaSALT = convert2gcmfaces(nan .* ones(nx,ny,numLevels,1));
    deltaALK = convert2gcmfaces(nan .* ones(nx,ny,numLevels,1));
    deltaAPCO2 = convert2gcmfaces(nan .* ones(nx,ny,1,1));
    deltaPO4 = convert2gcmfaces(nan .* ones(nx,ny,numLevels,1));
    
    deltaSnapTime = nan;
    deltaAverageTime = nan;
    
    if timeStep == 1
        
        %no pre-initial snapshot
        ETAN_SNAP(:,:,2) = rdmds([diagDir filename7],tt(1),'rec',1);
        SALT_SNAP(:,:,:,2) = rdmds([diagDir filename8],tt(1),'rec',1);
        DIC_SNAP(:,:,:,2) = rdmds([diagDir filename8],tt(1),'rec',2) .* mmol_to_mol;  %mol m^-3 s^-1
        
    elseif timeStep == numFiles %no final snapshot
        
        ETAN_SNAP(:,:,1) = rdmds([diagDir filename7],tt(timeStep-1),'rec',1);
        SALT_SNAP(:,:,:,1) = rdmds([diagDir filename8],tt(timeStep-1),'rec',1);
        DIC_SNAP(:,:,:,1) = rdmds([diagDir filename8],tt(timeStep-1),'rec',2) .* mmol_to_mol;  %mol m^-3 s^-1
        
    else %timeStep~=1 & timeStep~=numFiles
        
        ETAN_SNAP = rdmds([diagDir filename7],ttSnap,'rec',1);
        SALT_SNAP = rdmds([diagDir filename8],ttSnap,'rec',1);
        DIC_SNAP = rdmds([diagDir filename8],ttSnap,'rec',2) .* mmol_to_mol;  %mol m^-3 s^-1
    
        deltaAPCO2 = mygrid.mskC(:,:,1) .* convert2gcmfaces(rdmds([diagDir filename5],tt(timeStep),'rec',3) - ...
            rdmds([diagDir filename5],tt(timeStep-1),'rec',3)) .* atm_2_uatm; %uatm
        
        deltaTHETA = mygrid.mskC .* convert2gcmfaces(rdmds([diagDir filename6],tt(timeStep),'rec',2) - ...
            rdmds([diagDir filename6],tt(timeStep-1),'rec',2));
        
        deltaSALT = mygrid.mskC .* convert2gcmfaces(rdmds([diagDir filename6],tt(timeStep),'rec',3) - ...
            rdmds([diagDir filename6],tt(timeStep-1),'rec',3));
        
        deltaALK = mygrid.mskC .* convert2gcmfaces(rdmds([diagDir filename6],tt(timeStep),'rec',4) - ...
            rdmds([diagDir filename6],tt(timeStep-1),'rec',4)) .* mmol_to_mol;  %mol m^-3 s^-1

        deltaPO4 = mygrid.mskC .* convert2gcmfaces(rdmds([diagDir filename6],tt(timeStep),'rec',5) - ...
            rdmds([diagDir filename6],tt(timeStep),'rec',5)) .* mmol_to_mol;  %mol m^-3 s^-1
        
        deltaSIO2 = mygrid.mskC .* convert2gcmfaces(rdmds([diagDir filename6],tt(timeStep),'rec',6) - ...
            rdmds([diagDir filename6],tt(timeStep),'rec',6)) .* mmol_to_mol;  %mol m^-3 s^-1
        
        deltaSnapTime = dt(timeStep) .* deltaT;
        deltaAverageTime = (tt(timeStep) - tt(timeStep-1)) .* deltaT;
        
    end
    
    ETAN_SNAP = convert2gcmfaces(ETAN_SNAP);
    SALT_SNAP = convert2gcmfaces(SALT_SNAP);
    DIC_SNAP = convert2gcmfaces(DIC_SNAP);
    
    %% 
    
    %surface area-weighted means
    meanSST =  nansum(nansum((THETA(:,:,1) .* mygrid.RAC))) ./ nansum(nansum(nansum(mygrid.RAC)));
    meanSSS =  nansum(nansum((SALT(:,:,1) .* mygrid.RAC))) ./ nansum(nansum(nansum(mygrid.RAC)));
    meanSSALK =  nansum(nansum((ALK(:,:,1) .* mygrid.RAC))) ./ nansum(nansum(nansum(mygrid.RAC)));
    meanAPCO2 =  nansum(nansum((APCO2 .* mygrid.RAC))) ./ nansum(nansum(nansum(mygrid.RAC)));
    meanSSPO4 =  nansum(nansum((PO4(:,:,1) .* mygrid.RAC))) ./ nansum(nansum(nansum(mygrid.RAC)));
    meanSSSIO2 =  nansum(nansum((SIO2(:,:,1) .* mygrid.RAC))) ./ nansum(nansum(nansum(mygrid.RAC)));
    
    %convert to umol kg^-1 for CO2SYS, output is mol m^-3
    [deltaDIC_deltaTHETA, deltaDIC_deltaSALT, deltaDIC_deltaALK, deltaDIC_deltaAPCO2] = compute_dic_coeff( ...
        meanSST, ...
        meanSSS, ...
        meanSSALK .* molPerM3_2_umolPerKg, ....
        meanSSPO4 .* molPerM3_2_umolPerKg, ...
        meanSSSIO2 .* molPerM3_2_umolPerKg, ....
        meanAPCO2);
    
    %disp(['deltaDIC/detlaSST = ' num2str(deltaDIC_deltaTHETA)]);
    %disp(['deltaDIC/deltaSSS = ' num2str(deltaDIC_deltaSALT)]);
    %disp(['deltaDIC/deltaSSALK = ' num2str(deltaDIC_deltaALK)]);
    %disp(['deltaDIC/deltaAPCO2 = ' num2str(deltaDIC_deltaAPCO2)]);
    
    %%
    %volume budget, s^-1
    
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
        
        S_snap(:,:,:,nt) = (SALT_SNAP(:,:,:,nt) .* (1+mk3D(ETAN_SNAP(:,:,nt) ./ mygrid.Depth,dzMat)));
        
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
    
    surfS = mygrid.mskC(:,:,1) .* SALT(:,:,1);
    
    %%
    %salinity budget
    
    rstarfac = (mygrid.Depth + ETAN) ./ mygrid.Depth;
    
    %tendency
    tendSal = mygrid.mskC .* (SALT_SNAP(:,:,:,2) - SALT_SNAP(:,:,:,1)) ...
        ./ (secPerHour .* dt(timeStep));
    
    %advection
    adv_hConvSal = (-SALT .* adv_hConvV + adv_hConvS) ./ rstarfac;
    adv_vConvSal = (-SALT .* adv_vConvV + adv_vConvS) ./ rstarfac;
    
    %diffusion
    dif_hConvSal = dif_hConvS ./ rstarfac;
    dif_vConvSal = dif_vConvS ./ rstarfac;
    
    %forcing
    forcSal = (-SALT .* forcV + forcS) ./ rstarfac;
    
    %%
    %DIC budget
    
    D_snap = 0 .* DIC_SNAP;
    
    for nt = 1:2
        
        D_snap(:,:,:,nt) = (DIC_SNAP(:,:,:,nt) .* (1+mk3D(ETAN_SNAP(:,:,nt) ./ mygrid.Depth,dzMat)));
        
    end
    
    %tendency, mol m^-3 s^-1
    tendDIC = (D_snap(:,:,:,2) - D_snap(:,:,:,1)) ./ (secPerHour .* dt(timeStep));
    
    %horizontal divergences, mol m^-3 s^-1
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
    
    %air-sea CO2 flux tendency, mol m^-3 s^-1
    forcDIC = 0 .* oceSPtnd;
    
    for nz = 1:numLevels
        
        if nz == 1
            
            forcDIC(:,:,1) = DICTFLX;
            
        else
            
            forcDIC(:,:,nz) = 0;
            
        end
        
    end
    
    BIO_DIC = BIO_DIC - forcDIC; %remove air-sea CO2 flux from gDAR
    
    forcDIC = mygrid.mskC .* forcDIC;
    
    %biology, mol m^-3 s^-1
    bioDIC = mygrid.mskC .* BIO_DIC;
        
    surfDIC = mygrid.mskC(:,:,1) .* DIC(:,:,1);
    meanSurf_DIC = nansum(surfDIC(:) .* mygrid.RAC(:)) ./ nansum(mygrid.RAC(:));
    
    virtualFluxDIC = mygrid.mskC .* (forcSal - forcS) .* (surfDIC ./ surfS) .* mygrid.hFacC(:,:,1);
    
    %z* correction, as done in salinity budget
    
    tendDICzs = mygrid.mskC .* (DIC_SNAP(:,:,:,2) - DIC_SNAP(:,:,:,1)) ...
        ./ (secPerHour .* dt(timeStep));
    
    %advection
    adv_hConvDICzs = (-DIC .* adv_hConvV + adv_hConvDIC) ./ rstarfac;
    adv_vConvDICzs = (-DIC .* adv_vConvV + adv_vConvDIC) ./ rstarfac;
    
    %diffusion
    dif_vConvDICzs = dif_vConvDIC ./ rstarfac;
    dif_hConvDICzs = dif_hConvDIC ./ rstarfac;
    
    %forcing
    forcDICzs = forcDIC ./ rstarfac;
    
    %biology
    bioDICzs =  bioDIC ./ rstarfac;
    
    %virtual flux
    virtualFluxDICzs = virtualFluxDIC ./ rstarfac;
    
    %% 
    
    %solubility terms
    dDIC_dTHETA = deltaTHETA .* deltaDIC_deltaTHETA;
    dDIC_dSALT = deltaSALT .* deltaDIC_deltaSALT;
    dDIC_dALK = deltaALK .* deltaDIC_deltaALK;
    dDIC_dAPCO2 = deltaAPCO2 .* deltaDIC_deltaAPCO2;
    
    %%
    %convert gcmfaces objects to matrices
    
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
    
    tendDIC = convert2gcmfaces(tendDICzs);
    adv_hConvDIC = convert2gcmfaces(adv_hConvDICzs);
    adv_vConvDIC = convert2gcmfaces(adv_vConvDICzs);
    dif_hConvDIC = convert2gcmfaces(dif_hConvDICzs);
    dif_vConvDIC = convert2gcmfaces(dif_vConvDICzs);
    forcDIC = convert2gcmfaces(forcDICzs);
    bioDIC = convert2gcmfaces(bioDICzs);
    virtualFluxDIC = convert2gcmfaces(virtualFluxDICzs);
    
    dDIC_dTHETA = convert2gcmfaces(dDIC_dTHETA);
    dDIC_dSALT = convert2gcmfaces(dDIC_dSALT);
    dDIC_dALK = convert2gcmfaces(dDIC_dALK);
    dDIC_dAPCO2 = convert2gcmfaces(dDIC_dAPCO2);
    
    %%
    %vertically integrate
    
    intTendV = nansum(tendV(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel) .* convert2gcmfaces(mygrid.hFacC(:,:,1:intLevel))),3);
    intHConvV = nansum(adv_hConvV(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel) .* convert2gcmfaces(mygrid.hFacC(:,:,1:intLevel))),3);
    intVConvV = nansum(adv_vConvV(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel) .* convert2gcmfaces(mygrid.hFacC(:,:,1:intLevel))),3);
    intForcV = nansum(forcV(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel) .* convert2gcmfaces(mygrid.hFacC(:,:,1:intLevel))),3);
    
    intTendS = nansum(tendS(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel) .* convert2gcmfaces(mygrid.hFacC(:,:,1:intLevel))),3);
    intHAdvS = nansum(adv_hConvS(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel) .* convert2gcmfaces(mygrid.hFacC(:,:,1:intLevel))),3);
    intVAdvS = nansum(adv_vConvS(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel) .* convert2gcmfaces(mygrid.hFacC(:,:,1:intLevel))),3);
    intHDifS = nansum(dif_hConvS(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel) .* convert2gcmfaces(mygrid.hFacC(:,:,1:intLevel))),3);
    intVDifS = nansum(dif_vConvS(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel) .* convert2gcmfaces(mygrid.hFacC(:,:,1:intLevel))),3);
    intForcS = nansum(forcS(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel) .* convert2gcmfaces(mygrid.hFacC(:,:,1:intLevel))),3);
    
    intTendSal = nansum(tendSal(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel) .* convert2gcmfaces(mygrid.hFacC(:,:,1:intLevel))),3);
    intHAdvSal = nansum(adv_hConvSal(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel) .* convert2gcmfaces(mygrid.hFacC(:,:,1:intLevel))),3);
    intVAdvSal = nansum(adv_vConvSal(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel) .* convert2gcmfaces(mygrid.hFacC(:,:,1:intLevel))),3);
    intHDifSal = nansum(dif_hConvSal(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel) .* convert2gcmfaces(mygrid.hFacC(:,:,1:intLevel))),3);
    intVDifSal = nansum(dif_vConvSal(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel) .* convert2gcmfaces(mygrid.hFacC(:,:,1:intLevel))),3);
    intForcSal = nansum(forcSal(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel) .* convert2gcmfaces(mygrid.hFacC(:,:,1:intLevel))),3);
    
    intHAdvDIC = nansum(adv_hConvDIC(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel) .* convert2gcmfaces(mygrid.hFacC(:,:,1:intLevel))),3);
    intVAdvDIC = nansum(adv_vConvDIC(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel) .* convert2gcmfaces(mygrid.hFacC(:,:,1:intLevel))),3);
    intHDifDIC = nansum(dif_hConvDIC(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel) .* convert2gcmfaces(mygrid.hFacC(:,:,1:intLevel))),3);
    intVDifDIC = nansum(dif_vConvDIC(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel) .* convert2gcmfaces(mygrid.hFacC(:,:,1:intLevel))),3);
    intVirtualFluxDIC = nansum(virtualFluxDIC(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel) .* convert2gcmfaces(mygrid.hFacC(:,:,1:intLevel))),3);
        
    %CO2 flux budget
    intTendDIC(:,:,timeStep) = nansum(tendDIC(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel) .* convert2gcmfaces(mygrid.hFacC(:,:,1:intLevel))),3) .* deltaSnapTime; 

    intForcDIC(:,:,timeStep) = nansum(forcDIC(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel) .* convert2gcmfaces(mygrid.hFacC(:,:,1:intLevel))),3) .* deltaAverageTime; 
    
    int_dDIC_dTHETA(:,:,timeStep) = nansum(dDIC_dTHETA(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel) ...
        .* convert2gcmfaces(mygrid.hFacC(:,:,1:intLevel))),3);
    
    int_dDIC_dSALT(:,:,timeStep) = nansum(dDIC_dSALT(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel) ...
        .* convert2gcmfaces(mygrid.hFacC(:,:,1:intLevel))),3);
    
    int_dDIC_dALK(:,:,timeStep) = nansum(dDIC_dALK(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel) ...
        .* convert2gcmfaces(mygrid.hFacC(:,:,1:intLevel))),3);
    
    int_dDIC_dAPCO2(:,:,timeStep) = nansum(dDIC_dAPCO2(:,:,1) .* convert2gcmfaces(dzMat(:,:,1:intLevel) ...
        .* convert2gcmfaces(mygrid.hFacC(:,:,1:intLevel))),3);
    
    intSolDIC(:,:,timeStep) =  (int_dDIC_dTHETA(:,:,timeStep) +  int_dDIC_dSALT(:,:,timeStep) ...
        +  int_dDIC_dALK(:,:,timeStep) +  int_dDIC_dAPCO2(:,:,timeStep));
    
    intBioDIC(:,:,timeStep) = nansum(bioDIC(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel) .* convert2gcmfaces(mygrid.hFacC(:,:,1:intLevel))),3) * deltaAverageTime; 

    intPhysicsDIC(:,:,timeStep) = (intHAdvDIC + intVAdvDIC + intHDifDIC + intVDifDIC + intVirtualFluxDIC) .*  deltaAverageTime;
    
    %% 
    
    intSumDIC1(:,:,timeStep) = intTendDIC(:,:,timeStep) ...
        - intBioDIC(:,:,timeStep) - intPhysicsDIC(:,:,timeStep);
    
    intSumDIC2(:,:,timeStep) = intSolDIC(:,:,timeStep) ...
        - intBioDIC(:,:,timeStep) - intPhysicsDIC(:,:,timeStep);

    %% 
    
    if plotCO2Budget
        
        hFig1 = figure(1);
        set(hFig1,'units','normalized','outerposition',[0 0 1 1]);
        set(gcf,'color',[1 1 1]);
        set(gca,'color',[0.5 0.5 0.5]);
        
        cMin = -2;
        cMax = 2;
        
        subplot(2,6,1);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,intTendDIC(:,:,timeStep));
        shading flat
        
        caxis([cMin cMax]);
        
        colormap(colors1);
        
        %colorbar
        
        axis tight
        
        xlim([0 360]);
        ylim([-90 90]);
        
        set(gca,'xtick',[0:180:360]);
        set(gca,'ytick',[-90:30:90]);
        
        ylabel('Latitude (deg)');
        
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'DIC Tendency (mol m^-^2), ';'timeStep: ';num2str(tt(timeStep))},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,6,2);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
       
        pcolor(mygrid.XC.f1,mygrid.YC.f1,intForcDIC(:,:,timeStep));
        
        shading flat
        
        caxis([cMin cMax]);
        
        colormap(colors1);
        
        %colorbar
        
        axis tight
        
        xlim([0 360]);
        ylim([-90 90]);
        
        set(gca,'xtick',[0:180:360]);
        set(gca,'ytick',[-90:30:90]);
        
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Air-sea CO_2 Flux';'(mol m^-^2)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,6,3);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,int_dDIC_dTHETA(:,:,timeStep));
        shading flat
        
        caxis([cMin cMax]);
        
        colormap(colors1);
        
        %colorbar
        
        axis tight
        
        xlim([0 360]);
        ylim([-90 90]);
        
        set(gca,'xtick',[0:180:360]);
        set(gca,'ytick',[-90:30:90]);
        
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'dDIC/dTHETA';'(mol m^-^2)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,6,4);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,int_dDIC_dSALT(:,:,timeStep));
        shading flat
        
        caxis([cMin cMax]);
        
        colormap(colors1);
        
        %colorbar
        
        axis tight
        
        xlim([0 360]);
        ylim([-90 90]);
        
        set(gca,'xtick',[0:180:360]);
        set(gca,'ytick',[-90:30:90]);
        
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'dDIC/dSALT';'(mol m^-^2)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,6,5);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,int_dDIC_dALK(:,:,timeStep));
        shading flat
        
        caxis([cMin cMax]);
        
        colormap(colors1);
        
        %colorbar
        
        axis tight
        
        xlim([0 360]);
        ylim([-90 90]);
        
        set(gca,'xtick',[0:180:360]);
        set(gca,'ytick',[-90:30:90]);
        
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'dDIC/dALK';'(mol m^-^2)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,6,6);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,int_dDIC_dAPCO2(:,:,timeStep));
        shading flat
        
        caxis([cMin cMax]);
        
        colormap(colors1);
        
        %colorbar
        
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
        
        title({'dDIC/dAPCO2';'(mol m^-^2)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,6,7);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,intSolDIC(:,:,timeStep));
        shading flat
        
        caxis([cMin cMax]);
        
        colormap(colors1);
        
        %colorbar
        
        axis tight
        
        xlim([0 360]);
        ylim([-90 90]);
        
        set(gca,'xtick',[0:180:360]);
        set(gca,'ytick',[-90:30:90]);
        
        xlabel('Longitude (deg)');
        
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Solubility';'(mol m^-^2)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,6,8);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,intBioDIC(:,:,timeStep));
        shading flat
        
        caxis([cMin cMax]);
        
        colormap(colors1);
        
        %colorbar
        
        axis tight
        
        xlim([0 360]);
        ylim([-90 90]);
        
        set(gca,'xtick',[0:180:360]);
        set(gca,'ytick',[-90:30:90]);
        
        xlabel('Longitude (deg)');
        
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Biology';'(mol m^-^2)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,6,9);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,intPhysicsDIC(:,:,timeStep));
        shading flat
        
        caxis([cMin cMax]);
        
        colormap(colors1);
        
        %colorbar
        
        axis tight
        
        xlim([0 360]);
        ylim([-90 90]);
        
        set(gca,'xtick',[0:180:360]);
        set(gca,'ytick',[-90:30:90]);
        
        xlabel('Longitude (deg)');
        
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Physics';'(mol m^-^2)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,6,10);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,intSumDIC1(:,:,timeStep));
        
        shading flat
        
        caxis([cMin cMax]);
        
        colormap(colors1);
        
        %colorbar
        
        axis tight
        
        xlim([0 360]);
        ylim([-90 90]);
        
        set(gca,'xtick',[0:180:360]);
        set(gca,'ytick',[-90:30:90]);
        
        xlabel('Longitude (deg)');
        
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Total Budget';'(mol m^-^2)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,6,11);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,intForcDIC(:,:,timeStep) - intSumDIC1(:,:,timeStep));
        
        shading flat
        caxis([cMin cMax]);
        
        colormap(colors1);
        
        %colorbar
        
        axis tight
        
        xlim([0 360]);
        ylim([-90 90]);
        
        set(gca,'xtick',[0:180:360]);
        set(gca,'ytick',[-90:30:90]);
        
        xlabel('Longitude (deg)');
        
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Total Budget - Tendency';'(mol m^-^2)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,6,12);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,intForcDIC(:,:,timeStep) - intSumDIC2(:,:,timeStep));
        
        shading flat
        caxis([cMin cMax]);
        
        colormap(colors1);
        
        %colorbar
        
        axis tight
        
        xlim([0 360]);
        ylim([-90 90]);
        
        set(gca,'xtick',[0:180:360]);
        set(gca,'ytick',[-90:30:90]);
        
        xlabel('Longitude (deg)');
        
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Total Budget w/ Solubility - Tendency';'(mol m^-^2)'},'FontWeight','Bold','FontSize',fs);
        
        drawnow
        pause
        
    end
    
    %clear gcmfaces objects to avoid memory leaks
    
    clear ETAN oceFWflx SFLUX  oceSPflx UVELMASS VVELMASS WVELMASS ...
    SALT ADVr_SLT ADVx_SLT ADVy_SLT DFrI_SLT DFrE_SLT DFxE_SLT DFyE_SLT oceSPtnd ...
    DIC ADVx_DIC ADVy_DIC ADVr_DIC DFxE_DIC DFyE_DIC DFrE_DIC DFrI_DIC DICTFLX ...
    CONSUMP_DIC CONSUMP_DIC_PIC REMIN_DOC REMIN_POC DISSC_PIC ...
    ETAN_SNAP SALT_SNAP DIC_SNAP ...
    tendV adv_hConvV adv_vConvV forcV ...
    S_snap tendS adv_hConvS dif_hConvS adv_vConvS dif_vConvS forcS surfS ...
    rstarfac tendSal adv_hConvSal adv_vConvSal dif_hConvSal dif_vConvSal forcSal ...
    tendDIC adv_hConvDIC dif_hConvDIC adv_vConvDIC dif_vConvDIC forcDIC bioDIC surfDIC virtualFluxDIC ...
    tendDICzs adv_hConvDICzs adv_vConvDICzs dif_vConvDICzs dif_hConvDICzs forcDICzs bioDICzs virtualFluxDICzs ...
    PCO2 APCO2 SSDIC THETA ALK PO4 SIO2 
    
    close all
    
end

%% 

timeStepStart = (1 .* 12) + 1;
timeStepEnd = (4 .* 12) + 1;

cum_intTendDIC = nansum(intTendDIC(:,:,timeStepStart:timeStepEnd),3);
cum_intForcDIC = nansum(intForcDIC(:,:,timeStepStart:timeStepEnd),3);
cum_int_dDIC_dTHETA = nansum(int_dDIC_dTHETA(:,:,timeStepStart:timeStepEnd),3);
cum_int_dDIC_dSALT = nansum(int_dDIC_dSALT(:,:,timeStepStart:timeStepEnd),3);
cum_int_dDIC_dALK = nansum(int_dDIC_dALK(:,:,timeStepStart:timeStepEnd),3);
cum_int_dDIC_dAPCO2 = nansum(int_dDIC_dAPCO2(:,:,timeStepStart:timeStepEnd),3);

cum_intSolDIC = nansum(intSolDIC(:,:,timeStepStart:timeStepEnd),3);
cum_intBioDIC = nansum(intBioDIC(:,:,timeStepStart:timeStepEnd),3);
cum_intPhysicsDIC = nansum(intPhysicsDIC(:,:,timeStepStart:timeStepEnd),3);
    
cum_intSumDIC1 = nansum(intSumDIC1(:,:,timeStepStart:timeStepEnd),3);
cum_intSumDIC2 = nansum(intSumDIC2(:,:,timeStepStart:timeStepEnd),3);

%global
intTendDIC_all = squeeze(squeeze(nansum(nansum(intTendDIC .* RACMat.f1(:,:,1)))));
intForcDIC_all = squeeze(squeeze(nansum(nansum(intForcDIC .* RACMat.f1(:,:,1)))));
int_dDIC_dTHETA_all = squeeze(squeeze(nansum(nansum(int_dDIC_dTHETA .* RACMat.f1(:,:,1)))));
int_dDIC_dSALT_all = squeeze(squeeze(nansum(nansum(int_dDIC_dSALT .* RACMat.f1(:,:,1)))));
int_dDIC_dALK_all = squeeze(squeeze(nansum(nansum(int_dDIC_dALK .* RACMat.f1(:,:,1)))));
int_dDIC_dAPCO2_all = squeeze(squeeze(nansum(nansum(int_dDIC_dAPCO2 .* RACMat.f1(:,:,1)))));
intSolDIC_all = squeeze(squeeze(nansum(nansum(intSolDIC .* RACMat.f1(:,:,1)))));
intBioDIC_all = squeeze(squeeze(nansum(nansum(intBioDIC .* RACMat.f1(:,:,1)))));
intPhysicsDIC_all = squeeze(squeeze(nansum(nansum(intPhysicsDIC .* RACMat.f1(:,:,1)))));
intSumDIC1_all = squeeze(squeeze(nansum(nansum(intSumDIC1 .* RACMat.f1(:,:,1)))));
intSumDIC2_all = squeeze(squeeze(nansum(nansum(intSumDIC2 .* RACMat.f1(:,:,1)))));

hFig2 = figure(2);
set(hFig2,'units','normalized','outerposition',[0 0 1 1]);
set(gcf,'color',[1 1 1]);
set(gca,'color',[0.5 0.5 0.5]);

cMin = -20;
cMax = 20;

colors1 = flipud(cbrewer('div','RdBu',500));

subplot(2,6,1);

hold on

set(gca,'color',[0.5 0.5 0.5]);

pcolor(mygrid.XC.f1,mygrid.YC.f1,cum_intTendDIC);
shading flat

caxis([cMin cMax]);

colormap(colors1);

%colorbar

axis tight

xlim([0 360]);
ylim([-90 90]);

set(gca,'xtick',[0:180:360]);
set(gca,'ytick',[-90:30:90]);

ylabel('Latitude (deg)');

box on
set(gca,'LineWidth',lw);
set(gca,'FontSize',fs);

title({'DIC Tendency (mol m^-^2), ';'timeStep: ';num2str(tt(timeStep))},'FontWeight','Bold','FontSize',fs);

subplot(2,6,2);

hold on

set(gca,'color',[0.5 0.5 0.5]);

pcolor(mygrid.XC.f1,mygrid.YC.f1,cum_intForcDIC);

shading flat

caxis([cMin cMax]);

colormap(colors1);

%colorbar

axis tight

xlim([0 360]);
ylim([-90 90]);

set(gca,'xtick',[0:180:360]);
set(gca,'ytick',[-90:30:90]);

box on
set(gca,'LineWidth',lw);
set(gca,'FontSize',fs);

title({'Air-sea CO_2 Flux';'(mol m^-^2)'},'FontWeight','Bold','FontSize',fs);

subplot(2,6,3);

hold on

set(gca,'color',[0.5 0.5 0.5]);

pcolor(mygrid.XC.f1,mygrid.YC.f1,cum_int_dDIC_dTHETA);
shading flat

caxis([cMin cMax]);

colormap(colors1);

%colorbar

axis tight

xlim([0 360]);
ylim([-90 90]);

set(gca,'xtick',[0:180:360]);
set(gca,'ytick',[-90:30:90]);

box on
set(gca,'LineWidth',lw);
set(gca,'FontSize',fs);

title({'dDIC/dTHETA';'(mol m^-^2)'},'FontWeight','Bold','FontSize',fs);

subplot(2,6,4);

hold on

set(gca,'color',[0.5 0.5 0.5]);

pcolor(mygrid.XC.f1,mygrid.YC.f1,cum_int_dDIC_dSALT);
shading flat

caxis([cMin cMax]);

colormap(colors1);

%colorbar

axis tight

xlim([0 360]);
ylim([-90 90]);

set(gca,'xtick',[0:180:360]);
set(gca,'ytick',[-90:30:90]);

box on
set(gca,'LineWidth',lw);
set(gca,'FontSize',fs);

title({'dDIC/dSALT';'(mol m^-^2)'},'FontWeight','Bold','FontSize',fs);

subplot(2,6,5);

hold on

set(gca,'color',[0.5 0.5 0.5]);

pcolor(mygrid.XC.f1,mygrid.YC.f1,cum_int_dDIC_dALK);
shading flat

caxis([cMin cMax]);

colormap(colors1);

%colorbar

axis tight

xlim([0 360]);
ylim([-90 90]);

set(gca,'xtick',[0:180:360]);
set(gca,'ytick',[-90:30:90]);

box on
set(gca,'LineWidth',lw);
set(gca,'FontSize',fs);

title({'dDIC/dALK';'(mol m^-^2)'},'FontWeight','Bold','FontSize',fs);

subplot(2,6,6);

hold on

set(gca,'color',[0.5 0.5 0.5]);

pcolor(mygrid.XC.f1,mygrid.YC.f1,cum_int_dDIC_dAPCO2);
shading flat

caxis([cMin cMax]);

colormap(colors1);

%colorbar

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

title({'dDIC/dAPCO2';'(mol m^-^2)'},'FontWeight','Bold','FontSize',fs);

subplot(2,6,7);

hold on

set(gca,'color',[0.5 0.5 0.5]);

pcolor(mygrid.XC.f1,mygrid.YC.f1,cum_intSolDIC);
shading flat

caxis([cMin cMax]);

colormap(colors1);

%colorbar

axis tight

xlim([0 360]);
ylim([-90 90]);

set(gca,'xtick',[0:180:360]);
set(gca,'ytick',[-90:30:90]);

xlabel('Longitude (deg)');

box on
set(gca,'LineWidth',lw);
set(gca,'FontSize',fs);

title({'Solubility';'(mol m^-^2)'},'FontWeight','Bold','FontSize',fs);

subplot(2,6,8);

hold on

set(gca,'color',[0.5 0.5 0.5]);

pcolor(mygrid.XC.f1,mygrid.YC.f1,cum_intBioDIC);
shading flat

caxis([cMin cMax]);

colormap(colors1);

%colorbar

axis tight

xlim([0 360]);
ylim([-90 90]);

set(gca,'xtick',[0:180:360]);
set(gca,'ytick',[-90:30:90]);

xlabel('Longitude (deg)');

box on
set(gca,'LineWidth',lw);
set(gca,'FontSize',fs);

title({'Biology';'(mol m^-^2)'},'FontWeight','Bold','FontSize',fs);

subplot(2,6,9);

hold on

set(gca,'color',[0.5 0.5 0.5]);

pcolor(mygrid.XC.f1,mygrid.YC.f1,cum_intPhysicsDIC);
shading flat

caxis([cMin cMax]);

colormap(colors1);

%colorbar

axis tight

xlim([0 360]);
ylim([-90 90]);

set(gca,'xtick',[0:180:360]);
set(gca,'ytick',[-90:30:90]);

xlabel('Longitude (deg)');

box on
set(gca,'LineWidth',lw);
set(gca,'FontSize',fs);

title({'Physics';'(mol m^-^2)'},'FontWeight','Bold','FontSize',fs);

subplot(2,6,10);

hold on

set(gca,'color',[0.5 0.5 0.5]);

pcolor(mygrid.XC.f1,mygrid.YC.f1,cum_intSumDIC2);

shading flat

caxis([cMin cMax]);

colormap(colors1);

%colorbar

axis tight

xlim([0 360]);
ylim([-90 90]);

set(gca,'xtick',[0:180:360]);
set(gca,'ytick',[-90:30:90]);

xlabel('Longitude (deg)');

box on
set(gca,'LineWidth',lw);
set(gca,'FontSize',fs);

title({'Total Budget';'(mol m^-^2)'},'FontWeight','Bold','FontSize',fs);

subplot(2,6,11);

hold on

set(gca,'color',[0.5 0.5 0.5]);

pcolor(mygrid.XC.f1,mygrid.YC.f1,cum_intForcDIC - cum_intSumDIC2);

shading flat
caxis([cMin cMax]);

colormap(colors1);

%colorbar

axis tight

xlim([0 360]);
ylim([-90 90]);

set(gca,'xtick',[0:180:360]);
set(gca,'ytick',[-90:30:90]);

xlabel('Longitude (deg)');

box on
set(gca,'LineWidth',lw);
set(gca,'FontSize',fs);

title({'Budget Total - Tend';'(mol m^-^2)'},'FontWeight','Bold','FontSize',fs);

drawnow

hFig3 = figure(3);
set(hFig3,'units','normalized','outerposition',[0 0 1 1]);
set(gcf,'color',[1 1 1]);

lw = 3;

colors2 = parula(10);

hold on

%p1 = plot(intTendDIC_all,'Color','k','LineWidth',lw);
p2 = plot(intForcDIC_all,'Color',colors2(1,:),'LineWidth',lw)
%p3 = plot(int_dDIC_dTHETA_all,'Color',colors2(2,:),'LineWidth',lw)
%p4 = plot(int_dDIC_dSALT_all,'Color',colors2(3,:),'LineWidth',lw) 
%p5 = plot(int_dDIC_dALK_all,'Color',colors2(4,:),'LineWidth',lw) 
%p6 = plot(int_dDIC_dAPCO2_all,'Color',colors2(5,:),'LineWidth',lw)

%p7 = plot(intSolDIC_all,'Color',colors2(6,:),'LineWidth',lw) 
%p8 = plot(intBioDIC_all,'Color',colors2(7,:),'LineWidth',lw) 
%p9 = plot(intPhysicsDIC_all,'Color',colors2(8,:),'LineWidth',lw) 

p10 = plot(intSumDIC1_all,'Color',colors2(9,:),'LineWidth',lw)
p11 = plot(intSumDIC2_all,'Color',colors2(10,:),'LineWidth',lw) 

%legend([p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11],{'TendDIC','ForcDIC', ...
%    'dDIC_dTHETA','dDIC_dSALT','dDIC_dALK','dDIC_dAPCO2', ...
%    'Solubility','Biology','Physics','Sum1','Sum2'});
    
axis tight

box on
set(gca,'LineWidth',lw);
set(gca,'FontSize',fs);

drawnow