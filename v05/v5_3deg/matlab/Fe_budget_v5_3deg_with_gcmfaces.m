clear
close all

%%
%settings, modify as needed

intLevel = 15; %integration k level

%set to 1 to plot budget terms
plotVolumeBudget = 0;
plotSalinityBudget = 0;
plotFeBudget = 1;

gridDir = '../../../../darwin3/run/';
modelDir = '../../../../darwin3/run/';

%%
%plotting params.

fs = 16;
lw = 2;

colors = parula(500);
%colors = flipud(cbrewer('div','RdBu',500));

%%
%constants

secPerDay = 86400;
secPerHour = 3600;
hoursPerDay = 24;
deltaT = 3600;

rhoConst = 1029;
mmol_to_mol = 10^-3;
umol_to_mol = 10^-6;

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

diagDir = [modelDir 'diags/'];

filename1 = 'budg2d_zflux_set1';
filename2 = 'trsp_3d_set3';
filename3 = 'trsp_3d_set2';
filename4 = 'state_3d_set1';
filename5 = 'state_3d_set2';
filename6 = 'state_2d_set1';
filename7 = 'state_3d_snap_set1';
filename8 = 'budg2d_snap_set1';
filename9 = 'trsp_3d_set1';
filename10 = 'fe_physics';
filename11 = 'fe_darwin_3d';
filename12 = 'fe_darwin_2d';

%%

fn = dir([diagDir filename6 '.*.data']);
tt = zeros(length(fn),1);

for i = 1:length(tt)
    
    nme = fn(i).name;
    tt(i) = str2num(nme(end-10:end-5));
    
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
    ETAN = convert2gcmfaces(rdmds([diagDir filename6],ttAverage,'rec',1));
    oceFWflx = convert2gcmfaces(rdmds([diagDir filename1],ttAverage,'rec',1));
    SFLUX = convert2gcmfaces(rdmds([diagDir filename1],ttAverage,'rec',5));
    oceSPflx = convert2gcmfaces(rdmds([diagDir filename1],ttAverage,'rec',7));
    
    %load three-dimensional time-averaged fields
    UVELMASS = convert2gcmfaces(rdmds([diagDir filename9],ttAverage,'rec',1));
    VVELMASS = convert2gcmfaces(rdmds([diagDir filename9],ttAverage,'rec',2));
    WVELMASS = convert2gcmfaces(rdmds([diagDir filename9],ttAverage,'rec',3));
    
    SALT = convert2gcmfaces(rdmds([diagDir filename4],ttAverage,'rec',2));
    ADVr_SLT = convert2gcmfaces(rdmds([diagDir filename2],ttAverage,'rec',4));
    ADVx_SLT = convert2gcmfaces(rdmds([diagDir filename3],ttAverage,'rec',7));
    ADVy_SLT = convert2gcmfaces(rdmds([diagDir filename3],ttAverage,'rec',8));
    DFrI_SLT = convert2gcmfaces(rdmds([diagDir filename2],ttAverage,'rec',6));
    DFrE_SLT = convert2gcmfaces(rdmds([diagDir filename2],ttAverage,'rec',5));
    DFxE_SLT = convert2gcmfaces(rdmds([diagDir filename3],ttAverage,'rec',5));
    DFyE_SLT = convert2gcmfaces(rdmds([diagDir filename3],ttAverage,'rec',6));
    oceSPtnd = convert2gcmfaces(rdmds([diagDir filename5],ttAverage,'rec',3));
    
    Fe = convert2gcmfaces(rdmds([diagDir filename10],ttAverage,'rec',1)) .* mmol_to_mol; %mol m^-3
    ADVx_Fe = convert2gcmfaces(rdmds([diagDir filename10],ttAverage,'rec',2)) .* mmol_to_mol; %mol s^-1
    ADVy_Fe = convert2gcmfaces(rdmds([diagDir filename10],ttAverage,'rec',3)) .* mmol_to_mol; %mol s^-1
    ADVr_Fe = convert2gcmfaces(rdmds([diagDir filename10],ttAverage,'rec',4)) .* mmol_to_mol; %mol s^-1
    DFxE_Fe = convert2gcmfaces(rdmds([diagDir filename10],ttAverage,'rec',5)) .* mmol_to_mol; %mol s^-1
    DFyE_Fe = convert2gcmfaces(rdmds([diagDir filename10],ttAverage,'rec',6)) .* mmol_to_mol; %mol s^-1
    DFrE_Fe = convert2gcmfaces(rdmds([diagDir filename10],ttAverage,'rec',7)) .* mmol_to_mol; %mol s^-1
    DFrI_Fe = convert2gcmfaces(rdmds([diagDir filename10],ttAverage,'rec',8)) .* mmol_to_mol; %mol s^-1
   
    %3D fe
    C_Fe = convert2gcmfaces(rdmds([diagDir filename11],ttAverage,'rec',1)) .* mmol_to_mol;  %mol m^-3 s^-1
    S_Fe = convert2gcmfaces(rdmds([diagDir filename11],ttAverage,'rec',2)) .* mmol_to_mol;  %mol m^-3 s^-1
    DARWIN_Fe = convert2gcmfaces(rdmds([diagDir filename11],ttAverage,'rec',3)) .* mmol_to_mol;  %mol m^-3 s^-1
    SEDFe = convert2gcmfaces(rdmds([diagDir filename11],ttAverage,'rec',4)) .* mmol_to_mol;  %mol m^-3 s^-1
    FREEFe = convert2gcmfaces(rdmds([diagDir filename11],ttAverage,'rec',5)) .* mmol_to_mol;  %mol m^-3 s^-1

    %2D fe 
    SFCSOLFe = convert2gcmfaces(rdmds([diagDir filename12],ttAverage,'rec',1)) .* mmol_to_mol;  %mol m^-2 s^-1
    
    %%
    %snapshots
    
    ETAN_SNAP = nan*ones(nx,ny,2);
    SALT_SNAP = nan .* ones(nx,ny,numLevels,2);
    Fe_SNAP = nan .* ones(nx,ny,numLevels,2);
    
    if timeStep == 1
        
        %no pre-initial snapshot
        ETAN_SNAP(:,:,2) = rdmds([diagDir filename8],tt(1),'rec',1);
        SALT_SNAP(:,:,:,2) = rdmds([diagDir filename7],tt(1),'rec',1);
        Fe_SNAP(:,:,:,2) = rdmds([diagDir filename7],tt(1),'rec',3) .* mmol_to_mol;  %mol m^-3 s^-1
        
    elseif timeStep == numFiles %no final snapshot
        
        ETAN_SNAP(:,:,1) = rdmds([diagDir filename8],tt(timeStep-1),'rec',1);
        SALT_SNAP(:,:,:,1) = rdmds([diagDir filename7],tt(timeStep-1),'rec',1);
        Fe_SNAP(:,:,:,1) = rdmds([diagDir filename7],tt(timeStep-1),'rec',3) .* mmol_to_mol;  %mol m^-3 s^-1
        
    else %timeStep~=1 & timeStep~=numFiles
        
        ETAN_SNAP = rdmds([diagDir filename8],ttSnap,'rec',1);
        SALT_SNAP = rdmds([diagDir filename7],ttSnap,'rec',1);
        Fe_SNAP = rdmds([diagDir filename7],ttSnap,'rec',3) .* mmol_to_mol;  %mol m^-3 s^-1
        
    end
    
    ETAN_SNAP = convert2gcmfaces(ETAN_SNAP);
    SALT_SNAP = convert2gcmfaces(SALT_SNAP);
    Fe_SNAP = convert2gcmfaces(Fe_SNAP);
    
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
    %Fe budget
    
    F_snap = 0 .* Fe_SNAP;
    
    for nt = 1:2
        
        F_snap(:,:,:,nt) = (Fe_SNAP(:,:,:,nt) .* (1+mk3D(ETAN_SNAP(:,:,nt) ./ mygrid.Depth,dzMat)));
        
    end
    
    %tendency, mol m^-3 s^-1
    tendFe = (F_snap(:,:,:,2) - F_snap(:,:,:,1)) ./ (secPerHour .* dt(timeStep));
    
    %horizontal divergences, mol m^-3 s^-1
    adv_hConvFe = calc_UV_conv(ADVx_Fe,ADVy_Fe) ./ VVV;
    dif_hConvFe = calc_UV_conv(DFxE_Fe,DFyE_Fe) ./ VVV;
    
    %vertical divergences, mol m^-3 s^-1
    adv_vConvFe = 0 .* ADVx_Fe;
    dif_vConvFe = 0 .* ADVx_Fe;
    
    for nz = 1:numLevels
        
        nzp1 = min([nz+1,numLevels]);
        
        adv_vConvFe(:,:,nz) = squeeze(ADVr_Fe(:,:,nzp1) .* double(nz~=numLevels) - ADVr_Fe(:,:,nz));
        dif_vConvFe(:,:,nz) = squeeze(DFrI_Fe(:,:,nzp1) .* double(nz~=numLevels) - DFrI_Fe(:,:,nz) ...
            + DFrE_Fe(:,:,nzp1) .* double(nz~=numLevels) - DFrE_Fe(:,:,nz));
        
    end
    
    adv_vConvFe = adv_vConvFe ./ VVV;
    dif_vConvFe = dif_vConvFe ./ VVV;
    
    %air-sea CO2 flux tendency, mol m^-3 s^-1
    forcFe = 0 .* oceSPtnd;
    
    for nz = 1:numLevels
        
        if nz == 1
            
            forcFe(:,:,1) =  SFCSOLFe ./ dzMat(:,:,1);
            
        else
            
            forcFe(:,:,nz) = 0;
            
        end
        
    end
    
    DARWIN_Fe = DARWIN_Fe - forcFe - SEDFe; %remove iron dust and sediment from gDAR
    
    forcFe = mygrid.mskC .* forcFe;
    sedFe = mygrid.mskC .* SEDFe;
    freeFe = mygrid.mskC .* FREEFe;
    
    %darwin, mol m^-3 s^-1
    darwinFe = mygrid.mskC .* DARWIN_Fe;
    
    surfFe = mygrid.mskC(:,:,1) .* Fe(:,:,1);
    meanSurf_Fe = nansum(surfFe(:) .* mygrid.RAC(:)) ./ nansum(mygrid.RAC(:));
    
    %z* correction, as done in salinity budget
    
    tendFeZs = mygrid.mskC .* (Fe_SNAP(:,:,:,2) - Fe_SNAP(:,:,:,1)) ...
        ./ (secPerHour .* dt(timeStep));
    
    %advection
    adv_hConvFeZs = (-Fe .* adv_hConvV + adv_hConvFe) ./ rstarfac;
    adv_vConvFeZs = (-Fe .* adv_vConvV + adv_vConvFe) ./ rstarfac;
    
    %diffusion
    dif_vConvFeZs = dif_vConvFe ./ rstarfac;
    dif_hConvFeZs = dif_hConvFe ./ rstarfac;
    
    %forcing
    forcFeZs = forcFe ./ rstarfac;
    
    %sed
    sedFeZs = sedFe ./ rstarfac;
    
    %free
    freeFeZs = freeFe ./ rstarfac;
    
    %biology
    darwinFeZs =  darwinFe ./ rstarfac;
    
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
    
    tendFe = convert2gcmfaces(tendFeZs);
    adv_hConvFe = convert2gcmfaces(adv_hConvFeZs);
    adv_vConvFe = convert2gcmfaces(adv_vConvFeZs);
    dif_hConvFe = convert2gcmfaces(dif_hConvFeZs);
    dif_vConvFe = convert2gcmfaces(dif_vConvFeZs);
    forcFe = convert2gcmfaces(forcFeZs);
    sedFe = convert2gcmfaces(sedFeZs);
    freeFe = convert2gcmfaces(freeFeZs);
    darwinFe = convert2gcmfaces(darwinFeZs);
    
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
    
    intTendFe = nansum(tendFe(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel) .* convert2gcmfaces(mygrid.hFacC(:,:,1:intLevel))),3);
    intHAdvFe = nansum(adv_hConvFe(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel) .* convert2gcmfaces(mygrid.hFacC(:,:,1:intLevel))),3);
    intVAdvFe = nansum(adv_vConvFe(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel) .* convert2gcmfaces(mygrid.hFacC(:,:,1:intLevel))),3);
    intHDifFe = nansum(dif_hConvFe(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel) .* convert2gcmfaces(mygrid.hFacC(:,:,1:intLevel))),3);
    intVDifFe = nansum(dif_vConvFe(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel) .* convert2gcmfaces(mygrid.hFacC(:,:,1:intLevel))),3);
    intForcFe = nansum(forcFe(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel) .* convert2gcmfaces(mygrid.hFacC(:,:,1:intLevel))),3);
    intSedFe = nansum(sedFe(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel) .* convert2gcmfaces(mygrid.hFacC(:,:,1:intLevel))),3);
    intFreeFe = nansum(freeFe(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel) .* convert2gcmfaces(mygrid.hFacC(:,:,1:intLevel))),3);
    intBioFe = nansum(darwinFe(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel) .* convert2gcmfaces(mygrid.hFacC(:,:,1:intLevel))),3);
    
    %%
    %plot budgets
    
    if plotVolumeBudget
        
        hFig1 = figure(1);
        set(hFig1,'units','normalized','outerposition',[0 0 1 1]);
        set(gcf,'color',[1 1 1]);
        set(gca,'color',[0.5 0.5 0.5]);
        
        cMin = -10^-8;
        cMax = 10^-8;
        
        subplot(161);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,intTendV);
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
        
        title({'Tendency (m s^-^1), ';'Timestep: ';num2str(tt(timeStep))},'FontWeight','Bold','FontSize',fs);
        
        subplot(162);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,intHConvV);
        shading flat
        
        caxis([cMin .* 10^3 cMax .* 10^3]);
        
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
        
        title({'Horizontal';'Advection (m s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(163);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,intVConvV);
        shading flat
        
        caxis([cMin .* 10^3 cMax .* 10^3]);
        
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
        
        title({'Vertical';'Advection (m s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(164);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,intForcV);
        shading flat
        
        caxis([cMin .* 10 cMax .* 10]);
        
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
        
        title({'Surface Volume';'Forcing (m s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(165);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,intHConvV + intVConvV + intForcV);
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
        
        title({'Total';'Budget (m s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(166);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,intTendV - (intForcV + intHConvV + intVConvV));
        shading flat
        
        caxis([cMin .* 10^-3 cMax .* 10^-3]);
        
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
        
        title('Residual (m s^-^1)','FontWeight','Bold','FontSize',fs);
        
        drawnow
        
    end
    
    if plotSalinityBudget
        
        hFig2 = figure(2);
        set(hFig2,'units','normalized','outerposition',[0 0 1 1]);
        set(gcf,'color',[1 1 1]);
        
        cMin = -10^-5;
        cMax = 10^-5;
        
        subplot(281);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,intTendS);
        shading flat
        
        caxis([cMin cMax]);
        
        colormap(colors);
        
        colorbar
        
        axis tight
        
        xlim([0 360]);
        ylim([-90 90]);
        
        set(gca,'xtick',[0:180:360]);
        set(gca,'ytick',[-90:30:90]);
        
        ylabel('Latitude (deg)');
        
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Salt Tendency (psu m s^-^1)';'Timestep: ';num2str(tt(timeStep))},'FontWeight','Bold','FontSize',fs);
        
        subplot(282);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,intHAdvS);
        shading flat
        
        caxis([cMin .* 10 cMax .* 10]);
        
        colormap(colors);
        
        colorbar
        
        axis tight
        
        xlim([0 360]);
        ylim([-90 90]);
        
        set(gca,'xtick',[0:180:360]);
        set(gca,'ytick',[-90:30:90]);
        
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Horizontal';'Advection';'(psu m s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(283);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,intVAdvS);
        shading flat
        
        caxis([cMin .* 10 cMax .* 10]);
        
        colormap(colors);
        
        colorbar
        
        axis tight
        
        xlim([0 360]);
        ylim([-90 90]);
        
        set(gca,'xtick',[0:180:360]);
        set(gca,'ytick',[-90:30:90]);
        
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Vertical';'Advection';'(psu m s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(284);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,intHDifS);
        shading flat
        
        caxis([cMin cMax]);
        
        colormap(colors);
        
        colorbar
        
        axis tight
        
        xlim([0 360]);
        ylim([-90 90]);
        
        set(gca,'xtick',[0:180:360]);
        set(gca,'ytick',[-90:30:90]);
        
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Horizontal';'Diffusion';'(psu m s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(285);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,intVDifS);
        shading flat
        
        caxis([cMin cMax]);
        
        colormap(colors);
        
        colorbar
        
        axis tight
        
        xlim([0 360]);
        ylim([-90 90]);
        
        set(gca,'xtick',[0:180:360]);
        set(gca,'ytick',[-90:30:90]);
        
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Vertical';'Diffusion';'(psu m s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(286);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,intForcS);
        shading flat
        
        caxis([cMin .* 10^-2 cMax .* 10^-2]);
        
        colormap(colors);
        
        colorbar
        
        axis tight
        
        xlim([0 360]);
        ylim([-90 90]);
        
        set(gca,'xtick',[0:180:360]);
        set(gca,'ytick',[-90:30:90]);
        
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Surface';'Volume';'Forcing';'(psu m s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(287);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,intHAdvS + intVAdvS + intHDifS + intVDifS + intForcS);
        shading flat
        
        caxis([cMin cMax]);
        
        colormap(colors);
        
        colorbar
        
        axis tight
        
        xlim([0 360]);
        ylim([-90 90]);
        
        set(gca,'xtick',[0:180:360]);
        set(gca,'ytick',[-90:30:90]);
        
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Budget';'Total';'(psu m s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(288);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,intTendS - (intHAdvS + intVAdvS + intHDifS + intVDifS + intForcS));
        shading flat
        
        caxis([cMin .* 10^-3 cMax .* 10^-3]);
        
        colormap(colors);
        
        colorbar
        
        axis tight
        
        xlim([0 360]);
        ylim([-90 90]);
        
        set(gca,'xtick',[0:180:360]);
        set(gca,'ytick',[-90:30:90]);
        
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Budget Total - Tend';'(psu m s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(289);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,intTendSal);
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
        
        title({'Salinity Tendency (psu m s^-^1),';'Timestep: ';num2str(tt(timeStep))},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,8,10);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,intHAdvSal);
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
        
        subplot(2,8,11);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,intVAdvSal);
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
        
        subplot(2,8,12);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,intHDifSal);
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
        
        subplot(2,8,13);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,intVDifSal);
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
        
        subplot(2,8,14);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,intForcSal);
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
        
        subplot(2,8,15);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,(intHAdvSal + intVAdvSal + intHDifSal + intVDifSal + intForcSal));
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
        
        subplot(2,8,16);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,intTendSal - (intHAdvSal + intVAdvSal + intHDifSal + intVDifSal + intForcSal));
        shading flat
        
        caxis([cMin .* 10^-3 cMax .* 10^-3]);
        
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
        
    end
    
    if plotFeBudget
        
        hFig3 = figure(3);
        set(hFig3,'units','normalized','outerposition',[0 0 1 1]);
        set(gcf,'color',[1 1 1]);
        set(gca,'color',[0.5 0.5 0.5]);
        
        cMin = -10^-11;
        cMax = 10^-11;
        
        subplot(2,5,1);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,intTendFe);
        shading flat
        
        caxis([cMin cMax]);
        
        colormap(colors);
        
        colorbar
        
        axis tight
        
        xlim([0 360]);
        ylim([-90 90]);
        
        set(gca,'xtick',[0:180:360]);
        set(gca,'ytick',[-90:30:90]);
        
        ylabel('Latitude (deg)');
        
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'dFe/dt (mol m^-^2 s^-^1), ';'timeStep: ';num2str(tt(timeStep))},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,2);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,intHAdvFe);
        shading flat
        
        caxis([cMin cMax]);
        
        colormap(colors);
        
        colorbar
        
        axis tight
        
        xlim([0 360]);
        ylim([-90 90]);
        
        set(gca,'xtick',[0:180:360]);
        set(gca,'ytick',[-90:30:90]);
        
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Horizontal Advection';'(mol m^-^2 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,3);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,intVAdvFe);
        shading flat
        
        caxis([cMin cMax]);
        
        colormap(colors);
        
        colorbar
        
        axis tight
        
        xlim([0 360]);
        ylim([-90 90]);
        
        set(gca,'xtick',[0:180:360]);
        set(gca,'ytick',[-90:30:90]);
        
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Vertical Advection';'(mol m^-^2 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,4);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,intHDifFe);
        shading flat
        
        caxis([cMin cMax]);
        
        colormap(colors);
        
        colorbar
        
        axis tight
        
        xlim([0 360]);
        ylim([-90 90]);
        
        set(gca,'xtick',[0:180:360]);
        set(gca,'ytick',[-90:30:90]);
        
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Horizontal Diffusion';'(mol m^-^2 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,5);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,intVDifFe);
        shading flat
        
        caxis([cMin cMax]);
        
        colormap(colors);
        
        colorbar
        
        axis tight
        
        xlim([0 360]);
        ylim([-90 90]);
        
        set(gca,'xtick',[0:180:360]);
        set(gca,'ytick',[-90:30:90]);
        
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Vertical Diffusion';'(mol m^-^2 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,6);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
       
        pcolor(mygrid.XC.f1,mygrid.YC.f1,intForcFe);
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
        
        title({'Surface Dust Flux';'(mol m^-^2 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,7);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,intSedFe);
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
        
        title({'Sediment Flux';'(mol m^-^2 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,8);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,intBioFe);
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
        
        title({'Scav + Bio';'(mol m^-^2 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,9);
        
          hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,intFreeFe);
        
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
        
        title({'Free Iron';'(mol m^-^2 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,10);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolor(mygrid.XC.f1,mygrid.YC.f1,intTendFe - (intHAdvFe + intVAdvFe + intHDifFe + intVDifFe ...
            + intForcFe + intSedFe + intFreeFe + intBioFe));
        
        shading flat
        
        caxis([cMin.*10^-2 cMax.*10^-2]);
        
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
        
        title({'Budget Residual';'(mol m^-^2 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        drawnow
        
    end
    
    %clear gcmfaces objects to avoid memory leaks
    
    clear ETAN oceFWflx SFLUX  oceSPflx UVELMASS VVELMASS WVELMASS ...
    SALT ADVr_SLT ADVx_SLT ADVy_SLT DFrI_SLT DFrE_SLT DFxE_SLT DFyE_SLT oceSPtnd ...
    Fe ADVx_Fe ADVy_Fe ADVr_Fe DFxE_Fe DFyE_Fe DFrE_Fe DFrI_Fe FeTFLX ...
    CONSUMP_Fe CONSUMP_Fe_PIC REMIN_DOC REMIN_POC DISSC_PIC ...
    ETAN_SNAP SALT_SNAP Fe_SNAP ...
    tendV adv_hConvV adv_vConvV forcV ...
    S_snap tendS adv_hConvS dif_hConvS adv_vConvS dif_vConvS forcS surfS ...
    rstarfac tendSal adv_hConvSal adv_vConvSal dif_hConvSal dif_vConvSal forcSal ...
    tendFe adv_hConvFe dif_hConvFe adv_vConvFe dif_vConvFe forcFe darwinFe surfFe virtualFluxFe ...
    tendFeZs adv_hConvFeZs adv_vConvFeZs dif_vConvFeZs dif_hConvFeZs forcFeZs darwinFeZs virtualFluxFeZs   
    
end
