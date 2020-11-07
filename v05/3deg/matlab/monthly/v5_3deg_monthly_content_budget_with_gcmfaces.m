clear
close all

%%
%settings, modify as needed

savePlot = 0;
intLevel = 2; %integration k level

%set to 1 to plot budget terms
plotVolumeBudget = 0;
plotSalinityBudget = 0;

plotDICBudget = 0;
plotNO3Budget = 0;
plotNO2Budget = 0;
plotNH4Budget = 0;
plotPO4Budget = 0;
plotFeBudget = 1;
plotSiO2Budget = 0;

gridDir = '../../../../../darwin3/run/';
modelDir = '../../../../../darwin3/run/';

figureDir = '/Users/carrolld/Documents/research/v05_budget/figures/v5_3deg_budget/';

%%
%plotting params.

fs = 12;
lw = 2;

colors = parula(500);
colors = flipud(cbrewer('div','RdBu',500));
colors = cmocean('balance',500);

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

mygrid.domainPeriodicity=[1 0];

nx = mygrid.ioSize(1);
ny = mygrid.ioSize(2);
numLevels = numel(mygrid.RC);

dzMatF = mk3D(mygrid.DRF, mygrid.hFacC);
dzMat = dzMatF .* mygrid.hFacC;
RACMat = mk3D(mygrid.RAC, mygrid.hFacC);
VVV = mygrid.mskC .* mygrid.hFacC .* mk3D(mygrid.RAC,mygrid.mskC) .* mk3D(mygrid.DRF,mygrid.mskC);
VVV2 = mygrid.mskC .* mk3D(mygrid.RAC,mygrid.mskC) .* mk3D(mygrid.DRF,mygrid.mskC);

%%

diagDir = [modelDir 'diags/budget/'];

filename1 = 'average_2d';
filename2 = 'average_velmass_3d';
filename3 = 'average_salt_3d';

filename4 = 'average_DIC_3d';
filename5 = 'average_NO3_3d';
filename6 = 'average_NO2_3d';
filename7 = 'average_NH4_3d';
filename8 = 'average_PO4_3d';
filename9 = 'average_Fe_3d';
filename10 = 'average_Fe_darwin_2d';
filename11 = 'average_SiO2_3d';

filename12 = 'snap_2d';
filename13 = 'snap_3d';

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
    NO3_Epr = convert2gcmfaces(rdmds([diagDir filename1],ttAverage,'rec',8)) .* mmol_to_mol; %mol m^-3 s^-1
    NO2_Epr = convert2gcmfaces(rdmds([diagDir filename1],ttAverage,'rec',9)) .* mmol_to_mol; %mol m^-3 s^-1
    NH4_Epr = convert2gcmfaces(rdmds([diagDir filename1],ttAverage,'rec',10)) .* mmol_to_mol; %mol m^-3 s^-1
    PO4_Epr = convert2gcmfaces(rdmds([diagDir filename1],ttAverage,'rec',11)) .* mmol_to_mol; %mol m^-3 s^-1
    Fe_Epr = convert2gcmfaces(rdmds([diagDir filename1],ttAverage,'rec',12)) .* mmol_to_mol; %mol m^-3 s^-1
    SiO2_Epr = convert2gcmfaces(rdmds([diagDir filename1],ttAverage,'rec',13)) .* mmol_to_mol; %mol m^-3 s^-1
    
    %load three-dimensional time-averaged fields
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
    gDAR_DIC = convert2gcmfaces(rdmds([diagDir filename4],ttAverage,'rec',9)) .* mmol_to_mol;  %mol m^-3 s^-1
    
    %NO3 content
    NO3 = convert2gcmfaces(rdmds([diagDir filename5],ttAverage,'rec',1)) .* mmol_to_mol; %mol m^-3
    ADVx_NO3 = convert2gcmfaces(rdmds([diagDir filename5],ttAverage,'rec',2)) .* mmol_to_mol; %mol s^-1
    ADVy_NO3 = convert2gcmfaces(rdmds([diagDir filename5],ttAverage,'rec',3)) .* mmol_to_mol; %mol s^-1
    ADVr_NO3 = convert2gcmfaces(rdmds([diagDir filename5],ttAverage,'rec',4)) .* mmol_to_mol; %mol s^-1
    DFxE_NO3 = convert2gcmfaces(rdmds([diagDir filename5],ttAverage,'rec',5)) .* mmol_to_mol; %mol s^-1
    DFyE_NO3 = convert2gcmfaces(rdmds([diagDir filename5],ttAverage,'rec',6)) .* mmol_to_mol; %mol s^-1
    DFrE_NO3 = convert2gcmfaces(rdmds([diagDir filename5],ttAverage,'rec',7)) .* mmol_to_mol; %mol s^-1
    DFrI_NO3 = convert2gcmfaces(rdmds([diagDir filename5],ttAverage,'rec',8)) .* mmol_to_mol; %mol s^-1
    gDAR_NO3 = convert2gcmfaces(rdmds([diagDir filename5],ttAverage,'rec',9)) .* mmol_to_mol; %mol s^-1
    C_NO3 = convert2gcmfaces(rdmds([diagDir filename5],ttAverage,'rec',10)) .* mmol_to_mol;  %mol m^-3 s^-1
    S_NO3 = convert2gcmfaces(rdmds([diagDir filename5],ttAverage,'rec',11)) .* mmol_to_mol;  %mol m^-3 s^-1
    
    %NO2 content
    NO2 = convert2gcmfaces(rdmds([diagDir filename6],ttAverage,'rec',1)) .* mmol_to_mol; %mol m^-3
    ADVx_NO2 = convert2gcmfaces(rdmds([diagDir filename6],ttAverage,'rec',2)) .* mmol_to_mol; %mol s^-1
    ADVy_NO2 = convert2gcmfaces(rdmds([diagDir filename6],ttAverage,'rec',3)) .* mmol_to_mol; %mol s^-1
    ADVr_NO2 = convert2gcmfaces(rdmds([diagDir filename6],ttAverage,'rec',4)) .* mmol_to_mol; %mol s^-1
    DFxE_NO2 = convert2gcmfaces(rdmds([diagDir filename6],ttAverage,'rec',5)) .* mmol_to_mol; %mol s^-1
    DFyE_NO2 = convert2gcmfaces(rdmds([diagDir filename6],ttAverage,'rec',6)) .* mmol_to_mol; %mol s^-1
    DFrE_NO2 = convert2gcmfaces(rdmds([diagDir filename6],ttAverage,'rec',7)) .* mmol_to_mol; %mol s^-1
    DFrI_NO2 = convert2gcmfaces(rdmds([diagDir filename6],ttAverage,'rec',8)) .* mmol_to_mol; %mol s^-1
    gDAR_NO2 = convert2gcmfaces(rdmds([diagDir filename6],ttAverage,'rec',9)) .* mmol_to_mol; %mol s^-1
    C_NO2 = convert2gcmfaces(rdmds([diagDir filename6],ttAverage,'rec',10)) .* mmol_to_mol;  %mol m^-3 s^-1
    S_NO2 = convert2gcmfaces(rdmds([diagDir filename6],ttAverage,'rec',11)) .* mmol_to_mol;  %mol m^-3 s^-1
    
    %NH4 content
    NH4 = convert2gcmfaces(rdmds([diagDir filename7],ttAverage,'rec',1)) .* mmol_to_mol; %mol m^-3
    ADVx_NH4 = convert2gcmfaces(rdmds([diagDir filename7],ttAverage,'rec',2)) .* mmol_to_mol; %mol s^-1
    ADVy_NH4 = convert2gcmfaces(rdmds([diagDir filename7],ttAverage,'rec',3)) .* mmol_to_mol; %mol s^-1
    ADVr_NH4 = convert2gcmfaces(rdmds([diagDir filename7],ttAverage,'rec',4)) .* mmol_to_mol; %mol s^-1
    DFxE_NH4 = convert2gcmfaces(rdmds([diagDir filename7],ttAverage,'rec',5)) .* mmol_to_mol; %mol s^-1
    DFyE_NH4 = convert2gcmfaces(rdmds([diagDir filename7],ttAverage,'rec',6)) .* mmol_to_mol; %mol s^-1
    DFrE_NH4 = convert2gcmfaces(rdmds([diagDir filename7],ttAverage,'rec',7)) .* mmol_to_mol; %mol s^-1
    DFrI_NH4 = convert2gcmfaces(rdmds([diagDir filename7],ttAverage,'rec',8)) .* mmol_to_mol; %mol s^-1
    gDAR_NH4 = convert2gcmfaces(rdmds([diagDir filename7],ttAverage,'rec',9)) .* mmol_to_mol; %mol s^-1
    C_NH4 = convert2gcmfaces(rdmds([diagDir filename7],ttAverage,'rec',10)) .* mmol_to_mol;  %mol m^-3 s^-1
    S_NH4 = convert2gcmfaces(rdmds([diagDir filename7],ttAverage,'rec',11)) .* mmol_to_mol;  %mol m^-3 s^-1
    
    %PO4 content
    PO4 = convert2gcmfaces(rdmds([diagDir filename8],ttAverage,'rec',1)) .* mmol_to_mol; %mol m^-3
    ADVx_PO4 = convert2gcmfaces(rdmds([diagDir filename8],ttAverage,'rec',2)) .* mmol_to_mol; %mol s^-1
    ADVy_PO4 = convert2gcmfaces(rdmds([diagDir filename8],ttAverage,'rec',3)) .* mmol_to_mol; %mol s^-1
    ADVr_PO4 = convert2gcmfaces(rdmds([diagDir filename8],ttAverage,'rec',4)) .* mmol_to_mol; %mol s^-1
    DFxE_PO4 = convert2gcmfaces(rdmds([diagDir filename8],ttAverage,'rec',5)) .* mmol_to_mol; %mol s^-1
    DFyE_PO4 = convert2gcmfaces(rdmds([diagDir filename8],ttAverage,'rec',6)) .* mmol_to_mol; %mol s^-1
    DFrE_PO4 = convert2gcmfaces(rdmds([diagDir filename8],ttAverage,'rec',7)) .* mmol_to_mol; %mol s^-1
    DFrI_PO4 = convert2gcmfaces(rdmds([diagDir filename8],ttAverage,'rec',8)) .* mmol_to_mol; %mol s^-1
    gDAR_PO4 = convert2gcmfaces(rdmds([diagDir filename8],ttAverage,'rec',9)) .* mmol_to_mol; %mol s^-1
    C_PO4 = convert2gcmfaces(rdmds([diagDir filename8],ttAverage,'rec',10)) .* mmol_to_mol;  %mol m^-3 s^-1
    S_PO4 = convert2gcmfaces(rdmds([diagDir filename8],ttAverage,'rec',11)) .* mmol_to_mol;  %mol m^-3 s^-1
    
    %Fe content
    Fe = convert2gcmfaces(rdmds([diagDir filename9],ttAverage,'rec',1)) .* mmol_to_mol; %mol m^-3
    ADVx_Fe = convert2gcmfaces(rdmds([diagDir filename9],ttAverage,'rec',2)) .* mmol_to_mol; %mol s^-1
    ADVy_Fe = convert2gcmfaces(rdmds([diagDir filename9],ttAverage,'rec',3)) .* mmol_to_mol; %mol s^-1
    ADVr_Fe = convert2gcmfaces(rdmds([diagDir filename9],ttAverage,'rec',4)) .* mmol_to_mol; %mol s^-1
    DFxE_Fe = convert2gcmfaces(rdmds([diagDir filename9],ttAverage,'rec',5)) .* mmol_to_mol; %mol s^-1
    DFyE_Fe = convert2gcmfaces(rdmds([diagDir filename9],ttAverage,'rec',6)) .* mmol_to_mol; %mol s^-1
    DFrE_Fe = convert2gcmfaces(rdmds([diagDir filename9],ttAverage,'rec',7)) .* mmol_to_mol; %mol s^-1
    DFrI_Fe = convert2gcmfaces(rdmds([diagDir filename9],ttAverage,'rec',8)) .* mmol_to_mol; %mol s^-1
    gDAR_Fe = convert2gcmfaces(rdmds([diagDir filename9],ttAverage,'rec',9)) .* mmol_to_mol;  %mol m^-3 s^-1
    C_Fe = convert2gcmfaces(rdmds([diagDir filename9],ttAverage,'rec',10)) .* mmol_to_mol;  %mol m^-3 s^-1
    S_Fe = convert2gcmfaces(rdmds([diagDir filename9],ttAverage,'rec',11)) .* mmol_to_mol;  %mol m^-3 s^-1
    SEDFe = convert2gcmfaces(rdmds([diagDir filename9],ttAverage,'rec',12)) .* mmol_to_mol;  %mol m^-3 s^-1
    FREEFe = convert2gcmfaces(rdmds([diagDir filename9],ttAverage,'rec',13)) .* mmol_to_mol;  %mol m^-3 s^-1
    SFCSOLFe = convert2gcmfaces(rdmds([diagDir filename10],ttAverage,'rec',1)) .* mmol_to_mol;  %mol m^-2 s^-1
    
    %SiO2 content
    SiO2 = convert2gcmfaces(rdmds([diagDir filename11],ttAverage,'rec',1)) .* mmol_to_mol; %mol m^-3
    ADVx_SiO2 = convert2gcmfaces(rdmds([diagDir filename11],ttAverage,'rec',2)) .* mmol_to_mol; %mol s^-1
    ADVy_SiO2 = convert2gcmfaces(rdmds([diagDir filename11],ttAverage,'rec',3)) .* mmol_to_mol; %mol s^-1
    ADVr_SiO2 = convert2gcmfaces(rdmds([diagDir filename11],ttAverage,'rec',4)) .* mmol_to_mol; %mol s^-1
    DFxE_SiO2 = convert2gcmfaces(rdmds([diagDir filename11],ttAverage,'rec',5)) .* mmol_to_mol; %mol s^-1
    DFyE_SiO2 = convert2gcmfaces(rdmds([diagDir filename11],ttAverage,'rec',6)) .* mmol_to_mol; %mol s^-1
    DFrE_SiO2 = convert2gcmfaces(rdmds([diagDir filename11],ttAverage,'rec',7)) .* mmol_to_mol; %mol s^-1
    DFrI_SiO2 = convert2gcmfaces(rdmds([diagDir filename11],ttAverage,'rec',8)) .* mmol_to_mol; %mol s^-1
    gDAR_SiO2 = convert2gcmfaces(rdmds([diagDir filename11],ttAverage,'rec',9)) .* mmol_to_mol;  %mol m^-3 s^-1
    C_SiO2 = convert2gcmfaces(rdmds([diagDir filename11],ttAverage,'rec',10)) .* mmol_to_mol;  %mol m^-3 s^-1
    S_SiO2 = convert2gcmfaces(rdmds([diagDir filename11],ttAverage,'rec',11)) .* mmol_to_mol;  %mol m^-3 s^-1
    
    %%
    %load snapshots
    
    ETAN_SNAP = nan*ones(nx,ny,2);
    SALT_SNAP = nan .* ones(nx,ny,numLevels,2);
    DIC_SNAP = nan .* ones(nx,ny,numLevels,2);
    NO3_SNAP = nan .* ones(nx,ny,numLevels,2);
    NO2_SNAP = nan .* ones(nx,ny,numLevels,2);
    NH4_SNAP = nan .* ones(nx,ny,numLevels,2);
    PO4_SNAP = nan .* ones(nx,ny,numLevels,2);
    Fe_SNAP = nan .* ones(nx,ny,numLevels,2);
    SiO2_SNAP = nan .* ones(nx,ny,numLevels,2);
    
    if timeStep == 1
        
        %no pre-initial snapshot
        ETAN_SNAP(:,:,2) = rdmds([diagDir filename12],tt(1),'rec',1);
        
        SALT_SNAP(:,:,:,2) = rdmds([diagDir filename13],tt(1),'rec',1);
        DIC_SNAP(:,:,:,2) = rdmds([diagDir filename13],tt(1),'rec',2) .* mmol_to_mol;  %mol m^-3 s^-1
        NO3_SNAP(:,:,:,2) = rdmds([diagDir filename13],tt(1),'rec',3) .* mmol_to_mol;  %mol m^-3 s^-1
        NO2_SNAP(:,:,:,2) = rdmds([diagDir filename13],tt(1),'rec',4) .* mmol_to_mol;  %mol m^-3 s^-1
        NH4_SNAP(:,:,:,2) = rdmds([diagDir filename13],tt(1),'rec',5) .* mmol_to_mol;  %mol m^-3 s^-1
        PO4_SNAP(:,:,:,2) = rdmds([diagDir filename13],tt(1),'rec',6) .* mmol_to_mol;  %mol m^-3 s^-1
        Fe_SNAP(:,:,:,2) = rdmds([diagDir filename13],tt(1),'rec',7) .* mmol_to_mol;  %mol m^-3 s^-1
        SiO2_SNAP(:,:,:,2) = rdmds([diagDir filename13],tt(1),'rec',8) .* mmol_to_mol;  %mol m^-3 s^-1
        
    elseif timeStep == numFiles %no final snapshot
        
        ETAN_SNAP(:,:,1) = rdmds([diagDir filename12],tt(timeStep-1),'rec',1);
        
        SALT_SNAP(:,:,:,1) = rdmds([diagDir filename13],tt(timeStep-1),'rec',1);
        DIC_SNAP(:,:,:,1) = rdmds([diagDir filename13],tt(timeStep-1),'rec',2) .* mmol_to_mol;  %mol m^-3 s^-1
        NO3_SNAP(:,:,:,1) = rdmds([diagDir filename13],tt(timeStep-1),'rec',3) .* mmol_to_mol;  %mol m^-3 s^-1
        NO2_SNAP(:,:,:,1) = rdmds([diagDir filename13],tt(timeStep-1),'rec',4) .* mmol_to_mol;  %mol m^-3 s^-1
        NH4_SNAP(:,:,:,1) = rdmds([diagDir filename13],tt(timeStep-1),'rec',5) .* mmol_to_mol;  %mol m^-3 s^-1
        PO4_SNAP(:,:,:,1) = rdmds([diagDir filename13],tt(timeStep-1),'rec',6) .* mmol_to_mol;  %mol m^-3 s^-1
        Fe_SNAP(:,:,:,1) = rdmds([diagDir filename13],tt(timeStep-1),'rec',7) .* mmol_to_mol;  %mol m^-3 s^-1
        SiO2_SNAP(:,:,:,1) = rdmds([diagDir filename13],tt(timeStep-1),'rec',8) .* mmol_to_mol;  %mol m^-3 s^-1
        
    else %timeStep~=1 & timeStep~=numFiles
        
        ETAN_SNAP = rdmds([diagDir filename12],ttSnap,'rec',1);
        
        SALT_SNAP = rdmds([diagDir filename13],ttSnap,'rec',1);
        DIC_SNAP = rdmds([diagDir filename13],ttSnap,'rec',2) .* mmol_to_mol;  %mol m^-3 s^-1
        NO3_SNAP = rdmds([diagDir filename13],ttSnap,'rec',3) .* mmol_to_mol;  %mol m^-3 s^-1
        NO2_SNAP = rdmds([diagDir filename13],ttSnap,'rec',4) .* mmol_to_mol;  %mol m^-3 s^-1
        NH4_SNAP = rdmds([diagDir filename13],ttSnap,'rec',5) .* mmol_to_mol;  %mol m^-3 s^-1
        PO4_SNAP = rdmds([diagDir filename13],ttSnap,'rec',6) .* mmol_to_mol;  %mol m^-3 s^-1
        Fe_SNAP = rdmds([diagDir filename13],ttSnap,'rec',7) .* mmol_to_mol;  %mol m^-3 s^-1
        SiO2_SNAP = rdmds([diagDir filename13],ttSnap,'rec',8) .* mmol_to_mol;  %mol m^-3 s^-1
        
    end
    
    ETAN_SNAP = convert2gcmfaces(ETAN_SNAP);
    SALT_SNAP = convert2gcmfaces(SALT_SNAP);
    DIC_SNAP = convert2gcmfaces(DIC_SNAP);
    NO3_SNAP = convert2gcmfaces(NO3_SNAP);
    NO2_SNAP = convert2gcmfaces(NO2_SNAP);
    NH4_SNAP = convert2gcmfaces(NH4_SNAP);
    PO4_SNAP = convert2gcmfaces(PO4_SNAP);
    Fe_SNAP = convert2gcmfaces(Fe_SNAP);
    SiO2_SNAP = convert2gcmfaces(SiO2_SNAP);
    
    %%
    %volume budget, s^-1
    
    sStarMean = (mygrid.Depth + ETAN) ./ mygrid.Depth;
    
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
    
    %%
    %salinity budget
    
    rstarfac = (mygrid.Depth + ETAN) ./ mygrid.Depth;
    
    %tendency
    tendSal = mygrid.mskC .* (SALT_SNAP(:,:,:,2) - SALT_SNAP(:,:,:,1)) ...
        ./ (secPerHour .* dt(timeStep));
    
    %advection
    adv_hConvSal = ((-SALT .* adv_hConvV) + adv_hConvS) ./ rstarfac;
    adv_vConvSal = ((-SALT .* adv_vConvV) + adv_vConvS) ./ rstarfac;
    
    %diffusion
    dif_hConvSal = dif_hConvS ./ rstarfac;
    dif_vConvSal = dif_vConvS ./ rstarfac;
    
    %forcing
    forcSal = ((-SALT .* forcV) + forcS) ./ rstarfac;
    
    %%
    %DIC budget
    
    D_snap = 0 .* DIC_SNAP;
    sStar = 0 .* DIC_SNAP;
    
    for nt = 1:2
        
        sStar(:,:,:,nt) = (1+mk3D(ETAN_SNAP(:,:,nt) ./ mygrid.Depth,dzMat));
        
        D_snap(:,:,:,nt) = sStar(:,:,:,nt) .* DIC_SNAP(:,:,:,nt);
        
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
    
    %DIC tendency due to air-sea CO2 flux, mol m^-3 s^-1
    forcDIC1 = 0 .* oceSPtnd;
    forcDIC2 = 0 .* oceSPtnd;
    
    %DIC tendency due to E/P/runoff, mol m^-3 s^-1
    virtualFluxDIC = 0 .* oceSPtnd;
    
    for nz = 1:numLevels
        
        if nz == 1
            
            forcDIC1(:,:,1) = DICTFLX;
            forcDIC2(:,:,1) = fluxCO2;
            
            virtualFluxDIC(:,:,1) =  DIC_Epr;
            
        else
            
            forcDIC1(:,:,nz) = 0;
            forcDIC2(:,:,nz) = 0;
            
            virtualFluxDIC(:,:,nz) =  0;
            
        end
        
    end
    
    gDAR_DIC = (gDAR_DIC - (forcDIC2 ./ dzMat)) ./ mygrid.hFacC; %remove air-sea CO2 flux from gDAR, so it is just biology
    %gDAR_DIC = mygrid.mskC .* (gDAR_DIC ./ mygrid.hFacC); %remove air-sea CO2 flux from gDAR, so it is just biology
    
    forcDIC = forcDIC2 ./ dzMat;
    forcDIC = mygrid.mskC .* forcDIC;
    virtualFluxDIC = mygrid.mskC .* virtualFluxDIC .* 0;
    
    %biology, mol m^-3 s^-1
    bioDIC = mygrid.mskC .* gDAR_DIC;
    
    %%
    %NO3 budget
    
    D_snap = 0 .* NO3_SNAP;
    sStar = 0 .* NO3_SNAP;
    
    for nt = 1:2
        
        sStar(:,:,:,nt) = (1+mk3D(ETAN_SNAP(:,:,nt) ./ mygrid.Depth,dzMat));
        
        D_snap(:,:,:,nt) = sStar(:,:,:,nt) .* NO3_SNAP(:,:,:,nt);
        
    end
    
    %tendency, mol m^-3 s^-1
    tendNO3 = (D_snap(:,:,:,2) - D_snap(:,:,:,1)) ./ (secPerHour .* dt(timeStep));
    
    %horizontal divergences, mol m^-3 s^-1
    adv_hConvNO3 = calc_UV_conv(ADVx_NO3,ADVy_NO3) ./ VVV;
    dif_hConvNO3 = calc_UV_conv(DFxE_NO3,DFyE_NO3) ./ VVV;
    
    %vertical divergences, mol m^-3 s^-1
    adv_vConvNO3 = 0 .* ADVx_NO3;
    dif_vConvNO3 = 0 .* ADVx_NO3;
    
    for nz = 1:numLevels
        
        nzp1 = min([nz+1,numLevels]);
        
        adv_vConvNO3(:,:,nz) = squeeze(ADVr_NO3(:,:,nzp1) .* double(nz~=numLevels) - ADVr_NO3(:,:,nz));
        dif_vConvNO3(:,:,nz) = squeeze(DFrI_NO3(:,:,nzp1) .* double(nz~=numLevels) - DFrI_NO3(:,:,nz) ...
            + DFrE_NO3(:,:,nzp1) .* double(nz~=numLevels) - DFrE_NO3(:,:,nz));
        
    end
    
    adv_vConvNO3 = adv_vConvNO3 ./ VVV;
    dif_vConvNO3 = dif_vConvNO3 ./ VVV;
    
    %NO3 tendency due to E/P/runoff, mol m^-3 s^-1
    virtualFluxNO3 = 0 .* oceSPtnd;
    
    for nz = 1:numLevels
        
        if nz == 1
            
            virtualFluxNO3(:,:,1) =  NO3_Epr;
            
        else
            
            virtualFluxNO3(:,:,nz) =  0;
            
        end
        
    end
    
    virtualFluxNO3 = mygrid.mskC .* virtualFluxNO3;
    
    gDAR_NO3 = mygrid.mskC .* (gDAR_NO3 ./ mygrid.hFacC);
    
    %bio consumption
    bioC_NO3 = mygrid.mskC .* (-C_NO3 ./ mygrid.hFacC);
    
    %bio source
    bioS_NO3 = mygrid.mskC .* (S_NO3 ./ mygrid.hFacC);
    
    clear sStar
    
    %% 
    %NO2 budget
    
    D_snap = 0 .* NO2_SNAP;
    sStar = 0 .* NO2_SNAP;
    
    for nt = 1:2
        
        sStar(:,:,:,nt) = (1+mk3D(ETAN_SNAP(:,:,nt) ./ mygrid.Depth,dzMat));
        
        D_snap(:,:,:,nt) = sStar(:,:,:,nt) .* NO2_SNAP(:,:,:,nt);
        
    end
    
    %tendency, mol m^-3 s^-1
    tendNO2 = (D_snap(:,:,:,2) - D_snap(:,:,:,1)) ./ (secPerHour .* dt(timeStep));
    
    %horizontal divergences, mol m^-3 s^-1
    adv_hConvNO2 = calc_UV_conv(ADVx_NO2,ADVy_NO2) ./ VVV;
    dif_hConvNO2 = calc_UV_conv(DFxE_NO2,DFyE_NO2) ./ VVV;
    
    %vertical divergences, mol m^-3 s^-1
    adv_vConvNO2 = 0 .* ADVx_NO2;
    dif_vConvNO2 = 0 .* ADVx_NO2;
    
    for nz = 1:numLevels
        
        nzp1 = min([nz+1,numLevels]);
        
        adv_vConvNO2(:,:,nz) = squeeze(ADVr_NO2(:,:,nzp1) .* double(nz~=numLevels) - ADVr_NO2(:,:,nz));
        dif_vConvNO2(:,:,nz) = squeeze(DFrI_NO2(:,:,nzp1) .* double(nz~=numLevels) - DFrI_NO2(:,:,nz) ...
            + DFrE_NO2(:,:,nzp1) .* double(nz~=numLevels) - DFrE_NO2(:,:,nz));
        
    end
    
    adv_vConvNO2 = adv_vConvNO2 ./ VVV;
    dif_vConvNO2 = dif_vConvNO2 ./ VVV;
    
    %NO2 tendency due to E/P/runoff, mol m^-3 s^-1
    virtualFluxNO2 = 0 .* oceSPtnd;
    
    for nz = 1:numLevels
        
        if nz == 1
            
            virtualFluxNO2(:,:,1) =  NO2_Epr;
            
        else
            
            virtualFluxNO2(:,:,nz) =  0;
            
        end
        
    end
    
    virtualFluxNO2 = mygrid.mskC .* virtualFluxNO2;
    
    gDAR_NO2 = mygrid.mskC .* (gDAR_NO2 ./ mygrid.hFacC);
    
    %bio consumption
    bioC_NO2 = mygrid.mskC .* (-C_NO2 ./ mygrid.hFacC);
    
    %bio source
    bioS_NO2 = mygrid.mskC .* (S_NO2 ./ mygrid.hFacC);
       
    clear sStar
    
    %% 
    %NH4 budget
    
    D_snap = 0 .* NH4_SNAP;
    sStar = 0 .* NH4_SNAP;
    
    for nt = 1:2
        
        sStar(:,:,:,nt) = (1+mk3D(ETAN_SNAP(:,:,nt) ./ mygrid.Depth,dzMat));
        
        D_snap(:,:,:,nt) = sStar(:,:,:,nt) .* NH4_SNAP(:,:,:,nt);
        
    end
    
    %tendency, mol m^-3 s^-1
    tendNH4 = (D_snap(:,:,:,2) - D_snap(:,:,:,1)) ./ (secPerHour .* dt(timeStep));
    
    %horizontal divergences, mol m^-3 s^-1
    adv_hConvNH4 = calc_UV_conv(ADVx_NH4,ADVy_NH4) ./ VVV;
    dif_hConvNH4 = calc_UV_conv(DFxE_NH4,DFyE_NH4) ./ VVV;
    
    %vertical divergences, mol m^-3 s^-1
    adv_vConvNH4 = 0 .* ADVx_NH4;
    dif_vConvNH4 = 0 .* ADVx_NH4;
    
    for nz = 1:numLevels
        
        nzp1 = min([nz+1,numLevels]);
        
        adv_vConvNH4(:,:,nz) = squeeze(ADVr_NH4(:,:,nzp1) .* double(nz~=numLevels) - ADVr_NH4(:,:,nz));
        dif_vConvNH4(:,:,nz) = squeeze(DFrI_NH4(:,:,nzp1) .* double(nz~=numLevels) - DFrI_NH4(:,:,nz) ...
            + DFrE_NH4(:,:,nzp1) .* double(nz~=numLevels) - DFrE_NH4(:,:,nz));
        
    end
    
    adv_vConvNH4 = adv_vConvNH4 ./ VVV;
    dif_vConvNH4 = dif_vConvNH4 ./ VVV;
    
    %NH4 tendency due to E/P/runoff, mol m^-3 s^-1
    virtualFluxNH4 = 0 .* oceSPtnd;
    
    for nz = 1:numLevels
        
        if nz == 1
            
            virtualFluxNH4(:,:,1) =  NH4_Epr;
            
        else
            
            virtualFluxNH4(:,:,nz) =  0;
            
        end
        
    end
    
    virtualFluxNH4 = mygrid.mskC .* virtualFluxNH4;
    
    gDAR_NH4 = mygrid.mskC .* (gDAR_NH4 ./ mygrid.hFacC);
    
    %bio consumption
    bioC_NH4 = mygrid.mskC .* (-C_NH4 ./ mygrid.hFacC);
    
    %bio source
    bioS_NH4 = mygrid.mskC .* (S_NH4 ./ mygrid.hFacC);
    
    clear sStar
    
    %%
    %PO4 budget
    
    D_snap = 0 .* PO4_SNAP;
    sStar = 0 .* PO4_SNAP;
    
    for nt = 1:2
        
        sStar(:,:,:,nt) = (1+mk3D(ETAN_SNAP(:,:,nt) ./ mygrid.Depth,dzMat));
        
        D_snap(:,:,:,nt) = sStar(:,:,:,nt) .* PO4_SNAP(:,:,:,nt);
        
    end
    
    %tendency, mol m^-3 s^-1
    tendPO4 = (D_snap(:,:,:,2) - D_snap(:,:,:,1)) ./ (secPerHour .* dt(timeStep));
    
    %horizontal divergences, mol m^-3 s^-1
    adv_hConvPO4 = calc_UV_conv(ADVx_PO4,ADVy_PO4) ./ VVV;
    dif_hConvPO4 = calc_UV_conv(DFxE_PO4,DFyE_PO4) ./ VVV;
    
    %vertical divergences, mol m^-3 s^-1
    adv_vConvPO4 = 0 .* ADVx_PO4;
    dif_vConvPO4 = 0 .* ADVx_PO4;
    
    for nz = 1:numLevels
        
        nzp1 = min([nz+1,numLevels]);
        
        adv_vConvPO4(:,:,nz) = squeeze(ADVr_PO4(:,:,nzp1) .* double(nz~=numLevels) - ADVr_PO4(:,:,nz));
        dif_vConvPO4(:,:,nz) = squeeze(DFrI_PO4(:,:,nzp1) .* double(nz~=numLevels) - DFrI_PO4(:,:,nz) ...
            + DFrE_PO4(:,:,nzp1) .* double(nz~=numLevels) - DFrE_PO4(:,:,nz));
        
    end
    
    adv_vConvPO4 = adv_vConvPO4 ./ VVV;
    dif_vConvPO4 = dif_vConvPO4 ./ VVV;
    
    %PO4 tendency due to E/P/runoff, mol m^-3 s^-1
    virtualFluxPO4 = 0 .* oceSPtnd;
    
    for nz = 1:numLevels
        
        if nz == 1
            
            virtualFluxPO4(:,:,1) =  PO4_Epr;
            
        else
            
            virtualFluxPO4(:,:,nz) =  0;
            
        end
        
    end
    
    virtualFluxPO4 = mygrid.mskC .* virtualFluxPO4;
    
    %darwin tendency
    gDAR_PO4 = mygrid.mskC .* (gDAR_PO4 ./ mygrid.hFacC);
    
    %bio consumption
    bioC_PO4 = mygrid.mskC .* (-C_PO4 ./ mygrid.hFacC);
    
    %bio source
    bioS_PO4 = mygrid.mskC .* (S_PO4 ./ mygrid.hFacC);
    
    clear sStar
    
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
    
    %Fe tendency from iron dust, mol m^-3 s^-1
    forcFe = 0 .* oceSPtnd;
    
    %Fe tendency from E/P/runoff, mol m^-3 s^-1
    virtualFluxFe = 0 .* oceSPtnd;
    
    for nz = 1:numLevels
        
        if nz == 1
            
            forcFe(:,:,1) =  SFCSOLFe ./ dzMat(:,:,1);
            virtualFluxFe(:,:,1) =  Fe_Epr;
            
        else
            
            forcFe(:,:,nz) = 0;
            virtualFluxFe(:,:,nz) =  0;
            
        end
        
    end
    
    virtualFluxFe = mygrid.mskC .* virtualFluxFe;
    
    gDar_Fe = mygrid.mskC .* (gDAR_Fe ./ mygrid.hFacC);
    
    bioFe = ((gDAR_Fe./ mygrid.hFacC) - forcFe - SEDFe); %remove iron dust and sediment from gDAR
    
    forcFe = mygrid.mskC .* forcFe ./ mygrid.hFacC;
    sedFe = mygrid.mskC .* SEDFe ./ mygrid.hFacC;
    freeFe = mygrid.mskC .* (FREEFe  ./ mygrid.hFacC);
    
    clear sStar
    
    %%
    %SiO2 budget
    
    D_snap = 0 .* SiO2_SNAP;
    sStar = 0 .* SiO2_SNAP;
    
    for nt = 1:2
        
        sStar(:,:,:,nt) = (1+mk3D(ETAN_SNAP(:,:,nt) ./ mygrid.Depth,dzMat));
        
        D_snap(:,:,:,nt) = sStar(:,:,:,nt) .* SiO2_SNAP(:,:,:,nt);
        
    end
    
    %tendency, mol m^-3 s^-1
    tendSiO2 = (D_snap(:,:,:,2) - D_snap(:,:,:,1)) ./ (secPerHour .* dt(timeStep));
    
    %horizontal divergences, mol m^-3 s^-1
    adv_hConvSiO2 = calc_UV_conv(ADVx_SiO2,ADVy_SiO2) ./ VVV;
    dif_hConvSiO2 = calc_UV_conv(DFxE_SiO2,DFyE_SiO2) ./ VVV;
    
    %vertical divergences, mol m^-3 s^-1
    adv_vConvSiO2 = 0 .* ADVx_SiO2;
    dif_vConvSiO2 = 0 .* ADVx_SiO2;
    
    for nz = 1:numLevels
        
        nzp1 = min([nz+1,numLevels]);
        
        adv_vConvSiO2(:,:,nz) = squeeze(ADVr_SiO2(:,:,nzp1) .* double(nz~=numLevels) - ADVr_SiO2(:,:,nz));
        dif_vConvSiO2(:,:,nz) = squeeze(DFrI_SiO2(:,:,nzp1) .* double(nz~=numLevels) - DFrI_SiO2(:,:,nz) ...
            + DFrE_SiO2(:,:,nzp1) .* double(nz~=numLevels) - DFrE_SiO2(:,:,nz));
        
    end
    
    adv_vConvSiO2 = adv_vConvSiO2 ./ VVV;
    dif_vConvSiO2 = dif_vConvSiO2 ./ VVV;
    
    %SiO2 tendency due to E/P/runoff, mol m^-3 s^-1
    virtualFluxSiO2 = 0 .* oceSPtnd;
    
    for nz = 1:numLevels
        
        if nz == 1
            
            virtualFluxSiO2(:,:,1) =  SiO2_Epr;
            
        else
            
            virtualFluxSiO2(:,:,nz) =  0;
            
        end
        
    end
    
    virtualFluxSiO2 = mygrid.mskC .* virtualFluxSiO2;
    
    gDAR_SiO2 = mygrid.mskC .* (gDAR_SiO2 ./ mygrid.hFacC);
    
    %bio consumption
    bioC_SiO2 = mygrid.mskC .* (-C_SiO2 ./ mygrid.hFacC);
    
    %bio source
    bioS_SiO2 = mygrid.mskC .* (S_SiO2 ./ mygrid.hFacC);
    
    clear sStar
    
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
    
    tendDIC = convert2gcmfaces(tendDIC);
    adv_hConvDIC = convert2gcmfaces(adv_hConvDIC);
    adv_vConvDIC = convert2gcmfaces(adv_vConvDIC);
    dif_hConvDIC = convert2gcmfaces(dif_hConvDIC);
    dif_vConvDIC = convert2gcmfaces(dif_vConvDIC);
    forcDIC = convert2gcmfaces(forcDIC);
    virtualFluxDIC = convert2gcmfaces(virtualFluxDIC);
    gDARDIC = convert2gcmfaces(gDAR_DIC);
    bioDIC = convert2gcmfaces(bioDIC);
    
    tendNO3 = convert2gcmfaces(tendNO3);
    adv_hConvNO3 = convert2gcmfaces(adv_hConvNO3);
    adv_vConvNO3 = convert2gcmfaces(adv_vConvNO3);
    dif_hConvNO3 = convert2gcmfaces(dif_hConvNO3);
    dif_vConvNO3 = convert2gcmfaces(dif_vConvNO3);
    virtualFluxNO3 = convert2gcmfaces(virtualFluxNO3);
    gDARNO3 = convert2gcmfaces(gDAR_NO3);
    bioCNO3 = convert2gcmfaces(bioC_NO3);
    bioSNO3 = convert2gcmfaces(bioS_NO3);
    
    tendNO2 = convert2gcmfaces(tendNO2);
    adv_hConvNO2 = convert2gcmfaces(adv_hConvNO2);
    adv_vConvNO2 = convert2gcmfaces(adv_vConvNO2);
    dif_hConvNO2 = convert2gcmfaces(dif_hConvNO2);
    dif_vConvNO2 = convert2gcmfaces(dif_vConvNO2);
    virtualFluxNO2 = convert2gcmfaces(virtualFluxNO2);
    gDARNO2 = convert2gcmfaces(gDAR_NO2);
    bioCNO2 = convert2gcmfaces(bioC_NO2);
    bioSNO2 = convert2gcmfaces(bioS_NO2);
    
    tendNH4 = convert2gcmfaces(tendNH4);
    adv_hConvNH4 = convert2gcmfaces(adv_hConvNH4);
    adv_vConvNH4 = convert2gcmfaces(adv_vConvNH4);
    dif_hConvNH4 = convert2gcmfaces(dif_hConvNH4);
    dif_vConvNH4 = convert2gcmfaces(dif_vConvNH4);
    virtualFluxNH4 = convert2gcmfaces(virtualFluxNH4);
    gDARNH4 = convert2gcmfaces(gDAR_NH4);
    bioCNH4 = convert2gcmfaces(bioC_NH4);
    bioSNH4 = convert2gcmfaces(bioS_NH4);
    
    tendPO4 = convert2gcmfaces(tendPO4);
    adv_hConvPO4 = convert2gcmfaces(adv_hConvPO4);
    adv_vConvPO4 = convert2gcmfaces(adv_vConvPO4);
    dif_hConvPO4 = convert2gcmfaces(dif_hConvPO4);
    dif_vConvPO4 = convert2gcmfaces(dif_vConvPO4);
    virtualFluxPO4 = convert2gcmfaces(virtualFluxPO4);
    gDARPO4 = convert2gcmfaces(gDAR_PO4);
    bioCPO4 = convert2gcmfaces(bioC_PO4);
    bioSPO4 = convert2gcmfaces(bioS_PO4);
    
    tendFe = convert2gcmfaces(tendFe);
    adv_hConvFe = convert2gcmfaces(adv_hConvFe);
    adv_vConvFe = convert2gcmfaces(adv_vConvFe);
    dif_hConvFe = convert2gcmfaces(dif_hConvFe);
    dif_vConvFe = convert2gcmfaces(dif_vConvFe);
    forcFe = convert2gcmfaces(forcFe);
    virtualFluxFe = convert2gcmfaces(virtualFluxFe);
    gDARFe = convert2gcmfaces(gDAR_Fe);
    bioFe = convert2gcmfaces(bioFe);
    sedFe = convert2gcmfaces(sedFe);
    freeFe = convert2gcmfaces(freeFe);
    
    tendSiO2 = convert2gcmfaces(tendSiO2);
    adv_hConvSiO2 = convert2gcmfaces(adv_hConvSiO2);
    adv_vConvSiO2 = convert2gcmfaces(adv_vConvSiO2);
    dif_hConvSiO2 = convert2gcmfaces(dif_hConvSiO2);
    dif_vConvSiO2 = convert2gcmfaces(dif_vConvSiO2);
    virtualFluxSiO2 = convert2gcmfaces(virtualFluxSiO2);
    gDARSiO2 = convert2gcmfaces(gDAR_SiO2);
    bioCSiO2 = convert2gcmfaces(bioC_SiO2);
    bioSSiO2 = convert2gcmfaces(bioS_SiO2);
    
    %%
    
    intTendV = tendV(:,:,intLevel);
    intHConvV = adv_hConvV(:,:,intLevel);
    intVConvV = adv_vConvV(:,:,intLevel);
    intForcV = forcV(:,:,intLevel);
    
    intTendS = tendS(:,:,intLevel);
    intHAdvS = adv_hConvS(:,:,intLevel);
    intVAdvS = adv_vConvS(:,:,intLevel);
    intHDifS = dif_hConvS(:,:,intLevel);
    intVDifS = dif_vConvS(:,:,intLevel);
    intForcS = forcS(:,:,intLevel);
    
    intTendSal = tendSal(:,:,intLevel);
    intHAdvSal = adv_hConvSal(:,:,intLevel);
    intVAdvSal = adv_vConvSal(:,:,intLevel);
    intHDifSal = dif_hConvSal(:,:,intLevel);
    intVDifSal = dif_vConvSal(:,:,intLevel);
    intForcSal = forcSal(:,:,intLevel);
    
    intTendDIC = tendDIC(:,:,intLevel);
    intHAdvDIC = adv_hConvDIC(:,:,intLevel);
    intVAdvDIC = adv_vConvDIC(:,:,intLevel);
    intHDifDIC = dif_hConvDIC(:,:,intLevel);
    intVDifDIC = dif_vConvDIC(:,:,intLevel);
    intForcDIC = forcDIC(:,:,intLevel);
    intVirtualFluxDIC = virtualFluxDIC(:,:,intLevel);
    intGDARDIC = gDARDIC(:,:,intLevel);
    intBioDIC = bioDIC(:,:,intLevel);
    
    intDICResidual = intTendDIC - (intHAdvDIC + intVAdvDIC ...
        + intHDifDIC + intVDifDIC + intForcDIC + intVirtualFluxDIC.*0 + intBioDIC);

    intTendNO3 = tendNO3(:,:,intLevel);
    intHAdvNO3 = adv_hConvNO3(:,:,intLevel);
    intVAdvNO3 = adv_vConvNO3(:,:,intLevel);
    intHDifNO3 = dif_hConvNO3(:,:,intLevel);
    intVDifNO3 = dif_vConvNO3(:,:,intLevel);
    intVirtualFluxNO3 = virtualFluxNO3(:,:,intLevel);
    intGDARNO3 = gDARNO3(:,:,intLevel);
    intBioCNO3 = bioCNO3(:,:,intLevel);
    intBioSNO3 = bioSNO3(:,:,intLevel);
    
    intNO3Residual = intTendNO3 - (intHAdvNO3 + intVAdvNO3 ...
        + intHDifNO3 + intVDifNO3 + intVirtualFluxNO3.*0 ...
        + intBioCNO3 + intBioSNO3 + intGDARNO3.*0);
    
    intTendNO2 = tendNO2(:,:,intLevel);
    intHAdvNO2 = adv_hConvNO2(:,:,intLevel);
    intVAdvNO2 = adv_vConvNO2(:,:,intLevel);
    intHDifNO2 = dif_hConvNO2(:,:,intLevel);
    intVDifNO2 = dif_vConvNO2(:,:,intLevel);
    intVirtualFluxNO2 = virtualFluxNO2(:,:,intLevel);
    intGDARNO2 = gDARNO2(:,:,intLevel);
    intBioCNO2 = bioCNO2(:,:,intLevel);
    intBioSNO2 = bioSNO2(:,:,intLevel);
    
    intNO2Residual = intTendNO2 - (intHAdvNO2 + intVAdvNO2 ...
        + intHDifNO2 + intVDifNO2 + intVirtualFluxNO2.*0 ...
        + intBioCNO2.*0 + intBioSNO2.*0 + intGDARNO2);
    
    intTendNH4 = tendNH4(:,:,intLevel);
    intHAdvNH4 = adv_hConvNH4(:,:,intLevel);
    intVAdvNH4 = adv_vConvNH4(:,:,intLevel);
    intHDifNH4 = dif_hConvNH4(:,:,intLevel);
    intVDifNH4 = dif_vConvNH4(:,:,intLevel);
    intVirtualFluxNH4 = virtualFluxNH4(:,:,intLevel);
    intGDARNH4 = gDARNH4(:,:,intLevel);
    intBioCNH4 = bioCNH4(:,:,intLevel);
    intBioSNH4 = bioSNH4(:,:,intLevel);
    
    intNH4Residual = intTendNH4 - (intHAdvNH4 + intVAdvNH4 ...
        + intHDifNH4 + intVDifNH4 + intVirtualFluxNH4.*0 ...
        + intBioCNH4 + intBioSNH4 + intGDARNH4.*0);
    
    intTendPO4 = tendPO4(:,:,intLevel);
    intHAdvPO4 = adv_hConvPO4(:,:,intLevel);
    intVAdvPO4 = adv_vConvPO4(:,:,intLevel);
    intHDifPO4 = dif_hConvPO4(:,:,intLevel);
    intVDifPO4 = dif_vConvPO4(:,:,intLevel);
    intVirtualFluxPO4 = virtualFluxPO4(:,:,intLevel);
    intGDARPO4 = gDARPO4(:,:,intLevel);
    intBioCPO4 = bioCPO4(:,:,intLevel);
    intBioSPO4 = bioSPO4(:,:,intLevel);
    
    intPO4Residual = intTendPO4 - (intHAdvPO4 + intVAdvPO4 ...
        + intHDifPO4 + intVDifPO4 + intVirtualFluxPO4.*0 ...
        + intBioCPO4 + intBioSPO4);
    
    intTendFe = tendFe(:,:,intLevel);
    intHAdvFe = adv_hConvFe(:,:,intLevel);
    intVAdvFe = adv_vConvFe(:,:,intLevel);
    intHDifFe = dif_hConvFe(:,:,intLevel);
    intVDifFe = dif_vConvFe(:,:,intLevel);
    intForcFe = forcFe(:,:,intLevel);
    intVirtualFluxFe = virtualFluxFe(:,:,intLevel);
    intGDARFe = gDARFe(:,:,intLevel);
    
    intBioFe = bioFe(:,:,intLevel);
    intSedFe = sedFe(:,:,intLevel);
    intFreeFe = freeFe(:,:,intLevel);
    
    intFeResidual = intTendFe - (intHAdvFe + intVAdvFe + intHDifFe + intVDifFe ...
        + intGDARFe + intFreeFe.*0);
    
    intTendSiO2 = tendSiO2(:,:,intLevel);
    intHAdvSiO2 = adv_hConvSiO2(:,:,intLevel);
    intVAdvSiO2 = adv_vConvSiO2(:,:,intLevel);
    intHDifSiO2 = dif_hConvSiO2(:,:,intLevel);
    intVDifSiO2 = dif_vConvSiO2(:,:,intLevel);
    intVirtualFluxSiO2 = virtualFluxSiO2(:,:,intLevel);
    intGDARSiO2 = gDARSiO2(:,:,intLevel);
    intBioCSiO2 = bioCSiO2(:,:,intLevel);
    intBioSSiO2 = bioSSiO2(:,:,intLevel);
    
    intSiO2Residual = intTendSiO2 - (intHAdvSiO2 + intVAdvSiO2 ...
        + intHDifSiO2 + intVDifSiO2 + intVirtualFluxSiO2.*0 ...
        + intBioCSiO2 + intBioSSiO2);
    
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
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intTendV);
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
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intHConvV);
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
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intVConvV);
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
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intForcV);
        shading flat
        
        caxis([cMin .* 50 cMax .* 50]);
        
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
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intHConvV + intVConvV + intForcV);
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
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intTendV - (intForcV + intHConvV + intVConvV));
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
        
        if savePlot
            
            print('-dpng',[figureDir 'volume_budget' num2str(timeStep) '.png']);
            
        end
        
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
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intTendS);
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
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intHAdvS);
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
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intVAdvS);
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
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intHDifS);
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
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intVDifS);
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
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intForcS);
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
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intHAdvS + intVAdvS + intHDifS + intVDifS + intForcS);
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
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intTendS - (intHAdvS + intVAdvS + intHDifS + intVDifS + intForcS));
        shading flat
        
        caxis([cMin .* 10^-7 cMax .* 10^-7]);
        
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
        
        title({'Budget Total - Tendency';'(psu m s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(289);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intTendSal);
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
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intHAdvSal);
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
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intVAdvSal);
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
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intHDifSal);
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
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intVDifSal);
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
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intForcSal);
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
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,(intHAdvSal + intVAdvSal + intHDifSal + intVDifSal + intForcSal));
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
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intTendSal - (intHAdvSal + intVAdvSal + intHDifSal + intVDifSal + intForcSal));
        shading flat
        
        caxis([cMin .* 10^-7 cMax .* 10^-7]);
        
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
        
        title({'Budget Total - Tendency';'(psu m s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        if savePlot
            
            print('-dpng',[figureDir 'salinity_budget' num2str(timeStep) '.png']);
            
        end
        
        drawnow
        
    end
    
    if plotDICBudget
        
        hFig3 = figure(3);
        set(hFig3,'units','normalized','outerposition',[0 0 1 1]);
        set(gcf,'color',[1 1 1]);
        set(gca,'color',[0.5 0.5 0.5]);
        
        cMin = -10^-8;
        cMax = 10^-8;
        
        subplot(2,5,1);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intTendDIC);
        shading flat
        
        caxis([cMin*10^2 cMax*10^2]);
        
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
        
        title({'DIC Tendency (mol m^-^3 s^-^1), ';'timeStep: ';num2str(tt(timeStep))},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,2);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intHAdvDIC);
        shading flat
        
        caxis([cMin*10^2 cMax*10^2]);
        
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
        
        title({'Horizontal Advection';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,3);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intVAdvDIC);
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
        
        title({'Vertical Advection';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,4);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intHDifDIC);
        shading flat
        
        caxis([cMin*10^2 cMax*10^2]);
        
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
        
        title({'Horizontal Diffusion';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,5);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intVDifDIC);
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
        
        title({'Vertical Diffusion';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,6);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intForcDIC);
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
        
        title({'Air-sea CO_2 Flux';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,7);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intVirtualFluxDIC);
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
        
        title({'P/E/runoff';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,8);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intBioDIC);
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
        
        title({'Biology';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,9);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intHAdvDIC + intVAdvDIC + intHDifDIC + intVDifDIC ...
            + intForcDIC + intVirtualFluxDIC + intBioDIC);
        
        shading flat
        
        caxis([cMin*10^2 cMax*10^2]);
        
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
        
        title({'Total Budget';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,10);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intDICResidual);
        
        shading flat
        %caxis([cMin cMax]);
        caxis([cMin*10^-5 cMax*10^-5]);
        
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
        
        title({'Budget Total - Tendency';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        drawnow
        
        if savePlot
            
            print('-dpng',[figureDir 'DIC_budget' num2str(timeStep) '.png']);
            
        end
        
    end
    
    if plotNO3Budget
        
        hFig4 = figure(4);
        set(hFig4,'units','normalized','outerposition',[0 0 1 1]);
        set(gcf,'color',[1 1 1]);
        set(gca,'color',[0.5 0.5 0.5]);
        
        cMin = -10^-9;
        cMax = 10^-9;
        
        subplot(2,5,1);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intTendNO3);
        shading flat
        
        caxis([cMin*10^2 cMax*10^2]);
        
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
        
        title({'NO3 Tendency (mol m^-^3 s^-^1), ';'timeStep: ';num2str(tt(timeStep))},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,2);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intHAdvNO3);
        shading flat
        
        caxis([cMin*10^2 cMax*10^2]);
        
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
        
        title({'Horizontal Advection';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,3);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intVAdvNO3);
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
        
        title({'Vertical Advection';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,4);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intHDifNO3);
        shading flat
        
        caxis([cMin*10^2 cMax*10^2]);
        
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
        
        title({'Horizontal Diffusion';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,5);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intVDifNO3);
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
        
        title({'Vertical Diffusion';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,6);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intBioCNO3);
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
        
        title({'Bio C';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,7);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intBioSNO3);
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
        
        title({'Bio S';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,8);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intGDARNO3);
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
        
        title({'gDAR NO3';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,9);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intHAdvNO3 + intVAdvNO3 + intHDifNO3 + intVDifNO3 ...
            + intGDARNO3);
        
        shading flat
        
        caxis([cMin*10^2 cMax*10^2]);
        
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
        
        title({'Total Budget';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,10);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intNO3Residual);
        
        shading flat
        %caxis([cMin*10^-2 cMax*10^-2]);
        caxis([cMin*10^-5 cMax*10^-5]);
        
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
        
        title({'Budget Total - Tendency';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        drawnow
        
        if savePlot
            
            print('-dpng',[figureDir 'NO3_budget' num2str(timeStep) '.png']);
            
        end
        
    end
    
    if plotNO2Budget
        
        hFig5 = figure(5);
        set(hFig5,'units','normalized','outerposition',[0 0 1 1]);
        set(gcf,'color',[1 1 1]);
        set(gca,'color',[0.5 0.5 0.5]);
        
        cMin = -10^-9;
        cMax = 10^-9;
        
        subplot(2,5,1);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intTendNO2);
        shading flat
        
        caxis([cMin*10^2 cMax*10^2]);
        
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
        
        title({'NO2 Tendency (mol m^-^3 s^-^1), ';'timeStep: ';num2str(tt(timeStep))},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,2);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intHAdvNO2);
        shading flat
        
        caxis([cMin*10^2 cMax*10^2]);
        
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
        
        title({'Horizontal Advection';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,3);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intVAdvNO2);
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
        
        title({'Vertical Advection';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,4);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intHDifNO2);
        shading flat
        
        caxis([cMin*10^2 cMax*10^2]);
        
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
        
        title({'Horizontal Diffusion';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,5);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intVDifNO2);
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
        
        title({'Vertical Diffusion';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,6);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intBioCNO2);
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
        
        title({'Bio C';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,7);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intBioSNO2);
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
        
        title({'Bio S';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,8);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intGDARNO2);
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
        
        title({'gDAR NO2';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,9);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intHAdvNO2 + intVAdvNO2 + intHDifNO2 + intVDifNO2 ...
            + intGDARNO2);
        
        shading flat
        
        caxis([cMin*10^2 cMax*10^2]);
        
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
        
        title({'Total Budget';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,10);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intNO2Residual);
        
        shading flat
        %caxis([cMin*10^-2 cMax*10^-2]);
        caxis([cMin*10^-5 cMax*10^-5]);
        
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
        
        title({'Budget Total - Tendency';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        drawnow
        
        if savePlot
            
            print('-dpng',[figureDir 'NO2_budget' num2str(timeStep) '.png']);
            
        end
        
    end
    
    if plotNH4Budget
        
        hFig6 = figure(6);
        set(hFig6,'units','normalized','outerposition',[0 0 1 1]);
        set(gcf,'color',[1 1 1]);
        set(gca,'color',[0.5 0.5 0.5]);
        
        cMin = -10^-9;
        cMax = 10^-9;
        
        subplot(2,5,1);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intTendNH4);
        shading flat
        
        caxis([cMin*10^2 cMax*10^2]);
        
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
        
        title({'NH4 Tendency (mol m^-^3 s^-^1), ';'timeStep: ';num2str(tt(timeStep))},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,2);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intHAdvNH4);
        shading flat
        
        caxis([cMin*10^2 cMax*10^2]);
        
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
        
        title({'Horizontal Advection';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,3);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intVAdvNH4);
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
        
        title({'Vertical Advection';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,4);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intHDifNH4);
        shading flat
        
        caxis([cMin*10^2 cMax*10^2]);
        
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
        
        title({'Horizontal Diffusion';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,5);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intVDifNH4);
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
        
        title({'Vertical Diffusion';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,6);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intBioCNH4);
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
        
        title({'Bio C';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,7);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intBioSNH4);
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
        
        title({'Bio S';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,8);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intGDARNH4);
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
        
        title({'gDAR NH4';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,9);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intHAdvNH4 + intVAdvNH4 + intHDifNH4 + intVDifNH4 ...
            + intGDARNH4);
        
        shading flat
        
        caxis([cMin*10^2 cMax*10^2]);
        
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
        
        title({'Total Budget';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,10);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intNH4Residual);
        
        shading flat
        %caxis([cMin*10^-2 cMax*10^-2]);
        caxis([cMin*10^-5 cMax*10^-5]);
        
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
        
        title({'Budget Total - Tendency';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        drawnow
        
        if savePlot
            
            print('-dpng',[figureDir 'NH4_budget' num2str(timeStep) '.png']);
            
        end
        
    end
    
    if plotPO4Budget
        
        hFig7 = figure(7);
        set(hFig7,'units','normalized','outerposition',[0 0 1 1]);
        set(gcf,'color',[1 1 1]);
        set(gca,'color',[0.5 0.5 0.5]);
        
        cMin = -10^-10;
        cMax = 10^-10;
        
        subplot(2,5,1);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intTendPO4);
        shading flat
        
        caxis([cMin*10^2 cMax*10^2]);
        
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
        
        title({'PO4 Tendency (mol m^-^3 s^-^1), ';'timeStep: ';num2str(tt(timeStep))},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,2);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intHAdvPO4);
        shading flat
        
        caxis([cMin*10^2 cMax*10^2]);
        
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
        
        title({'Horizontal Advection';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,3);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intVAdvPO4);
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
        
        title({'Vertical Advection';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,4);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intHDifPO4);
        shading flat
        
        caxis([cMin*10^2 cMax*10^2]);
        
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
        
        title({'Horizontal Diffusion';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,5);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intVDifPO4);
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
        
        title({'Vertical Diffusion';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,6);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intBioCPO4);
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
        
        title({'Bio C';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,7);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intBioSPO4);
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
        
        title({'Bio S';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,8);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intGDARPO4);
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
        
        title({'gDAR';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,9);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intHAdvPO4 + intVAdvPO4 + intHDifPO4 + intVDifPO4 ...
            + intGDARPO4);
        
        shading flat
        
        caxis([cMin*10^2 cMax*10^2]);
        
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
        
        title({'Total Budget';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,10);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intPO4Residual);
        
        shading flat
        %caxis([cMin*10^-2 cMax*10^-2]);
        caxis([cMin*10^-5 cMax*10^-5]);
        
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
        
        title({'Budget Total - Tendency';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        drawnow
        
        if savePlot
            
            print('-dpng',[figureDir 'PO4_budget' num2str(timeStep) '.png']);
            
        end
        
    end
    
    if plotFeBudget
        
        hFig8 = figure(8);
        set(hFig8,'units','normalized','outerposition',[0 0 1 1]);
        set(gcf,'color',[1 1 1]);
        set(gca,'color',[0.5 0.5 0.5]);
        
        cMin = -10^-14;
        cMax = 10^-14;
        
        subplot(2,5,1);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intTendFe);
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
        
        title({'dFe/dt (mol m^-^3 s^-^1), ';'timeStep: ';num2str(tt(timeStep))},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,2);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intHAdvFe);
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
        
        title({'Horizontal Advection';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,3);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intVAdvFe);
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
        
        title({'Vertical Advection';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,4);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intHDifFe);
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
        
        title({'Horizontal Diffusion';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,5);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intVDifFe);
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
        
        title({'Vertical Diffusion';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,6);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intForcFe);
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
        
        title({'Surface Dust Flux';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,7);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intSedFe);
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
        
        title({'Sediment Flux';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,8);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intBioFe);
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
        
        title({'Scav + Bio';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,9);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intFreeFe);
        
        shading flat
        caxis([cMin*10^-5 cMax*10^-5]);
        
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
        
        title({'Free Iron';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,10);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intFeResidual);
        
        shading flat
        
        %caxis([cMin.*10^-2 cMax.*10^-2]);
        caxis([cMin.*10^-5 cMax.*10^-5]);
        
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
        
        title({'Budget Residual';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        drawnow
        
    end
    
    if plotSiO2Budget
        
        hFig9 = figure(9);
        set(hFig9,'units','normalized','outerposition',[0 0 1 1]);
        set(gcf,'color',[1 1 1]);
        set(gca,'color',[0.5 0.5 0.5]);
        
        cMin = -10^-9;
        cMax = 10^-9;
        
        subplot(2,5,1);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intTendSiO2);
        shading flat
        
        caxis([cMin*10^2 cMax*10^2]);
        
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
        
        title({'SiO2 Tendency (mol m^-^3 s^-^1), ';'timeStep: ';num2str(tt(timeStep))},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,2);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intHAdvSiO2);
        shading flat
        
        caxis([cMin*10^2 cMax*10^2]);
        
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
        
        title({'Horizontal Advection';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,3);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intVAdvSiO2);
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
        
        title({'Vertical Advection';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,4);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intHDifSiO2);
        shading flat
        
        caxis([cMin*10^2 cMax*10^2]);
        
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
        
        title({'Horizontal Diffusion';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,5);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intVDifSiO2);
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
        
        title({'Vertical Diffusion';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,6);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intBioCSiO2);
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
        
        title({'Bio C';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,7);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intBioSSiO2);
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
        
        title({'Bio S';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,8);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intGDARSiO2);
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
        
        title({'gDAR SiO2';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,9);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intHAdvSiO2 + intVAdvSiO2 + intHDifSiO2 + intVDifSiO2 ...
            + intGDARSiO2);
        
        shading flat
        
        caxis([cMin*10^2 cMax*10^2]);
        
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
        
        title({'Total Budget';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,10);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,intSiO2Residual);
        
        shading flat
        %caxis([cMin*10^-1 cMax*10^-1]);
        caxis([cMin*10^-5 cMax*10^-5]);
        
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
        
        title({'Budget Total - Tendency';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        drawnow
        
        if savePlot
            
            print('-dpng',[figureDir 'SiO2_budget' num2str(timeStep) '.png']);
            
        end
        
    end
    
    %clear gcmfaces objects to avoid memory leaks
    
    test.DIC = nansum(nansum(intDICResidual));
    test.NO3 = nansum(nansum(intNO3Residual));
    test.NO2 = nansum(nansum(intNO2Residual));
    test.NH4 = nansum(nansum(intNH4Residual));
    test.PO4 = nansum(nansum(intPO4Residual));
    test.Fe = nansum(nansum(intFeResidual));
    test.SiO2 = nansum(nansum(intSiO2Residual));

    clear ETAN oceFWflx SFLUX  oceSPflx UVELMASS VVELMASS WVELMASS ...
        SALT ADVr_SLT ADVx_SLT ADVy_SLT DFrI_SLT DFrE_SLT DFxE_SLT DFyE_SLT oceSPtnd ...
        DIC ADVx_DIC ADVy_DIC ADVr_DIC DFxE_DIC DFyE_DIC DFrE_DIC DFrI_DIC DICTFLX ...
        CONSUMP_DIC CONSUMP_DIC_PIC REMIN_DOC REMIN_POC DISSC_PIC ...
        ETAN_SNAP SALT_SNAP DIC_SNAP NO3_SNAP NO2_SNAP NH4_SNAP PO4_SNAP Fe_SNAP SiO2_SNAP ...
        tendV adv_hConvV adv_vConvV forcV ...
        S_snap tendS adv_hConvS dif_hConvS adv_vConvS dif_vConvS forcS surfS ...
        rstarfac tendSal adv_hConvSal adv_vConvSal dif_hConvSal dif_vConvSal forcSal ...
        
    clear tendDIC adv_hConvDIC dif_hConvDIC adv_vConvDIC dif_vConvDIC forcDIC gDARDIC bioDIC virtualFluxDIC ...
        
    clear intTendDIC intHAdvDIC intVAdvDIC intHDifDIC intVDifDIC intForcDIC intVirtualFluxDIC intGDARDIC intBioDIC intDICResidual

    %pause

%close all

end
