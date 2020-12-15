clear
close all

%%
%settings, modify as needed

%addpath(genpath('/nobackup/dcarrol2/MATLAB'));

savePlot = 1;
saveMat = 1; %save budget .mat files 

useVol = 0; %integrate w/ volume
useLLC270 = 1; %use LLC 270 grid

startIntLevel = 1; %vertical integration start k level
endIntLevel = 15;

%set to 1 to plot budget terms
plotBudgetVolume = 1;
plotBudgetSalinity = 1;
plotBudgetDIC = 1;
plotBudgetNO3 = 1;
plotBudgetNO2 = 1;
plotBudgetNH4 = 1;
plotBudgetPO4 = 1;
plotBudgetFe = 1;
plotBudgetSiO2 = 1;

nanString = 'omitnan';

%%

gridDir = '/nobackup/dcarrol2/grid/LLC_270/';

modelDir = '/nobackup/dcarrol2/v05_DIC_budget/darwin3/run/';
figureDir = 'figures/';
saveDir = 'mat/';

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
    
    flipString = '';
    
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
    gDAR_DIC = convert2gcmfaces(rdmds([diagDir filename4],ttAverage,'rec',9)) .* mmol_to_mol;  %mol m^-3 s^-1
    
    %bio decomposition
    cDIC = convert2gcmfaces(rdmds([diagDir filename4],ttAverage,'rec',10)) .* mmol_to_mol;  %mol m^-3 s^-1
    cDIC_PIC = convert2gcmfaces(rdmds([diagDir filename4],ttAverage,'rec',11)) .* mmol_to_mol;  %mol m^-3 s^-1
    respDIC = convert2gcmfaces(rdmds([diagDir filename4],ttAverage,'rec',12)) .* mmol_to_mol;  %mol m^-3 s^-1
    rDIC_DOC = convert2gcmfaces(rdmds([diagDir filename4],ttAverage,'rec',13)) .* mmol_to_mol;  %mol m^-3 s^-1
    rDIC_POC = convert2gcmfaces(rdmds([diagDir filename4],ttAverage,'rec',14)) .* mmol_to_mol;  %mol m^-3 s^-1
    dDIC_PIC = convert2gcmfaces(rdmds([diagDir filename4],ttAverage,'rec',15)) .* mmol_to_mol;  %mol m^-3 s^-1
    
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
    forcDIC = 0 .* oceSPtnd;
    
    %DIC tendency due to E/P/runoff, mol m^-3 s^-1
    virtualFluxDIC = 0 .* oceSPtnd;
    
    for nz = 1:numLevels
        
        if nz == 1
            
            forcDIC(:,:,1) = fluxCO2;
            virtualFluxDIC(:,:,1) =  DIC_Epr;
            
        else
            
            forcDIC(:,:,nz) = 0;
            virtualFluxDIC(:,:,nz) =  0;
            
        end
        
    end
    
    gDAR_DIC = (gDAR_DIC - (forcDIC ./ dzMat)) ./ mygrid.hFacC; %remove air-sea CO2 flux from gDAR, so it is just biology
    
    forcDIC = mygrid.mskC .* (forcDIC ./ dzMat);
    
    virtualFluxDIC = mygrid.mskC .* virtualFluxDIC .* 0;
    
    %net biology, mol m^-3 s^-1
    bioDIC = mygrid.mskC .* gDAR_DIC;
    
    %individual biology terms
    bioCons_DIC = mygrid.mskC .* (-cDIC ./ mygrid.hFacC);
    bioCons_DIC_PIC = mygrid.mskC .* (-cDIC_PIC ./ mygrid.hFacC);
    bioResp_DIC = mygrid.mskC .* (respDIC ./ mygrid.hFacC);
    bioRemin_DIC_DOC = mygrid.mskC .* (rDIC_DOC ./ mygrid.hFacC);
    bioRemin_DIC_POC = mygrid.mskC .* (rDIC_POC ./ mygrid.hFacC);
    bioDissc_DIC_PIC = mygrid.mskC .* (dDIC_PIC ./ mygrid.hFacC);
    
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
    
    forcFe = mygrid.mskC .* (forcFe ./ mygrid.hFacC);
    sedFe = mygrid.mskC .* (SEDFe ./ mygrid.hFacC);
    freeFe = mygrid.mskC .* (FREEFe  ./ mygrid.hFacC);
    
    bioFe = gDAR_Fe - forcFe - SEDFe; %remove iron dust and sediment flu from gDAR
    
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
    
    bioConsDIC = convert2gcmfaces(bioCons_DIC);
    bioConsDIC_PIC = convert2gcmfaces(bioCons_DIC_PIC);
    bioRespDIC = convert2gcmfaces(bioResp_DIC);
    bioReminDIC_DOC = convert2gcmfaces(bioRemin_DIC_DOC);
    bioReminDIC_POC = convert2gcmfaces(bioRemin_DIC_POC);
    bioDisscDIC_PIC = convert2gcmfaces(bioDissc_DIC_PIC);
    
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
    
    %NO3
    intTendNO3 = sum((tendNO3(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intHAdvNO3 = sum((adv_hConvNO3(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intVAdvNO3 = sum((adv_vConvNO3(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intHDifNO3 = sum((dif_hConvNO3(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intVDifNO3 = sum((dif_vConvNO3(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intVirtualFluxNO3 = sum((virtualFluxNO3(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intGDARNO3 = sum((gDARNO3(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intBioCNO3 = sum((bioCNO3(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intBioSNO3 = sum((bioSNO3(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    
    intTotalNO3 = intHAdvNO3 + intVAdvNO3 ...
        + intHDifNO3 + intVDifNO3 + intVirtualFluxNO3.*0 ...
        + intBioCNO3 + intBioSNO3;
    
    intResidualNO3 = intTendNO3 - intTotalNO3;
    
    %NO2
    intTendNO2 = sum((tendNO2(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intHAdvNO2 = sum((adv_hConvNO2(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intVAdvNO2 = sum((adv_vConvNO2(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intHDifNO2 = sum((dif_hConvNO2(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intVDifNO2 = sum((dif_vConvNO2(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intVirtualFluxNO2 = sum((virtualFluxNO2(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intGDARNO2 = sum((gDARNO2(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intBioCNO2 = sum((bioCNO2(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intBioSNO2 = sum((bioSNO2(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    
    intTotalNO2 = intHAdvNO2 + intVAdvNO2 ...
        + intHDifNO2 + intVDifNO2 + intVirtualFluxNO2.*0 ...
        + intBioCNO2 + intBioSNO2;
    
    intResidualNO2 = intTendNO2 - intTotalNO2;
    
    %NH4
    intTendNH4 = sum((tendNH4(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intHAdvNH4 = sum((adv_hConvNH4(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intVAdvNH4 = sum((adv_vConvNH4(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intHDifNH4 = sum((dif_hConvNH4(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intVDifNH4 = sum((dif_vConvNH4(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intVirtualFluxNH4 = sum((virtualFluxNH4(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intGDARNH4 = sum((gDARNH4(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intBioCNH4 = sum((bioCNH4(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intBioSNH4 = sum((bioSNH4(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    
    intTotalNH4 = intHAdvNH4 + intVAdvNH4 ...
        + intHDifNH4 + intVDifNH4 + intVirtualFluxNH4.*0 ...
        + intBioCNH4 + intBioSNH4;
    
    intResidualNH4 = intTendNH4 - intTotalNH4;
    
    %PO4
    intTendPO4 = sum((tendPO4(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intHAdvPO4 = sum((adv_hConvPO4(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intVAdvPO4 = sum((adv_vConvPO4(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intHDifPO4 = sum((dif_hConvPO4(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intVDifPO4 = sum((dif_vConvPO4(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intVirtualFluxPO4 = sum((virtualFluxPO4(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intGDARPO4 = sum((gDARPO4(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intBioCPO4 = sum((bioCPO4(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intBioSPO4 = sum((bioSPO4(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    
    intTotalPO4 = intHAdvPO4 + intVAdvPO4 ...
        + intHDifPO4 + intVDifPO4 + intVirtualFluxPO4.*0 ...
        + intBioCPO4 + intBioSPO4;
    
    intResidualPO4 = intTendPO4 - intTotalPO4;
    
    %Fe
    intTendFe = sum((tendFe(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intHAdvFe = sum((adv_hConvFe(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intVAdvFe = sum((adv_vConvFe(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intHDifFe = sum((dif_hConvFe(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intVDifFe = sum((dif_vConvFe(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intForcFe = sum((forcFe(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intVirtualFluxFe = sum((virtualFluxFe(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intGDARFe = sum((gDARFe(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intBioFe = sum((bioFe(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intSedFe = sum((sedFe(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intFreeFe = sum((freeFe(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    
    intTotalFe = intHAdvFe + intVAdvFe + intHDifFe + intVDifFe ...
        + intGDARFe + intFreeFe;
    
    intResidualFe = intTendFe - intTotalFe;
    
    %SiO2
    intTendSiO2 = sum((tendSiO2(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intHAdvSiO2 = sum((adv_hConvSiO2(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intVAdvSiO2 = sum((adv_vConvSiO2(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intHDifSiO2 = sum((dif_hConvSiO2(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intVDifSiO2 = sum((dif_vConvSiO2(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intVirtualFluxSiO2 = sum((virtualFluxSiO2(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intGDARSiO2 = sum((gDARSiO2(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intBioCSiO2 = sum((bioCSiO2(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    intBioSSiO2 = sum((bioSSiO2(:,:,startIntLevel:endIntLevel) .* vol(:,:,startIntLevel:endIntLevel)),3,nanString);
    
    intTotalSiO2 = intHAdvSiO2 + intVAdvSiO2 ...
        + intHDifSiO2 + intVDifSiO2 + intVirtualFluxSiO2.*0 ...
        + intBioCSiO2 + intBioSSiO2;
    
    intResidualSiO2 = intTendSiO2 - intTotalSiO2;
    
    %%
    %plot budgets
    
    if plotBudgetVolume
        
        hFig1 = figure(1);
        set(hFig1,'units','normalized','outerposition',[0 0 1 0.5]);
        set(gcf,'color',[1 1 1]);
        
        if useVol
            
            cMin = -10^3;
            cMax = 10^3;
            
        else
            
            cMin = -10^-10;
            cMax = 10^-10;
            
        end
        
        subplot(161);
        
        eval([plotString '(intTendV' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Volume Tendency, ';'Timestep: ';num2str(tt(timeStep))},'FontWeight','Bold','FontSize',fs);
        
        subplot(162);
        
        eval([plotString '(intHConvV' flipString ')']);
        
        shading flat
        caxis([cMin.*10^3 cMax.*10^3]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Horizontal';'Advection'},'FontWeight','Bold','FontSize',fs);
        
        subplot(163);
        
        eval([plotString '(intVConvV' flipString ')']);
        
        shading flat
        caxis([cMin.*10^3 cMax.*10^3]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Vertical';'Advection'},'FontWeight','Bold','FontSize',fs);
        
        subplot(164);
        
        eval([plotString '(intForcV' flipString ')']);
        
        shading flat
        caxis([cMin.*10^1 cMax.*10^1]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Surface Volume';'Forcing'},'FontWeight','Bold','FontSize',fs);
        
        subplot(165);
        
        eval([plotString '(intTotalV' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Total'},'FontWeight','Bold','FontSize',fs);
        
        subplot(166);
        
        eval([plotString '(intResidualV' flipString ')']);
        
        shading flat
        caxis([cMin.*10^-6 cMax.*10^-6]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title('Residual','FontWeight','Bold','FontSize',fs);
        
        drawnow
        
        if savePlot
            
            print('-dpng',[figureDir 'volume_budget' num2str(tt(timeStep)) '.png']);
            
        end
        
    end
    
    if plotBudgetSalinity
        
        hFig2 = figure(2);
        set(hFig2,'units','normalized','outerposition',[0 0 1 1]);
        set(gcf,'color',[1 1 1]);
        
        if useVol
            
            cMin = -10^6;
            cMax = 10^6;
            
        else
            
            cMin = -10^-7;
            cMax = 10^-7;
            
        end
        
        subplot(281);
        
        eval([plotString '(intTendS' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Salt Tendency,';'Timestep: ';num2str(tt(timeStep))},'FontWeight','Bold','FontSize',fs);
        
        subplot(282);
        
        eval([plotString '(intHAdvS' flipString ')']);
        
        shading flat
        caxis([cMin.*10^1 cMax.*10^1]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Horizontal';'Advection'},'FontWeight','Bold','FontSize',fs);
        
        subplot(283);
        
        pcolor(intVAdvS);
        
        eval([plotString '(intVAdvS' flipString ')']);
        
        shading flat
        caxis([cMin.*10^1 cMax.*10^1]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Vertical';'Advection'},'FontWeight','Bold','FontSize',fs);
        
        subplot(284);
        
        eval([plotString '(intHDifS' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Horizontal';'Diffusion'},'FontWeight','Bold','FontSize',fs);
        
        subplot(285);
        
        eval([plotString '(intVDifS' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Vertical';'Diffusion'},'FontWeight','Bold','FontSize',fs);
        
        subplot(286);
        
        eval([plotString '(intForcS' flipString ')']);
        
        shading flat
        caxis([cMin.*10^-2 cMax.*10^-2]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Surface';'Volume';'Forcing'},'FontWeight','Bold','FontSize',fs);
        
        subplot(287);
        
        eval([plotString '(intTotalS' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Total'},'FontWeight','Bold','FontSize',fs);
        
        subplot(288);
        
        eval([plotString '(intResidualS' flipString ')']);
        
        shading flat
        caxis([cMin.*10^-7 cMax.*10^-7]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Residual'},'FontWeight','Bold','FontSize',fs);
        
        subplot(289);
        
        eval([plotString '(intTendSal' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Salinity Tendency,';'Timestep: ';num2str(tt(timeStep))},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,8,10);
        
        eval([plotString '(intHAdvSal' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Horizontal';'Advection'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,8,11);
        
        eval([plotString '(intVAdvSal' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Vertical';'Advection'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,8,12);
        
        eval([plotString '(intHDifSal' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Horizontal';'Diffusion'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,8,13);
        
        eval([plotString '(intVDifSal' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Vertical';'Diffusion'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,8,14);
        
        eval([plotString '(intForcSal' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Surface';'Volume';'Forcing'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,8,15);
        
        pcolor(intTotalSal);
        
        eval([plotString '(intTotalSal' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Total'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,8,16);
        
        pcolor(intResidualSal);
        
        eval([plotString '(intResidualSal' flipString ')']);
        
        shading flat
        caxis([cMin.*10^-7 cMax.*10^-7]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Residual'},'FontWeight','Bold','FontSize',fs);
        
        if savePlot
            
            print('-dpng',[figureDir 'salinity_budget' num2str(tt(timeStep)) '.png']);
            
        end
        
        drawnow
        
    end
    
    if plotBudgetDIC
        
        hFig3 = figure(3);
        set(hFig3,'units','normalized','outerposition',[0 0 1 1]);
        set(gcf,'color',[1 1 1]);
        
        if useVol
            
            cMin = -10^6;
            cMax = 10^6;
            
        else
            
            cMin = -10^-7;
            cMax = 10^-7;
            
        end
        
        subplot(2,5,1);
        
        eval([plotString '(intTendDIC' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'DIC Tendency,';'Timestep: ';num2str(tt(timeStep))},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,2);
        
        eval([plotString '(intHAdvDIC' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Horizontal Advection'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,3);
        
        eval([plotString '(intVAdvDIC' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Vertical Advection'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,4);
        
        eval([plotString '(intHDifDIC' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Horizontal Diffusion'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,5);
        
        eval([plotString '(intVDifDIC' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Vertical Diffusion'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,6);
        
        eval([plotString '(intForcDIC' flipString ')']);
        
        shading flat
        caxis([cMin.*10^-1 cMax.*10^-1]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Air-sea CO_2 Flux'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,7);
        
        eval([plotString '(intVirtualFluxDIC' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'P/E/runoff'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,8);
        
        eval([plotString '(intBioDIC' flipString ')']);
        
        shading flat
        caxis([cMin.*10^-1 cMax.*10^-1]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Biology'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,9);
        
        eval([plotString '(intTotalDIC' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Total'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,10);
        
        eval([plotString '(intResidualDIC' flipString ')']);
        
        shading flat
        caxis([cMin.*10^-6 cMax.*10^-6]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Residual'},'FontWeight','Bold','FontSize',fs);
        
        drawnow
        
        if savePlot
            
            print('-dpng',[figureDir 'DIC_budget' num2str(tt(timeStep)) '.png']);
            
        end
        
    end
    
    if plotBudgetNO3
        
        hFig4 = figure(4);
        set(hFig4,'units','normalized','outerposition',[0 0 1 1]);
        set(gcf,'color',[1 1 1]);
        
        if useVol
            
            cMin = -10^4;
            cMax = 10^4;
            
        else
            
            cMin = -10^-9;
            cMax = 10^-9;
            
        end
        
        subplot(2,5,1);
        
        eval([plotString '(intTendNO3' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'NO3 Tendency, ';'Timestep: ';num2str(tt(timeStep))},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,2);
        
        eval([plotString '(intHAdvNO3' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Horizontal Advection'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,3);
        
        eval([plotString '(intVAdvNO3' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Vertical Advection'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,4);
        
        eval([plotString '(intHDifNO3' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Horizontal Diffusion'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,5);
        
        eval([plotString '(intVDifNO3' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Vertical Diffusion'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,6);
        
        eval([plotString '(intBioCNO3' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Bio Consum'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,7);
        
        eval([plotString '(intBioSNO3' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Bio Source'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,8);
        
        eval([plotString '(intGDARNO3' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'gDAR'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,9);
        
        eval([plotString '(intTotalNO3' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Total'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,10);
        
        eval([plotString '(intResidualNO3' flipString ')']);
        
        shading flat
        caxis([cMin.*10^-6 cMax.*10^-6]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Residual'},'FontWeight','Bold','FontSize',fs);
        
        drawnow
        
        if savePlot
            
            print('-dpng',[figureDir 'NO3_budget' num2str(tt(timeStep)) '.png']);
            
        end
        
    end
    
    if plotBudgetNO2
        
        hFig5 = figure(5);
        set(hFig5,'units','normalized','outerposition',[0 0 1 1]);
        set(gcf,'color',[1 1 1]);
        
        if useVol
            
            cMin = -10^3;
            cMax = 10^3;
            
        else
            
            cMin = -10^-10;
            cMax = 10^-10;
            
        end
        
        subplot(2,5,1);
        
        eval([plotString '(intTendNO2' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'NO2 Tendency, ';'Timestep: ';num2str(tt(timeStep))},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,2);
        
        eval([plotString '(intHAdvNO2' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Horizontal Advection'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,3);
        
        eval([plotString '(intVAdvNO2' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Vertical Advection'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,4);
        
        eval([plotString '(intHDifNO2' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Horizontal Diffusion'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,5);
        
        eval([plotString '(intVDifNO2' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Vertical Diffusion'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,6);
        
        eval([plotString '(intBioCNO2' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Bio Consum'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,7);
        
        eval([plotString '(intBioSNO2' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Bio Source'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,8);
        
        eval([plotString '(intGDARNO2' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'gDAR'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,9);
        
        eval([plotString '(intTotalNO2' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Total'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,10);
        
        eval([plotString '(intResidualNO2' flipString ')']);
        
        shading flat
        caxis([cMin.*10^-6 cMax.*10^-6]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Residual'},'FontWeight','Bold','FontSize',fs);
        
        drawnow
        
        if savePlot
            
            print('-dpng',[figureDir 'NO2_budget' num2str(tt(timeStep)) '.png']);
            
        end
        
    end
    
    if plotBudgetNH4
        
        hFig6 = figure(6);
        set(hFig6,'units','normalized','outerposition',[0 0 1 1]);
        set(gcf,'color',[1 1 1]);
        
        if useVol
            
            cMin = -10^2;
            cMax = 10^2;
            
        else
            
            cMin = -10^-10;
            cMax = 10^-10;
            
        end
        
        subplot(2,5,1);
        
        eval([plotString '(intTendNH4' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'NH4 Tendency, ';'Timestep: ';num2str(tt(timeStep))},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,2);
        
        eval([plotString '(intHAdvNH4' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Horizontal Advection'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,3);
        
        eval([plotString '(intVAdvNH4' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Vertical Advection'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,4);
        
        eval([plotString '(intHDifNH4' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Horizontal Diffusion'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,5);
        
        eval([plotString '(intVDifNH4' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Vertical Diffusion'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,6);
        
        eval([plotString '(intBioCNH4' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Bio Consum'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,7);
        
        eval([plotString '(intBioSNH4' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Bio Source'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,8);
        
        eval([plotString '(intGDARNH4' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'gDAR'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,9);
        
        eval([plotString '(intTotalNH4' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Total'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,10);
        
        eval([plotString '(intResidualNH4' flipString ')']);
        
        shading flat
        caxis([cMin.*10^-6 cMax.*10^-6]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Residual'},'FontWeight','Bold','FontSize',fs);
        
        drawnow
        
        if savePlot
            
            print('-dpng',[figureDir 'NH4_budget' num2str(tt(timeStep)) '.png']);
            
        end
        
    end
    
    if plotBudgetPO4
        
        hFig7 = figure(7);
        set(hFig7,'units','normalized','outerposition',[0 0 1 1]);
        set(gcf,'color',[1 1 1]);
        
        if useVol
            
            cMin = -10^2;
            cMax = 10^2;
            
        else
            
            cMin = -10^-10;
            cMax = 10^-10;
            
        end
        
        subplot(2,5,1);
        
        eval([plotString '(intTendPO4' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'PO4 Tendency, ';'Timestep: ';num2str(tt(timeStep))},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,2);
        
        eval([plotString '(intHAdvPO4' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Horizontal Advection';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,3);
        
        eval([plotString '(intVAdvPO4' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Vertical Advection'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,4);
        
        eval([plotString '(intHDifPO4' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Horizontal Diffusion'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,5);
        
        eval([plotString '(intVDifPO4' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Vertical Diffusion'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,6);
        
        eval([plotString '(intBioCPO4' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Bio Consum'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,7);
        
        eval([plotString '(intBioSPO4' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Bio Source'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,8);
        
        eval([plotString '(intGDARPO4' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'gDAR'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,9);
        
        eval([plotString '(intTotalPO4' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Total'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,10);
        
        eval([plotString '(intResidualPO4' flipString ')']);
        
        shading flat
        caxis([cMin.*10^-6 cMax.*10^-6]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Residual'},'FontWeight','Bold','FontSize',fs);
        
        drawnow
        
        if savePlot
            
            print('-dpng',[figureDir 'PO4_budget' num2str(tt(timeStep)) '.png']);
            
        end
        
    end
    
    if plotBudgetFe
        
        hFig8 = figure(8);
        set(hFig8,'units','normalized','outerposition',[0 0 1 1]);
        set(gcf,'color',[1 1 1]);
        
        if useVol
            
            cMin = -10^-1;
            cMax = 10^-1;
            
        else
            
            cMin = -10^-13;
            cMax = 10^-13;
            
        end
        
        subplot(2,5,1);
        
        eval([plotString '(intTendFe' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Fe Tendency, ';'Timestep: ';num2str(tt(timeStep))},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,2);
        
        eval([plotString '(intHAdvFe' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Horizontal Advection'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,3);
        
        eval([plotString '(intVAdvFe' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Vertical Advection'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,4);
        
        eval([plotString '(intHDifFe' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Horizontal Diffusion'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,5);
        
        eval([plotString '(intVDifFe' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Vertical Diffusion'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,6);
        
        eval([plotString '(intForcFe' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Surface Dust Flux'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,7);
        
        eval([plotString '(intSedFe' flipString ')']);
        
        shading flat
        caxis([cMin.*10^-1 cMax.*10^-1]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Sediment Flux'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,8);
        
        eval([plotString '(intBioFe' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Scav + Bio'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,9);
        
        eval([plotString '(intFreeFe' flipString ')']);
        
        shading flat
        caxis([cMin.*10^-5 cMax.*10^-5]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Free Iron'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,10);
        
        eval([plotString '(intResidualFe' flipString ')']);
        
        shading flat
        caxis([cMin.*10^-3 cMax.*10^-3]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Residual'},'FontWeight','Bold','FontSize',fs);
        
        drawnow
        
        if savePlot
            
            print('-dpng',[figureDir 'Fe_budget' num2str(tt(timeStep)) '.png']);
            
        end
    end
    
    if plotBudgetSiO2
        
        hFig9 = figure(9);
        set(hFig9,'units','normalized','outerposition',[0 0 1 1]);
        set(gcf,'color',[1 1 1]);
        
        if useVol
            
            cMin = -10^4;
            cMax = 10^4;
            
        else
            
            cMin = -10^-9;
            cMax = 10^-9;
            
        end
        
        subplot(2,5,1);
        
        eval([plotString '(intTendSiO2' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'SiO2 Tendency, ';'Timestep: ';num2str(tt(timeStep))},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,2);
        
        eval([plotString '(intHAdvSiO2' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Horizontal Advection'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,3);
        
        eval([plotString '(intVAdvSiO2' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Vertical Advection'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,4);
        
        eval([plotString '(intHDifSiO2' flipString ')']);
        
        shading flat
        caxis([cMin.*10^2 cMax.*10^2]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Horizontal Diffusion'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,5);
        
        eval([plotString '(intVDifSiO2' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Vertical Diffusion'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,6);
        
        eval([plotString '(intBioCSiO2' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Bio Consum'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,7);
        
        eval([plotString '(intBioSSiO2' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Bio Source'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,8);
        
        eval([plotString '(intGDARSiO2' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'gDAR'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,9);
        
        eval([plotString '(intTotalSiO2' flipString ')']);
        
        shading flat
        caxis([cMin cMax]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Total Budget'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,5,10);
        
        eval([plotString '(intResidualSiO2' flipString ')']);
        
        shading flat
        caxis([cMin.*10^-6 cMax.*10^-6]);
        colormap(colors);
        colorbar
        axis tight
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Residual'},'FontWeight','Bold','FontSize',fs);
        
        drawnow
        
        if savePlot
            
            print('-dpng',[figureDir 'SiO2_budget' num2str(tt(timeStep)) '.png']);
            
        end
        
    end
    
    %%
    %save budget
    
    if saveMat
        
        save([saveDir 'volume_' num2str(tt(timeStep)) '_budget.mat'],'intTendV','intHConvV','intVConvV','intForcV','intTotalV','intResidualV','-v7.3');
        
        save([saveDir 'salt_' num2str(tt(timeStep)) '_budget.mat'],'intTendS','intHAdvS','intVAdvS','intHDifS','intVDifS','intForcS','intTotalS','intResidualS','-v7.3');
        
        save([saveDir 'salinity_' num2str(tt(timeStep)) '_budget.mat'],'intTendSal','intHAdvSal','intVAdvSal','intHDifSal','intVDifSal','intForcSal','intTotalSal','intResidualSal','-v7.3');
        
        save([saveDir 'DIC_' num2str(tt(timeStep)) '_budget.mat'],'intTendDIC','intHAdvDIC','intVAdvDIC','intHDifDIC','intVDifDIC','intForcDIC', ...
            'intBioConsDIC','intBioConsDIC_PIC','intBioRespDIC','intBioReminDIC_DOC','intBioReminDIC_POC','intBioDisscDIC_PIC','intTotalDIC','intResidualDIC','-v7.3');
        
        save([saveDir 'NO3_' num2str(tt(timeStep)) '_budget.mat'],'intTendNO3','intHAdvNO3','intVAdvNO3','intHDifNO3','intVDifNO3','intBioCNO3','intBioSNO3','intTotalNO3','intResidualNO3','-v7.3');
        
        save([saveDir 'NO2_' num2str(tt(timeStep)) '_budget.mat'],'intTendNO2','intHAdvNO2','intVAdvNO2','intHDifNO2','intVDifNO2','intBioCNO2','intBioSNO2','intTotalNO2','intResidualNO2','-v7.3');
        
        save([saveDir 'NH4_' num2str(tt(timeStep)) '_budget.mat'],'intTendNH4','intHAdvNH4','intVAdvNH4','intHDifNH4','intVDifNH4','intBioCNH4','intBioSNH4','intTotalNH4','intResidualNH4','-v7.3');
        
        save([saveDir 'PO4_' num2str(tt(timeStep)) '_budget.mat'],'intTendPO4','intHAdvPO4','intVAdvPO4','intHDifPO4','intVDifPO4','intBioCPO4','intBioSPO4','intTotalPO4','intResidualPO4','-v7.3');
        
        save([saveDir 'Fe_' num2str(tt(timeStep)) '_budget.mat'],'intTendFe','intHAdvFe','intVAdvFe','intHDifFe','intVDifFe','intForcFe','intBioFe','intSedFe','intFreeFe','intTotalFe','intResidualFe','-v7.3');
        
        save([saveDir 'SiO2_' num2str(tt(timeStep)) '_budget.mat'],'intTendSiO2','intHAdvSiO2','intVAdvSiO2','intHDifSiO2','intVDifSiO2','intBioCSiO2','intBioSSiO2','intTotalSiO2','intResidualSiO2','-v7.3');
        
    end
    
    %%
    %clear gcmfaces objects to avoid memory leaks
    
    clear ETAN oceFWflx SFLUX  oceSPflx UVELMASS VVELMASS WVELMASS ...
        SALT ADVr_SLT ADVx_SLT ADVy_SLT DFrI_SLT DFrE_SLT DFxE_SLT DFyE_SLT oceSPtnd ...
        DIC ADVx_DIC ADVy_DIC ADVr_DIC DFxE_DIC DFyE_DIC DFrE_DIC DFrI_DIC DICTFLX ...
        CONSUMP_DIC CONSUMP_DIC_PIC REMIN_DOC REMIN_POC DISSC_PIC ...
        ETAN_SNAP SALT_SNAP DIC_SNAP NO3_SNAP NO2_SNAP NH4_SNAP PO4_SNAP Fe_SNAP SiO2_SNAP ...
        
    clear tendV adv_hConvV adv_vConvV forcV
    clear tendS adv_hConvS adv_vConvS dif_hConvS dif_vConvS forcS
    clear tendSal adv_hConvSal adv_vConvSal dif_hConvSal dif_vConvSal forcSal
    clear tendDIC adv_hConvDIC adv_vConvDIC dif_hConvDIC dif_vConvDIC forcDIC virtualFluxDIC gDAR_DIC bioDIC
    clear bioCons_DIC bioCons_DIC_PIC bioResp_DIC bioRemin_DIC_DOC bioRemin_DIC_POC bioDissc_DIC_PIC
    clear tendNO3 adv_hConvNO3 adv_vConvNO3 dif_hConvNO3 dif_vConvNO3 virtualFluxNO3 gDAR_NO3 bioC_NO3 bioS_NO3
    clear tendNO2 adv_hConvNO2 adv_vConvNO2 dif_hConvNO2 dif_vConvNO2 virtualFluxNO2 gDAR_NO2 bioC_NO2 bioS_NO2
    clear tendNH4 adv_hConvNH4 adv_vConvNH4 dif_hConvNH4 dif_vConvNH4 virtualFluxNH4 gDAR_NH4 bioC_NH4 bioS_NH4
    clear tendPO4 adv_hConvPO4 adv_vConvPO4 dif_hConvPO4 dif_vConvPO4 virtualFluxPO4 gDAR_PO4 bioC_PO4 bioS_PO4
    clear tendFe adv_hConvFe adv_vConvFe dif_hConvFe dif_vConvFe forcFe virtualFluxFe gDAR_Fe bioFe sedFe freeFe
    clear tendSiO2 dv_hConvSiO2 adv_vConvSiO2 dif_hConvSiO2 dif_vConvSiO2 virtualFluxSiO2 gDAR_SiO2 bioC_SiO2 bioS_SiO2
    
    if(savePlot)
        
        close all
        
    end
    
end
