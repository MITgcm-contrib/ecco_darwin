clear
close all

%%
%settings, modify as needed

kLevel = 1; %depth level for budget computation

%set to 1 to plot budget
plotPO4Budget = 0;
plotFeBudget = 1;

gridDir = '../../../../../darwin3/run/';
modelDir = '../../../../../darwin3/run/';

%%
%plotting params.

fs = 12;
lw = 2;

colors = parula(500);

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
filename2 = 'average_PO4_3d';
filename3 = 'average_Fe_3d';
filename4 = 'average_Fe_darwin_2d';

filename5 = 'snap_2d';
filename6 = 'snap_3d';

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
    
    %PO4 content
    PO4 = convert2gcmfaces(rdmds([diagDir filename2],ttAverage,'rec',1)) .* mmol_to_mol; %mol m^-3
    ADVx_PO4 = convert2gcmfaces(rdmds([diagDir filename2],ttAverage,'rec',2)) .* mmol_to_mol; %mol s^-1
    ADVy_PO4 = convert2gcmfaces(rdmds([diagDir filename2],ttAverage,'rec',3)) .* mmol_to_mol; %mol s^-1
    ADVr_PO4 = convert2gcmfaces(rdmds([diagDir filename2],ttAverage,'rec',4)) .* mmol_to_mol; %mol s^-1
    DFxE_PO4 = convert2gcmfaces(rdmds([diagDir filename2],ttAverage,'rec',5)) .* mmol_to_mol; %mol s^-1
    DFyE_PO4 = convert2gcmfaces(rdmds([diagDir filename2],ttAverage,'rec',6)) .* mmol_to_mol; %mol s^-1
    DFrE_PO4 = convert2gcmfaces(rdmds([diagDir filename2],ttAverage,'rec',7)) .* mmol_to_mol; %mol s^-1
    DFrI_PO4 = convert2gcmfaces(rdmds([diagDir filename2],ttAverage,'rec',8)) .* mmol_to_mol; %mol s^-1
    gDAR_PO4 = convert2gcmfaces(rdmds([diagDir filename2],ttAverage,'rec',9)) .* mmol_to_mol; %mol s^-1
    C_PO4 = convert2gcmfaces(rdmds([diagDir filename2],ttAverage,'rec',10)) .* mmol_to_mol;  %mol m^-3 s^-1
    S_PO4 = convert2gcmfaces(rdmds([diagDir filename2],ttAverage,'rec',11)) .* mmol_to_mol;  %mol m^-3 s^-1
    
    %Fe content
    Fe = convert2gcmfaces(rdmds([diagDir filename3],ttAverage,'rec',1)) .* mmol_to_mol; %mol m^-3
    ADVx_Fe = convert2gcmfaces(rdmds([diagDir filename3],ttAverage,'rec',2)) .* mmol_to_mol; %mol s^-1
    ADVy_Fe = convert2gcmfaces(rdmds([diagDir filename3],ttAverage,'rec',3)) .* mmol_to_mol; %mol s^-1
    ADVr_Fe = convert2gcmfaces(rdmds([diagDir filename3],ttAverage,'rec',4)) .* mmol_to_mol; %mol s^-1
    DFxE_Fe = convert2gcmfaces(rdmds([diagDir filename3],ttAverage,'rec',5)) .* mmol_to_mol; %mol s^-1
    DFyE_Fe = convert2gcmfaces(rdmds([diagDir filename3],ttAverage,'rec',6)) .* mmol_to_mol; %mol s^-1
    DFrE_Fe = convert2gcmfaces(rdmds([diagDir filename3],ttAverage,'rec',7)) .* mmol_to_mol; %mol s^-1
    DFrI_Fe = convert2gcmfaces(rdmds([diagDir filename3],ttAverage,'rec',8)) .* mmol_to_mol; %mol s^-1
    gDAR_Fe = convert2gcmfaces(rdmds([diagDir filename3],ttAverage,'rec',9)) .* mmol_to_mol;  %mol m^-3 s^-1
    C_Fe = convert2gcmfaces(rdmds([diagDir filename3],ttAverage,'rec',10)) .* mmol_to_mol;  %mol m^-3 s^-1
    S_Fe = convert2gcmfaces(rdmds([diagDir filename3],ttAverage,'rec',11)) .* mmol_to_mol;  %mol m^-3 s^-1
    SEDFe = convert2gcmfaces(rdmds([diagDir filename3],ttAverage,'rec',12)) .* mmol_to_mol;  %mol m^-3 s^-1
    FREEFe = convert2gcmfaces(rdmds([diagDir filename3],ttAverage,'rec',13)) .* mmol_to_mol;  %mol m^-3 s^-1
    SFCSOLFe = convert2gcmfaces(rdmds([diagDir filename4],ttAverage,'rec',1)) .* mmol_to_mol;  %mol m^-2 s^-1
    
    %%
    %load snapshots
    
    ETAN_SNAP = nan*ones(nx,ny,2);
    PO4_SNAP = nan .* ones(nx,ny,numLevels,2);
    Fe_SNAP = nan .* ones(nx,ny,numLevels,2);
    
    if timeStep == 1
        
        %no pre-initial snapshot
        ETAN_SNAP(:,:,2) = rdmds([diagDir filename5],tt(1),'rec',1);
        PO4_SNAP(:,:,:,2) = rdmds([diagDir filename6],tt(1),'rec',6) .* mmol_to_mol;  %mol m^-3 s^-1
        Fe_SNAP(:,:,:,2) = rdmds([diagDir filename6],tt(1),'rec',7) .* mmol_to_mol;  %mol m^-3 s^-1
        
    elseif timeStep == numFiles %no final snapshot
        
        ETAN_SNAP(:,:,1) = rdmds([diagDir filename5],tt(timeStep-1),'rec',1);
        PO4_SNAP(:,:,:,1) = rdmds([diagDir filename6],tt(timeStep-1),'rec',6) .* mmol_to_mol;  %mol m^-3 s^-1
        Fe_SNAP(:,:,:,1) = rdmds([diagDir filename6],tt(timeStep-1),'rec',7) .* mmol_to_mol;  %mol m^-3 s^-1
        
    else %timeStep~=1 & timeStep~=numFiles
        
        ETAN_SNAP = rdmds([diagDir filename5],ttSnap,'rec',1);
        PO4_SNAP = rdmds([diagDir filename6],ttSnap,'rec',6) .* mmol_to_mol;  %mol m^-3 s^-1
        Fe_SNAP = rdmds([diagDir filename6],ttSnap,'rec',7) .* mmol_to_mol;  %mol m^-3 s^-1
        
    end
    
    PO4_SNAP = convert2gcmfaces(PO4_SNAP);
    Fe_SNAP = convert2gcmfaces(Fe_SNAP);
    
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
    
    %darwin tendency
    gDAR_PO4 = mygrid.mskC .* (gDAR_PO4 ./ mygrid.hFacC);
    
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
    
    gDAR_Fe = mygrid.mskC .* (gDAR_Fe ./ mygrid.hFacC);
    
    %%
    %convert gcmfaces objects to matrices
 
    tendPO4 = convert2gcmfaces(tendPO4);
    adv_hConvPO4 = convert2gcmfaces(adv_hConvPO4);
    adv_vConvPO4 = convert2gcmfaces(adv_vConvPO4);
    dif_hConvPO4 = convert2gcmfaces(dif_hConvPO4);
    dif_vConvPO4 = convert2gcmfaces(dif_vConvPO4);
    gDARPO4 = convert2gcmfaces(gDAR_PO4);

    tendFe = convert2gcmfaces(tendFe);
    adv_hConvFe = convert2gcmfaces(adv_hConvFe);
    adv_vConvFe = convert2gcmfaces(adv_vConvFe);
    dif_hConvFe = convert2gcmfaces(dif_hConvFe);
    dif_vConvFe = convert2gcmfaces(dif_vConvFe);
    gDARFe = convert2gcmfaces(gDAR_Fe);
    
    %%
     
    budget.TendPO4 = tendPO4(:,:,kLevel);
    budget.HAdvPO4 = adv_hConvPO4(:,:,kLevel);
    budget.VAdvPO4 = adv_vConvPO4(:,:,kLevel);
    budget.HDifPO4 = dif_hConvPO4(:,:,kLevel);
    budget.VDifPO4 = dif_vConvPO4(:,:,kLevel);
    budget.GDARPO4 = gDARPO4(:,:,kLevel);
    
    budget.PO4Total = budget.HAdvPO4 + budget.VAdvPO4 + budget.HDifPO4 + budget.VDifPO4 + budget.GDARPO4;
    budget.PO4Residual = budget.TendPO4 - budget.PO4Total;
    
    budget.TendFe = tendFe(:,:,kLevel);
    budget.HAdvFe = adv_hConvFe(:,:,kLevel);
    budget.VAdvFe = adv_vConvFe(:,:,kLevel);
    budget.HDifFe = dif_hConvFe(:,:,kLevel);
    budget.VDifFe = dif_vConvFe(:,:,kLevel);
    budget.GDARFe = gDARFe(:,:,kLevel);
    
    budget.FeTotal = budget.HAdvFe + budget.VAdvFe + budget.HDifFe + budget.VDifFe + budget.GDARFe;
    budget.FeResidual = budget.TendFe - budget.FeTotal;
    
    %%
    %plot budgets
    
    if plotPO4Budget
        
        hFig1 = figure(1);
        set(hFig1,'units','normalized','outerposition',[0 0 1 1]);
        set(gcf,'color',[1 1 1]);
        set(gca,'color',[0.5 0.5 0.5]);
        
        cMin = -10^-10;
        cMax = 10^-10;
        
        subplot(2,4,1);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,budget.TendPO4);
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
        
        title({'PO4 Tendency (mol m^-^3 s^-^1), ';'timeStep: ';num2str(tt(timeStep))},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,4,2);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,budget.HAdvPO4);
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
        
        subplot(2,4,3);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,budget.VAdvPO4);
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
        
        subplot(2,4,4);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,budget.HDifPO4);
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
        
        subplot(2,4,5);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,budget.VDifPO4);
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
        
        subplot(2,4,6);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,budget.GDARPO4);
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
        
        subplot(2,4,7);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,budget.PO4Total);
        
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
        
        title({'Total Budget';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,4,8);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,budget.PO4Residual);
        
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
        
        title({'Budget Total - Tendency';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        drawnow
         
    end
    
    if plotFeBudget
        
        hFig2 = figure(2);
        set(hFig2,'units','normalized','outerposition',[0 0 1 1]);
        set(gcf,'color',[1 1 1]);
        set(gca,'color',[0.5 0.5 0.5]);
        
        cMin = -10^-14;
        cMax = 10^-14;
        
        subplot(2,4,1);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,budget.TendFe);
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
        
        subplot(2,4,2);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,budget.HAdvFe);
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
        
        subplot(2,4,3);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,budget.VAdvFe);
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
        
        subplot(2,4,4);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,budget.HDifFe);
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
        
        subplot(2,4,5);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,budget.VDifFe);
        shading flat
        
        caxis([cMin cMax]);
        
        colormap(colors);
        
        colorbar
        
        axis tight
        
        ylabel('Latitude (deg)');
        
        xlim([0 360]);
        ylim([-90 90]);
        
        set(gca,'xtick',[0:180:360]);
        set(gca,'ytick',[-90:30:90]);
        
        box on
        set(gca,'LineWidth',lw);
        set(gca,'FontSize',fs);
        
        title({'Vertical Diffusion';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
      
        subplot(2,4,6);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,budget.GDARFe);
        
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
        
        title({'gDAR Fe';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);
        
        subplot(2,4,7);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,budget.FeTotal);
        
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
        
        title({'Budget Total';'(mol m^-^3 s^-^1)'},'FontWeight','Bold','FontSize',fs);

        subplot(2,4,8);
        
        hold on
        
        set(gca,'color',[0.5 0.5 0.5]);
        
        pcolorcen(mygrid.XC.f1,mygrid.YC.f1,budget.FeResidual);
        
        shading flat
        
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
    
end
