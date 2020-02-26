clear
close all

saveMovie = 1;

gridDir = '/Users/carrolld/Documents/research/v4_3deg/MITgcm/run_zstar_months_mpi/';
modelDir = '/Users/carrolld/Documents/research/v4_3deg/MITgcm/run_zstar_months_mpi/';

movieDir = '/Users/carrolld/Documents/research/carbon/movies/CO2_flux_budget/from_hong/v4_3deg/';

%%

nF = 1;
fileFormat='straight';

%mygrid = [];
global mygrid
grid_load(gridDir,nF,fileFormat);

nx = mygrid.ioSize(1);
ny = mygrid.ioSize(2);
nL = numel(mygrid.RC) ;

dzMatF = mk3D(mygrid.DRF, mygrid.hFacC);
dzMat = dzMatF .* mygrid.hFacC;
RACMat = mk3D(mygrid.RAC, mygrid.hFacC);
VVV = mygrid.mskC .* mygrid.hFacC .* mk3D(mygrid.RAC,mygrid.mskC) .* mk3D(mygrid.DRF,mygrid.mskC);
nLevels = numel(mygrid.RC);

%%

dirIn = [modelDir 'diags/'];
nameFld1 = 'budg2d_zflux_set1';
nameFld6 = 'state_2d_set1';
nameFld8 = 'budg2d_snap_set1';
nameFld9 = 'trsp_3d_set1';

%%

fn = dir([dirIn nameFld6 '.*.data']);
tt = zeros(length(fn),1);

for i = 1:length(tt)
    
    nme = fn(i).name;
    tt(i) = str2num(nme(end-10:end-5));
    
end

TotalNumOfFiles = numel(tt);
dt = diff([0 tt']);
dt = dt ./72 .*24; %hours
numTimeSteps = numel(tt);

secPerHour = 3600;
rhoconst = 1029;

%%

for NT = 1:12
    
    disp(num2str(NT));
    
    tic
    
    %set file numbers to load in
    if NT==1
        
        tt_avrg = tt(NT);
        tt_snap = [1 tt(NT)];
        
    else %NT~=1
        
        tt_avrg = tt(NT);
        tt_snap = [tt(NT-1) tt(NT)];
        
    end
    
    %load two-dimensional monthly averaged fields
    oceFWflx = convert2gcmfaces(rdmds([dirIn nameFld1],tt_avrg,'rec',1));
    ETAN = convert2gcmfaces(rdmds([dirIn nameFld6],tt_avrg,'rec',1));
    
    %load three-dimensional monthly averaged fields
    UVELMASS = convert2gcmfaces(rdmds([dirIn nameFld9],tt_avrg,'rec',1));
    VVELMASS = convert2gcmfaces(rdmds([dirIn nameFld9],tt_avrg,'rec',2));
    WVELMASS = convert2gcmfaces(rdmds([dirIn nameFld9],tt_avrg,'rec',3));
    
    %%
    
    ETAN_SNAP_raw = nan*ones(nx,ny,2);
    
    if NT == 1
        
        %no pre-initial snapshot
        ETAN_SNAP_raw(:,:,2) = rdmds([dirIn nameFld8],tt(1),'rec',1);
        
    elseif NT == TotalNumOfFiles %no final snapshot
        
        ETAN_SNAP_raw(:,:,1) = rdmds([dirIn nameFld8],tt(NT-1),'rec',1);
        
    else %NT~=1 & NT~=TotalNumOfFiles
        
        ETAN_SNAP_raw  = rdmds([dirIn nameFld8],tt_snap,'rec',1);
        
    end
    
    ETAN_SNAP = convert2gcmfaces(ETAN_SNAP_raw);
    
    %%
    
    %total tendency
    tendM = (1 ./ mk3D(mygrid.Depth,mygrid.mskC)) .* mk3D((ETAN_SNAP(:,:,2) - ETAN_SNAP(:,:,1)) ...
        ./ (secPerHour .* dt(NT)),mygrid.mskC);
    
    %horizontal convergence
    hConvM = mygrid.mskC .* calc_UV_conv(UVELMASS,VVELMASS,{'dh'}) ./ (RACMat.*mygrid.hFacC);
    
    %vertical divergence
    vConvM = 0*hConvM;
    
    for nz = 1:nLevels
        
        nzp1 = min([nz+1,nLevels]);
        
        vConvM(:,:,nz) = squeeze(WVELMASS(:,:,nzp1)*double(nz~=nLevels) - WVELMASS(:,:,nz)*double(nz~=1)) ...
            ./ (dzMat(:,:,nz));
        
    end
    
    %forcing
    forcM = mygrid.mskC .* mk3D(oceFWflx,mygrid.mskC) ./ (dzMat .* rhoconst);
    forcM(:,:,2:nL) = 0 .* mygrid.mskC(:,:,2:nL);
    
    TOT_V = tendM;
    FRC_V = forcM;
    hCON_V = hConvM;
    vCON_V = vConvM;
    
    tendV = convert2gcmfaces(TOT_V);
    hConvV = convert2gcmfaces(hCON_V);
    vConvV = convert2gcmfaces(vCON_V);
    forcV = convert2gcmfaces(FRC_V);
    
    %%
    
    intLevel = 1;
    
    intTendV = nansum(tendV(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel)),3);
    intHConvV = nansum(hConvV(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel)),3);
    intVConvV = nansum(vConvV(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel)),3);
    intForcV = nansum(forcV(:,:,1:intLevel) .* convert2gcmfaces(dzMat(:,:,1:intLevel)),3);
    
    hFig1 = figure(1);
    set(hFig1,'units','normalized','outerposition',[0 0 1 1]);
    set(gcf,'color',[1 1 1]);
    set(gca,'color',[0.5 0.5 0.5]);
    
    fs = 16;
    lw = 2;
    
    colors = cbrewer('div','RdBu',100);
    cMin = -10^-8;
    cMax = 10^-8;
    
    subplot(161);
    
    hold on
    
    set(gca,'color',[0.5 0.5 0.5]);
    
    pcolorcen(intTendV');
    
    caxis([cMin cMax]);
    
    colormap(colors);
    
    colorbar
    
    axis tight
    
    box on
    set(gca,'LineWidth',lw);
    set(gca,'FontSize',fs);
    
    title(['Tendency, NT: ' num2str(NT)]);
    
    subplot(162);
    
    hold on
    
    set(gca,'color',[0.5 0.5 0.5]);
    
    pcolorcen(intHConvV');
    
    caxis([cMin*1000 cMax*1000]);
    
    colormap(colors);
    
    colorbar
    
    axis tight
    
    box on
    set(gca,'LineWidth',lw);
    set(gca,'FontSize',fs);
    
    title('Hoz Conv');
    
    subplot(163);
    
    hold on
    
    set(gca,'color',[0.5 0.5 0.5]);
    
    pcolorcen(intVConvV');
    
    caxis([cMin*1000 cMax*1000]);
    
    colormap(colors);
    
    colorbar
    
    axis tight
    
    box on
    set(gca,'LineWidth',lw);
    set(gca,'FontSize',fs);
    
    title('Vert Conv');
    
    subplot(164);
    
    hold on
    
    set(gca,'color',[0.5 0.5 0.5]);
    
    pcolorcen(intForcV');
    
    caxis([cMin*10 cMax*10]);
    
    colormap(colors);
    
    colorbar
    
    axis tight
    
    box on
    set(gca,'LineWidth',lw);
    set(gca,'FontSize',fs);
    
    title({'Surface Volume';'Forcing'});
    
    subplot(165);
    
    hold on
    
    set(gca,'color',[0.5 0.5 0.5]);

    pcolorcen((intHConvV + intVConvV + intForcV)');
    
    caxis([cMin cMax]);
    
    colormap(colors);
    
    colorbar
    
    axis tight
    
    box on
    set(gca,'LineWidth',lw);
    set(gca,'FontSize',fs);
    
    title('Total');
    
    subplot(166);
    
    hold on
    
    set(gca,'color',[0.5 0.5 0.5]);
  
    pcolorcen((intTendV - (intForcV + intHConvV + intVConvV))');
    
    caxis([cMin cMax]);
    
    colormap(colors);
    
    colorbar
    
    axis tight
    
    box on
    set(gca,'LineWidth',lw);
    set(gca,'FontSize',fs);
    
    title('Residual');
    
    drawnow
    
    if saveMovie
        
        export_fig([movieDir 'volume_budget_' num2str(NT) '.png'],'-png');
        
    end
    
    toc
    
end
