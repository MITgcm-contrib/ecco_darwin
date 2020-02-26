clear
close all

saveMovie = 1;

gridDir = '/Users/carrolld/Documents/research/v4_3deg/MITgcm/run_zstar_months_mpi/';
modelDir = '/Users/carrolld/Documents/research/v4_3deg/MITgcm/run_zstar_months_mpi/';
saveDir = '/Users/carrolld/Documents/research/carbon/mat/hong_LLC270/';
movieDir = '/Users/carrolld/Documents/research/carbon/figures/CO2_flux_budget/testing/v4_3deg/';

%% 

nx = 128;
ny = 64;
nz = 15;

XC = readbin([modelDir 'XC.data'],[nx ny],1,'real*4');
YC = readbin([modelDir 'YC.data'],[nx ny],1,'real*4');
hFacC = readbin([modelDir 'hFacC.data'],[nx ny nz],1,'real*4');
RAC = readbin([modelDir 'RAC.data'],[nx ny],1,'real*4');
DXG = readbin([modelDir 'DXG.data'],[nx ny],1,'real*4');
DYG = readbin([modelDir 'DYG.data'],[nx ny],1,'real*4');
DRF = readbin([modelDir 'dRF.data'],[nz],1,'real*4');
RC = readbin([modelDir 'RC.data'],[nz],1,'real*4');

depth = readbin([modelDir 'Depth.data'],[nx ny],1,'real*4');

dzMatF = repmat(DRF,1,nx,ny);
dzMatF = permute(dzMatF,[2 3 1]);
dzMat = dzMatF .* hFacC;

RACMat = repmat(RAC,1,1,nz);

mskC = hFacC .* 1;
mskC(hFacC == 0) = nan;

VVV = mskC .* hFacC .* RACMat .* dzMatF;
nLevels = numel(RC);
%secPerHour = 3600;

rhoconst=1029;

%% 

modelDir = [modelDir 'diags/'];
nameFld1 = 'budg2d_zflux_set1';
nameFld6 = 'state_2d_set1';
nameFld8 = 'budg2d_snap_set1';
nameFld9 = 'trsp_3d_set1';

files = dir([modelDir nameFld6 '.*.data']);
timeStep = zeros(length(files),1);

for i=1:length(timeStep)

    filename = files(i).name;
    timeStep(i) = str2num(filename(end-10:end-5));

end

numFiles = numel(timeStep);

dt = diff([0 timeStep']);
dt = dt ./72 .* 24; %timeStepsteps in hours, 72 timeStepsteps/day * 24 hours/day = hours

numTimeSteps = numel(timeStep);
secPerHour = 3600;

%% 

for nTimeSteps = 1:numTimeSteps, disp(num2str(nTimeSteps))
    
    tic
    
    if nTimeSteps == 1
        
        timeStep_average = timeStep(nTimeSteps);
        timeStep_snap = [1 timeStep(nTimeSteps)];
        
    else % nTimeSteps~=1
        
        timeStep_average = timeStep(nTimeSteps);
        timeStep_snap = [timeStep(nTimeSteps-1) timeStep(nTimeSteps)];
        
    end
    
    % load two-dimensional monthly averaged fields
    oceFWflx = rdmds([modelDir nameFld1],timeStep_average,'rec',1);
    ETAN = rdmds([modelDir nameFld6],timeStep_average,'rec',1);
    
    % load three-dimensional monthly averaged fields
    UVELMASS = rdmds([modelDir nameFld9],timeStep_average,'rec',1);
    VVELMASS = rdmds([modelDir nameFld9],timeStep_average,'rec',2);
    WVELMASS = rdmds([modelDir nameFld9],timeStep_average,'rec',3);
    
    % load "snapshot" fields (for computing total tendencies)
    ETAN_SNAP_raw = nan .* ones(nx,ny,2);
    
    if nTimeSteps == 1
        
        % no pre-initial snapshot
        ETAN_SNAP_raw(:,:,2) = rdmds([modelDir nameFld8],timeStep(1),'rec',1);
        
    elseif nTimeSteps == numFiles %no final snapshot
        
        ETAN_SNAP_raw(:,:,1) = rdmds([modelDir nameFld8],timeStep(nTimeSteps-1),'rec',1);
        
    else % nTimeSteps~=1 & nTimeSteps~=numFiles
        
        ETAN_SNAP_raw = rdmds([modelDir nameFld8],timeStep_snap,'rec',1);
        
    end
    
    ETAN_SNAP  = ETAN_SNAP_raw;
    
   %% 
 
    %total tendency
    tendM = 1 ./ repmat(depth,1,1,nz) .* repmat((ETAN_SNAP(:,:,2)-ETAN_SNAP(:,:,1)) ./ (secPerHour .* dt(nTimeSteps)),1,1,nz);

    facW = repmat(DYG,1,1,nz);
    facS = repmat(DXG,1,1,nz);
    
    FLDU = UVELMASS .* facW;
    FLDV = VVELMASS .* facS;
    
    FLDU(nx+1,:,:) = FLDU(1,:,:);
    FLDV(:,ny+1,:) = FLDV(:,end,:);
    
    fldDIV = (FLDU(1:end-1,:,:) - FLDU(2:end,:,:)) + ...
              (FLDV(:,1:end-1,:) - FLDV(:,2:end,:));

    hConvM = mskC .* fldDIV ./ (RACMat .* hFacC);
          
    %vertical divergence
    vConvM = 0 .* hConvM;
    
    for k = 1:nLevels
    
        kp1 = min([k+1,nLevels]);
        
        vConvM(:,:,k) = squeeze((WVELMASS(:,:,kp1) .* double(k~=nLevels)) - (WVELMASS(:,:,k) .* double(k~=1))) ...
            ./ (dzMat(:,:,k));   
        
    end
    
    %forcing
    forcM = mskC .* repmat(oceFWflx,1,1,nz) ./ (dzMat .* rhoconst);   
    forcM(:,:,2:nz) = 0 .* mskC(:,:,2:nz);
    
    intLevel = 1;
    
    intTendV = nansum(tendM(:,:,1:intLevel) .* dzMat(:,:,1:intLevel),3);
    intHConvV = nansum(hConvM(:,:,1:intLevel) .* dzMat(:,:,1:intLevel),3);
    intVConvV = nansum(vConvM(:,:,1:intLevel) .* dzMat(:,:,1:intLevel),3);
    intForcV = nansum(forcM(:,:,1:intLevel) .* dzMat(:,:,1:intLevel),3);
 
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
    
    pcolorcen(XC,YC,intTendV);
    
    caxis([cMin cMax]);
    
    colormap(colors);
    
    colorbar
    
    axis tight
    
    box on
    set(gca,'LineWidth',lw);
    set(gca,'FontSize',fs);
    
    title(['Tendency, NT: ' num2str(nTimeSteps)]);
    
    subplot(162);
    
    hold on
    
    set(gca,'color',[0.5 0.5 0.5]);
    
    pcolorcen(XC,YC,intHConvV);
    
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
    
    pcolorcen(XC,YC,intVConvV);
    
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
    
    pcolorcen(XC,YC,intForcV);
    
    caxis([cMin*10 cMax*10]);
    
    colormap(colors);
    
    colorbar
    
    axis tight
        
    box on
    set(gca,'LineWidth',lw);
    set(gca,'FontSize',fs);

    title('Surface Volume Forcing');

    subplot(165);
    
    hold on
    
    set(gca,'color',[0.5 0.5 0.5]);
    
    pcolorcen(XC,YC,intHConvV + intVConvV + intForcV);
    
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
    
    pcolorcen(XC,YC,intTendV - (intForcV + intHConvV + intVConvV));
    
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
        
        export_fig([movieDir 'volume_budget_' num2str(nTimeSteps) '.png'],'-png');
        
    end
    
end
