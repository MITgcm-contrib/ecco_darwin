clear 
close all

%compare snapshots of (SALT * WVELMASS) / s*
%to Dimitris' new WVELSALT snapshot

gridDir = '../../../../../darwin3/run/';
dataDir1 = '../../../../../darwin3/run/diags/budget/';
dataDir2 = '../../../../../darwin3/run_equation_12/diags/budget/';

%% 

nx = 128;
ny = 64;
nz = 15;

XC = readbin([gridDir 'XC.data'],[nx ny],1,'real*8');
YC = readbin([gridDir 'YC.data'],[nx ny],1,'real*8');
hFacC = readbin([gridDir 'hFacC.data'],[nx ny nz],1,'real*8');
RAC = readbin([gridDir 'RAC.data'],[nx ny],1,'real*8');
DXG = readbin([gridDir 'DXG.data'],[nx ny],1,'real*8');
DYG = readbin([gridDir 'DYG.data'],[nx ny],1,'real*8');
DRF = readbin([gridDir 'dRF.data'],[nz],1,'real*8');
RC = readbin([gridDir 'RC.data'],[nz],1,'real*8');

depth = readbin([gridDir 'Depth.data'],[nx ny],1,'real*8');

%% 

%old offline diagnostics
fileName1 = 'snap_2d';
fileName2 = 'snap_3d';
fileName3 = 'snap_velmass_3d';

%new online diagnostic
fileName4 = 'snap_velmass_3d'; 

%% 

intLevel = 2;

for i = 1:23

    ttSnap = [i-1 i];

    ETAN_SNAP = rdmds([dataDir1 fileName1],ttSnap,'rec',1);
    SALT_SNAP = rdmds([dataDir1 fileName2],ttSnap,'rec',1);
    WVELMASS_SNAP = rdmds([dataDir1 fileName3],ttSnap,'rec',3);
    WVELSLT_SNAP = rdmds([dataDir2 fileName4],ttSnap,'rec',3); 

    %old offline diagnostics
    rstarfac = (depth + ETAN_SNAP(:,:,1)) ./ depth;
    SALT = SALT_SNAP(:,:,intLevel,1);
    WVELMASS = WVELMASS_SNAP(:,:,intLevel,1);
 
    test1 = (-SALT .* WVELMASS) ./ rstarfac;
    
    %new online diagnostic
    test2 = -WVELSLT_SNAP(:,:,intLevel,1); 
    
    hold on
    
    set(gca,'color',[0.5 0.5 0.5]);
    
    pcolorcen(test1' - test2');
    
    colorbar
    caxis([-10^-5 10^-5]);
    
    drawnow
    
    clear temp
    
    disp(num2str(i));
    
end
