%Test code to read in GoM files and upload...
% to github for practice 

addpath('/nobackup/jzaissbo/matlabcode/lib');

nx=20;
ny=15;
nz=47;


DataDir='/nobackup/dmenemen/ecco_darwin/GoM/';

theta=readbin([DataDir 'THETA/THETA_' num2str(nx) ...
    'x' num2str(ny) 'x' num2str(nz) '.' ...
    '20020821T120000'],[nx ny nz]);

oxy=readbin([DataDir 'O2/O2_' num2str(nx) ...
    'x' num2str(ny) 'x' num2str(nz) '.' ...
    '20020821T120000'],[nx ny nz]);

mask=readbin([DataDir 'grid/hFacC_' num2str(nx) ...
    'x' num2str(ny) 'x' num2str(nz)],[nx ny nz]);

tmp=theta(:,:,1);
tmp(mask(:,:,1)==0) = nan;

figure()
pcolorcen(tmp')
c=colorbar();
colormap(cmocean('haline'))
title('Theta [^oC]')

tmp=oxy(:,:,40);
tmp(mask(:,:,40)==0)=nan;

figure()
pcolorcen(tmp')
c=colorbar();
colormap(cmocean('matter'))
title('Oxygen [mmol/m^3]')

