clear
close all;

gridDir = '/Users/carrolld/Documents/research/GoM/grid/';
writeDir = '/Users/carrolld/Documents/research/GoM/model_setup/diffusivity_mask/';

cd(writeDir);

delete *

%%

nx = 960;
ny = nx;
nz = 90;

XC = readbin([gridDir 'XC.data'],[nx ny 1],1,'real*4');
YC = readbin([gridDir 'YC.data'],[nx ny 1],1,'real*4');

hFacC = readbin([gridDir 'hFacC.data'],[nx ny 1],1,'real*4');
depth = readbin([gridDir 'depth.data'],[nx ny],1,'real*4');

%%

temp = ones(nx,ny);

diffKrFile = 'diff100x_GoM960x960x90_runoff';
diffKrT = 5.44e-7 * 100;

temp = temp.*diffKrT;
temp(depth < 50) = 10^-3;

pcolorcen(temp);

%%

temp = repmat(temp,[1 1 nz]);

writebin([writeDir diffKrFile],temp);

%%

unique(temp)

%%

