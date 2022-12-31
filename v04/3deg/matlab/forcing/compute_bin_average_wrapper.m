clear
close all

%% 

codeDir = '/Users/carrolld/Documents/research/carbon/m_files/bin_average/standard/';
gridDir = '/Users/carrolld/Documents/research/carbon/simulations/grid/LLC_270/';
saveDir = '/Users/carrolld/Documents/research/carbon/mat/bin_average/standard/';

saveFilename = 'bin_average_test.mat'

numFaces = 13;
nx = 270;
ny = 270 .* numFaces;

XC = readbin([gridDir 'XC.data'],[nx ny],1,'real*4');
YC = readbin([gridDir 'YC.data'],[nx ny],1,'real*4');

RAC = readbin([gridDir 'RAC.data'],[nx ny],1,'real*4');

cd(codeDir);
compute_bin_average(XC,YC,RAC,3,3,saveDir,saveFilename);

%% 

clear
close all

nx = 270;

gridDir = '/Users/carrolld/Documents/research/carbon/simulations/grid/LLC_270/';
dataDir = '/Users/carrolld/Documents/research/carbon/m_files/CMS/abhishek/offline/'
saveDir = '/Users/carrolld/Documents/research/carbon/mat/bin_average/standard/';

load([saveDir 'bin_average_test.mat']);

CFLX = quikread_llc([dataDir 'DICCFLX3hrly.0000655272.data'],270);
cflx = reshape(bin_average*double(CFLX(:)),length(lon),length(lat));

figure(1);
quikplot_llc(CFLX);
colorbar;

figure(2);
pcolorcen(cflx');
colorbar;

% verify conservation
RAC = quikread_llc([gridDir 'RAC.data'],nx);

disp((sum(sum(RAC.*CFLX)) - sum(sum(AREA.*cflx))) / sum(sum(AREA.*cflx)));
