clear 
close all

dataDir = '/Users/carrolld/Documents/research/carbon/mat/bin_average/exf/';

cd(dataDir);

%% 

load([dataDir 'atemp.mat']);

[xx yy] = meshgrid(blon,blat);
xx = xx';
yy = yy';

for i = 1:12
    
    subplot(121);
    
    pcolorcen(atemp.lon,atemp.lat,squeeze(atemp.meanField(:,:,i)));
    
    colorbar;
    
    caxis([0 30]);
    
    subplot(122);
    
    pcolorcen(xx,yy,squeeze(atemp.binMeanField(:,:,i)));
    
    colorbar;
    
    caxis([0 30]);

    drawnow
    
    pause(0.1)
    
end

%% 

load([dataDir 'aqh.mat']);

for i = 1:12
    
    subplot(121);
    
    pcolorcen(aqh.lon,aqh.lat,squeeze(aqh.meanField(:,:,i)));
    
    colorbar;
    
    caxis([0 0.05]);
    
    subplot(122);
    
    pcolorcen(xx,yy,squeeze(aqh.binMeanField(:,:,i)));
    
    colorbar;
    
    caxis([0 0.05]);
    
    drawnow
    
    pause(0.1)
    
end

%% 

load([dataDir 'precip.mat']);

for i = 1:12
    
    subplot(121);
    
    pcolorcen(precip.lon,precip.lat,squeeze(precip.meanField(:,:,i)));
    
    colorbar
    
    caxis([0 10^-7]);
    
    subplot(122);
    
    pcolorcen(xx,yy,squeeze(precip.binMeanField(:,:,i)));
    
    colorbar
    
    caxis([0 10^-7]);
    
    drawnow
    
    pause(0.1)
    
end

%% 

load([dataDir 'uwind.mat']);

for i = 1:12
    
    subplot(121);
    
    pcolorcen(uwind.lon,uwind.lat,squeeze(uwind.meanField(:,:,i)));
    
    colorbar;
    
    caxis([0 10]);
    
    subplot(122);
    
    pcolorcen(xx,yy,squeeze(uwind.binMeanField(:,:,i)));
    
    colorbar;
    
    caxis([0 10]);
    
    drawnow
    
    pause(0.1)
    
end

%% 

load([dataDir 'vwind.mat']);

for i = 1:12
    
    subplot(121);
    
    pcolorcen(vwind.lon,vwind.lat,squeeze(vwind.meanField(:,:,i)));
    
    colorbar;
    
    caxis([0 10]);
    
    subplot(122);
    
    pcolorcen(xx,yy,squeeze(vwind.binMeanField(:,:,i)));
    
    colorbar;
    
    caxis([0 10]);
    
    drawnow
    
    pause(0.1)
    
end

%% 

load([dataDir 'swdown.mat']);

for i = 1:12
    
    subplot(121);
    
    pcolorcen(swdown.lon,swdown.lat,squeeze(swdown.meanField(:,:,i)));
    
    colorbar;
    
    caxis([-400 50]);
    
    subplot(122);
    
    pcolorcen(xx,yy,squeeze(swdown.binMeanField(:,:,i)));
    
    colorbar;
    
    caxis([-400 50]);
    
    drawnow
    
    pause(0.1)
    
end

%% 


load([dataDir 'lwdown.mat']);

for i = 1:12
    
    subplot(121);
    
    pcolorcen(lwdown.lon,lwdown.lat,squeeze(lwdown.meanField(:,:,i)));
    
    colorbar;
    
    caxis([-400 50]);
    
    subplot(122);
    
    pcolorcen(xx,yy,squeeze(lwdown.binMeanField(:,:,i)));
    
    colorbar;
    
    caxis([-400 50]);
    
    drawnow
    
    pause(0.1)
    
end

%% 

load([dataDir 'apco2.mat']);

for i = 1:12
    
    subplot(121);
    
    pcolorcen(apco2.lon,apco2.lat,squeeze(apco2.meanField(:,:,i)));
    
    colorbar;
    
    subplot(122);
    
    pcolorcen(xx,yy,squeeze(apco2.binMeanField(:,:,i)));
    
    colorbar;
    
    drawnow
    
    pause(0.1)
    
end

%% 


load([dataDir 'runoff.mat']);

for i = 1:12
    
    subplot(121);
    
    quikplot_llc(log10(squeeze(runoff.meanField(:,:,i))));
    
    colorbar;
    
    caxis([-50 0]);
    
    subplot(122);
    
    pcolorcen(xx,yy,log10(squeeze(runoff.binMeanField(:,:,i))));
    
    colorbar;
    
    caxis([-50 0]);
    
    drawnow
    
    pause(0.1)
    
end

%% 


load([dataDir 'irondust.mat']);

for i = 1:12
    
    subplot(121);
    
    quikplot_llc(log10(squeeze(irondust.meanField(:,:,i))));
    
    colorbar;
    
    %caxis([-400 50]);
    
    subplot(122);
    
    pcolorcen(xx,yy,log10(squeeze(irondust.binMeanField(:,:,i))));
    
    colorbar;
    
    %caxis([-400 50]);
    
    drawnow
    
    pause(0.1)
    
end

