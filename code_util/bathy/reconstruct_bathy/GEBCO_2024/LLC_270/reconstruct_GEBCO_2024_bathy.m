clear
close all;

tic

saveBathy = 1;
maskDryCells = 0;

useShallowDepthCrit = 0;
useDeepDepthCrit = 0;

gridDir = '/Users/carrolld/Documents/research/bathy/bin/grid/LLC_270/';

dataDir1 = '/Users/carrolld/Documents/research/bathy/mat/indices/GEBCO_2024/';
dataDir2 = '/Users/carrolld/Documents/research/bathy/grid/LLC_270/';
saveDir = '/Users/carrolld/Documents/research/bathy/mat/bathy/LLC_270/GEBCO_2024/';

%%

if maskDryCells
    
    suffix = 'wet_dustin';
    
else
    
    suffix = 'all_dustin';
    
end

%%

numFacets = 5;
numFaces = 13;

nx = 4320;
ny = nx .* numFaces;

%%

for i = 1:numFacets
    
    eval(['load([dataDir1 ''GEBCO_2024_LLC_270_indices_facet_' num2str(i) '_' suffix '.mat'']);']);
    
    field =  bathy.medianDepth;
    
    field(isnan(field)) = 0;
    
    if useShallowDepthCrit
        
        field(field <= 5) = 0;
        %field(field > -5 & field < 0) = -5; %mackenzie setup
        
        saveSuffix = [suffix '_5m_crit'];
        
    elseif useDeepDepthCrit
        
        field(field >= 1 & field <= 10) = 10;
        
        saveSuffix = [suffix '_10m_crit'];
        
    else
        
        saveSuffix = [suffix '_method'];
        
    end
    
    facet{i}.bathy = field;
    
    clear bathy field
    
    disp(num2str(i));
    
end

toc

%%

nx = 270;
ny = 3510;

depth = zeros(nx,ny);

depth(1:nx*nx*3) = facet{1}.bathy;
depth(nx*nx*3+1:nx*nx*6) = facet{2}.bathy;
depth(nx*nx*6+1:nx*nx*7) = facet{3}.bathy;
depth(nx*nx*7+1:nx*nx*10) = facet{4}.bathy;
depth(nx*nx*10+1:nx*nx*13) = facet{5}.bathy;

%%

%b = depth;
%b2=1+0*b;
%b2(find(b))=0;
%b3=imfill(b2,'holes');
%bf=b;
%bf(find(b3))=0;
%depth = bf;

%%

hFig1 = figure(1);
set(hFig1,'units','normalized','outerposition',[0 0 1 1]);
set(gcf,'color',[1 1 1]);

quikplot_llc(depth);

axis tight

%%

if saveBathy
    
    writebin([saveDir  'LLC_270_bathy_' saveSuffix '.bin'],depth,1,'real*4');
    save([saveDir  'LLC_270_bathy_' saveSuffix '.mat'],'depth','-v7.3');
    
    cd(saveDir);
    
    close all
    
    testDepth = readbin([saveDir 'LLC_270_bathy_' saveSuffix '.bin'],[nx ny],1,'real*4');
    testDepth(testDepth == 0) = nan;
    
    hFig1 = figure(1);
    set(hFig1,'units','normalized','outerposition',[0 0 1 1]);
    set(gcf,'color',[1 1 1]);
    
    hold on
    
    set(gca,'Color',[0.65 0.65 0.65]);
    
    quikplot_llc(testDepth);
    
    caxis([0 5000]);
    
end

%%
