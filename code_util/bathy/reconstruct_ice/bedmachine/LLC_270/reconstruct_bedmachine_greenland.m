clear
close all;

tic

saveBathy = 1;

gridDir = '/Users/carrolld/Documents/research/carbon/grid/LLC_270/';

dataDir1 = '/Users/carrolld/Documents/research/bathy/mat/ice/bedmachine/greenland/LLC_270/';
dataDir2 = '/Users/carrolld/Documents/research/bathy/grid/LLC_270/';
saveDir = '/Users/carrolld/Documents/research/bathy/mat/ice/bedmachine/greenland/LLC_270/';

%%

numFacets = 5;
numFaces = 13;

nx = 270;
ny = nx .* numFaces;

%%

for i = 1:numFacets
    
    eval(['load([dataDir1 ''bedmachine_greenland_LLC_270_indices_facet_' num2str(i) '_all_dustin.mat'']);']);
    
    field1 =  bathy.medianDepth;
    field2 =  bathy.medianMask-1;

    field1(field1 == -1) = nan;
    field2(field2 == -1) = nan;
    
    facet{i}.ice = field1;
    facet{i}.mask = field2;
    
    clear bathy field1 field2
    
    disp(num2str(i));
    
end

toc

%%

ice = zeros(nx,ny);

ice(1:nx*nx*3) = facet{1}.ice;
ice(nx*nx*3+1:nx*nx*6) = facet{2}.ice;
ice(nx*nx*6+1:nx*nx*7) = facet{3}.ice;
ice(nx*nx*7+1:nx*nx*10) = facet{4}.ice;
ice(nx*nx*10+1:nx*nx*13) = facet{5}.ice;

mask = zeros(nx,ny);

mask(1:nx*nx*3) = facet{1}.mask;
mask(nx*nx*3+1:nx*nx*6) = facet{2}.mask;
mask(nx*nx*6+1:nx*nx*7) = facet{3}.mask;
mask(nx*nx*7+1:nx*nx*10) = facet{4}.mask;
mask(nx*nx*10+1:nx*nx*13) = facet{5}.mask;

%%

% b = depth;
% b2=1+0*b;
% b2(find(b))=0;
% b3=imfill(b2,'holes');
% bf=b;
% bf(find(b3))=0;

%depth = bf;

%%

if saveBathy
    
    save([saveDir  'LLC_270_bedmachine_greenland_ice_mask_dustin_method.mat'],'ice','mask','-v7.3');
        
end

%%

quikplot_llc(mask);

colorbar

%%

