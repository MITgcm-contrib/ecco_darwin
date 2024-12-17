clear
close all;

savePlot = 1;
plotIce = 1;
saveIce = 1;

gridDir = '/Users/carrolld/Documents/research/carbon/grid/LLC_270/';

dataDir1 = '/Users/carrolld/Documents/research/bathy/mat/ice/bedmachine/greenland/LLC_270/';
dataDir2 = '/Users/carrolld/Documents/research/bathy/mat/ice/bedmachine/antarctica/LLC_270/';
figureDir = '/Users/carrolld/Documents/research/bathy/figures/ice/';
saveDir = '/Users/carrolld/Documents/research/bathy/mat/ice/bedmachine/merged/LLC_270/';

%%

numFacets = 5;
numFaces = 13;

nx = 270;
ny = nx .* numFaces;

hFacC = readbin([gridDir 'hFacC.data'],[nx ny],1,'real*4');

%%

B1 = load([dataDir1 'LLC_270_bedmachine_greenland_ice_mask_dustin_method.mat']);
B2 = load([dataDir2 'LLC_270_bedmachine_antarctica_ice_mask_dustin_method.mat']);

%%

ice = cat(3,B1.ice,B2.ice);
ice = nansum(ice,3);

mask = cat(3,B1.mask,B2.mask);
mask = nansum(mask,3);

%%

clear B1 B2

%%

if plotIce

    hFig1 = figure(1);
    set(hFig1,'units','normalized','outerposition',[0 0 1 1]);
    set(gcf,'color',[1 1 1]);

    colors1 = cmocean('ice',1000);
    colors2 = cmocean('phase',7);

    colors2 = colors2(1:5,:);

    lw = 2;
    fs = 22;

    cc1 = subplot(121);

    hold on

    set(gca,'Color',[0.65 0.65 0.65]);

    quikplot_llc(ice);

    colormap(cc1,colors1);

    caxis([-1000 0]);

    hcb1 = colorbar;
    set(get(hcb1,'ylabel'),'String','Ice Base (m)');

    axis tight

    xlabel('nx');
    ylabel('ny');

    grid off
    box on
    set(gca,'GridLineStyle','--');
    set(gca,'LineWidth',lw);
    set(gca,'FontSize',fs);

    title('BedMachine Ice Base Depth');

    cc2 = subplot(122);

    hold on

    set(gca,'Color',[0.65 0.65 0.65]);

    quikplot_llc(mask);

    caxis([0 5]);

    colormap(cc2,colors2);

    colorbar

    axis tight

    xlabel('nx');
    ylabel('ny');

    grid off
    box on
    set(gca,'GridLineStyle','--');
    set(gca,'LineWidth',lw);
    set(gca,'FontSize',fs);

    title('Mask (0 = Ocean, 1 = Ice-free Land, 2 = Grounded Ice, 3 = Floating Ice, 4 = Non-Greenland Land/Lake Vostok)');
    drawnow

    if savePlot

        export_fig(hFig1,[figureDir 'BedMachine_Greenland_Antarctica_merged.png'],'-png','-r150');

    end

end

%%

if saveIce

    save([saveDir  'LLC_270_bedmachine_merged_ice_mask.mat'],'ice','mask','-v7.3');

end

%%
