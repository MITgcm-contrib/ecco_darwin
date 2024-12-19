clear
close all;

savePlot = 1;
plotBathy = 1;
saveIce = 1;

gridDir = '/Users/carrolld/Documents/research/carbon/grid/LLC_270/';

dataDir = '/Users/carrolld/Documents/research/bathy/mat/bathy/LLC_270/Schaffer_2019/';
figureDir = '/Users/carrolld/Documents/research/bathy/figures/bathy/LLC_270/Schaffer_2019/';
saveDir = '/Users/carrolld/Documents/research/bathy/mat/bathy/LLC_270/Schaffer_2019/';

%%

numFacets = 5;
numFaces = 13;

nx = 270;
ny = nx .* numFaces;

hFacC = readbin([gridDir 'hFacC.data'],[nx ny],1,'real*4');

global mygrid

mygrid = [];

grid_load(gridDir,5,'compact');

%%

B1 = load([dataDir 'LLC_270_bathy_all_dustin_method.mat']);
B2 = load([dataDir 'LLC_270_ice_all_dustin_method.mat']);

%%

B1.depth = -B1.depth;

depth = B1.depth + B2.depth;

B1.depth(B1.depth == 0) = nan;
B2.depth(B2.depth == 0) = nan;

%depth(depth == 0) = nan;

tempDepth = convert2gcmfaces(depth);

for i = 1:5

    eval(['temp = tempDepth.f' num2str(i) ';']);

    % remove lakes
    b2 = 1+0.*temp;
    b2(find(temp))=0;
    b3 = imfill(b2,'holes');
    bf = temp;
    bf(find(b3)) = 0;

    bf(bf == 0) = nan;

    %pcolorcen(bf)
    %drawnow
    %pause

    eval(['tempDepth.f' num2str(i) ' = bf;']);

end

depth = convert2gcmfaces(tempDepth);

%%

if plotBathy

    hFig1 = figure(1);
    set(hFig1,'units','normalized','outerposition',[0 0 1 1]);
    set(gcf,'color',[1 1 1]);

    colors1 = cmocean('deep',1000);
    colors2 = flipud(cmocean('ice',1000));

    lw = 2;
    fs = 26;

    cc1 = subplot(131);

    hold on

    set(gca,'Color',[0.65 0.65 0.65]);

    quikplot_llc(B1.depth);

    colormap(cc1,colors1);

    caxis([-5000 0]);

    hcb1 = colorbar;
    set(get(hcb1,'ylabel'),'String','Depth (m)');

    axis tight

    xlabel('nx');
    ylabel('ny');

    grid off
    box on
    set(gca,'GridLineStyle','--');
    set(gca,'LineWidth',lw);
    set(gca,'FontSize',fs);

    title('Schaffer et al. (2016) Bathy');

    cc2 = subplot(132);

    hold on

    set(gca,'Color',[0.65 0.65 0.65]);

    quikplot_llc(B2.depth);

    caxis([0 500]);

    colormap(cc2,colors2);

    hcb2 = colorbar;
    set(get(hcb2,'ylabel'),'String','Ice Base (m)');

    axis tight

    xlabel('nx');
    ylabel('ny');

    grid off
    box on
    set(gca,'GridLineStyle','--');
    set(gca,'LineWidth',lw);
    set(gca,'FontSize',fs);

    title('Schaffer et al. (2016) Ice Base');

    cc3 = subplot(133);

    hold on

    set(gca,'Color',[0.65 0.65 0.65]);

    quikplot_llc(depth);

    caxis([-5000 0]);

    colormap(cc3,colors1);

    hcb3 = colorbar;
    set(get(hcb3,'ylabel'),'String','Depth (m)');

    axis tight

    xlabel('nx');
    ylabel('ny');

    grid off
    box on
    set(gca,'GridLineStyle','--');
    set(gca,'LineWidth',lw);
    set(gca,'FontSize',fs);

    title('Schaffer et al. (2016) Bathy + Ice Base');

    drawnow

    if savePlot

        export_fig(hFig1,[figureDir 'Schaffer_2016_merged.png'],'-png','-r150');

    end

end

%%

if saveIce

    save([saveDir  'LLC_270_Schaffer_2019_merged_bathy.mat'],'depth','-v7.3');

end

%%
