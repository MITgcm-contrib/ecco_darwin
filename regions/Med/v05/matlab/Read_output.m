% Dynamics of Limiting nutriets in Mediterranean locations


load('cm_dblue_white_orange2.mat')

%A=rdmds('daily_test1/SST.0000000576');

%A=rdmds('daily_test1/pickup.0000002880');

%   [NO3,ITS]=rdmds('daily/NO3',NaN);
%   [PO4,ITS]=rdmds('daily/PO4',NaN);
% [FeT,ITS]=rdmds('daily/FeT',NaN);
% [c1 ,ITS]=rdmds('daily/c1',NaN);
% [c2 ,ITS]=rdmds('daily/c2',NaN);
% [c3 ,ITS]=rdmds('daily/c3',NaN);
% [c4 ,ITS]=rdmds('daily/c4',NaN);
% [c5 ,ITS]=rdmds('daily/c5',NaN);
% [c6 ,ITS]=rdmds('daily/c6',NaN);
% [c7 ,IdTS]=rdmds('daily/c7',NaN);

%[SST ,ITS]=rdmds('daily_test1/SST',NaN);

load('Mediterranean_grid','Med_grid') %Made by SST file (when SST=0-->land)

loc=[35 70 95; 45 20 15];

%plot_time=500-360:500;
plot_time=1:length(NO3);

figure(1)
subplot(2,2,1)
pcolor(1:144,1:72,Med_grid')
colormap(cm)
shading flat
%colorbar
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
title('Mediterranean grid')
%set(gca,'Color',grey)
hold on
%plot(loc(1,:),loc(2,:),'ko','LineWidth',2)

labelsM = {'L1','L2','L3'};
pointAm=loc(1,:);
pointBm=loc(2,:);
for ma=1:length(pointAm)
    plot(pointAm(ma),pointBm(ma),'k.','MarkerSize',15,'LineWidth',2)
    text(pointAm(ma),pointBm(ma),labelsM(ma),'VerticalAlignment','bottom','HorizontalAlignment','left','Color','k','FontSize',12,'FontWeight','bold')
    
end
hold off

for i=1:3
    subplot(2,2,i+1)
    plot(plot_time,squeeze(NO3(loc(1,i),loc(2,i),5,plot_time)),'LineWidth',1.5)
    hold on
    plot(plot_time,squeeze(PO4(loc(1,i),loc(2,i),5,plot_time)*16),'LineWidth',1.5)
    %plot(plot_time,squeeze(FeT(loc(1,i),loc(2,i),5,plot_time)*1000),'LineWidth',1)
    hold off
    xlim([plot_time(1), plot_time(end)])
    xlabel('days')
    ylabel('μmol L^-^1')
    
    if i==1
        title('L1')
        legend('NO3','PO4 *16','Location','Best')
    elseif i==2
        title('L2')
    elseif i==3
        title('L3')
    end
end


