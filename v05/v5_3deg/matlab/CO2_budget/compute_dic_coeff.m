function [dDIC_dSST, dDIC_dSSS, dDIC_dSSAlk, dDIC_dPCO2] = compute_dic_coeff(meanSST,meanSSSl, ...
    meanSSAlk,meanSSPO4,meanSSSiO2,meanAPCO2)

rhoConst = 1029;

thetaVec = [-2:40]';
saltVec = [0:40]';
alkVec = [0:1:3500]';
aPCO2Vec = [1:1:1000]';

%%

hFig1 = figure(1);
set(hFig1,'units','normalized','outerposition',[0 0 1 1]);
set(gcf,'color',[1 1 1]);

fs = 26;
fsTitle = 20;
lw = 2;
mw = 100;

%compute temperature dependence

A = CO2SYS(meanAPCO2,meanSSAlk,4,1,meanSSSl,thetaVec,thetaVec,...
    0,0,meanSSSiO2,meanSSPO4,1,4,1);

tcsat = A(:,2); %umol kg^-1

[p S] = polyfit(thetaVec,tcsat,1);
f = polyval(p,thetaVec);

subplot(221)

hold on

scatter(thetaVec,tcsat,mw,'MarkerFaceColor','b','MarkerEdgeColor','none')
plot([thetaVec(1) thetaVec(end)],[(p(2)+p(1)*thetaVec(1)) (p(2)+p(1)*thetaVec(end))],'r--','LineWidth',lw)

ylabel('Csat (umol kg^-^1)','FontSize',fs)
xlabel('SST (\circC)','FontSize',fs)

axis tight

box on
grid on
set(gca,'gridlinestyle','--');
set(gca,'FontSize',fs);
set(gca,'LineWidth',lw);

title(['Mean SST = ' num2str(meanSST) ', dDIC/dSST = ',num2str(p(1)),' umol kg^-^1 K^-^1'], ...
    'FontSize',fsTitle,'FontWeight','Normal')

dDIC_dSST = p(1) .* rhoConst .* 10^-6; %umol kg^-1 psu^-1 to mol m^-3 psu^-1

%%
%calculate salinity dependence

A = CO2SYS(meanAPCO2,meanSSAlk,4,1,saltVec,meanSST,meanSST,...
    0,0,meanSSSiO2,meanSSPO4,1,4,1);

scsat = A(:,2);

[p S] = polyfit(saltVec,scsat,1);
f = polyval(p,saltVec);

subplot(222)

hold on

scatter(saltVec,scsat,mw,'MarkerFaceColor','b','MarkerEdgeColor','none')
plot([saltVec(1) saltVec(end)],[(p(2)+p(1)*saltVec(1)) (p(2)+p(1)*saltVec(end))],'r--','LineWidth',lw)

ylabel('DIC (\mumol kg^-^1)','FontSize',fs)
xlabel('SSS','FontSize',fs)

axis tight

box on
grid on
set(gca,'gridlinestyle','--');
set(gca,'FontSize',fs);
set(gca,'LineWidth',lw);

title(['Mean SSS = ' num2str(meanSSSl) ', dDIC/dSSS = ',num2str(p(1)),' \mumol kg^-^1 psu^-^1'], ...
    'FontSize',fsTitle,'FontWeight','Normal')

dDIC_dSSS = p(1) .* rhoConst .* 10^-6; %umol kg^-1 psu^-1 to mol m^-3 psu^-1

%%
%calculate ALK dependence

A = CO2SYS(meanAPCO2,alkVec,4,1,meanSSSl,meanSST,meanSST,...
    0,0,meanSSSiO2,meanSSPO4,1,4,1);

acsat = A(:,2); % umol kg^-^1

[p S] = polyfit(alkVec,acsat,1);
f = polyval(p,alkVec);

subplot(223)

hold on

scatter(alkVec,acsat,mw,'MarkerFaceColor','b','MarkerEdgeColor','none')
plot([alkVec(1) alkVec(end)],[(p(2)+p(1)*alkVec(1)) (p(2)+p(1)*alkVec(end))],'r--','LineWidth',lw)

ylabel('DIC (\mumol kg^-^1)','FontSize',fs);
xlabel('SSAlk (\mumol kg^-^1)','FontSize',fs);

axis tight

box on
grid on
set(gca,'gridlinestyle','--');
set(gca,'FontSize',fs);
set(gca,'LineWidth',lw);

title(['Mean SSAlk = ' num2str(meanSSAlk) ', dDIC/dSSAlk = ' num2str(p(1))], ...
    'FontSize',fsTitle,'FontWeight','Normal');

dDIC_dSSAlk = p(1); %umol/kg/(umol/kg)

%%
%calculate pCO2 dependence

A = CO2SYS(aPCO2Vec,meanSSAlk,4,1,meanSSSl,meanSST,meanSST,...
    0,0,meanSSSiO2,meanSSPO4,1,4,1);

pco2csat = A(:,2); % umol kg^-1
revelle = A(:,29);

dDIC_dPCO2 = (ones(size(aPCO2Vec)) .* pco2csat) ./ (revelle .* aPCO2Vec);

[p S] = polyfit(aPCO2Vec(15:end),pco2csat(15:end),1);
f = polyval(p,aPCO2Vec(15:end));

dDIC_dPCO22 = dDIC_dPCO2 .* rhoConst .* 10^-6;  %umol kg^-1 uatm^-1 to mol m^-3 uatm^-1
dDIC_dPCO22 = [aPCO2Vec, dDIC_dPCO22];

xi = find(dDIC_dPCO22(:,1) >= meanAPCO2,1);

subplot(224)

hold on

scatter(aPCO2Vec,pco2csat,mw,'MarkerFaceColor','b','MarkerEdgeColor','none')
plot(aPCO2Vec(2:end),(pco2csat(1:end-1)+dDIC_dPCO2(2:end).*unique(diff(aPCO2Vec))),'r--','LineWidth',lw)
scatter(aPCO2Vec(xi),pco2csat(xi),mw*2,'MarkerFaceColor','g','MarkerEdgeColor','none');

ylabel('DIC (\mumol kg^-^1)','FontSize',fs);
xlabel('pCO2 (\muAtm)','FontSize',fs);

xlim([0 3000]);

axis tight

box on
grid on
set(gca,'gridlinestyle','--');
set(gca,'FontSize',fs);
set(gca,'LineWidth',lw);

title(['Mean apCO2 = ' num2str(meanAPCO2) ', dDIC/dPCO2 = ' num2str(dDIC_dPCO22(xi,2)) '(\mumol kg^-^1 \muAtm^-^1)'], ...
    'FontSize',fsTitle,'FontWeight','Normal')

dDIC_dPCO2 = dDIC_dPCO22(xi,2);

end