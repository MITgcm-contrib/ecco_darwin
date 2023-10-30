% Compare ECCO-Darwin model output to Gonzalez et al. (2016)
% https://doi.org/10.1016/j.jmarsys.2016.03.007

% Replace following line with location of ECCO_Darwin/Med directory
cd ~/Links/Public/ECCO_Darwin/Med

% model grid dimensions
nx=144;
ny=72;
nz=47;

% file name suffices
s1=['_' int2str(nz)];
s2=['_' int2str(nx) 'x' int2str(ny)];
s3=[s2 'x' int2str(nz)];

% Find model indices near Saronikos Seawatch buoy (37.61°N 23.56°E)
XC=readbin(['grid/XC' s2],[nx ny]);
YC=readbin(['grid/YC' s2],[nx ny]);
i=closest(23.56,XC(:,1));
j=closest(37.61,YC(1,:));

% Plot Fig. 1, model depth and buoy location
Depth=readbin(['grid/Depth' s2],[nx ny]);
clf
pcolorcen(XC',YC',log10(Depth)');
axis([13 28 33 41])
caxis([0 3])
hold on
plot(XC(i,1),YC(1,j),'k*','linewidth',2)
title('Location of Saronikos Seawatch buoy (37.61°N 23.56°E)')
xlabel('Longitude East (degrees)')
ylabel('Latitude North (degrees)')
h=colorbar('h');
set(h,'Ticks',[0 1 2 3],'TickLabels',{'0 m','10 m','100 m','1 km'})
print -djpeg figures/Gonzalez2016_fig1

% Plot Fig. 3, SST and SSS from 9/2013 to 11/2014
dte=datenum(2013,10:23,1);
SST=0*dte;
SSS=0*dte;
for m=1:length(dte)
    tmp=readbin(['THETA/THETA' s3 '.' datestr(dte(m),30)],[nx ny]);
    SST(m)=tmp(i,j);
    tmp=readbin(['SALT/SALT' s3 '.' datestr(dte(m),30)],[nx ny]);
    SSS(m)=tmp(i,j);
end
clf
yyaxis right
plot(dte-15,SST,'.-')
xlim([datenum(2013,9,1) datenum(2014,11,1)])
ylim([12 30])
ylabel('SST (deg C)')
datetick('x',12,'keeplimits')
yyaxis left
plot(dte-15,SSS,'.-')
ylim([37.5 39.5])
ylabel('SSS (g/kg)')
grid
title('SSS and SST at Saronikos Seawatch buoy')
print -djpeg figures/Gonzalez2016_fig3

% Plot Fig. 5, surface pH and CO2 fugacity from 9/2013 to 11/2014
pH=0*dte;
fCO2=0*dte;
for m=1:length(dte)
    tmp=readbin(['pH/pH' s3 '.' datestr(dte(m),30)],[nx ny]);
    pH(m)=tmp(i,j);
    tmp=readbin(['fugCO2/fugCO2' s2 '.' datestr(dte(m),30)],[nx ny]);
    fCO2(m)=tmp(i,j);
end
clf
colororder({'g','b'})
yyaxis right
plot(dte-15,pH,'.-')
xlim([datenum(2013,9,1) datenum(2014,11,1)])
ylim([7.95 8.2])
ylabel('Acidity (Ph)')
datetick('x',12,'keeplimits')
yyaxis left
plot(dte-15,fCO2,'.-')
ylabel('Fugacity of CO2 (atm)')
grid
title('CO2 fugacity and acidity')
print -djpeg figures/Gonzalez2016_fig5

% Plot Fig. 7, air-sea CO2 flux from 9/2013 to 11/2014
CO2flux=0*dte;
for m=1:length(dte)
    tmp=readbin(['CO2_flux/CO2_flux' s2 '.' datestr(dte(m),30)],[nx ny]);
    CO2flux(m)=tmp(i,j);
end
clf
plot(dte-15,-CO2flux*60*60*24,'k.-')
xlim([datenum(2013,9,1) datenum(2014,11,1)])
ylim([-30 30])
ylabel('Flux of CO_2 (mmol C m^{-2} d^{-1})')
datetick('x',12,'keeplimits')
grid
title('Air-sea CO_2 flux')
print -djpeg figures/Gonzalez2016_fig7
