% Compare ECCO-Darwin model output to Solidoro et al. (2022)
% https://doi.org/10.3389/fmars.2021.781522

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

% Plot Fig. 1, model bathymetry
XC=readbin(['grid/XC' s2],[nx ny]);
YC=readbin(['grid/YC' s2],[nx ny]);
Depth=readbin(['grid/Depth' s2],[nx ny]);
clf
pcolorcen(XC',YC',log10(Depth)');
caxis([0 3.6])
title('Model bathymetry')
xlabel('Longitude East (degrees)')
ylabel('Latitude North (degrees)')
h=colorbar('h');
set(h,'Ticks',[0 1 2 3],'TickLabels',{'0 m','10 m','100 m','1 km'})
print -djpeg figures/Solidoro2022_fig1

% Plot Fig.2, time series of apCO2, SST, and surface DIC
hFacC=readbin(['grid/hFacC' s3],[nx ny]);
hFacC(1:16,54:end)=0;
hFacC(100:end,44:end)=0;
dte=datenum(1992,2:379,1);
apCO2=0*dte;
SST=0*dte;
DIC=0*dte;
for m=1:length(dte)
    tmp=readbin(['apCO2/apCO2' s2 '.' datestr(dte(m),30)],[nx ny]);
    apCO2(m)=sum(hFacC(:).*tmp(:))/sum(hFacC(:));
    tmp=readbin(['THETA/THETA' s3 '.' datestr(dte(m),30)],[nx ny]);
    SST(m)=sum(hFacC(:).*tmp(:))/sum(hFacC(:));
    tmp=readbin(['DIC/DIC' s3 '.' datestr(dte(m),30)],[nx ny]);
    DIC(m)=sum(hFacC(:).*tmp(:))/sum(hFacC(:));
end
clf
subplot(211)
plot(dte-15,apCO2*1e6,'k.-')
xlim([datenum(1992,1,1) datenum(2023,7,1)])
ylim([349 424])
grid
datetick('x',10,'keeplimits')
ylabel('apCO_2 (\muatm)')
title('Mediterranean basin-averaged atmospheric CO2 concentration')
subplot(212)
yyaxis right
plot(dte-15,SST,'.-')
hold on
plot(dte-15,movmean(SST,24),'-','linewidth',2)
xlim([datenum(1992,1,1) datenum(2023,7,1)])
ylim([14 27])
ylabel('SST (^\circC)')
datetick('x',10,'keeplimits')
yyaxis left
plot(dte-15,DIC,'.-')
hold on
plot(dte-15,movmean(DIC,24),'-','linewidth',2)
ylim([2210 2390])
ylabel('DIC (mmol C m^{-3})')
grid
title('Mediterranean basin-averaged SST and DIC')
print -djpeg figures/Solidoro2022_fig2

% Plot Fig.4a, East and West Basin
hFacC=readbin(['grid/hFacC' s3],[nx ny length(iz)]);
hFacC(1:16,54:end,:)=0;
hFacC(100:end,44:end,:)=0;
WestBasin=hFacC;
WestBasin(56:end,47:end,:)=0;
WestBasin(67:end,36:end,:)=0;
WestBasin(68:end,:,:)=0;
WestBasin(56:end,1:31,:)=0;
WestBasin(49:end,1:28,:)=0;
WestBasin(52:end,1:29,:)=0;
WestBasin(54:end,1:30,:)=0;
EastBasin=hFacC;
EastBasin(find(WestBasin))=0;
clf
subplot(211)
pcolorcen(XC',YC',WestBasin(:,:,1)'+hFacC(:,:,1)');
colormap(cmap)
xlabel('Longitude East (\circ)')
ylabel('Latitude North (\circ)')
title('West Basin (white)')
subplot(212)
pcolorcen(XC',YC',EastBasin(:,:,1)'+hFacC(:,:,1)');
colormap(cmap)
xlabel('Longitude East (\circ)')
ylabel('Latitude North (\circ)')
title('East Basin (white)')
print -djpeg figures/Solidoro2022_fig4a

% Plot Fig.4b, time-mean 2000-2010 surface Chlorophyll-a
RC=readbin(['grid/RC' s1],nz);
iz=1:19;
dte=datenum(2000,2:133,1);
Chl=zeros(nx,ny,length(iz),12);
for m=1:length(dte)
    tmp=zeros(nx,ny,length(iz));
    for c=1:5
        fnm=['Chl' int2str(c)];
        fnm=[fnm '/' fnm s3 '.' datestr(dte(m),30)];
        tmp=tmp+readbin(fnm,[nx ny length(iz)]);
    end
    Chl(:,:,:,mod(m-1,12)+1)=Chl(:,:,:,mod(m-1,12)+1)+tmp;
end
Chl=Chl/11;
clf
pcolorcen(XC',YC',Chl(:,:,1)');
caxis([0 .4])
colorbar('h')
colormap(cmap)
xlabel('Longitude East (\circ)')
ylabel('Latitude North (\circ)')
title('Mean 2000-2010 surface Chlorophyll-a concentration (mg Chl-a m^{-3})')
print -djpeg figures/Solidoro2022_fig4b

% Plot Fig.4c, 2000-2010 Chl-a seasonal cycle in West and East Basins
Chl_West=zeros(12,length(iz));
Chl_East=zeros(12,length(iz));
for mo=1:12
    for k=1:length(iz)
        tmp1=WestBasin(:,:,k).*Chl(:,:,k,mo);
        tmp2=WestBasin(:,:,k);
        Chl_West(mo,k)=sum(tmp1)/sum(tmp2);
        tmp1=EastBasin(:,:,k).*Chl(:,:,k,mo);
        tmp2=EastBasin(:,:,k);
        Chl_East(mo,k)=sum(tmp1)/sum(tmp2);
    end
end
clf
subplot(211)
pcolorcen(1:12,RC(iz),Chl_West');
axis([.5 12.5 -250 0]) 
caxis([0 .45])
colorbar
xlabel('month')
ylabel('depth (m)')
title('2000-2010 Chl-a seasonal cycle in West Med (mg Chl-a m^{-3})')
subplot(212)
pcolorcen(1:12,RC(iz),Chl_East');
axis([.5 12.5 -250 0]) 
caxis([0 .26])
colorbar
xlabel('month')
ylabel('depth (m)')
title('2000-2010 Chl-a seasonal cycle in East Med (mg Chl-a m^{-3})')
print -djpeg figures/Solidoro2022_fig4c
