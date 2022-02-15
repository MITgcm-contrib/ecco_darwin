% Snap GlobalNEWS2 2000 time-mean Qact runoff to jra55_do grid

%%%%%%%%%%%%%%%%%%%%% on pleiades
clear, close all, cd ~dmenemen/forcing/jra55_do

% Compute jra55_do time-mean year-2000 runoff
yr=2000; loy=365; if mod(yr,4)==0, loy=366; end
lon=0.125:0.25:360; lat=-89.875:0.25:90;
nx=length(lon); ny=length(lat); [jlat jlon]=meshgrid(lat,lon);
cellarea=readbin('cellarea.bin',[nx ny]);  % area of grid cell in m^2
pin='/nobackup/dcarrol2/LOAC/bin/jra55_do/v1.4.0/';
fin=[pin 'jra55_do_runoff_' int2str(yr)];
jra=readbin(fin,[nx ny loy]);              % load whole year, units are m/s
jra=sum(jra,3).*cellarea*24*60*60/1e9;     % cumulate & convert to km^3/yr
jlon=jlon(:); jlat=jlat(:); jra=jra(:);
ix=find(jra==0); jlon(ix)=0; jlat(ix)=0;
jra=sparse(jra); jlon=sparse(jlon); jlat=sparse(jlat);
save jra55_2000 j*

%%%%%%%%%%%%%%%%%%%%% on local workstation
clear, close all
cd ~dmenemen/projects/LOAC/runoff_products/GlobalNEWS
load globalnews

% Load GlobalNEWS2 mouth_lon, mouth_lat, Qact (km3/yr)
gns=xlsread('globalnews');
glon=gns(:,1); glat=gns(:,2); gns=gns(:,3);
ix=find(glon<0); glon(ix)=glon(ix)+360;
clear c* f* i* l* n* p* y*
save globalnews g*

% Plot JRA55 and GlobalNEWS runoff
figure(1), clf, plotland(.8*[1 1 1],12), hold on
for i=1:75
    ix=find(jra>(i-1)*100&jra<=i*100);
    if ~isempty(ix)
        plot(jlon(ix),jlat(ix),'ro','markersize',i+3)
    end
    ix=find(gns>(i-1)*100&gns<=i*100);
    if ~isempty(ix)
        plot(glon(ix),glat(ix),'bo','markersize',i+3)
    end
end
title('All runoff locations')
text(185,-50,['jra ' int2str(sum(jra)) ' km^3/yr'],'color','r')
text(185,-60,['gns ' int2str(sum(gns)) ' km^3/yr'],'color','b')
print -dpdf AllRunoff

% Remove jra55_do Antarctic locations
ix=find(jlat<-60); jra(ix)=0;
figure(2), clf, plotland(.8*[1 1 1],12), hold on
for i=1:75
    ix=find(jra>(i-1)*100&jra<=i*100);
    if ~isempty(ix)
        plot(jlon(ix),jlat(ix),'ro','markersize',i+3)
    end
    ix=find(gns>(i-1)*100&gns<=i*100);
    if ~isempty(ix)
        plot(glon(ix),glat(ix),'bo','markersize',i+3)
    end
end
title('No Antarctica')
text(185,-50,['jra ' int2str(sum(jra)) ' km^3/yr'],'color','r')
text(185,-60,['gns ' int2str(sum(gns)) ' km^3/yr'],'color','b')
print -dpdf NoAntacrctic

% Remove landlocked GlobalNEWS locations
gIX=0*glon; gD=gIX;
for i=1:length(gIX);
    d=sqrt((glon(i)-jlon).^2+(glat(i)-jlat).^2);
    [Y,I]=min(d); gIX(i)=I; gD(i)=sqrt(Y);
end
ix=find(gD>1.6); gns(ix)=0;
figure(3), clf, plotland(.8*[1 1 1],12), hold on
for i=1:75
    ix=find(jra>(i-1)*100&jra<=i*100);
    if ~isempty(ix)
        plot(jlon(ix),jlat(ix),'ro','markersize',i+3)
    end
    ix=find(gns>(i-1)*100&gns<=i*100);
    if ~isempty(ix)
        plot(glon(ix),glat(ix),'bo','markersize',i+3)
    end
end
clear I Y d gD gIX i*
title('No landlocked locations')
text(185,-50,['jra ' int2str(sum(jra)) ' km^3/yr'],'color','r')
text(185,-60,['gns ' int2str(sum(gns)) ' km^3/yr'],'color','b')
print -dpdf NoLandlocked

% Associate each GlobalNEWS location with jra55_do runoff
maxsep=2.5;                                % maximum separation between GlobalNEWS and JRA55 in degrees
gns2jra=sparse(length(jra),length(gns)); % GlobalNEWS and JRA55 tranformation matrix
gns_r=gns; jra_r=jra;                    % residual non-associated runoff for GlobalNEWS and JRA55
[tmp gx]=sort(gns,'descend');
for g=1:length(gx)
    jx=find(jra_r);
    d=sqrt((glon(gx(g))-jlon(jx)).^2+(glat(gx(g))-jlat(jx)).^2);
    while gns_r(gx(g))>0 & min(d)<=maxsep
        [tmp j]=min(d);
        if tmp>maxsep
            break
        else
            tmp=min([gns_r(gx(g)) jra_r(jx(j))]);
            gns2jra(gx(g),jx(j))=tmp/gns_r(gx(g));
            gns_r(gx(g))=gns_r(gx(g))-tmp;
            jra_r(jx(j))=jra_r(jx(j))-tmp;
            d(j)=100;
end,end,end

% Print residuals
figure(4), clf, plotland(.8*[1 1 1],12), hold on
for i=1:75
    ix=find(jra_r>(i-1)*100&jra_r<=i*100);
    if ~isempty(ix)
        plot(jlon(ix),jlat(ix),'ro','markersize',i+3)
    end
    ix=find(gns_r>(i-1)*100&gns_r<=i*100);
    if ~isempty(ix)
        plot(glon(ix),glat(ix),'bo','markersize',i+3)
    end
end
clear I Y d gD gIX i*
title('Residuals')
text(185,-50,['jra ' int2str(sum(jra_r)) ' km^3/yr'],'color','r')
text(185,-60,['gns ' int2str(sum(gns_r)) ' km^3/yr'],'color','b')
print -dpdf Residuals


disp(full([maxsep max(jra_r) max(gns_r) sum(jra_r) sum(gns_r) length(find(jra>0&jra==jra_r)) length(find(gns>0&gns==gns_r))]))
[tmp gx]=sort(gns,'descend');
2.50       7395.46       6380.42      12766.59      12795.27       3364.00       1493.00
3.00       1013.98        278.38       6215.61       6244.29       3282.00       1528.00
4.00       1008.29        276.70       5933.44       5962.12       3118.00       1644.00
5.00       1008.29        276.10       5620.04       5648.72       2940.00       1724.00

[tmp gx]=sort(gns,'ascend');
2.50       7434.36       6429.47      13183.29      13211.97       3667.00        403.00
3.00       1004.89        278.47       6482.96       6511.64       3563.00        382.00
4.00       1004.89        278.54       6241.80       6270.47       3417.00        348.00
5.00       1004.89        300.92       5985.80       6014.47       3296.00        319.00
