% Snap GlobalNEWS2 2000 time-mean Qact runoff to jra55_do grid
clear, close all

% Go to wherever you have GlobalNEWS and output from mk_jra55_2000.m
cd ~dmenemen/Documents/projects/LOAC/runoff_products/GlobalNEWS

% Load JRA55-do time-mean year-2000 runoff
% jlat/jlon: latitude/longitude of jra55_do
% jns: jra55_do runoff km^3/yr
load jra55_2000

% globalnews.xlsx is derived from
% GlobalNEWS2__RH2000Dataset-version1.0.xls
% obtained from ...
% the columns are, respectively:
% basins/mouth_lon
% basins/mouth_lat
% hydrology/Qact
% river export/Ld_DIN
% river export/Ld_DIP
% river export/Ld_DON
% river export/Ld_DOP
% river export/Ld_DOC
% river export/Ld_DSi
% river export/Ld_PN
% river export/Ld_PP
% river export/Ld_POC
% river export/Ld_TSS


% Load GlobalNEWS2 mouth_lon, mouth_lat, Qact (km3/yr)
gns=xlsread('globalnews');
glon=gns(:,1);  % GlobalNEWS2 longitude E (deg)
ix=find(glon<0); glon(ix)=glon(ix)+360;
glat=gns(:,2);  % GlobalNEWS2 latitude N (deg)
gQact=gns(:,3); % GlobalNEWS2 actual discharge (km^3/yr)
gDIN=gns(:,4);  % GlobalNEWS2 load DIN (Mg/yr)
gDIP=gns(:,5);  % GlobalNEWS2 load DIP (Mg/yr)
gDON=gns(:,6);  % GlobalNEWS2 load DON (Mg/yr)
gDOP=gns(:,7);  % GlobalNEWS2 load DON (Mg/yr)
gDOC=gns(:,8);  % GlobalNEWS2 load DOC (Mg/yr)
gDSi=gns(:,9);  % GlobalNEWS2 load DSi (Mg/yr)
gPN=gns(:,10);  % GlobalNEWS2 load PN (Mg/yr)
gPP=gns(:,11);  % GlobalNEWS2 load PP (Mg/yr)
gPOC=gns(:,12); % GlobalNEWS2 load POC (Mg/yr)
gTSS=gns(:,13); % GlobalNEWS2 load TSS (Mg/yr)
clear gns ix

% Plot JRA55 and GlobalNEWS runoff
figure(1), clf, plotland(.8*[1 1 1],12), hold on
for i=1:75
    ix=find(jra>(i-1)*100&jra<=i*100);
    if ~isempty(ix)
        plot(jlon(ix),jlat(ix),'ro','markersize',i+3)
    end
    ix=find(gQact>(i-1)*100&gQact<=i*100);
    if ~isempty(ix)
        plot(glon(ix),glat(ix),'bo','markersize',i+3)
    end
end
title('All runoff locations')
text(185,-50,['jra ' int2str(sum(jra)) ' km^3/yr'],'color','r')
text(185,-60,['gQact ' int2str(sum(gQact)) ' km^3/yr'],'color','b')
print -dpdf AllRunoff

% Remove JRA55 Antarctic locations
ix=find(jlat<-60); jra(ix)=0;
figure(2), clf, plotland(.8*[1 1 1],12), hold on
for i=1:75
    ix=find(jra>(i-1)*100&jra<=i*100);
    if ~isempty(ix)
        plot(jlon(ix),jlat(ix),'ro','markersize',i+3)
    end
    ix=find(gQact>(i-1)*100&gQact<=i*100);
    if ~isempty(ix)
        plot(glon(ix),glat(ix),'bo','markersize',i+3)
    end
end
title('No Antarctica')
text(185,-50,['jra ' int2str(sum(jra)) ' km^3/yr'],'color','r')
text(185,-60,['gQact ' int2str(sum(gQact)) ' km^3/yr'],'color','b')
print -dpdf NoAntacrctic

% Associate each JRA55 runoff with a GlobalNEWS location
maxsep=9;                        % max separation between GlobalNEWS and JRA55 in deg
gQact2jra=sparse(length(jra),1); % GlobalNEWS index for each non-zero JRA55 location
jx=find(jra);                    % indices of non-zero JRA55 locations
gx=find(gQact);                  % indices of non-zero GlobalNEWS locations
for j=1:length(jx)               % loop through all non-zero JRA55 locations
    d=sqrt((glon(gx)-jlon(jx(j))).^2+(glat(gx)-jlat(jx(j))).^2); % distance in deg
    dx=find(d<maxsep);           % search within maxsep degrees    
    % select GlobalNEWS location with closest runoff volume to JRA55 within maxsep    
    I=closest(jra(jx(j)),gQact(gx(dx)));
    gQact2jra(jx(j))=gx(dx(I));
end
clear I d* gx i* j jx m*
save GlobalNews_to_JRA55
