% Snap GlobalNEWS2 2000 time-mean Qact runoff to jra55_do grid
clear, close all

% Load JRA55-do time-mean year-2000 runoff
% jlat/jlon: latitude/longitude of jra55_do
% jns: jra55_do runoff km^3/yr
load jra55_2000

% File globalnews.xlsx is derived from
% GlobalNEWS2__RH2000Dataset-version1.0.xls
% The columns of globalnews.xlsx are, respectively:
%       basins / mouth_lon
%       basins / mouth_lat
%    hydrology / Qact
% river export / Ld_DIN
% river export / Ld_DIP
% river export / Ld_DON
% river export / Ld_DOP
% river export / Ld_DOC
% river export / Ld_DSi
% river export / Ld_PN
% river export / Ld_PP
% river export / Ld_POC
% river export / Ld_TSS

% Load GlobalNEWS2 mouth_lon, mouth_lat, Qact (km3/yr) and end of the basin
% (land or ocean)
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

load GlobalNEWS_basins_end.mat
gbasins_end = GlobalNEWSbasinsend;

load GlobalNEWS_basinsID.mat
gbasins = GlobalNEWS_basinsID;

% File DIC_final_globalnews.xlsx is derived from 
% computation of DIC fluxes from BDIC_GlobalNEWS.xlsx
gns=xlsread('DIC_final_globalnews.xlsx');
gDIC=gns(:,1); % GlobalNEWS2 load DIC (Mg/yr)
clear gns ix

% Correct DIC inputs from the Amazon
% relation from Li et al., 2017 overestimates DIC export three times than
% literature range (da Cuha et al., 2013, Probst et al, 1994)
gDIC(1) = gDIC(1)/3;

ix = find(strcmp(cellstr(gbasins_end),"Land"));
for f={'lat','lon','Qact','DIN','DIP','DON','DOP','DOC','DSi','PN','PP','POC','TSS','DIC'}
    eval(['g' f{1} '(ix) = [];'])
end

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
ix=find(jlat>=-60);
jra=jra(ix);
jlat=jlat(ix);
jlon=jlon(ix);
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
% Identification of JRA point according volume.  This may lead to certain
% high volume location of GlobalNEWS not being used.
[gQact,idx]= sort(gQact,'ascend');
for f={'lat','lon','DIN','DIP','DON','DOP','DOC','DSi','PN','PP','POC','TSS','DIC','basins'}
    eval(['g' f{1} ' = g' f{1} '(idx);'])
end

maxsep=5;                        % max separation between GlobalNEWS and JRA55 in deg
minvol=0;                        % minimum volume for matching volume algorithm
gQact2jra=zeros(length(jra),1);  % GlobalNEWS index for each non-zero JRA55 location
gx=find(gQact);                  % indices of non-zero GlobalNEWS locations
for j=1:length(jra)              % loop through all non-zero JRA55 locations
    d=sqrt((glon(gx)-jlon(j)).^2+(glat(gx)-jlat(j)).^2); % distance in deg
    % select closest non-zero GlobalNEWS location
    [M,I]=min(d);
    %gQact2jra(j)=gx(I);
    if jra(j) > minvol
        % select GlobalNEWS location with closest runoff
        % volume to JRA55 within maxsep degrees
        dx=find(d<maxsep);
        if ~isempty(dx)
            I=closest(jra(j),gQact(gx(dx)));
            gQact2jra(j)=gx(dx(I));
        end
    end
end

%set zero values to smallest discharge possible
gQact2jra(gQact2jra == 0) = 1;
%% PREVIOUS CORRECTION FOR DISTANCE SNAPPING
% WARNING: snap on mismatched volume for the Amazon

sort(unique(gQact(gQact2jra)),'descend');
sort(unique(gbasins(gQact2jra)),'ascend');

%%%%% iterations until the first 100 basins are captured
% find non-attributed basins because of non extension to the coastline of
% the basin
missing_basins = setdiff(gbasins(gx), unique(gbasins(gQact2jra)));

while min(missing_basins)<100
% ix = ismember(missing_basins,gbasins);
% missing_basins = missing_basins(ix);
missing_basins = sort(missing_basins,'descend');
%find the closest jra point of river mouth of these basins    
for i = 1 :length(missing_basins)
    ix = find(gbasins(gx) == missing_basins(i));
    if isempty(ix)
        continue
    end
    d = sqrt((glon(gx(ix))-jlon).^2+(glat(gx(ix))-jlat).^2); % distance in deg
    [M,I_jra]=min(d);
%     ix = find(gbasins == missing_basins(i));
%     d = sqrt((glon(ix)-jlon(I_jra)).^2+(glat(ix)-jlat(I_jra)).^2); % distance in deg
%     [M,I_gn]=min(d);
    gQact2jra(I_jra) = gx(ix);
end
%%%%%
missing_basins = setdiff(gbasins(gx), unique(gbasins(gQact2jra)));
end

disp("first 100 basins captured")
disp([num2str(length(unique(gbasins(gQact2jra)))) "basins captured"]);

figure(3), clf, plotland(.8*[1 1 1],12), hold on
for i=1:75
    ix=find(jra>(i-1)*100&jra<=i*100);
    if ~isempty(ix)
        plot(jlon(ix),jlat(ix),'ro','markersize',i+3)
    end
    ix=find(gQact(gQact2jra)>(i-1)*100&gQact(gQact2jra)<=i*100);
    if ~isempty(ix)
        plot(jlon(ix),jlat(ix),'b+','markersize',i+3)
    end
end
title('No Antarctica')
text(185,-50,['jra ' int2str(sum(jra)) ' km^3/yr'],'color','r')
text(185,-60,['gQact ' int2str(sum(gQact(gQact2jra))) ' km^3/yr'],'color','b')


%correct Amazon point to avoid high concentrations
% WARNING: snap on mismatched volume for the Amazon
tmp = gQact2jra(3418);
gQact2jra(3418) = gQact2jra(3443);
gQact2jra(3443) = 1;
%%% END CORRECTION

figure(3), clf, plotland(.8*[1 1 1],12), hold on
for i=1:75
    ix=find(jra>(i-1)*100&jra<=i*100);
    if ~isempty(ix)
        plot(jlon(ix),jlat(ix),'ro','markersize',i+3)
    end
    ix=find(gQact(gQact2jra)>(i-1)*100&gQact(gQact2jra)<=i*100);
    if ~isempty(ix)
        plot(jlon(ix),jlat(ix),'b+','markersize',i+3)
    end
end
title('No Antarctica')
text(185,-50,['jra ' int2str(sum(jra)) ' km^3/yr'],'color','r')
text(185,-60,['gQact ' int2str(sum(gQact(gQact2jra))) ' km^3/yr'],'color','b')

save GlobalNews_to_JRA55
