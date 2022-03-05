% Example snapping of GlobalNEWS2 to JRA55 locations
clear, close all
cd ~dmenemen/forcing/jra55_do/GlobalNEWS
load GlobalNews_to_JRA55

% glon  GlobalNEWS2 longitude E (deg)
% glat  GlobalNEWS2 latitude N (deg)
% gQact GlobalNEWS2 actual discharge (km^3/yr)
% gDIN  GlobalNEWS2 load DIN (Mg/yr)
% gDIP  GlobalNEWS2 load DIP (Mg/yr)
% gDON  GlobalNEWS2 load DON (Mg/yr)
% gDOP  GlobalNEWS2 load DON (Mg/yr)
% gDOC  GlobalNEWS2 load DOC (Mg/yr)
% gDSi  GlobalNEWS2 load DSi (Mg/yr)
% gPN   GlobalNEWS2 load PN (Mg/yr)
% gPP   GlobalNEWS2 load PP (Mg/yr)
% gPOC  GlobalNEWS2 load POC (Mg/yr)
% gTSS  GlobalNEWS2 load TSS (Mg/yr)

% jlat/jlon : latitude/longitude of jra55_do
% jra       : jra55_do year-2000 runoff km^3/yr

% gQact2jra : index of GlobalNEWS2 location
%             for each jra55_do location that

% Compute weights, that is, the ratio of jra55_do runoff volume relative to
% the GlobalNEWS2 location associated with each jra55_do location 
jraWeights=jra*0;
jraWeights=jra./gQact(gQact2jra);

% As a sanity check, snap year-2000 runoff (km^3/yr)
jQact=jra*0;
jQact=gQact(gQact2jra).*jraWeights;
disp([max(jra) max(gQact) max(abs(jra-jQact))])
disp([sum(jra) sum(jQact)])

% Projecct GlobalNEWS2 nutrients to JRA55 locations
for f={'DIN','DIP','DON','DOP','DOC','DSi','PN','PP','POC','TSS'}
    eval(['fld=g' f{1} ';'])
    tmp=fld(gQact2jra).*jraWeights;
    disp(['percent difference of sum(JRA) and sum(GNews2) for ' f{1} ': ' int2str(100*(sum(tmp)-sum(fld))/sum(tmp))])
end

% put runoff, etc, back on jra55 grid
% jra55_do original grid

% find indices IX the jlat/jlon location on the jra55_do grid
lon=0.125:0.25:360; nx=length(lon);
lat=-89.875:0.25:90; ny=length(lat); 
[LAT LON]=meshgrid(lat,lon);
IX=jlat;
for i=1:length(jlat)
    IX(i)=find(LAT==jlat(i)&LON==jlon(i)); 
end

% project gQact to jra55_do locations
JRA=0*LAT;
JRA(IX)=gQact(gQact2jra).*jraWeights;
clf, plotland(.8*[1 1 1],12), hold on
for i=1:75
    ix=find(JRA>(i-1)*100&JRA<=i*100);
    if ~isempty(ix)
        plot(LON(ix),LAT(ix),'ro','markersize',i+3)
    end
    ix=find(gQact>(i-1)*100&gQact<=i*100);
    if ~isempty(ix)
        plot(glon(ix),glat(ix),'bo','markersize',i+3)
    end
end
title('No Antarctica')
text(185,-50,['jra ' int2str(sum(jra)) ' km^3/yr'],'color','r')
text(185,-60,['gQact ' int2str(sum(gQact)) ' km^3/yr'],'color','b')
