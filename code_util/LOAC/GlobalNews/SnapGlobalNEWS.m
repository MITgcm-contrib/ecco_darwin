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

% find indices IX the jlat/jlon location on the jra55_do grid
lon=0.125:0.25:360; nx=length(lon);
lat=-89.875:0.25:90; ny=length(lat); 
[LAT LON]=meshgrid(lat,lon);
IX=jlat;
for i=1:length(jlat)
    IX(i)=find(LAT==jlat(i)&LON==jlon(i)); 
end

% Compute weights, that is, the ratio of jra55_do runoff volume relative to
% the GlobalNEWS2 location associated with each jra55_do location 
jraWeights=jra*0;
jraWeights=jra./gQact(gQact2jra);

% Projecct GlobalNEWS2 nutrients to JRA55 locations
pin='/nobackup/dcarrol2/LOAC/bin/jra55_do/v1.4.0/';
pout='~dmenemen/forcing/jra55_do/GlobalNEWS/GlobalNEWS2_on_jra55v1.4.0/';
for yr=1991:2021
    fin=[pin 'jra55_do_runoff_' int2str(yr)];
    loy=365;
    if mod(yr,4)==0, loy=366; end
    for dy=1:loy, disp([yr dy])
        Jravol=readbin(fin,[nx ny],1,'real*4',dy-1);
        for f={'DIN','DIP','DON','DOP','DOC','DSi','PN','PP','POC','TSS'}
            fout=[pout f{1} '_' int2str(yr)];
            eval(['fld=g' f{1} ';'])
            FLD=0*LAT;
            FLD(IX)=fld(gQact2jra).*jraWeights;
            writebin(fout,FLD,1,'real*4',dy-1);
        end
    end
end
