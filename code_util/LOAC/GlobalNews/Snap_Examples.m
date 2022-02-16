% Example snapping of GlobalNEWS2 to JRA55 locations
clear, close all
cd ~dmenemen/projects/LOAC/runoff_products/GlobalNEWS
load GlobalNews_to_JRA55

% Compute weights
ix=find(jra);
jraWeights=jra*0;
jraWeights(ix)=jra(ix)./gQact(gQact2jra(ix));

% As a sanity check, snap year-2000 runoff (km^3/yr)
ix=find(jra);
jQact=jra*0;
jQact(ix)=gQact(gQact2jra(ix)).*jraWeights(ix);
disp(full([max(jra) max(gQact) max(abs(jra-jQact))]))

% jra55_do original grid
lon=0.125:0.25:360; nx=length(lon);
lat=-89.875:0.25:90; ny=length(lat); 

% Projecct GlobalNEWS2 nutrients to JRA55 locations
ix=find(jra);
for f={'DIN','DIP','DON','DOP','DOC','DSi','PN','PP','POC','TSS'}
    tmp=jra*0;
    eval(['fld=g' f{1} ';'])
    tmp(ix)=fld(gQact2jra(ix)).*jraWeights(ix);
    disp(['percent difference of sum(JRA) and sum(GNews2) for ' f{1} ': ' int2str(100*(sum(tmp)-sum(fld))/sum(tmp))])
end
