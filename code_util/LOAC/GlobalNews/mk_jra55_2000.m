clear, close all
clear, close all

% pathnames below are on pleiades
pn1='~dmenemen/forcing/jra55_do/';
pn2='/nobackup/dcarrol2/LOAC/bin/jra55_do/v1.4.0/';

% Compute jra55_do time-mean year-2000 runoff
yr=2000;
loy=365;
if mod(yr,4)==0, loy=366; end

lon=0.125:0.25:360;
lat=-89.875:0.25:90;
nx=length(lon);
ny=length(lat);
[jlat jlon]=meshgrid(lat,lon);

fin=[pn1 'cellarea.bin'];
cellarea=readbin(fin,[nx ny]);  % area of grid cell in m^2

fin=[pn1 'jra55_do_runoff_' int2str(yr)];
jra=readbin(fin,[nx ny loy]);              % load whole year, units are m/s
jra=sum(jra,3).*cellarea*24*60*60/1e9;     % cumulate & convert to km^3/yr

jlon=jlon(:);
jlat=jlat(:);
jra=jra(:);
ix=find(jra==0); jlon(ix)=0; jlat(ix)=0;
jra=sparse(jra); jlon=sparse(jlon); jlat=sparse(jlat);

save jra55_2000 j*
