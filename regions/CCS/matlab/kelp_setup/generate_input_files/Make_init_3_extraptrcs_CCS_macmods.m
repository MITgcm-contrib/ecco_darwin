
% DOMAIN X Y Z DIMENSIONS
nx=80;
ny=100;
nz=48;
varo=ones(nx,ny,nz);

% In Seaweeds we have C:N:P = 550:35:1
% then C ~ 15.7 * N
% approx 15.7 to 16 and use N=10
% then C = N * 16 = 160 mmmol/m3


varC=varo.*160; % mmol/m3
varN=varo.*10; % mmol/m3
varB=varo.*50; % g/m2 (dry weight?)
%varMaFra=varo.*30*1e3;

%Chl5_80x100x48.16-Jan-1992

fnamoB='ini_CCS_80x100x48_Bmag.bin';

fnamoC='ini_CCS_80x100x48_Cmag.bin';

fnamoN='ini_CCS_80x100x48_Nmag.bin';

typp=1;

writebin(fnamoB,varB,typp);
writebin(fnamoC,varC,typp);
writebin(fnamoN,varN,typp);




