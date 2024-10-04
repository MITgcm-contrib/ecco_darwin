
nx=80;
ny=100;
nz=48;
varo=ones(nx,ny,nz);
varC=varo.*1e-1;
varN=varo.*1e-2;
varMaNum=varo.*1e7;
varMaFra=varo.*30*1e3;

%Chl5_80x100x48.16-Jan-1992

fnamo1='ini_CCS_80x100x48_MaNum.bin';

fnamo2='ini_CCS_80x100x48_MaFra.bin';

fnamo3='ini_CCS_80x100x48_MaC.bin';

fnamo4='ini_CCS_80x100x48_MaN.bin';

typp=1;

writebin(fnamo1,varMaNum,typp);
writebin(fnamo2,varMaFra,typp);
writebin(fnamo3,varC,typp);
writebin(fnamo4,varN,typp);




