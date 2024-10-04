

dim1 = 80; % points along X
dim2 = 48; % points along Z
dim3=100;  % points along Y
dimt=374; % time dimension (monthly)

Cval=1e-2;
Nval=Cval;

varMaNum=1e3;
varMaFra=30*1e3;

varN=ones(dim1,dim2,dimt);
%varN=varN*const_val;
varS=varN;

varW =ones(dim3,dim2,dimt); 
%varW=varW*const_val;

typp=1;


% Files North

fnamo1='CCS_kelp_80x100x48_MaNum_North.bin';

fnamo2='CCS_kelp_80x100x48_MaFra_North.bin';

fnamo3='CCS_kelp_80x100x48_MaC_North.bin';

fnamo4='CCS_kelp_80x100x48_MaN_North.bin';


% Files South

fnamo5='CCS_kelp_80x100x48_MaNum_South.bin';

fnamo6='CCS_kelp_80x100x48_MaFra_South.bin';

fnamo7='CCS_kelp_80x100x48_MaC_South.bin';

fnamo8='CCS_kelp_80x100x48_MaN_South.bin';

% Files West

fnamo9='CCS_kelp_80x100x48_MaNum_West.bin';

fnamo10='CCS_kelp_80x100x48_MaFra_West.bin';

fnamo11='CCS_kelp_80x100x48_MaC_West.bin';

fnamo12='CCS_kelp_80x100x48_MaN_West.bin';

writebin(fnamo1,varN.*varMaNum,typp);
writebin(fnamo2,varN.*varMaFra,typp);
writebin(fnamo3,varN.*Cval,typp);
writebin(fnamo4,varN.*Nval,typp);

writebin(fnamo5,varS.*varMaNum,typp);
writebin(fnamo6,varS.*varMaFra,typp);
writebin(fnamo7,varS.*Cval,typp);
writebin(fnamo8,varS.*Nval,typp);

writebin(fnamo9,varW.*varMaNum,typp);
writebin(fnamo10,varW.*varMaFra,typp);
writebin(fnamo11,varW.*Cval,typp);
writebin(fnamo12,varW.*Nval,typp);




