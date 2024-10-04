

dim1 = 80; % points along X
dim2 = 48; % points along Z
dim3=100;  % points along Y
dimt=374; % time dimension (monthly)


% Now generating values for the BC files
% for the new variables for MACMODS
% Bmag :: Biomass macroalgae (gC/m2)
% Cmag :: Carbon  macroalgae (mmol/m3)
% Nmag :: Nitrogen Macroalgae (mmol/m3)

% setting the values

Nval=10;  
Cval=16*Nval;
Bval=50;

Bmag=Bval;
Nmag=Nval;
Cmag=Cval;



varN=ones(dim1,dim2,dimt);
%varN=varN*const_val;
varS=varN;

varW =ones(dim3,dim2,dimt); 
%varW=varW*const_val;

typp=1;


% Files North

fnamo1='CCS_kelp_80x100x48_Bmag_North.bin';

fnamo2='CCS_kelp_80x100x48_Nmag_North.bin';

fnamo3='CCS_kelp_80x100x48_Cmag_North.bin';

writebin(fnamo1,varN.*Bmag,typp);
writebin(fnamo2,varN.*Nmag,typp);
writebin(fnamo3,varN.*Cmag,typp);


% Files South

fnamo4='CCS_kelp_80x100x48_Bmag_South.bin';

fnamo5='CCS_kelp_80x100x48_Nmag_South.bin';

fnamo6='CCS_kelp_80x100x48_Cmag_South.bin';

writebin(fnamo4,varS.*Bmag,typp);
writebin(fnamo5,varS.*Nmag,typp);
writebin(fnamo6,varS.*Cmag,typp);


% Files West

fnamo7='CCS_kelp_80x100x48_Bmag_West.bin';

fnamo8='CCS_kelp_80x100x48_Nmag_West.bin';

fnamo9='CCS_kelp_80x100x48_Cmag_West.bin';

writebin(fnamo7,varW.*Bmag,typp);
writebin(fnamo8,varW.*Nmag,typp);
writebin(fnamo9,varW.*Cmag,typp);








