% Check annual of bgc runoff spread on Fekete
clear, close all

addpath(genpath('/nobackup/dcarrol2/MATLAB'));

% Load LLC90 grid files
gridDir = '/nobackup/dcarrol2/LOAC/grid/ECCO_V4r5_raw/';
nx = 90;
ny = 1170;
nz = 50;
RAC = readbin([gridDir 'RAC.data'],[nx ny],1,'real*4');

load River_Fekete_Mapping_Continuous_RadiusLimited_x1.mat
pin = '/nobackup/rsavelli/LOAC/Fekete/bgc_runoff/spread_x1/';

DIN = readbin([pin 'DIN_Fekete_ECCO_V4r5.bin'],[nx ny 12],1,'real*4');
DIP = readbin([pin 'DIP_Fekete_ECCO_V4r5.bin'],[nx ny 12],1,'real*4');
DON = readbin([pin 'DON_Fekete_ECCO_V4r5.bin'],[nx ny 12],1,'real*4');
DOP = readbin([pin 'DOP_Fekete_ECCO_V4r5.bin'],[nx ny 12],1,'real*4');
DOC = readbin([pin 'DOC_Fekete_ECCO_V4r5.bin'],[nx ny 12],1,'real*4');
DSi = readbin([pin 'DSi_Fekete_ECCO_V4r5.bin'],[nx ny 12],1,'real*4');
PN = readbin([pin 'PN_Fekete_ECCO_V4r5.bin'],[nx ny 12],1,'real*4');
PP = readbin([pin 'PP_Fekete_ECCO_V4r5.bin'],[nx ny 12],1,'real*4');
POC = readbin([pin 'POC_Fekete_ECCO_V4r5.bin'],[nx ny 12],1,'real*4');
DIC = readbin([pin 'DIC_Fekete_ECCO_V4r5.bin'],[nx ny 12],1,'real*4');

% conversion factors of gram to mol
gP_to_molP = 0.03228539149637;
gN_to_molN = 0.071394404106606;
gC_to_molC = 0.083259093974539;
gSi_to_molSi = 0.03560556158872;


%% check global loads

tmp = DOC .* RAC;
fprintf('\nGlobal DOC load %d Tg C yr-1\n',nansum(tmp(:))*...
1/gC_to_molC*1e-15*60*60*24*30.5);
tmp = POC.*RAC;
fprintf('\nGlobal POC load %d Tg C yr-1\n',nansum(tmp(:))*...
1/gC_to_molC*1e-15*60*60*24*30.5);
tmp = DIC.*RAC;
fprintf('\nGlobal DIC load %d Tg C yr-1\n',nansum(tmp(:))*...
1/gC_to_molC*1e-15*60*60*24*30.5);
tmp = DON.*RAC;
fprintf('\nGlobal DON load %d Tg N yr-1\n',nansum(tmp(:))*...
1/gN_to_molN*1e-15*60*60*24*30.5);
tmp = DIN.*RAC;
fprintf('\nGlobal DIN load %d Tg N yr-1\n',nansum(tmp(:))*...
1/gN_to_molN*1e-15*60*60*24*30.5);
tmp = PN.*RAC;
fprintf('\nGlobal PN load %d Tg N yr-1\n',nansum(tmp(:))*...
1/gN_to_molN*1e-15*60*60*24*30.5);
tmp = DIP.*RAC;
fprintf('\nGlobal DIP load %d Tg P yr-1\n',nansum(tmp(:))*...
1/gP_to_molP*1e-15*60*60*24*30.5);
tmp = DOP.*RAC;
fprintf('\nGlobal DOP load %d Tg P yr-1\n',nansum(tmp(:))*...
1/gP_to_molP*1e-15*60*60*24*30.5);
tmp = PP.*RAC;
fprintf('\nGlobal PP load %d Tg P yr-1\n',nansum(tmp(:))*...
1/gP_to_molP*1e-15*60*60*24*30.5);
tmp = DSi.*RAC;
fprintf('\nGlobal DSi load %d Tg Si yr-1\n',nansum(tmp(:))*...
1/gSi_to_molSi*1e-15*60*60*24*30.5);

pause()
%% check individual rivers
%1	Amazon
%2	Nile
%3	Zaire
%4	Mississippi
%5	Ob
%6	Parana
%7	Yenisei
%8	Lena
%9	Niger
%10	Tamanrasett
%11	Chang Jiang
%12	Amur
%13	Mackenzie
%14	Ganges
%15	CHARI BOUSSO/ Lake Chad
%16	Volga
%17	Zambezi
%18	GreatArtesian/CopperCR
%19	Tarim
%20	Indus
%21	Nelson
%22	Kerulen
%23	Syr-Darya
%24	SAINT LAWRENCE
%25	Orinoco
ID = 1;
river_IDidx = find([GN2fekete.riverID] == ID);
fprintf('\nRiver ID %d\n', ID)
DOC= nansum(DOC,3);
tmp = DOC(GN2fekete(river_IDidx).LLC90).*RAC(GN2fekete(river_IDidx).LLC90);
fprintf('\nDOC load %d Tg C yr-1\n',nansum(tmp(:))*...
(1/gC_to_molC)*1e-15*60*60*24*30.5);
DIC = nansum(DIC,3);
tmp = DIC(GN2fekete(river_IDidx).LLC90).*RAC(GN2fekete(river_IDidx).LLC90);
fprintf('\nDIC load %d Tg C yr-1\n',nansum(tmp(:))*...
(1/gC_to_molC)*1e-15*60*60*24*30.5);
DIN = nansum(DIN,3);
tmp = DIN(GN2fekete(river_IDidx).LLC90).*RAC(GN2fekete(river_IDidx).LLC90);
fprintf('\nDIN load %d Tg N yr-1\n',nansum(tmp(:))*...
(1/gN_to_molN)*1e-15*60*60*24*30.5);
DON = nansum(DON,3);
tmp = DON(GN2fekete(river_IDidx).LLC90).*RAC(GN2fekete(river_IDidx).LLC90);
fprintf('\nDON load %d Tg N yr-1\n',nansum(tmp(:))*...
(1/gN_to_molN)*1e-15*60*60*24*30.5);
DIP = nansum(DIP,3);
tmp = DIP(GN2fekete(river_IDidx).LLC90).*RAC(GN2fekete(river_IDidx).LLC90);
fprintf('\nDIP load %d Tg P yr-1\n',nansum(tmp(:))*...
(1/gP_to_molP)*1e-15*60*60*24*30.5);
DOP = nansum(DOP,3);
tmp = DOP(GN2fekete(river_IDidx).LLC90).*RAC(GN2fekete(river_IDidx).LLC90);
fprintf('\nDOP load %d Tg P yr-1\n',nansum(tmp(:))*...
(1/gP_to_molP)*1e-15*60*60*24*30.5);
DSi = nansum(DSi,3);
tmp = DSi(GN2fekete(river_IDidx).LLC90).*RAC(GN2fekete(river_IDidx).LLC90);
fprintf('\nDSi load %d Tg Si yr-1\n',nansum(tmp(:))*...
(1/gSi_to_molSi)*1e-15*60*60*24*30.5);

