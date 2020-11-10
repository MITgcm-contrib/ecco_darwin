clear all


%global llc270
nX=270;nY=nX*13;nZ=50;
dirGrid='/nobackup/hzhang1/llc_1080/MITgcm/DM_270/GRID_up/';
XC=readbin([dirGrid 'XC.data'],[nX nY]);
YC=readbin([dirGrid 'YC.data'],[nX nY]);
HC=readbin([dirGrid 'hFacC.data'],[nX nY]);
DC=readbin([dirGrid 'Depth.data'],[nX nY]);

rc=-readbin([dirGrid 'RC.data'],nZ);       % depths to center of cell
rf=-readbin([dirGrid 'RF.data'],nZ+1);     % depths to cell faces
thk=diff(rf);                              % thicknesses


% domain-specific preamble
LONLIMS = [-145 -126];
LATLIMS = [68.5 72];

i1=232:270; j1=203:270; %face 3
i2=1:7;     j2=j1;      %face 4
kx = 1:44;
nme='Mac';                        % domain name
nt=312;                           % number of obcs time steps
%face 3
I3=(nX*6+1):(nX*7); %relative to (nX,nY)
%face 4 /tile 8
f=8;
I2=(1:3:(nX*3))+7*nX+f-8; %relative to (nX,nY)
xc7=XC(:,I3);yc7=YC(:,I3); dc7=DC(:,I3);
xc8=XC(:,I2);yc8=YC(:,I2); dc8=DC(:,I2);
xcnew=[xc7(i1,j1); xc8(i2,j2)];
ycnew=[yc7(i1,j1); yc8(i2,j2)];
dcnew=[dc7(i1,j1); dc8(i2,j2)];

% derived quantities
nx=length(i1)+length(i2); ny=length(j1); nz=length(kx);
dim=[num2str(nx) 'x' num2str(ny)];


% directory names (may need to be created or modified)
pin='run_global/';
pi2=dirGrid;
pout='run_template/';    % output path name

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make pickup files
% requires that pickup*.meta files be constructed mannually
%%
pic='pickup_ptracers.0000289296.data';
fn=[pin pic];
disp(fn)
KK=31;
fld50=zeros([nx ny nZ*KK]);
fld44=zeros([nx ny nz*KK]);
for k=1:nZ*KK
 SST=readbin(fn,[nX nY],1,'real*8',k-1);
 st7=SST(:,I3);st8=SST(:,I2);
 fld50(:,:,k)=[st7(i1,j1); st8(i2,j2)];
end
for k=1:KK
 fld44(:,:,(1:nz)+(k-1)*nz)=fld50(:,:,(1:nz)+(k-1)*nZ);
end
 writebin([pout pic],fld44,1,'real*8');

pic='pickup_darwin.0000289296.data';
fn=[pin pic];
disp(fn)
KK=1;
fld50=zeros([nx ny nZ*KK]);
fld44=zeros([nx ny nz*KK]);
for k=1:nZ*KK
 SST=readbin(fn,[nX nY],1,'real*8',k-1);
 st7=SST(:,I3);st8=SST(:,I2);
 fld50(:,:,k)=[st7(i1,j1); st8(i2,j2)];
end
for k=1:KK
 fld44(:,:,(1:nz)+(k-1)*nz)=fld50(:,:,(1:nz)+(k-1)*nZ);
end
 writebin([pout pic],fld44,1,'real*8');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate lateral boundary conditions for "iron"
%%
ironfile = 'llc270_Mahowald_2009_soluble_iron_dust.bin';
ironperiod = -12;
fn=[pin ironfile];
disp(fn)
runoff=readbin(fn,[nX nY 12]);
runoff2=    zeros([nx ny 12]);
for i=1:12
SST=runoff(:,:,i);
st7=SST(:,I3);st8=SST(:,I2);
runoff2(:,:,i)=[st7(i1,j1); st8(i2,j2)];
end
writebin([pout 'iron_dust_' dim '_' nme],runoff2);
clear runoff*


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate lateral boundary conditions
%TS='from270/';
TS='/nobackup/dcarrol2/v05_release/darwin3/run/diags/monthly/';
tmp=zeros([nx ny nz]);
fla=zeros([nx ny]);

flds={'DIC', 'NO3','NO2','NH4','PO4','FeT','SiO2','DOC','DON', ...
      'DOP','DOFe','POC','PON','POP','POFe','POSi','PIC','ALK','O2', ...
      'c1','c2','c3','c4','c5','c6','c7','Chl1','Chl2','Chl3', ... 
      'Chl4','Chl5'};
flds2={'tr1', 'tr2', 'tr3', 'tr4', 'tr5', 'tr6', 'tr7', 'tr8', 'tr9', 'tr10', ...
       'tr11','tr12','tr13','tr14','tr15','tr16','tr17','tr18','tr19','tr20', ...
       'tr21','tr22','tr23','tr24','tr25','tr26','tr27','tr28','tr29','tr30', ...
       'tr31'};
      
for f=1:length(flds)
fld=flds{f};
fld2=flds2{f};
	fnm=dir([TS fld '.*data']);

for i=1:length(fnm), mydisp(i)
  tmp=zeros([nx ny nz]);
  fla=zeros([nx ny]);
  for k=1:nz
  	fni=[TS fnm(i).name];
	if k==1; disp(fni); end
  	SST=readbin(fni,[nX nY],1,'real*4',k -1);
	st7=SST(:,I3);st8=SST(:,I2);
	fla=[st7(i1,j1); st8(i2,j2)];
	tmp(:,:,k)=fla;
  end	

  OBN=tmp(:  ,end,:);
  OBS=tmp(:  ,1,  :);
  OBE=tmp(end,:,  :);
  OBW=tmp(1  ,:,  :);
  if i==1
    writebin([pout 'OBN' fld2 '_' nme '_' dim '.bin'],OBN)
    writebin([pout 'OBS' fld2 '_' nme '_' dim '.bin'],OBS)
    writebin([pout 'OBW' fld2 '_' nme '_' dim '.bin'],OBW)
    writebin([pout 'OBE' fld2 '_' nme '_' dim '.bin'],OBE)
  end
  writebin([pout 'OBN' fld2 '_' nme '_' dim '.bin'],OBN,1,'real*4',i)
  writebin([pout 'OBS' fld2 '_' nme '_' dim '.bin'],OBS,1,'real*4',i)
  writebin([pout 'OBW' fld2 '_' nme '_' dim '.bin'],OBW,1,'real*4',i)
  writebin([pout 'OBE' fld2 '_' nme '_' dim '.bin'],OBE,1,'real*4',i)
end %i
end %fld

