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
xcnew=[xc7(i1,j1); xc8(i2,j2)];
ycnew=[yc7(i1,j1); yc8(i2,j2)];
dcnew=[dc7(i1,j1); dc8(i2,j2)];


% directory names (may need to be created or modified)
pin='run_global/';
pi2=dirGrid;
pout='run_template/';    % output path name

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make bathymetry file
fn=[pin 'bathy270_filled_noCaspian_r4'];
disp(fn)
topog=readbin(fn,[nX nY]);
SST=topog;
st7=SST(:,I3);st8=SST(:,I2);
topog=[st7(i1,j1); st8(i2,j2)];
writebin([pout 'bathy_' dim '_' nme],topog);
clear topog

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%runoff
fn=[pin 'runoff-2d-Fekete-1deg-mon-V4-SMOOTH.bin'];
disp(fn)
runoff=readbin(fn,[nX nY 12]);
runoff2=    zeros([nx ny 12]);
for i=1:12
SST=runoff(:,:,i);
st7=SST(:,I3);st8=SST(:,I2);
runoff2(:,:,i)=[st7(i1,j1); st8(i2,j2)];
end
writebin([pout 'runoff_' dim '_' nme],runoff2);
clear runoff*


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make pickup files
% requires that pickup*.meta files be constructed mannually
%%
pic='pickup.0000000001.data';
fn=[pin pic];
disp(fn)
fld50=zeros([nx ny nZ*8+3]);
fld44=zeros([nx ny nz*8+3]);
for k=1:403
 SST=readbin(fn,[nX nY],1,'real*8',k-1);
 st7=SST(:,I3);st8=SST(:,I2);
 fld50(:,:,k)=[st7(i1,j1); st8(i2,j2)];
end
for k=1:8
 fld44(:,:,(1:nz)+(k-1)*nz)=fld50(:,:,(1:nz)+(k-1)*nZ);
end
for k=1:3
 fld44(:,:,nz*8+k)=fld50(:,:,nZ*8+k);
end
 writebin([pout pic],fld44,1,'real*8');

pic='pickup_seaice.0000000001.data';
fn=[pin pic];
disp(fn)
fld44=zeros([nx ny 6]);
for k=1:6
 SST=readbin(fn,[nX nY],1,'real*8',k-1);
 st7=SST(:,I3);st8=SST(:,I2);
 fld44(:,:,k)=[st7(i1,j1); st8(i2,j2)];
end
 writebin([pout pic],fld44,1,'real*8');

pic='pickup_ggl90.0000000001.data';
fn=[pin pic];
disp(fn)
fld44=zeros([nx ny nz]);
for k=1:nz
 SST=readbin(fn,[nX nY],1,'real*8',k-1);
 st7=SST(:,I3);st8=SST(:,I2);
 fld44(:,:,k)=[st7(i1,j1); st8(i2,j2)];
end
 writebin([pout pic],fld44,1,'real*8');

clear tmp* fld*

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make weight file (for ctrl)
flds={'wprecip','wlwdown','wswdown','waqh','watemp','wuwind','wvwind',...
      'wkapgmFld','wkaprediFld','wdiffkrFld'}; 
kk=[ones(1,7) ones(1,3)*nz];
for i=1:length(flds)
fni=[pin  flds{i} '.data'];
	disp(fni)
fno=[pout flds{i} '_' dim '_' nme];
for k=1:kk(i)
        SST=readbin(fni,[nX nY],1,'real*4',k-1);
        st7=SST(:,I3);st8=SST(:,I2);
        tmp3=[st7(i1,j1); st8(i2,j2)];
        writebin(fno,tmp3,1,'real*4',k-1);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make smooth file (for smooth)
%flds={'2Dnorm001','2Dscales001','3Dnorm001','3DscalesZ001','3DscalesH001'};
%kk=[1 2 nz nz nz*2];
 flds={'2Dnorm001','2Dscales001','3Dnorm001','3DscalesZ001'};
 kk=[1 2 nz nz];
for i=1:length(flds)
fni=[pin  'smooth' flds{i} '.data'];
	disp(fni)
fno=[pout 'smooth' flds{i} '_' dim '_' nme];
for k=1:kk(i)
        SST=readbin(fni,[nX nY],1,'real*4',k-1);
        st7=SST(:,I3);st8=SST(:,I2);
        tmp3=[st7(i1,j1); st8(i2,j2)];
        writebin(fno,tmp3,1,'real*4',k-1);
end
end
flds={'3DscalesH001'};
fld50=zeros([nx ny nZ*2]);
fld44=zeros([nx ny nz*2]);
for i=1:length(flds)
fni=[pin  'smooth' flds{i} '.data'];
	disp(fni)
fno=[pout 'smooth' flds{i} '_' dim '_' nme];
for k=1:nZ*2
 SST=readbin(fn,[nX nY],1,'real*4',k-1);
 st7=SST(:,I3);st8=SST(:,I2);
 fld50(:,:,k)=[st7(i1,j1); st8(i2,j2)];
end
for k=1:2
 fld44(:,:,(1:nz)+(k-1)*nz)=fld50(:,:,(1:nz)+(k-1)*nZ);
end
        writebin(fno,fld44);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make xx file (for ctrl)
flds={'aqh','atemp','uwind','vwind','precip','swdown','lwdown',...
      'diffkr','kapgm','kapredi'};
kk=[ones(1,7)*706 ones(1,5)*nz]
for i=1:length(flds)
fni=[pin  'xx_' flds{i} '.0000000042.data'];
	disp(fni)
fno=[pout 'xx_' flds{i} '.0000000042.data'];
for k=1:kk(i),mydisp(k)
        SST=readbin(fni,[nX nY],1,'real*4',k-1);
        st7=SST(:,I3);st8=SST(:,I2);
        tmp3=[st7(i1,j1); st8(i2,j2)];
        writebin(fno,tmp3,1,'real*4',k-1);
end
end

