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
nt=datenum(2021,1,1)-datenum(1992,1,1);
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
%face-view
i2_ext=   j2;
j2_ext= (-i2(end)+3*nX+1) : (-i2(1)+3*nX+1);
%or      -i2(end:-1:1)+3*nX+1



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

regional_grid='../grid/';
genBC={'W','S'};                    % generate genBC boundary conditions
balanceBC='W';                          % balance E, W, N, or S boundary condition


%step=1; %bathy + grid + IC
step=2; %BC
if step==1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make bathymetry file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make pickup files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make weight file (for ctrl)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make smooth file (for smooth)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make xx file (for ctrl)

else %step2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate lateral boundary conditions
%%

TS='diags/';
%SI
flds={'SIarea', 'SIheff',  'SIhsnow', 'SIuice', 'SIvice'};
flds2={ 'a',      'h',       'sn',     'uice',   'vice'};
months=nt;
for f=1:length(flds)
fld=flds{f};
fld2=flds2{f};

for i=1:months
  tmp=zeros([nx ny]);  fla=zeros([nx ny]);

	ts=i*72;
	dd=datestr(datenum(1992,i,1),'yyyymm');
  	fni=[TS fld '.' myint2str(ts,10) '.data'];
	disp(fni)
  	SST=readbin(fni,[nX nY]);
%Mac cutout
	st7=SST(:,I3);st8=SST(:,I2);
	fla=[st7(i1,j1); st8(i2,j2)];
  tmp=fla;

  OBN=tmp(:  ,end);
  if f==5
  OBS=tmp(:  ,2); %V @south
  else
  OBS=tmp(:  ,1);
  end
  OBE=tmp(end,:);
  if f==4
  OBW=tmp(2  ,:); %U @west
  else
  OBW=tmp(1  ,:);
  end
  if i==1
    writebin([pout 'OBN' fld2 '_' nme '_' dim '.bin_29yrs'],OBN)
    writebin([pout 'OBS' fld2 '_' nme '_' dim '.bin_29yrs'],OBS)
    writebin([pout 'OBW' fld2 '_' nme '_' dim '.bin_29yrs'],OBW)
    writebin([pout 'OBE' fld2 '_' nme '_' dim '.bin_29yrs'],OBE)
  end
  writebin([pout 'OBN' fld2 '_' nme '_' dim '.bin_29yrs'],OBN,1,'real*4',i)
  writebin([pout 'OBS' fld2 '_' nme '_' dim '.bin_29yrs'],OBS,1,'real*4',i)
  writebin([pout 'OBW' fld2 '_' nme '_' dim '.bin_29yrs'],OBW,1,'real*4',i)
  writebin([pout 'OBE' fld2 '_' nme '_' dim '.bin_29yrs'],OBE,1,'real*4',i)
end %i
  i=i+1;
  writebin([pout 'OBN' fld2 '_' nme '_' dim '.bin_29yrs'],OBN,1,'real*4',i)
  writebin([pout 'OBS' fld2 '_' nme '_' dim '.bin_29yrs'],OBS,1,'real*4',i)
  writebin([pout 'OBW' fld2 '_' nme '_' dim '.bin_29yrs'],OBW,1,'real*4',i)
  writebin([pout 'OBE' fld2 '_' nme '_' dim '.bin_29yrs'],OBE,1,'real*4',i)
end %fld

end %step2
