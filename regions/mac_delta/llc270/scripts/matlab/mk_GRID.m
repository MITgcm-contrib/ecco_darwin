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
% create grid input

flds1={'LONC','LATC','DXF','DYF','RA', 'LONG','LATG','DXV','DYU','RAZ','DXC','DYC','RAW','RAS','DXG','DYG'};
flds2={'XC',  'YC',  'DXF','DYF','RAC','XG'  ,'YG'  ,'DXV','DYU','RAZ','DXC','DYC','RAW','RAS','DXG','DYG'};
for f=1:length(flds1)
	fld1=flds1{f};
	fld2=flds2{f};

	SST=readbin([dirGrid fld2 '.data'],[nX nY]);
	st7=SST(:,I3);st8=SST(:,I2);
	tmp4=[st7(i1,j1); st8(i2,j2)];
        writebin([pout fld1 '.bin'],tmp4)
end	


