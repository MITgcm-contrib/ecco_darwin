clear all


%global llc270
nX=270;nY=nX*13;nZ=50;

%face 3
I3=(nX*6+1):(nX*7); %relative to (nX,nY)
%face 4 /tile 8
f=8;
I2=(1:3:(nX*3))+7*nX+f-8; %relative to (nX,nY)


i1=232:270; j1=203:270; %face 3
i2=1:7;     j2=j1;      %face 4
kx = 1:44;


nx=length(i1)+length(i2); ny=length(j1); nz=length(kx);

%Mac270 grid
pp='grid/';
siz270=[nx ny nz];
xc270=readbin([pp 'XC.data'],[nx ny]);
yc270=readbin([pp 'YC.data'],[nx ny]);
xg270=readbin([pp 'XG.data'],[nx ny]);
yg270=readbin([pp 'YG.data'],[nx ny]);
hc270=readbin([pp 'hFacC.data'],[nx ny nz]);
hw270=readbin([pp 'hFacW.data'],[nx ny nz]);
hs270=readbin([pp 'hFacS.data'],[nx ny nz]);
dpt25=-readbin([pp 'RC.data'],nz);

%Mac4320 grid
nx4320=720; ny4320=1080; nz4320=80;
siz4320=[nx4320 ny4320 nz4320];
pp='GRID/';
xc4320=readbin([pp 'XC.data'],[nx4320 ny4320]);
yc4320=readbin([pp 'YC.data'],[nx4320 ny4320]);
xg4320=readbin([pp 'XG.data'],[nx4320 ny4320]);
yg4320=readbin([pp 'YG.data'],[nx4320 ny4320]);
hc4320=readbin([pp 'hFacC.data'],[nx4320 ny4320 nz4320]);
hw4320=readbin([pp 'hFacW.data'],[nx4320 ny4320 nz4320]);
hs4320=readbin([pp 'hFacS.data'],[nx4320 ny4320 nz4320]);
dpt90=-readbin([pp 'RC.data'],nz4320);

pout='run_template/';    % output path name

% generate lateral boundary conditions
        nrecords = [   31 ];
flds={'DIC', 'NO3','NO2','NH4','PO4','FeT','SiO2','DOC','DON', ...
      'DOP','DOFe','POC','PON','POP','POFe','POSi','PIC','ALK','O2', ...
      'c1','c2','c3','c4','c5','c6','c7','Chl1','Chl2','Chl3', ... 
      'Chl4','Chl5'};
flds2={'tr1', 'tr2', 'tr3', 'tr4', 'tr5', 'tr6', 'tr7', 'tr8', 'tr9', 'tr10', ...
       'tr11','tr12','tr13','tr14','tr15','tr16','tr17','tr18','tr19','tr20', ...
       'tr21','tr22','tr23','tr24','tr25','tr26','tr27','tr28','tr29','tr30', ...
       'tr31'};

TS='from270/270_daily/';
tmp=zeros([nx4320 ny4320 nz4320]);
dim=[num2str(nx4320) 'x' num2str(ny4320)];

for f=1:length(flds)
fld=flds{f};
fld2=flds2{f};
fnm=dir([TS fld '.*']);

for i=1:length(fnm), mydisp(i)
  tmp=readbin([TS fnm(i).name],[nx4320 ny4320 nz4320]);

  OBN=tmp(:  ,end,:);
  OBS=tmp(:  ,1,  :);
  OBE=tmp(end,:,  :);
  OBW=tmp(1  ,:,  :);
  if i==1
    writebin([pout 'OBN' fld2 '_' nme '_' dim '_daily.bin'],OBN)
    writebin([pout 'OBS' fld2 '_' nme '_' dim '_daily.bin'],OBS)
    writebin([pout 'OBW' fld2 '_' nme '_' dim '_daily.bin'],OBW)
    writebin([pout 'OBE' fld2 '_' nme '_' dim '_daily.bin'],OBE)
  end
  writebin([pout 'OBN' fld2 '_' nme '_' dim '_daily.bin'],OBN,1,'real*4',i)
  writebin([pout 'OBS' fld2 '_' nme '_' dim '_daily.bin'],OBS,1,'real*4',i)
  writebin([pout 'OBW' fld2 '_' nme '_' dim '_daily.bin'],OBW,1,'real*4',i)
  writebin([pout 'OBE' fld2 '_' nme '_' dim '_daily.bin'],OBE,1,'real*4',i)
end %i
end %fld




