clear all

%global mygrid
%load_llc270

%global llc270
nx=270; ny=nx*13; nz=44;
dirGrid='/nobackup/hzhang1/pub/llc270/GRID/';
xc270=readbin([dirGrid 'XC.data'],[nx*ny]);
yc270=readbin([dirGrid 'YC.data'],[nx*ny]);
xg270=readbin([dirGrid 'XG.data'],[nx*ny]);
yg270=readbin([dirGrid 'YG.data'],[nx*ny]);
hc270=readbin([dirGrid 'hFacC.data'],[nx*ny nz]);
hs270=readbin([dirGrid 'hFacS.data'],[nx*ny nz]);
hw270=readbin([dirGrid 'hFacW.data'],[nx*ny nz]);
xc270(xc270<0)=xc270(xc270<0)+360;
ic270=find(hc270==0);

%rc=-readbin([dirGrid 'RC.data'],nZ);       % depths to center of cell
%rf=-readbin([dirGrid 'RF.data'],nZ+1);     % depths to cell faces
%thk=diff(rf);                              % thicknesses

%mac_latlon
nx4320=1224; ny4320=744; nz4320=81;
pp='grid/';
xc4320=readbin([pp 'XC.data'],[nx4320 ny4320]);
yc4320=readbin([pp 'YC.data'],[nx4320 ny4320]);
xg4320=readbin([pp 'XG.data'],[nx4320 ny4320]);
yg4320=readbin([pp 'YG.data'],[nx4320 ny4320]);
hc4320=readbin([pp 'hFacC.data'],[nx4320 ny4320 nz4320]);
hw4320=readbin([pp 'hFacW.data'],[nx4320 ny4320 nz4320]);
hs4320=readbin([pp 'hFacS.data'],[nx4320 ny4320 nz4320]);

%llc270 limit to mac
lon1=min(xc4320(:))-.2; lon2=max(xc4320(:))+.2;
lat1=min(yc4320(:))-.2; lat2=max(yc4320(:))+.2;
ix=find(xc270>lon1 & xc270<lon2 & yc270>lat1 & yc270<lat2);
IX=length(ix);
xc270=xc270(ix); yc270=yc270(ix);


%vertical
load thk90.mat
GRID_25
dpt90=dpt90(1:nz4320);
dpt25=dpt25(1:nz);
dpt52 = [dpt90(1) dpt25 dpt90(end)];


% directory names (may need to be created or modified)
pin='/nobackup/dcarrol2/v05_latest/darwin3/run/diags/monthly/';
pout='run_template/';    % output path name
nme='Mac';


%Ptracer
	nrecords = [   31 ];
flds={'DIC', 'NO3','NO2','NH4','PO4','FeT','SiO2','DOC','DON', ...
      'DOP','DOFe','POC','PON','POP','POFe','POSi','PIC','ALK','O2', ...
      'c1','c2','c3','c4','c5','c6','c7','Chl1','Chl2','Chl3', ... 
      'Chl4','Chl5'};
flds2={'tr1', 'tr2', 'tr3', 'tr4', 'tr5', 'tr6', 'tr7', 'tr8', 'tr9', 'tr10', ...
       'tr11','tr12','tr13','tr14','tr15','tr16','tr17','tr18','tr19','tr20', ...
       'tr21','tr22','tr23','tr24','tr25','tr26','tr27','tr28','tr29','tr30', ...
       'tr31'};

for f=1:nrecords
fld=flds{f};
fld2=flds2{f};
        fnm=dir([pin fld '.*data']);
nt=length(fnm);
for i=1:nt, mydisp(i)
  fn=fnm(i).name;
  dy=datestr(datenum(1992,1,0)+str2num(fn(end-14:end-5))/72,'yyyymmdd');

 LLC270_f0=nan([nx*ny nz  ]);
 LLC270_f1=nan([IX    nz+2]);
 LLC270_f2=nan([IX    nz4320]);
 LLC4320_fld=zeros([nx4320 ny4320 nz4320]);

%read in
 LLC270_f0=readbin([pin fn],[nx*ny nz]);
 LLC270_f0(ic270)=nan;

 LLC270_f1(:, 2:end-1) = LLC270_f0(ix,:);
 LLC270_f1(:, 1      ) = LLC270_f0(ix,1);
 LLC270_f1(:,   end  ) = LLC270_f0(ix,end);
%vertical first
TL = LLC270_f2*0; TN = LLC270_f2*0;
for j=1:IX
     TL(j,:)=interp1(dpt52,LLC270_f1(j,:)',dpt90,'linear');
%    TN(j,:)=interp1(dpt52,LLC270_f1(j,:)',dpt90,'nearest');
end
%in = find(isnan(TL)); TL(in) = TN(in);
LLC270_f2 = TL;

%horizontal second
for k=1:nz4320
%source
	tmp=LLC270_f2(:,k);
	ii=find(~isnan(tmp));
        xc270_1=xc270(ii); yc270_1=yc270(ii); tm270_1=tmp(ii);
%target:
%       ii=find(hc4320(:,:,k)==0);
        jj=find(hc4320(:,:,k)~=0);
        xc_1=xc4320(jj); yc_1=yc4320(jj);
       	val_whole=zeros([nx4320 ny4320]);
	F = scatteredInterpolant(xc270_1,yc270_1,tm270_1);
	val=F(xc_1,yc_1);
	val_whole(jj)=val(:);
	LLC4320_fld(:,:,k)=val_whole;
end %k

% generate lateral boundary conditions
  tmp=LLC4320_fld;

  OBN=tmp(:  ,end,:);
  OBS=tmp(:  ,1,  :);
  OBE=tmp(end,:,  :);
  OBW=tmp(1  ,:,  :);

  writebin([pout 'OBN' fld2 '_' nme '.bin'],OBN,1,'real*4',i-1)
  writebin([pout 'OBS' fld2 '_' nme '.bin'],OBS,1,'real*4',i-1)
  writebin([pout 'OBW' fld2 '_' nme '.bin'],OBW,1,'real*4',i-1)
  writebin([pout 'OBE' fld2 '_' nme '.bin'],OBE,1,'real*4',i-1)
end %nt
end %f

