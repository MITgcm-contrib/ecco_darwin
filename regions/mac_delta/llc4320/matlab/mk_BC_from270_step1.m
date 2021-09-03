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
dpt25=dpt25';


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
dpt90=dpt90';


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

TS='/nobackupp2/dmenemen/public/llc_270/ecco_darwin_v5/output/monthly/';
for f=1:nrecords
fld=flds{f};
fld2=flds2{f};
        fnm=dir([TS fld '/' fld '.*data']);
nt=length(fnm);
for i=1:nt, mydisp(i)
  fn=fnm(i).name;
  dy=datestr(datenum(1992,1,0)+str2num(fn(end-14:end-5))/72,'yyyymmdd');

LLC270_f1=nan([nx ny nz+2]);
LLC270_f2=nan([nx ny nz4320]);
LLC4320_fld=zeros([nx4320 ny4320 nz4320]);

%read in
for k=1:nz
	STATE_fld=readbin([TS fld '/' fn],[nX nY],1,'real*4',k -1);
	xc7=STATE_fld(:,I3); xc8=STATE_fld(:,I2); 
	tmp=[xc7(i1,j1); xc8(i2,j2)];
	ii=find(hc270(:,:,k)==0);
	tmp(ii)=nan;
	LLC270_f1(:,:,k+1)=tmp;
end
	LLC270_f1(:,:, 1)  =LLC270_f1(:,:, 2);
	LLC270_f1(:,:,nz+2)=LLC270_f1(:,:,nz+1);

%vertical first
LLC270_f1=reshape(LLC270_f1,[nx*ny nz+2]);
LLC270_f2=reshape(LLC270_f2,[nx*ny nz4320]);
TL = LLC270_f2*0;
TN = LLC270_f2*0;
dpt52 = [dpt90(1) dpt25 dpt90(end)];
for j=1:size(LLC270_f1,1)
     TL(j,:)=interp1(dpt52,LLC270_f1(j,:)',dpt90,'linear');
     TN(j,:)=interp1(dpt52,LLC270_f1(j,:)',dpt90,'nearest');
end
in = find(isnan(TL)); TL(in) = TN(in);
LLC270_f2 = reshape(TL,[nx ny nz4320]);

%horizontal second
for k=1:nz4320
%source
	tmp=LLC270_f2(:,:,k);
	ii=find(~isnan(tmp));
        xc270_1=xc270(ii);
        yc270_1=yc270(ii);
        tm270_1=tmp(ii);
%target:
%       ii=find(hc4320(:,:,k)==0);
        jj=find(hc4320(:,:,k)~=0);
        xc_1=xc4320(jj); yc_1=yc4320(jj);
       	val_whole=zeros([nx4320 ny4320]);
	if length(xc270_1)>0
	F = scatteredInterpolant(xc270_1,yc270_1,tm270_1);
	val=F(xc_1,yc_1);
	if length(val)>0
	val_whole(jj)=val(:);
	end
	end
	LLC4320_fld(:,:,k)=val_whole;
end %k
fno=['from270/270_daily/' fld '.' dy];
disp(fno)
writebin(fno,LLC4320_fld)

end %nt
end %f

