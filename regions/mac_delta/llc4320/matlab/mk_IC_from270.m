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
pp='grid/'; %from /nobackupp18/nobackupp8/hzhang1/pub/Mac_Delta270/grid/
siz270=[nx ny nz];
xc270=readbin([pp 'XC.data'],[nx ny]);
yc270=readbin([pp 'YC.data'],[nx ny]);
xg270=readbin([pp 'XG.data'],[nx ny]);
yg270=readbin([pp 'YG.data'],[nx ny]);
hc270=readbin([pp 'hFacC.data'],[nx ny nz]);
hw270=readbin([pp 'hFacW.data'],[nx ny nz]);
hs270=readbin([pp 'hFacS.data'],[nx ny nz]);

%Mac4320 grid
nx4320=720; ny4320=1080; nz4320=80;
siz4320=[nx4320 ny4320 nz4320];
pp='GRID/'; %from /nobackupp18/nobackupp8/hzhang1/pub/Mac_Delta/GRID/
xc4320=readbin([pp 'XC.data'],[nx4320 ny4320]);
yc4320=readbin([pp 'YC.data'],[nx4320 ny4320]);
xg4320=readbin([pp 'XG.data'],[nx4320 ny4320]);
yg4320=readbin([pp 'YG.data'],[nx4320 ny4320]);
hc4320=readbin([pp 'hFacC.data'],[nx4320 ny4320 nz4320]);
hw4320=readbin([pp 'hFacW.data'],[nx4320 ny4320 nz4320]);
hs4320=readbin([pp 'hFacS.data'],[nx4320 ny4320 nz4320]);


fn = '/nobackupp2/dmenemen/public/llc_270/ecco_darwin_v5/input/darwin_initial_conditions/pickup_ptracers.0000000001.data';
fout='pTr01_720x1080_Mac';

LLC4320_fld=zeros([nx4320 ny4320]);
k=1;
disp(fn)
	STATE_fld=readbin(fn,[nX nY],1,'real*8',k -1);
%source
	xc7=STATE_fld(:,I3); xc8=STATE_fld(:,I2); 
	tmp=[xc7(i1,j1); xc8(i2,j2)];
%	ii=find(tmp~=0);
	ii=find(hc270(:,:,1)>0);
        xc270_1=xc270(ii);
        yc270_1=yc270(ii);
        tm270_1=tmp(ii);
%target:
        ii=find(hc4320(:,:,k)==0);
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
	LLC4320_fld=val_whole;

disp(fout)
writebin(fout,LLC4320_fld);

