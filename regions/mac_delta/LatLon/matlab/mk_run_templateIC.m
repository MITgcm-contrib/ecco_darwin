clear all

global mygrid
load_llc270

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
pin='from270/';
pout='run_template/';    % output path name

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make pickup files
% requires that pickup*.meta files be constructed mannually
%STATE
pk='pickup.0000368208.data'; %2006/1/1
flds={'Uvel','Vvel','Theta','Salt','GuNm1','GuNm2','GvNm1','GvNm2','EtaN','dEtaHdt','EtaH'};
fn=[pin pk];
disp(fn)

f=1; %U
 disp(flds{f})
 LLC270_fu=readbin(fn,[nx ny 50],1, 'real*8', f -1);
 LLC270_fu=convert2gcmfaces(LLC270_fu);
f=2; %V
 disp(flds{f})
 LLC270_fv=readbin(fn,[nx ny 50],1, 'real*8', f -1);
 LLC270_fv=convert2gcmfaces(LLC270_fv);
 [UE,VN]=calc_UEVNfromUXVY(LLC270_fu,LLC270_fv);
 UE = convert2gcmfaces(UE); UE=reshape(UE,[nx*ny 50]);
 VN = convert2gcmfaces(VN); VN=reshape(VN,[nx*ny 50]);

for f=1:2 %UV
 LLC270_f0=nan([nx*ny nz  ]);
 LLC270_f1=nan([IX    nz+2]);
 LLC270_f2=nan([IX    nz4320]);
 LLC4320_fld=zeros([nx4320 ny4320 nz4320]);

 if f==1 %U
 LLC270_f0=UE(:,1:nz);
 else
 LLC270_f0=VN(:,1:nz);
 end
 LLC270_f0(ic270)=nan;

 LLC270_f1(:, 2:end-1) = LLC270_f0(ix,:);
 LLC270_f1(:, 1      ) = LLC270_f0(ix,1);
 LLC270_f1(:,   end  ) = LLC270_f0(ix,end);
%vertical first
TL = LLC270_f2*nan; TN = LLC270_f2*nan;
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
        if length(xc270_1)>0
        F = scatteredInterpolant(xc270_1,yc270_1,tm270_1);
        val=F(xc_1,yc_1);
        if length(val)>0
        val_whole(jj)=val(:);
        end
        end
        LLC4320_fld(:,:,k)=val_whole;
end %k
fnout=[pout pk];
disp(fnout)
writebin(fnout,LLC4320_fld,1,'real*8',f -1);
end %f

for f=3:4 %TS
 disp(flds{f})
 LLC270_ff=nan([nx*ny 50  ]);
 LLC270_f0=nan([nx*ny nz  ]);
 LLC270_f1=nan([IX    nz+2]);
 LLC270_f2=nan([IX    nz4320]);
 LLC4320_fld=zeros([nx4320 ny4320 nz4320]);

 LLC270_ff=readbin(fn,[nx*ny 50],1, 'real*8', f -1);
 LLC270_f0=LLC270_ff(:,1:nz);
 LLC270_f0(ic270)=nan;

 LLC270_f1(:, 2:end-1) = LLC270_f0(ix,:);
 LLC270_f1(:, 1      ) = LLC270_f0(ix,1);
 LLC270_f1(:,   end  ) = LLC270_f0(ix,end);
%vertical first
TL = LLC270_f2*nan; TN = LLC270_f2*nan;
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
        if length(xc270_1)>0
        F = scatteredInterpolant(xc270_1,yc270_1,tm270_1);
        val=F(xc_1,yc_1);
        if length(val)>0
        val_whole(jj)=val(:);
        end
        end
        LLC4320_fld(:,:,k)=val_whole;
end %k
fnout=[pout pk];
disp(fnout)
writebin(fnout,LLC4320_fld,1,'real*8',f -1);
end %f

for f=5:8 %ZEROS for GuNm1/GuNm2/GvNm1/GvNm2
 disp(flds{f})
 LLC4320_fld=zeros([nx4320 ny4320 nz4320]);
fnout=[pout pk];
disp(fnout)
writebin(fnout,LLC4320_fld,1,'real*8',f -1);
end %f

k=1;
for f=9:11 %EtaN/dEtaHdt/EtaH
 LLC270_f0=nan([nx*ny 1]);
 LLC270_f1=nan([IX 1   ]);
 LLC270_f2=nan([IX 1   ]);
 LLC4320_fld=zeros([nx4320 ny4320]);

 LLC270_f0=readbin(fn,[nx*ny 1],1, 'real*8', 50*8 +f-8 -1);
 LLC270_f0(hc270(:,k)==0)=nan;
 LLC270_f2 = LLC270_f0(ix);

%source
        tmp=LLC270_f2;
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
        LLC4320_fld=val_whole;
fnout=[pout pk];
disp(fnout)
writebin(fnout,LLC4320_fld,1,'real*8', nz4320*8 +f-8 -1);
end %f

%SEAICE
pk='pickup_seaice.0000368208.data'; %2006/1/1
flds={'siTICE', 'siAREA', 'siHEFF', 'siHSNOW', 'siUICE', 'siVICE'};
fn=[pin pk];
disp(fn)
k=1;
for f=1:4 %TICE/AREA/HEFF/HSNOW
 LLC270_f0=nan([nx*ny 1]);
 LLC270_f1=nan([IX 1   ]);
 LLC270_f2=nan([IX 1   ]);
 LLC4320_fld=zeros([nx4320 ny4320]);

 LLC270_f0=readbin(fn,[nx*ny 1],1, 'real*8', f -1);
 LLC270_f0(hc270(:,k)==0)=nan;
 LLC270_f2 = LLC270_f0(ix);

%source
        tmp=LLC270_f2;
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
        LLC4320_fld=val_whole;
fnout=[pout pk];
disp(fnout)
writebin(fnout,LLC4320_fld,1,'real*8', f -1);
end %f

f=5; %UICE
 disp(flds{f})
 LLC270_fu=readbin(fn,[nx ny],1, 'real*8', f -1);
 LLC270_fu=convert2gcmfaces(LLC270_fu);
f=6; %VICE
 disp(flds{f})
 LLC270_fv=readbin(fn,[nx ny],1, 'real*8', f -1);
 LLC270_fv=convert2gcmfaces(LLC270_fv);
 [UE,VN]=calc_UEVNfromUXVY(LLC270_fu,LLC270_fv);
 UE = convert2gcmfaces(UE); UE=reshape(UE,[nx*ny 1]);
 VN = convert2gcmfaces(VN); VN=reshape(VN,[nx*ny 1]);
for f=5:6 %UVice
 LLC270_f0=nan([nx*ny 1]);
 LLC270_f1=nan([IX 1   ]);
 LLC270_f2=nan([IX 1   ]);
 LLC4320_fld=zeros([nx4320 ny4320]);

 if f==5 %U
 LLC270_f0=UE;
 else
 LLC270_f0=VN;
 end
 LLC270_f0(hc270(:,k)==0)=nan;
 LLC270_f2 = LLC270_f0(ix);
%source
        tmp=LLC270_f2;
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
        LLC4320_fld=val_whole;
fnout=[pout pk];
disp(fnout)
writebin(fnout,LLC4320_fld,1,'real*8', f -1);
end %f


%PTRACER
pk='pickup_ptracers.0000368208.data'; %2006/1/1
for f=1:31
flds{f}=['pTr' myint2str(f)];
end
fn=[pin pk];
disp(fn)
for f=1:31 %PTRACER
 disp(flds{f})
 LLC270_ff=nan([nx*ny 50  ]);
 LLC270_f0=nan([nx*ny nz  ]);
 LLC270_f1=nan([IX    nz+2]);
 LLC270_f2=nan([IX    nz4320]);
 LLC4320_fld=zeros([nx4320 ny4320 nz4320]);

 LLC270_ff=readbin(fn,[nx*ny 50],1, 'real*8', f -1);
 LLC270_f0=LLC270_ff(:,1:nz);
 LLC270_f0(ic270)=nan;

 LLC270_f1(:, 2:end-1) = LLC270_f0(ix,:);
 LLC270_f1(:, 1      ) = LLC270_f0(ix,1);
 LLC270_f1(:,   end  ) = LLC270_f0(ix,end);
%vertical first
TL = LLC270_f2*nan; TN = LLC270_f2*nan;
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
fnout=[pout pk];
disp(fnout)
writebin(fnout,LLC4320_fld,1,'real*8',f -1);
end %f

%DARWIN
pk='pickup_darwin.0000368208.data'; %2006/1/1
flds={'pH'};
fn=[pin pk];
disp(fn)
for f=1 %DARWIN
 disp(flds{f})
 LLC270_ff=nan([nx*ny 50  ]);
 LLC270_f0=nan([nx*ny nz  ]);
 LLC270_f1=nan([IX    nz+2]);
 LLC270_f2=nan([IX    nz4320]);
 LLC4320_fld=zeros([nx4320 ny4320 nz4320]);

 LLC270_ff=readbin(fn,[nx*ny 50],1, 'real*8', f -1);
 LLC270_f0=LLC270_ff(:,1:nz);
 LLC270_f0(ic270)=nan;

 LLC270_f1(:, 2:end-1) = LLC270_f0(ix,:);
 LLC270_f1(:, 1      ) = LLC270_f0(ix,1);
 LLC270_f1(:,   end  ) = LLC270_f0(ix,end);
%vertical first
TL = LLC270_f2*nan; TN = LLC270_f2*nan;
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
fnout=[pout pk];
disp(fnout)
writebin(fnout,LLC4320_fld,1,'real*8',f -1);
end %f


%RUNOFF
pk='runoff-2d-Fekete-1deg-mon-V4-SMOOTH.bin';
fn=[pin pk];
disp(fn)
k=1;
for f=1:12 
 LLC270_f0=nan([nx*ny 1]);
 LLC270_f1=nan([IX 1   ]);
 LLC270_f2=nan([IX 1   ]);
 LLC4320_fld=zeros([nx4320 ny4320]);

 LLC270_f0=readbin(fn,[nx*ny 1],1, 'real*4', f -1);
 LLC270_f0(hc270(:,k)==0)=nan;
 LLC270_f2 = LLC270_f0(ix);

%source
        tmp=LLC270_f2;
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
        LLC4320_fld=val_whole;
fnout=[pout 'runoff-2d-Fekete_Mac.bin'];
disp(fnout)
writebin(fnout,LLC4320_fld,1,'real*4', f -1);
end %f
