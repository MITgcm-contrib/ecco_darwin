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
%{
TS='from270/';
tmp=zeros([nx ny nz]);
fla=zeros([nx ny]);

flds={'Theta','Salt','U','V'};
flds2={'t',   's',   'u','v'};
kk=[1 2 1 2]; %T/S/U/V record
for f=1:length(flds)
fld=flds{f};
fld2=flds2{f};
if f==1|f==2 %T/S
	fnm=dir([TS 'state_3d_set1*']);
else	     %U/V
	fnm=dir([TS 'trsp_3d_set1*']);
end

for i=1:length(fnm), mydisp(i)
  tmp=zeros([nx ny nz]);
  fla=zeros([nx ny]);
  for k=1:nz
  	fni=[TS fnm(i).name];
	if k==1; disp(fni); end
  	SST=readbin(fni,[nX nY],1,'real*4',k-1 +(kk(f)-1)*nZ);
	st7=SST(:,I3);st8=SST(:,I2);
	fla=[st7(i1,j1); st8(i2,j2)];
	tmp(:,:,k)=fla;
  end	

  OBN=tmp(:  ,end,:);
  if f==4
  OBS=tmp(:  ,2,  :); %V @south
  else
  OBS=tmp(:  ,1,  :);
  end
  OBE=tmp(end,:,  :);
  if f==3
  OBW=tmp(2  ,:,  :); %U @west
  else
  OBW=tmp(1  ,:,  :);
  end
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

%SI
flds={'SIarea', 'SIheff',  'SIhsnow', 'SIuice', 'SIvice'};
flds2={ 'a',      'h',       'sn',     'uice',   'vice'};
kk=[2 3 4 24 25]; %srea/hell/snow/uice/vice record
for f=1:length(flds)
fld=flds{f};
fld2=flds2{f};
fnm=dir([TS 'state_2d_set1*']);
for i=1:length(fnm), mydisp(i)
  tmp=zeros([nx ny]);
  fla=zeros([nx ny]);

  	fni=[TS fnm(i).name];
	disp(fni)
  	SST=readbin(fni,[nX nY],1,'real*4',kk(f) -1);
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
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stabilize T/S
%{
[Y DY]=meshgrid(1:ny,-1.5:-1:-nz);
[X DX]=meshgrid(1:nx,-1.5:-1:-nz);
fin=[regional_grid 'hFacC.data'];
tmp=zeros([nx ny nz]);
tmp=readbin(fin,[nx ny nz]);
maskW=squeeze(tmp(1,:,:));
maskW(find(maskW))=1;  maskW(find(~maskW))=nan;
maskE=squeeze(tmp(end,:,:));
maskE(find(maskE))=1;  maskE(find(~maskE))=nan;
maskN=squeeze(tmp(:,end,:));
maskN(find(maskN))=1;  maskN(find(~maskN))=nan;
maskS=squeeze(tmp(:,1,:));
maskS(find(maskS))=1;  maskS(find(~maskS))=nan;
for t=1:nt+1, mydisp(t)
	if ismember('E',genBC)
  T=readbin([pout 'OBEt_' nme '_' dim '.bin'],[ny nz],1,'real*4',t-1).*maskE;
  S=readbin([pout 'OBEs_' nme '_' dim '.bin'],[ny nz],1,'real*4',t-1).*maskE;
  R=rho(S,T,0);
%  clf, subplot(211), mypcolor(1:ny,-(1:nz),R'); thincolorbar
%  caxis([1024 1028]);thincolorbar
%  tmp=diff(R'); idx=find(tmp<0);
%  hold on, plot(Y(idx),DY(idx),'k.')
  for j=1:ny
    idx=find(diff(R(j,:))<0);
    while ~isempty(idx)
      T(j,min(idx)+1)=T(j,min(idx));
      S(j,min(idx)+1)=S(j,min(idx));
      r=rho(S(j,:),T(j,:),0); idx=find(diff(r)<0);
  end, end
%  R=rho(S,T,0);
%  subplot(212), mypcolor(1:ny,-(1:nz),R'); thincolorbar
%  caxis([1024 1028]);thincolorbar
%  tmp=diff(R'); idx=find(tmp<0);
%  hold on, plot(Y(idx),DY(idx),'k.')
  for k=1:nz
    if any(~isnan(T(:,k)))
      T(:,k)=xpolate(T(:,k)); S(:,k)=xpolate(S(:,k));
    else
      T(:,k)=T(:,k-1); S(:,k)=S(:,k-1);
  end, end
  writebin([pout 'OBEs_' nme '_' dim '.stable'],S,1,'real*4',t-1);
  writebin([pout 'OBEt_' nme '_' dim '.stable'],T,1,'real*4',t-1);
	end %if ismember('E',genBC)
	if ismember('W',genBC)
  T=readbin([pout 'OBWt_' nme '_' dim '.bin'],[ny nz],1,'real*4',t-1).*maskW;
  S=readbin([pout 'OBWs_' nme '_' dim '.bin'],[ny nz],1,'real*4',t-1).*maskW;
  R=rho(S,T,0);
  for j=1:ny
    idx=find(diff(R(j,:))<0);
    while ~isempty(idx)
      T(j,min(idx)+1)=T(j,min(idx));
      S(j,min(idx)+1)=S(j,min(idx));
      r=rho(S(j,:),T(j,:),0); idx=find(diff(r)<0);
  end, end
  for k=1:nz
    if any(~isnan(T(:,k)))
      T(:,k)=xpolate(T(:,k)); S(:,k)=xpolate(S(:,k));
    else
      T(:,k)=T(:,k-1); S(:,k)=S(:,k-1);
  end, end
  writebin([pout 'OBWs_' nme '_' dim '.stable'],S,1,'real*4',t-1);
  writebin([pout 'OBWt_' nme '_' dim '.stable'],T,1,'real*4',t-1);
	end %if ismember('W',genBC)
	if ismember('N',genBC)
  T=readbin([pout 'OBNt_' nme '_' dim '.bin'],[nx nz],1,'real*4',t-1).*maskN;
  S=readbin([pout 'OBNs_' nme '_' dim '.bin'],[nx nz],1,'real*4',t-1).*maskN;
  R=rho(S,T,0);
  for j=1:nx
    idx=find(diff(R(j,:))<0);
    while ~isempty(idx)
      T(j,min(idx)+1)=T(j,min(idx));
      S(j,min(idx)+1)=S(j,min(idx));
      r=rho(S(j,:),T(j,:),0); idx=find(diff(r)<0);
  end, end
  for k=1:nz
    if any(~isnan(T(:,k)))
      T(:,k)=xpolate(T(:,k)); S(:,k)=xpolate(S(:,k));
    else
      T(:,k)=T(:,k-1); S(:,k)=S(:,k-1);
  end, end
  writebin([pout 'OBNs_' nme '_' dim '.stable'],S,1,'real*4',t-1);
  writebin([pout 'OBNt_' nme '_' dim '.stable'],T,1,'real*4',t-1);
	end %if ismember('N',genBC)
	if ismember('S',genBC)
  T=readbin([pout 'OBSt_' nme '_' dim '.bin'],[nx nz],1,'real*4',t-1).*maskS;
  S=readbin([pout 'OBSs_' nme '_' dim '.bin'],[nx nz],1,'real*4',t-1).*maskS;
  R=rho(S,T,0);
  for j=1:nx
    idx=find(diff(R(j,:))<0);
    while ~isempty(idx)
      T(j,min(idx)+1)=T(j,min(idx));
      S(j,min(idx)+1)=S(j,min(idx));
      r=rho(S(j,:),T(j,:),0); idx=find(diff(r)<0);
  end, end
  for k=1:nz
    if any(~isnan(T(:,k)))
      T(:,k)=xpolate(T(:,k)); S(:,k)=xpolate(S(:,k));
    else
      T(:,k)=T(:,k-1); S(:,k)=S(:,k-1);
  end, end
  writebin([pout 'OBSs_' nme '_' dim '.stable'],S,1,'real*4',t-1);
  writebin([pout 'OBSt_' nme '_' dim '.stable'],T,1,'real*4',t-1);
	end %if ismember('S',genBC)
end
clear tmp*
%}

%else %step2
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% balance BCs -- requires that model be integrrated once so that
% grid and mask files be available.
DXG  =readbin([regional_grid 'DXG.data'  ],[nx ny]   );
DYG  =readbin([regional_grid 'DYG.data'  ],[nx ny]   );
HFACS=readbin([regional_grid 'hFacS.data'],[nx ny nz]);
HFACW=readbin([regional_grid 'hFacW.data'],[nx ny nz]);
OBWmask=squeeze(HFACW(2 ,: ,:)); OBWmask([1 ny],:)=0;
OBEmask=squeeze(HFACW(nx,: ,:)); OBEmask([1 ny],:)=0;
OBSmask=squeeze(HFACS(: ,2 ,:)); OBSmask([1 nx],:)=0;
OBNmask=squeeze(HFACS(: ,ny,:)); OBNmask([1 nx],:)=0;
for k=1:nz
  OBWmask(:,k)=thk(k)*OBWmask(:,k).*DYG(2 ,: )';
  OBEmask(:,k)=thk(k)*OBEmask(:,k).*DYG(nx,: )';
  OBSmask(:,k)=thk(k)*OBSmask(:,k).*DXG(: ,2 ) ;
  OBNmask(:,k)=thk(k)*OBNmask(:,k).*DXG(: ,ny) ;
end
OBW=1:nt+1; OBE=1:nt+1; OBS=1:nt+1; OBN=1:nt+1;
for t=1:nt+1, mydisp(t)
  tmp=readbin([pout 'OBEu_' nme '_' dim '.bin'],[ny nz],1,'real*4',t-1);
  OBE(t)=sum(sum(tmp.*OBEmask));
  tmp=readbin([pout 'OBWu_' nme '_' dim '.bin'],[ny nz],1,'real*4',t-1);
  OBW(t)=sum(sum(tmp.*OBWmask));
  tmp=readbin([pout 'OBNv_' nme '_' dim '.bin'],[nx nz],1,'real*4',t-1);
  OBN(t)=sum(sum(tmp.*OBNmask));
  tmp=readbin([pout 'OBSv_' nme '_' dim '.bin'],[nx nz],1,'real*4',t-1);
  OBS(t)=sum(sum(tmp.*OBSmask));
end

n=1;
for t=1:nt+1, mydisp(t)
    switch balanceBC
      case 'W'
        tmp=readbin([pout 'OBWu_' nme '_' dim '.bin'],[ny nz],1,'real*4',t-1);
        tmp(find(OBWmask))=tmp(find(OBWmask))+ ...
            (OBE(t)-OBS(t)+OBN(t)-OBW(t))/sum(sum(OBWmask));
        writebin([pout 'OBWu_' nme '_' dim '.balance'],tmp,1,'real*4',t-1);
      case 'E'
        tmp=readbin([pout 'OBEu_' nme '_' dim '.bin'],[ny nz],1,'real*4',t-1);
        tmp(find(OBEmask))=tmp(find(OBEmask))+ ...
            (OBW(t)-OBN(t)+OBS(t)-OBE(t))/sum(sum(OBEmask));
        writebin([pout 'OBEu_' nme '_' dim '.balance'],tmp,1,'real*4',t-1);
      case 'N'
        tmp=readbin([pout 'OBNv_' nme '_' dim '.bin'],[nx nz],1,'real*4',t-1);
        tmp(find(OBNmask))=tmp(find(OBNmask))+ ...
            (-OBE(t)+OBW(t)+OBS(t)-OBN(t))/sum(sum(OBNmask));
        writebin([pout 'OBNv_' nme '_' dim '.balance'],tmp,1,'real*4',t-1);
      case 'S'
        tmp=readbin([pout 'OBSv_' nme '_' dim '.bin'],[nx nz],1,'real*4',t-1);
        tmp(find(OBSmask))=tmp(find(OBSmask))+ ...
            (OBE(t)-OBW(t)+OBN(t)-OBS(t))/sum(sum(OBSmask));
        writebin([pout 'OBSv_' nme '_' dim '.balance'],tmp,1,'real*4',t-1);
    end
end

OB2=OBW;
for t=1:nt+1, mydisp(t)
  tmp=readbin([pout 'OBWu_' nme '_' dim '.bin'],[ny nz],1,'real*4',t-1);
  OBW(t)=sum(sum(tmp.*OBWmask));
  tmp=readbin([pout 'OBEu_' nme '_' dim '.bin'],[ny nz],1,'real*4',t-1);
  OBE(t)=sum(sum(tmp.*OBEmask));
  tmp=readbin([pout 'OBSv_' nme '_' dim '.bin'],[nx nz],1,'real*4',t-1);
  OBS(t)=sum(sum(tmp.*OBSmask));
  tmp=readbin([pout 'OBNv_' nme '_' dim '.bin'],[nx nz],1,'real*4',t-1);
  OBN(t)=sum(sum(tmp.*OBNmask));
  tmp=readbin([pout 'OBWu_' nme '_' dim '.balance'],[ny nz],1,'real*4',t-1);
  OB2(t)=sum(sum(tmp.*OBWmask));
end
%{
t=1:nt+1;
clf, plot(t,OBW,t,OB2,t,OBS,'linewidth',2), grid
     legend('OBW','OBW2','OBS')
clf, plot(t,cumsum(OBW+OBS),t,cumsum(OB2+OBS),'linewidth',2), grid
     legend('cumsum(OBW+OBS)','cumsum(OBW2+OBS)')
%}
%%

end %step2
