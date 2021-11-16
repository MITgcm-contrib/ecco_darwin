clear all



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

rc=-readbin([pp 'RC.data'],nz4320);       % depths to center of cell
rf=-readbin([pp 'RF.data'],nz4320+1);     % depths to cell faces
thk=diff(rf);       

genBC={'W','E', 'N'};                    % generate genBC boundary conditions
balanceBC='N';                          % balance E, W, N, or S boundary condition
regional_grid='grid/';
nx=1224; ny=744; nz=81;

nme='Mac';
pout='run_template/';    % output path name

% generate lateral boundary conditions

flds={'THETA', 'SALTanom','uVel_C','vVel_C'};
flds2={'t',   's',   'u','v'};
TS='from270/Mac/';

for f=1:length(flds)
fld=flds{f};
fld2=flds2{f};
fnm=dir([TS fld '.*']);

for i=1:length(fnm), mydisp(i)
  tmp=readbin([TS fnm(i).name],[nx4320 ny4320 nz4320]);

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
  writebin([pout 'OBN' fld2 '_' nme '.bin'],OBN,1,'real*4',i-1)
  writebin([pout 'OBS' fld2 '_' nme '.bin'],OBS,1,'real*4',i-1)
  writebin([pout 'OBW' fld2 '_' nme '.bin'],OBW,1,'real*4',i-1)
  writebin([pout 'OBE' fld2 '_' nme '.bin'],OBE,1,'real*4',i-1)
end %i
end %fld

%SI
flds={'SIarea', 'SIheff',  'SIhsnow', 'SIuice', 'SIvice'};
flds2={ 'a',      'h',       'sn',     'uice',   'vice'};
for f=1:length(flds)
fld=flds{f};
fld2=flds2{f};
fnm=dir([TS fld '.*']);
for i=1:length(fnm), mydisp(i)
  tmp=readbin([TS fnm(i).name],[nx4320 ny4320]);
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
  writebin([pout 'OBN' fld2 '_' nme '.bin'],OBN,1,'real*4',i-1)
  writebin([pout 'OBS' fld2 '_' nme '.bin'],OBS,1,'real*4',i-1)
  writebin([pout 'OBW' fld2 '_' nme '.bin'],OBW,1,'real*4',i-1)
  writebin([pout 'OBE' fld2 '_' nme '.bin'],OBE,1,'real*4',i-1)
end %i
end %fld

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stabilize T/S
%%
nt = length(fnm)-1;
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
  T=readbin([pout 'OBEt_' nme '.bin'],[ny nz],1,'real*4',t-1).*maskE;
  S=readbin([pout 'OBEs_' nme '.bin'],[ny nz],1,'real*4',t-1).*maskE;
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
  writebin([pout 'OBEs_' nme '.stable'],S,1,'real*4',t-1);
  writebin([pout 'OBEt_' nme '.stable'],T,1,'real*4',t-1);
	end %if ismember('E',genBC)
	if ismember('W',genBC)
  T=readbin([pout 'OBWt_' nme '.bin'],[ny nz],1,'real*4',t-1).*maskW;
  S=readbin([pout 'OBWs_' nme '.bin'],[ny nz],1,'real*4',t-1).*maskW;
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
  writebin([pout 'OBWs_' nme '.stable'],S,1,'real*4',t-1);
  writebin([pout 'OBWt_' nme '.stable'],T,1,'real*4',t-1);
	end %if ismember('W',genBC)
	if ismember('N',genBC)
  T=readbin([pout 'OBNt_' nme '.bin'],[nx nz],1,'real*4',t-1).*maskN;
  S=readbin([pout 'OBNs_' nme '.bin'],[nx nz],1,'real*4',t-1).*maskN;
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
  writebin([pout 'OBNs_' nme '.stable'],S,1,'real*4',t-1);
  writebin([pout 'OBNt_' nme '.stable'],T,1,'real*4',t-1);
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
  writebin([pout 'OBSs_' nme '.stable'],S,1,'real*4',t-1);
  writebin([pout 'OBSt_' nme '.stable'],T,1,'real*4',t-1);
	end %if ismember('S',genBC)
end
clear tmp*
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
  tmp=readbin([pout 'OBEu_' nme '.bin'],[ny nz],1,'real*4',t-1);
  OBE(t)=sum(sum(tmp.*OBEmask));
  tmp=readbin([pout 'OBWu_' nme '.bin'],[ny nz],1,'real*4',t-1);
  OBW(t)=sum(sum(tmp.*OBWmask));
  tmp=readbin([pout 'OBNv_' nme '.bin'],[nx nz],1,'real*4',t-1);
  OBN(t)=sum(sum(tmp.*OBNmask));
  tmp=readbin([pout 'OBSv_' nme '.bin'],[nx nz],1,'real*4',t-1);
  OBS(t)=sum(sum(tmp.*OBSmask));
end

n=1;
for t=1:nt+1, mydisp(t)
    switch balanceBC
      case 'W'
        tmp=readbin([pout 'OBWu_' nme '.bin'],[ny nz],1,'real*4',t-1);
        tmp(find(OBWmask))=tmp(find(OBWmask))+ ...
            (OBE(t)-OBS(t)+OBN(t)-OBW(t))/sum(sum(OBWmask));
        writebin([pout 'OBWu_' nme '.balance'],tmp,1,'real*4',t-1);
      case 'E'
        tmp=readbin([pout 'OBEu_' nme '.bin'],[ny nz],1,'real*4',t-1);
        tmp(find(OBEmask))=tmp(find(OBEmask))+ ...
            (OBW(t)-OBN(t)+OBS(t)-OBE(t))/sum(sum(OBEmask));
        writebin([pout 'OBEu_' nme '.balance'],tmp,1,'real*4',t-1);
      case 'N'
        tmp=readbin([pout 'OBNv_' nme '.bin'],[nx nz],1,'real*4',t-1);
        tmp(find(OBNmask))=tmp(find(OBNmask))+ ...
            (-OBE(t)+OBW(t)+OBS(t)-OBN(t))/sum(sum(OBNmask));
        writebin([pout 'OBNv_' nme '.balance'],tmp,1,'real*4',t-1);
      case 'S'
        tmp=readbin([pout 'OBSv_' nme '.bin'],[nx nz],1,'real*4',t-1);
        tmp(find(OBSmask))=tmp(find(OBSmask))+ ...
            (OBE(t)-OBW(t)+OBN(t)-OBS(t))/sum(sum(OBSmask));
        writebin([pout 'OBSv_' nme '.balance'],tmp,1,'real*4',t-1);
    end
end

OB2=OBN;
for t=1:nt+1, mydisp(t)
  tmp=readbin([pout 'OBWu_' nme '.bin'],[ny nz],1,'real*4',t-1);
  OBW(t)=sum(sum(tmp.*OBWmask));
  tmp=readbin([pout 'OBEu_' nme '.bin'],[ny nz],1,'real*4',t-1);
  OBE(t)=sum(sum(tmp.*OBEmask));
  tmp=readbin([pout 'OBSv_' nme '.bin'],[nx nz],1,'real*4',t-1);
  OBS(t)=sum(sum(tmp.*OBSmask));
  tmp=readbin([pout 'OBNv_' nme '.bin'],[nx nz],1,'real*4',t-1);
  OBN(t)=sum(sum(tmp.*OBNmask));
  tmp=readbin([pout 'OBNv_' nme '.balance'],[nx nz],1,'real*4',t-1);
  OB2(t)=sum(sum(tmp.*OBNmask));
end
%{
t=1:nt+1;
clf, plot(t,OBW,t,OBE,t,OBN,t,OB2,'linewidth',2), grid
     legend('OBW','OBE','OBN','OBN2')
clf, plot(t,cumsum(-OBE+OBW-OBN),t,cumsum(-OBE+OBW-OB2),'linewidth',2), grid
     legend('cumsum(-OBE+OBW-OBN)','cumsum(-OBE+OBW-OB2)')
%}
%%

