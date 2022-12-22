clear all

nx=46; ny=68; nz=44;
nt=12*29;                           % number of obcs time steps

nme='Mac';
dim='46x68';

% directory names (may need to be created or modified)
pout='/nobackup/hzhang1/pub/Mac_Delta270/run_template/';

regional_grid='../grid/';
genBC={'W','S'};                    % generate genBC boundary conditions
balanceBC='W';                          % balance E, W, N, or S boundary condition

rc=-readbin([regional_grid 'RC.data'],nz);       % depths to center of cell
rf=-readbin([regional_grid 'RF.data'],nz+1);     % depths to cell faces
thk=diff(rf);                              % thicknesses
rac  =readbin([regional_grid 'RAC.data'  ],[nx ny]);
HFACC=readbin([regional_grid 'hFacC.data'],[nx ny]);
AREA=sum(rac(:).*HFACC(:));

%step=1; %generate balance2
step=2; %generate balance3

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
  tmp=readbin([pout 'OBWu_' nme '_' dim '.balance'],[ny nz],1,'real*4',t-1);
  OBW(t)=sum(sum(tmp.*OBWmask));
  tmp=readbin([pout 'OBSv_' nme '_' dim '.bin'],[nx nz],1,'real*4',t-1);
  OBS(t)=sum(sum(tmp.*OBSmask));
end

fn4='stdout.0000_bal'; %STDOUT.0000 from run of "OBWuFile='OBWu_Mac_46x68.balance'," 
fld='dynstat_eta_mean';
val4=mitgcmhistory(fn4,'time_tsnumber',fld);
p=polyfit(val4(:,1)*1200,val4(:,2),1);
p4=p(1); %m/s
adj=p4*AREA/sum(sum(OBWmask));
BAL='balance2';

if step==2
fn4='stdout.0000_bal2'; %STDOUT.0000 from run of "OBWuFile='OBWu_Mac_46x68.balance2',"
fld='dynstat_eta_mean';
val4=mitgcmhistory(fn4,'time_tsnumber',fld);
p=polyfit(val4(:,1)*1200,val4(:,2),1);
p4=p(1); %m/s
adj2=p4*AREA/sum(sum(OBWmask));
BAL='balance3';
adj = adj + adj2;
end

for t=1:nt+1, mydisp(t)
        tmp=readbin([pout 'OBWu_' nme '_' dim '.balance'],[ny nz],1,'real*4',t-1);
        tmp(find(OBWmask))=tmp(find(OBWmask)) - adj;
        writebin([pout 'OBWu_' nme '_' dim '.' BAL],tmp,1,'real*4',t-1);
end


%{
OB2=OBW;
for t=1:nt+1, mydisp(t)
  tmp=readbin([pout 'OBWu_' nme '_' dim '.balance'],[ny nz],1,'real*4',t-1);
  OBW(t)=sum(sum(tmp.*OBWmask));
  tmp=readbin([pout 'OBSv_' nme '_' dim '.bin'],[nx nz],1,'real*4',t-1);
  OBS(t)=sum(sum(tmp.*OBSmask));
  tmp=readbin([pout 'OBWu_' nme '_' dim '.balance2'],[ny nz],1,'real*4',t-1);
  OB2(t)=sum(sum(tmp.*OBWmask));
end

t=1:nt+1;
clf, plot(t,OBW,t,OB2,t,OBS,'linewidth',2), grid
     legend('OBW','OBW2','OBS')
clf, plot(t,cumsum(OBW+OBS),t,cumsum(OB2+OBS),'linewidth',2), grid
     legend('cumsum(OBW+OBS)','cumsum(OBW2+OBS)')
%}

