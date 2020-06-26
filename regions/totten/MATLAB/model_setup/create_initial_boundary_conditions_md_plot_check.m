% Plus and Minus sign for velocity needs to be checked. 
% Current setting is probably worng...
cd /nobackupnfs2/ynakaya2/latlon_totten
addpath /nobackupnfs2/ynakaya2/MITgcm/utils/matlab/cs_grid/latloncap/
addpath /nobackupnfs2/ynakaya2/MITgcm/utils/matlab/cs_grid/read_cs/
addpath /nobackupnfs2/ynakaya2/MITgcm/utils/matlab/cs_grid/
addpath /nobackupnfs2/ynakaya2/MITgcm/utils/matlab/

nme='Tot'; 
timet=1; %select num of month

%grid info
prec='real*8';
gdir='/nobackup/dcarrol2/v05/darwin3/run/';
pout='/nobackupnfs2/ynakaya2/latlon_totten/BC/'
fc=[1 2 3 4 5];
nx=270; 
tz=624; % num of timestep

regx=720; regy=280; regbb=600;
XC=rdmds('/nobackup/ynakaya2/MITgcm/run/totten_high_res/totten_high_16/XC');
YC=rdmds('/nobackup/ynakaya2/MITgcm/run/totten_high_res/totten_high_16/YC');
RC=rdmds('/nobackup/ynakaya2/MITgcm/run/totten_high_res/totten_high_16/RC');

%read grid files
fnamx=['LLC_270_grid/XC.data'];
fnamy=['LLC_270_grid/YC.data'];
 f=2; fldx2(:,:)=read_llc_fkij(fnamx,nx,fc(f),1,1:nx,1:3*nx,'real*4');
 f=4; fldx4(:,:)=read_llc_fkij(fnamx,nx,fc(f),1,1:nx,1:3*nx,'real*4');
 f=2; fldy2(:,:)=read_llc_fkij(fnamy,nx,fc(f),1,1:nx,1:3*nx,'real*4');
 f=4; fldy4(:,:)=read_llc_fkij(fnamy,nx,fc(f),1,1:nx,1:3*nx,'real*4');

flds={'Theta','Salt','U','V'};
flds2={'t',  's',   'u','v'};
kk=[1 2 1 2]; %T/S/U/V record

fnm1=dir([TS 'STATE/state_3d_set1*data']);
fnm2=dir([TS 'TRSP/trsp_3d_set1*data']);

for ii=timet:timet %, mydisp(i)
    fni1=[TS 'STATE/' fnm1(ii).name];
    fni2=[TS 'TRSP/'  fnm2(ii).name];
disp(fni1); disp(fni2);

for k=1:100;
  ff=2; fld2_D1(:,:,k)=read_llc_fkij(fni1,nx,fc(ff),k,1:nx,1:3*nx,'real*4');
end
for k=1:100;
  ff=4; fld4_D1(:,:,k)=read_llc_fkij(fni1,nx,fc(ff),k,1:nx,1:3*nx,'real*4');
end

for k=1:100;
  ff=2; fld2_D2(:,:,k)=read_llc_fkij(fni1,nx,fc(ff),k,1:nx,1:3*nx,'real*4');
end
for k=1:100;
  ff=4; fld4_D2(:,:,k)=read_llc_fkij(fni1,nx,fc(ff),k,1:nx,1:3*nx,'real*4');
end

for i=1:regbb; for j=1:regy; for k=1:50
  T(i,j,k)=fld2_D1(nx2(i),ny2(j),k);
end; end; end;
for i=regbb+1:regx; for j=1:regy; for k=1:50
  T(i,j,k)=fld4_D1(nx4(i),ny4(j),k);
end; end; end;

for i=1:regbb; for j=1:regy; for k=51:100
  S(i,j,k-50)=fld2_D1(nx2(i),ny2(j),k);
end; end; end;
for i=regbb+1:regx; for j=1:regy; for k=51:100
  S(i,j,k-50)=fld4_D1(nx4(i),ny4(j),k);
end; end; end;

for i=1:regbb; for j=1:regy; for k=1:50
  U(i,j,k)=fld2_D(nx2(i),ny2(j),k);
end; end; end;
for i=regbb+1:regx; for j=1:regy; for k=51:100
  U(i,j,k-50)=fld4_D(nx4(i),ny4(j),k);
end; end; end;

for i=1:regbb; for j=1:regy; for k=51:100
  V(i,j,k-50)=fld2_D(nx2(i),ny2(j),k);
end; end; end;
for i=regbb+1:regx; for j=1:regy; for k=1:50
  V(i,j,k)=-fld4_D(nx4(i),ny4(j)-1,k);
end; end; end;

end

%SI
flds={'SIarea', 'SIheff',  'SIhsnow', 'SIuice', 'SIvice'};
flds2={ 'a',      'h',       'sn',     'uice',   'vice'};
kk=[2 3 4 24 25]; %srea/hell/snow/uice/vice record
for f=1:length(flds)
fld=flds{f};
fld2=flds2{f};
fnm=dir([TS 'STATE/state_2d_set1*data']);
for ii=timet:timet ;%, mydisp(i)

fni=[TS 'STATE/' fnm(ii).name];
disp(fni)

for k=1:25;
ff=2; fld2_D(:,:,k)=read_llc_fkij(fni,nx,fc(ff),k,1:nx,1:3*nx,'real*4');
end
for k=1:25;
ff=4; fld4_D(:,:,k)=read_llc_fkij(fni,nx,fc(ff),k,1:nx,1:3*nx,'real*4');
end

if(f==1)
  for i=1:regbb; for j=1:regy;
    data1(i,j)=fld2_D(nx2(i),ny2(j),kk(f));
  end; end;
  for i=regbb+1:regx; for j=1:regy;
    data1(i,j)=fld4_D(nx4(i),ny4(j),kk(f));
  end; end;
elseif(f==2)
  for i=1:regbb; for j=1:regy;
    data2(i,j)=fld2_D(nx2(i),ny2(j),kk(f));
  end; end;
  for i=regbb+1:regx; for j=1:regy;
    data2(i,j)=fld4_D(nx4(i),ny4(j),kk(f)); 
  end; end;
elseif(f==3);
  for i=1:regbb; for j=1:regy; 
     data3(i,j)=fld2_D(nx2(i),ny2(j),kk(f));
  end; end; 
  for i=regbb+1:regx; for j=1:regy; 
     data3(i,j)=fld4_D(nx4(i),ny4(j),kk(f));
  end; end; 
elseif(f==4);
  for i=1:regbb; for j=1:regy; 
     data4(i,j)=fld2_D(nx2(i),ny2(j),24);
  end; end; 
  for i=regbb+1:regx; for j=1:regy; 
     data4(i,j)=fld4_D(nx4(i),ny4(j)-1,25);
  end; end; 
elseif(f==5);
  for i=1:regbb; for j=1:regy;
     data5(i,j)=fld2_D(nx2(i),ny2(j),25);
  end; end;
  for i=regbb+1:regx; for j=1:regy;
     data5(i,j)=-fld4_D(nx4(i),ny4(j)-1,24);
  end; end;
end
end
end

%create boundary condition 
TS='LLC_270_data/';

flds={'DIC', 'NO3','NO2','NH4','PO4','FeT','SiO2','DOC','DON', ...
      'DOP','DOFe','POC','PON','POP','POFe','POSi','PIC','ALK','O2', ...
      'c1','c2','c3','c4','c5','c6','c7','Chl1','Chl2','Chl3', ... 
      'Chl4','Chl5'}
flds2={'tr1','tr2','tr3','tr4','tr5','tr6','tr7','tr8','tr9','tr10', ...
       'tr11','tr12','tr13','tr14','tr15','tr16','tr17','tr18','tr19','tr20', ...
       'tr21','tr22','tr23','tr24','tr25','tr26','tr27','tr28','tr29','tr30', ...
       'tr31'};

clear fld2; clear fld4;
for f=1:length(flds)
fld=flds{f};
fld2=flds2{f};
fnm=dir([gdir flds{f} '.*data']);

for ii=timet:timet;  %, mydisp(i)
    fni=[gdir fnm(ii).name];
    disp(fni);

for k=1:50;
  ff=2; fld2_D(:,:,k)=read_llc_fkij(fni,nx,fc(ff),k,1:nx,1:3*nx,'real*4');
end
for k=1:50;
  ff=4; fld4_D(:,:,k)=read_llc_fkij(fni,nx,fc(ff),k,1:nx,1:3*nx,'real*4');
end

for i=1:regbb; for j=1:regy; for k=1:50
    datatr(i,j,k,f)=fld2_D(nx2(i),ny2(j),k);
end; end; end;
for i=regbb+1:regx; for j=1:regy; for k=1:50
    datatr(i,j,k,f)=fld4_D(nx4(i),ny4(j),k);
end; end; end;

end %ii 
end

for i=1:720;
 for j=1:280
  for k=numz:numz
   if(S(i,j,k)<10); 
    T(i,j,k)=NaN; S(i,j,k)=NaN; U(i,j,k)=NaN; V(i,j,k)=NaN; W(i,j,k)=NaN; datatr(i,j,k,:)=NaN;
   end
  end
 end
end



numx=1; numy=5; 
figure
subplot(numx,numy,1); mypcolor(T(:,:,numz)'); caxis([-2 2]); xticks([]); yticks([]);
subplot(numx,numy,2); mypcolor(S(:,:,numz)'); caxis([34 35]); xticks([]); yticks([]);   
subplot(numx,numy,3); mypcolor(U(:,:,numz)'); caxis([-0.2 0.2]); xticks([]); yticks([]);   
subplot(numx,numy,4); mypcolor(V(:,:,numz)');caxis([-0.2 0.2]); xticks([]); yticks([]);   
subplot(numx,numy,5); mypcolor(W(:,:,numz)');caxis([-0.0002 0.0002]); xticks([]); yticks([]);   

numx=6; numy=6;
figure
for kkk=1:31
 subplot(numx,numy,kkk); mypcolor(datatr(:,:,numz,kkk)');  xticks([]); yticks([]); % caxis([]);
end 

