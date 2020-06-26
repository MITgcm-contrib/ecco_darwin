% Plus and Minus sign for velocity needs to be checked. 
% Current setting is probably worng...
cd /nobackupnfs2/ynakaya2/latlon_totten
addpath /nobackupnfs2/ynakaya2/MITgcm/utils/matlab/cs_grid/latloncap/
addpath /nobackupnfs2/ynakaya2/MITgcm/utils/matlab/cs_grid/read_cs/
addpath /nobackupnfs2/ynakaya2/MITgcm/utils/matlab/cs_grid/
addpath /nobackupnfs2/ynakaya2/MITgcm/utils/matlab/

nme='Tot'; 

%grid info
prec='real*8';
gdir='/nobackup/dcarrol2/v05/darwin3/run/';
pout='/nobackupnfs2/ynakaya2/latlon_totten/BC2/'
fc=[1 2 3 4 5];
nx=270; 
tz=624; % num of timestep

regx=720; regy=280; regbb=624;
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
 
%search for nearest grid %140 is the boundary between f2 and f4
%regbb is needed to be adjusted for different setup
 for i=1:regx
     if(XC(i,1)<141.97)
       [M(i), nx2(i)]=min((XC(i,1)-fldx2(:,180)).^2);
     end
     if(XC(i,1)>141.97)
       [M(i), nx4(i)]=min((XC(i,1)-fldx4(:,180)).^2);
     end
 end
for j=1:regy
	[M(j), ny2(j)]=min((YC(1,j)-fldy2(180,:)).^2);
end
for j=1:regy
    [M(j), ny4(j)]=min((YC(1,j)-fldy4(180,:)).^2);
end

%check the conversion
%{
 for i=1:regbb; for j=1:regy;
	 data(i,j)=fldx2(nx2(i),ny2(j));
 end; end;
 for i=regbb+1:regx; for j=1:regy;
	 data(i,j)=fldx4(nx4(i),ny4(j));
 end; end

 for i=1:regbb; for j=1:regy;
	 data(i,j)=fldy2(nx2(i),ny2(j));
 end; end;
 for i=regbb+1:regx; for j=1:regy;
	 data(i,j)=fldy4(nx4(i),ny4(j));
 end; end
%}
%{

%create initial condition
  fnami=[gdir 'pickup.0000000001.data'];  
  for k=1:403; 
    f=2; fld2(:,:,k)=read_llc_fkij(fnami,nx,fc(f),k,1:nx,1:3*nx,'real*8');
  end
  for k=1:403;
    f=4; fld4(:,:,k)=read_llc_fkij(fnami,nx,fc(f),k,1:nx,1:3*nx,'real*8');
  end

 clear data;
 for i=1:regbb; for j=1:regy; for k=101:150
	 data(i,j,k-100)=fld2(nx2(i),ny2(j),k);
 end; end; end;
 for i=regbb+1:regx; for j=1:regy; for k=101:150
	 data(i,j,k-100)=fld4(nx4(i),ny4(j),k);
 end; end; end;
writebin('T_LLC_270.bin',data)

 clear data;
 for i=1:regbb; for j=1:regy; for k=151:200
	 data(i,j,k-150)=fld2(nx2(i),ny2(j),k);
 end; end; end;
 for i=regbb+1:regx; for j=1:regy; for k=151:200
	 data(i,j,k-150)=fld4(nx4(i),ny4(j),k);
 end; end; end;
writebin('S_LLC_270.bin',data)

 clear data;
 for i=1:regbb; for j=1:regy; for k=1:50
	 data(i,j,k)=fld2(nx2(i),ny2(j),k);
 end; end; end;
 for i=regbb+1:regx; for j=1:regy; for k=51:100
	 data(i,j,k-50)=fld4(nx4(i),ny4(j),k);
 end; end; end;
writebin('U_LLC_270.bin',data)

 clear data;
 for i=1:regbb; for j=1:regy; for k=51:100
	 data(i,j,k-50)=fld2(nx2(i),ny2(j),k);
 end; end; end;
 for i=regbb+1:regx; for j=1:regy; for k=1:50
	 data(i,j,k)=-fld4(nx4(i),ny4(j)-1,k);
 end; end; end;
writebin('V_LLC_270.bin',data)
%}

%{
%create boundary condition 
TS='LLC_270_data/';

flds={'Theta','Salt','U','V'};
flds2={'t',   's',   'u','v'};
kk=[1 2 1 2]; %T/S/U/V record
for f=1:length(flds)
fld=flds{f};
fld2=flds2{f};
if f==1|f==2 %T/S
  fnm=dir([ 'STATE/state_3d_set1*data']);
else %U/V
  fnm=dir([TS 'TRSP/trsp_3d_set1*data']);
end
for ii=1:length(fnm) %, mydisp(i)
    
%    tmp=zeros([nx ny nz]);
%    fla=zeros([nx ny]);

if f==1|f==2
	fni=[TS 'STATE/' fnm(ii).name];
else
    fni=[TS 'TRSP/' fnm(ii).name];
end
disp(fni);


for k=1:100;
  ff=2; fld2_D(:,:,k)=read_llc_fkij(fni,nx,fc(ff),k,1:nx,1:3*nx,'real*4');
end
for k=1:100;
  ff=4; fld4_D(:,:,k)=read_llc_fkij(fni,nx,fc(ff),k,1:nx,1:3*nx,'real*4');
end


clear data;
if(f==1);
for i=1:regbb; for j=1:regy; for k=1:50
  data(i,j,k)=fld2_D(nx2(i),ny2(j),k);
end; end; end;
for i=regbb+1:regx; for j=1:regy; for k=1:50
  data(i,j,k)=fld4_D(nx4(i),ny4(j),k);
end; end; end;
elseif(f==2);
for i=1:regbb; for j=1:regy; for k=51:100
  data(i,j,k-50)=fld2_D(nx2(i),ny2(j),k);
end; end; end;
for i=regbb+1:regx; for j=1:regy; for k=51:100
  data(i,j,k-50)=fld4_D(nx4(i),ny4(j),k);
end; end; end;
elseif(f==3);
for i=1:regbb; for j=1:regy; for k=1:50
  data(i,j,k)=fld2_D(nx2(i),ny2(j),k);
end; end; end;
for i=regbb+1:regx; for j=1:regy; for k=51:100
  data(i,j,k-50)=fld4_D(nx4(i),ny4(j),k);
end; end; end;
elseif(f==4);
for i=1:regbb; for j=1:regy; for k=51:100
  data(i,j,k-50)=fld2_D(nx2(i),ny2(j),k);
end; end; end;
for i=regbb+1:regx; for j=1:regy; for k=1:50
  data(i,j,k)=-fld4_D(nx4(i),ny4(j)-1,k);
end; end; end;
end

if(f==1 || f==2)
    data_n=data;
for iii=1:regx; for jjj=1:regy; for kkk=1:50;
if(data(iii,jjj,kkk)==0); data_n(iii,jjj,kkk)=NaN; end;
end;end;end;

for kkk=1:50
data_n2(1:10,:,kkk)=inpaint_nans_2(data_n(1:10,:,kkk),5);
end
for kkk=1:50
data_n2(711:720,:,kkk)=inpaint_nans_2(data_n(711:720,:,kkk),5);
end
for kkk=1:50
data_n2(1:720,270:280,kkk)=inpaint_nans_2(data_n(1:720,270:280,kkk),5);
end

for iii=1:regx; for jjj=1:regy; for kkk=1:50;
if(data(iii,jjj,kkk)==0); data(iii,jjj,kkk)=data_n2(iii,jjj,kkk); end;
end;end;end;
end

  OBN=data(:  ,end,:);
  OBS=data(:  ,1,  :);
  OBE=data(end,:,  :);
  OBW=data(1  ,:,  :);

  if ii==1
    writebin([pout 'OBN' fld2 '_' nme '_' '.bin'],OBN)
    writebin([pout 'OBS' fld2 '_' nme '_' '.bin'],OBS)
    writebin([pout 'OBW' fld2 '_' nme '_' '.bin'],OBW)
    writebin([pout 'OBE' fld2 '_' nme '_' '.bin'],OBE)
  end
    writebin([pout 'OBN' fld2 '_' nme '_' '.bin'],OBN,1,'real*4',ii)
    writebin([pout 'OBS' fld2 '_' nme '_' '.bin'],OBS,1,'real*4',ii)
    writebin([pout 'OBW' fld2 '_' nme '_' '.bin'],OBW,1,'real*4',ii)
    writebin([pout 'OBE' fld2 '_' nme '_' '.bin'],OBE,1,'real*4',ii)
  end %ii
end %fld


    clear data;
clear OBN; clear OBS;
clear OBE; clear OBW;
%SI
flds={'SIarea', 'SIheff',  'SIhsnow', 'SIuice', 'SIvice'};
flds2={ 'a',      'h',       'sn',     'uice',   'vice'};
kk=[2 3 4 24 25]; %srea/hell/snow/uice/vice record
for f=1:length(flds)
fld=flds{f};
fld2=flds2{f};
fnm=dir([TS 'STATE/state_2d_set1*data']);
for ii=1:length(fnm) ;%, mydisp(i)

fni=[TS 'STATE/' fnm(ii).name];
disp(fni)

for k=1:25;
ff=2; fld2_D(:,:,k)=read_llc_fkij(fni,nx,fc(ff),k,1:nx,1:3*nx,'real*4');
end
for k=1:25;
ff=4; fld4_D(:,:,k)=read_llc_fkij(fni,nx,fc(ff),k,1:nx,1:3*nx,'real*4');
end

clear data;
if(f<3.5);
  for i=1:regbb; for j=1:regy; 
     data(i,j)=fld2_D(nx2(i),ny2(j),kk(f));
  end; end; 
  for i=regbb+1:regx; for j=1:regy; 
     data(i,j)=fld4_D(nx4(i),ny4(j),kk(f));
  end; end; 
elseif(f==4);
  for i=1:regbb; for j=1:regy; 
     data(i,j)=fld2_D(nx2(i),ny2(j),24);
  end; end; 
  for i=regbb+1:regx; for j=1:regy; 
     data(i,j)=fld4_D(nx4(i),ny4(j)-1,25);
  end; end; 
elseif(f==5);
  for i=1:regbb; for j=1:regy;
     data(i,j)=fld2_D(nx2(i),ny2(j),25);
  end; end;
  for i=regbb+1:regx; for j=1:regy;
     data(i,j)=-fld4_D(nx4(i),ny4(j)-1,24);
  end; end;
end


   OBN=data(:  ,end);
   OBS=data(:  ,1);
   OBE=data(end,:);
   OBW=data(1  ,:);

  if ii==1
    writebin([pout 'OBN' fld2 '_' nme '_' '.bin'],OBN)
    writebin([pout 'OBS' fld2 '_' nme '_' '.bin'],OBS)
    writebin([pout 'OBW' fld2 '_' nme '_' '.bin'],OBW)
    writebin([pout 'OBE' fld2 '_' nme '_' '.bin'],OBE)
  end
    writebin([pout 'OBN' fld2 '_' nme '_' '.bin'],OBN,1,'real*4',ii)
    writebin([pout 'OBS' fld2 '_' nme '_' '.bin'],OBS,1,'real*4',ii)
    writebin([pout 'OBW' fld2 '_' nme '_' '.bin'],OBW,1,'real*4',ii)
    writebin([pout 'OBE' fld2 '_' nme '_' '.bin'],OBE,1,'real*4',ii)
end %ii
end %fld
%}

%%next boundary conditions for biogeochemical tracer
  fnami=[gdir 'pickup_ptracers.0000000001.data'];  
  clear fld2; clear fld4;
  for k=1:31*50; 
    f=2; fld2(:,:,k)=read_llc_fkij(fnami,nx,fc(f),k,1:nx,1:3*nx,'real*8');
  end
  for k=1:31*50;
    f=4; fld4(:,:,k)=read_llc_fkij(fnami,nx,fc(f),k,1:nx,1:3*nx,'real*8');
  end
  
 clear data;
 for i=1:regbb; for j=1:regy; for k=1:31*50
	 data(i,j,k)=fld2(nx2(i),ny2(j),k);
 end; end; end;
 for i=regbb+1:regx; for j=1:regy; for k=1:31*50
	 data(i,j,k)=fld4(nx4(i),ny4(j),k);
 end; end; end;

data2=data;
%search for interpolation
for j=regy-1:-1:1;
for i=1:regx; 
for k=1:50*31
if(data(i,j,k) ==0); 
   data(i,j,k)=data(i,j+1,k);
end
end
end
end

for kk=1:31
for k=2:50
for j=regy-1:-1:1;
for i=1:regx; 
if(data(i,j,k+(kk-1)*50) ==0 ); 
   data(i,j,k+(kk-1)*50)=data(i,j,k+(kk-1)*50-1);
end
end
end
end
end    
writebin('pickup_ptracers.0000000001.data_reg_tot',data,1, 'real*8');

%Create pickup file from first timestep
flds={'DIC', 'NO3','NO2','NH4','PO4','FeT','SiO2','DOC','DON', ...
      'DOP','DOFe','POC','PON','POP','POFe','POSi','PIC','ALK','O2', ...
      'c1','c2','c3','c4','c5','c6','c7','Chl1','Chl2','Chl3', ...
      'Chl4','Chl5'}
%Do not know why but wrong order...  1-14 (ok) 15 (POFE), 16-31 (+1)

clear fld2; clear fld4;
for f=1:length(flds)
  fld=flds{f};
  fnm=dir([gdir flds{f} '.*data']);
ii=1;
   fni=[gdir fnm(ii).name];
   disp(fni);
    for k=1:50;
     ff=2; fld2_D(:,:,k)=read_llc_fkij(fni,nx,fc(ff),k,1:nx,1:3*nx,'real*4');
    end
    for k=1:50;
     ff=4; fld4_D(:,:,k)=read_llc_fkij(fni,nx,fc(ff),k,1:nx,1:3*nx,'real*4');
    end

    clear data;
    for i=1:regbb; for j=1:regy; for k=1:50
      data(i,j,k)=fld2_D(nx2(i),ny2(j),k);
    end; end; end;
    for i=regbb+1:regx; for j=1:regy; for k=1:50
      data(i,j,k)=fld4_D(nx4(i),ny4(j),k);
    end; end; end;

    if(f==1 || f==2)
      data_n=data;
     for iii=1:regx; for jjj=1:regy; for kkk=1:50;
       if(data(iii,jjj,kkk)==0); data_n(iii,jjj,kkk)=NaN; end;
     end;end;end;

     for kkk=1:50
	  data_n2(:,:,kkk)=inpaint_nans_2(data_n(:,:,kkk),5);
     end

     for iii=1:regx; for jjj=1:regy; for kkk=1:50;
     if(data(iii,jjj,kkk)==0); data(iii,jjj,kkk)=data_n2(iii,jjj,kkk); end;
     end;end;end;
    end
    datap(:,:,:,f)=data(:,:,:);

end
    writebin('pickup_ptracers.0000000001.data_reg_tot3',datap,1, 'real*8');


clear data; clear OBN; clear OBS; clear OBW; clear OBE; 
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

for ii=1: 6; %length(fnm) %, mydisp(i)
    
%    tmp=zeros([nx ny nz]);
%    fla=zeros([nx ny]);


    fni=[gdir fnm(ii).name];
    disp(fni);

for k=1:50;
  ff=2; fld2_D(:,:,k)=read_llc_fkij(fni,nx,fc(ff),k,1:nx,1:3*nx,'real*4');
end
for k=1:50;
  ff=4; fld4_D(:,:,k)=read_llc_fkij(fni,nx,fc(ff),k,1:nx,1:3*nx,'real*4');
end

clear data;

for i=1:regbb; for j=1:regy; for k=1:50
  data(i,j,k)=fld2_D(nx2(i),ny2(j),k);
end; end; end;
for i=regbb+1:regx; for j=1:regy; for k=1:50
  data(i,j,k)=fld4_D(nx4(i),ny4(j),k);
end; end; end;


if(f==1 || f==2)
    data_n=data;
for iii=1:regx; for jjj=1:regy; for kkk=1:50;
if(data(iii,jjj,kkk)==0); data_n(iii,jjj,kkk)=NaN; end;
end;end;end;

for kkk=1:50
data_n2(1:10,:,kkk)=inpaint_nans_2(data_n(1:10,:,kkk),5);
end
for kkk=1:50
data_n2(711:720,:,kkk)=inpaint_nans_2(data_n(711:720,:,kkk),5);
end
for kkk=1:50
data_n2(1:720,270:280,kkk)=inpaint_nans_2(data_n(1:720,270:280,kkk),5);
end

for iii=1:regx; for jjj=1:regy; for kkk=1:50;
if(data(iii,jjj,kkk)==0); data(iii,jjj,kkk)=data_n2(iii,jjj,kkk); end;
end;end;end;
end

  OBN=data(:  ,end,:);
  OBS=data(:  ,1,  :);
  OBE=data(end,:,  :);
  OBW=data(1  ,:,  :);

  if ii==1
    writebin([pout 'OBN' fld2 '_' nme '.bin'],OBN)
    writebin([pout 'OBS' fld2 '_' nme '.bin'],OBS)
    writebin([pout 'OBW' fld2 '_' nme '.bin'],OBW)
    writebin([pout 'OBE' fld2 '_' nme '.bin'],OBE)
  end
    writebin([pout 'OBN' fld2 '_' nme '.bin'],OBN,1,'real*4',ii)
    writebin([pout 'OBS' fld2 '_' nme '.bin'],OBS,1,'real*4',ii)
    writebin([pout 'OBW' fld2 '_' nme '.bin'],OBW,1,'real*4',ii)
    writebin([pout 'OBE' fld2 '_' nme '.bin'],OBE,1,'real*4',ii)
  end %ii 
end

%convert data

%%next boundary conditions for biogeochemical tracer
    fnami=[gdir 'llc270_Mahowald_2009_soluble_iron_dust.bin'];
clear fld2; clear fld4;
for k=1:12;
f=2; fld2(:,:,k)=read_llc_fkij(fnami,nx,fc(f),k,1:nx,1:3*nx,'real*4');
  end
  for k=1:12;
f=4; fld4(:,:,k)=read_llc_fkij(fnami,nx,fc(f),k,1:nx,1:3*nx,'real*4');
  end

  clear data;
for i=1:regbb; for j=1:regy; for k=1:12
				     data(i,j,k)=fld2(nx2(i),ny2(j),k);
end; end; end;
for i=regbb+1:regx; for j=1:regy; for k=1:12
					  data(i,j,k)=fld4(nx4(i),ny4(j),k);
end; end; end;

writebin('reg_tot_Mahowald_2009_soluble_iron_dust.bin',data);

