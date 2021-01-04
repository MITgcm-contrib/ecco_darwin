clear
%adapted from
%https://github.com/MITgcm/gcmfaces/blob/master/gcmfaces_diags/diags_set_MLD.m

dirModel='/nobackup/hzhang1/iter42'
dirMat='';
dirDiag='';

dirModel=[dirModel '/'];
dirMat=[dirModel 'mat/'];
dirDiag=[dirModel 'diags/STATE/']; 

gcmfaces_global; global mygrid
nx=270;ny=nx*13;nz=50;
dirGrid='/nobackup/hzhang1/llc_1080/MITgcm/DM_270/GRID_up/';
nF=5;fileFormat='compact';
grid_load(dirGrid,nF,fileFormat);


  listDiags='fldMldBoyer fldMldSuga fldMldKara';
  listFlds={    'THETA','SALT'};
  listFiles={  'state_3d_set1'};

TT=12*28; %1992-2019 month
fldMldKara=zeros(nx,ny,TT);
fldMldSuga=zeros(nx,ny,TT);
fldMldBoyer=zeros(nx,ny,TT);

for mn=1:TT
disp(mn)

	ts=(datenum(1992,mn+1,1)-datenum(1992,1,1) )*72;
	THETA1=readbin([dirDiag 'state_3d_set1.' myint2str(ts,10) '.data'],[nx ny nz],1,'real*4',1 -1);
	SALT1 =readbin([dirDiag 'state_3d_set1.' myint2str(ts,10) '.data'],[nx ny nz],1,'real*4',2 -1);

	THETA=convert2gcmfaces(THETA1);
	SALT =convert2gcmfaces(SALT1 );

        fldT=THETA.*mygrid.mskC; fldS=SALT.*mygrid.mskC;
        %
        %prepare to compute potential density:
        fldP=0*mygrid.mskC; for kk=1:length(mygrid.RC); fldP(:,:,kk)=-mygrid.RC(kk); end;

        T=convert2vector(fldT);
        S=convert2vector(fldS);
        msk=convert2vector(mygrid.mskC);
        P=convert2vector(fldP);
        %compute potential density:
        RHO=0*msk; alpha=0*msk;
        tmp1=find(~isnan(msk));
        RHO(tmp1) = density(T(tmp1),S(tmp1),P(tmp1));
        fldRhoPot=convert2vector(RHO);
        alpha(tmp1) = density(T(tmp1)+1e-4,S(tmp1),P(tmp1));
        fldAlpha=(convert2vector(alpha)-fldRhoPot)/1e-4;

        clear T S P msk RHO RHOis tmp1;

        %compute mld:
        tmp1=NaN*mygrid.mskC(:,:,1);
        for kk=1:50;
          tmp2=fldRhoPot(:,:,kk)-fldRhoPot(:,:,1);
          %if we pass RHO(1)+0.03 for the first time (or we reach the bottom)
          %then mld is the velocity point above RC(kk), which is RF(kk)
          jj=find((tmp2>0.03|isnan(tmp2))&isnan(tmp1));
         tmp1(jj)=-mygrid.RF(kk);
        end;
        fldMldBoyer(:,:,mn)=convert2gcmfaces(tmp1);

        %compute mld:
        tmp1=NaN*mygrid.mskC(:,:,1);
        for kk=1:50;
          tmp2=fldRhoPot(:,:,kk)-fldRhoPot(:,:,1);
          %if we pass RHO(1)+0.125 for the first time (or we reach the bottom)
          %then mld is the velocity point above RC(kk), which is RF(kk)
          jj=find((tmp2>0.125|isnan(tmp2))&isnan(tmp1));
         tmp1(jj)=-mygrid.RF(kk);
        end;
        fldMldSuga(:,:,mn)=convert2gcmfaces(tmp1);

        %compute mld:
        tmp1=NaN*mygrid.mskC(:,:,1);
        fldRhoPotMax=fldRhoPot(:,:,1)-0.8*fldAlpha(:,:,1);
        for kk=1:50;
          tmp2=fldRhoPot(:,:,kk)-fldRhoPotMax;
          %if we pass RHO(1)+0.8*alpha(1) for the first time (or we reach the bottom)
          %then mld is the velocity point above RC(kk), which is RF(kk)
          jj=find((tmp2>0|isnan(tmp2))&isnan(tmp1));
         tmp1(jj)=-mygrid.RF(kk);
        end;
        fldMldKara(:,:,mn)=convert2gcmfaces(tmp1);

end %mn
	
	fn='MldBoyer.bin';
	writebin([dirMat fn],fldMldBoyer)
	fn='MldSuga.bin';
	writebin([dirMat fn],fldMldSuga)
	fn='MldKara.bin';
	writebin([dirMat fn],fldMldKara)
