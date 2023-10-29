%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build runtime input files needed to run the California
% Coastal System for Carmen's kelp project (CCS_kelp)
% Requires ouput generated by ExtractFields.m
% Generates nominal January 16, 1992 initial conditions
% Boundary conditions span January 16, 1992 to February 13,
% 2023 with obcsperiod=2615438 (~30.27 days)

% This code is best viewed using a "folding" package with the opening
% and closing folds marked by, respectively, "% {{{" and "% }}}".
% For emacs, I use folding.el available here:
% https://github.com/jaalto/project-emacs--folding-mode/blob/master/folding.el

% {{{ define desired region and initialize some variables
region_name='CCS_kelp';
minlat=20;
maxlat=45.5;
minlon=-131.5;
maxlon=-105;
NX=270;
prec='real*4';
gdir= '/nobackup/dmenemen/llc/llc_270/grid/';
pin =['/nobackup/dmenemen/ecco_darwin/' region_name '/'];
pout=['/nobackup/dmenemen/ecco_darwin/' region_name '/run_template/'];
eval(['mkdir ' pout])
% }}}

% {{{ extract indices for desired region
fnm='/nobackupp19/dmenemen/public/llc_270/iter42/input/bathy270_filled_noCaspian_r4';
[tmp fc ix jx] = ...
    quikread_llc(fnm,NX,1,prec,gdir,minlat,maxlat,minlon,maxlon);
m(1)=0;
for f=1:length(fc)
    m(f+1)=length(ix{fc(f)});
end
n=length(jx{fc(1)});
fld=zeros(sum(m),n);
for f=1:length(fc)
    fld((sum(m(1:f))+1):sum(m(1:(f+1))),:)=tmp{fc(f)};
end
quikpcolor(fld');
[nx ny]=size(fld);
RF=-readbin([gdir 'RF.data'],51);
kx=1:min(find(RF(2:end)>max(abs(fld(:)))));
nz=length(kx);
suf1=['_' int2str(sum(m)) 'x' int2str(n)];
suf2=[suf1 'x' int2str(length(kx))];
% }}}

% {{{ Make bathymetry file
close all
writebin([pout 'BATHY' suf1  '_' region_name],fld);
% }}}

% {{{ Make delYFile
fnm=[pin '/grid/YG' suf1];
tmp=readbin(fnm,[nx ny]);
delY=zeros(ny,1);
delY(1:ny-1)=diff(tmp(1,:));
delY(ny)=2*delY(ny-1)-delY(ny-2);
writebin([pout 'delYFile'],delY);
% }}}

% {{{ Generate initial conditions
eval(['!mkdir ' pout 'init'])
eval(['cd ' pout 'init'])
fld={'ETAN'};
sufin=[suf1 '.19920201T000000'];
sufout=[suf1 '.16-Jan-1992'];
eval(['!cp ' pin fld{1} '/' fld{1} sufin ' ' fld{1} sufout])
sufin=[suf2 '.19920201T000000'];
sufout=[suf2 '.16-Jan-1992'];
for fld={'THETA', 'DIC', 'NO3', 'NO2', 'NH4', 'PO4', 'FeT', 'SiO2', ...
         'DOC', 'DON', 'DOP', 'DOFe', 'POC', 'PON', 'POP', 'POFe', ...
         'POSi', 'PIC', 'ALK', 'O2', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', ...
         'c7', 'Chl1', 'Chl2', 'Chl3', 'Chl4', 'Chl5', 'SALT', 'U', 'V'}
    eval(['!cp ' pin fld{1} '/' fld{1} sufin ' fld{1} sufout])
end
% }}}

% {{{ Generate lateral boundary conditions
eval(['!mkdir ' pout 'obcs'])
eval(['cd ' pout 'obcs'])

% {{{ Tracer fields
for fld={'THETA', 'DIC', 'NO3', 'NO2', 'NH4', 'PO4', 'FeT', 'SiO2', ...
         'DOC', 'DON', 'DOP', 'DOFe', 'POC', 'PON', 'POP', 'POFe', ...
         'POSi', 'PIC', 'ALK', 'O2', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', ...
         'c7', 'Chl1', 'Chl2', 'Chl3', 'Chl4', 'Chl5', 'SALT'}
    disp(fld{1})
    pnm=[pin fld{1} '/'];
    fnm=dir([pnm '*' fld{1} '*01T000000']);
    for t=1:length(fnm)
        fin=[pnm fnm(t).name];
        tmp=readbin(fin,[nx ny nz]);
        
        % western boundary condition
        fout=[region_name suf2 '_' fld{1} '_West'];
        writebin(fout,squeeze(tmp(1,:,:)),1,prec,t-1)
        
        % southern boundary condition
        fout=[region_name suf2 '_' fld{1} '_South'];
        writebin(fout,squeeze(tmp(:,1,:)),1,prec,t-1)
        
        % norththern boundary condition
        fout=[region_name suf2 '_' fld{1} '_North'];
        writebin(fout,squeeze(tmp(:,end,:)),1,prec,t-1)
    end
end
% }}}

% {{{ Horizontal velocity
DRF  =readbin([gdir 'DRF.data'],nz);
RAC  =readbin([pin '/grid/RAC'   suf1],[nx ny]);
RAS  =readbin([pin '/grid/RAS'   suf1],[nx ny]);
RAW  =readbin([pin '/grid/RAW'   suf1],[nx ny]);
Depth=readbin([pin '/grid/Depth' suf1],[nx ny]);
hFacS=readbin([pin '/grid/hFacS' suf2],[nx ny nz]);
hFacW=readbin([pin '/grid/hFacW' suf2],[nx ny nz]);
DepthS=0*Depth; DepthW=0*Depth;
for k=1:nz
    DepthS=DepthS+hFacS(:,:,k)*DRF(k);
    DepthW=DepthW+hFacW(:,:,k)*DRF(k);
end
fne=dir([pin 'ETAN/ETAN*01T000000']);

% {{{ U
disp('U')
fnm=dir([pin 'U/U*01T000000']);
for t=1:length(fnm)
    eta=readbin([fne(t).folder '/' fne(t).name],[nx ny]);
    tmp=readbin([fnm(t).folder '/' fnm(t).name],[nx ny nz]);
    
    % western boundary condition
    fout=[region_name suf2 '_U_West'];
    rStarFac=( ( eta(1,:).*RAC(1,:) + eta(2,:).*RAC(2,:) ) ...
               ./ RAW(2,:) / 2 + DepthW(2,:) ) ./ DepthW(2,:);
    hFac=squeeze(hFacW(2,:,:));
    for k=1:nz
        hFac(:,k)=hFac(:,k).*rStarFac';
    end
    obc=squeeze(tmp(2,:,:));
    in=find(hFac~=0&isfinite(hFac));
    obc(in)=obc(in)./hFac(in);
    writebin(fout,obc,1,prec,t-1)
       
    % southern boundary condition
    fout=[region_name suf2 '_U_South'];
    hFac=squeeze(hFacW(:,1,:));
    obc=squeeze(tmp(:,1,:));
    in=find(hFac~=0&isfinite(hFac));
    obc(in)=obc(in)./hFac(in);
    writebin(fout,obc,1,prec,t-1)
       
    % norththern boundary condition
    fout=[region_name suf2 '_U_North'];
    hFac=squeeze(hFacW(:,end,:));
    obc=squeeze(tmp(:,end,:));
    in=find(hFac~=0&isfinite(hFac));
    obc(in)=obc(in)./hFac(in);
    writebin(fout,obc,1,prec,t-1)
end
% }}}

% {{{ V
disp('V')
fnm=dir([pin 'V/V*01T000000']);
for t=1:length(fnm)
    eta=readbin([fne(t).folder '/' fne(t).name],[nx ny]);
    tmp=readbin([fnm(t).folder '/' fnm(t).name],[nx ny nz]);
    
    % western boundary condition
    fout=[region_name suf2 '_V_West'];
    hFac=squeeze(hFacS(1,:,:));
    obc=squeeze(tmp(1,:,:));
    in=find(hFac~=0&isfinite(hFac));
    obc(in)=obc(in)./hFac(in);
    writebin(fout,obc,1,prec,t-1)
       
    % southern boundary condition
    fout=[region_name suf2 '_V_South'];
    rStarFac=( ( eta(:,1).*RAC(:,1) + eta(:,2).*RAC(:,2) ) ...
               ./ RAS(:,2) / 2 + DepthS(:,2) ) ./ DepthS(:,2);
    hFac=squeeze(hFacS(:,2,:));
    for k=1:nz
        hFac(:,k)=hFac(:,k).*rStarFac;
    end
    obc=squeeze(tmp(:,2,:));
    in=find(hFac~=0&isfinite(hFac));
    obc(in)=obc(in)./hFac(in);
    writebin(fout,obc,1,prec,t-1)

    % norththern boundary condition
    fout=[region_name suf2 '_V_North'];
    rStarFac=( ( eta(:,end-1).*RAC(:,end-1) + eta(:,end).*RAC(:,end) ) ...
               ./ RAS(:,end) / 2 + DepthS(:,end) ) ./ DepthS(:,end);
    hFac=squeeze(hFacS(:,end,:));
    for k=1:nz
        hFac(:,k)=hFac(:,k).*rStarFac;
    end
    obc=squeeze(tmp(:,end,:));
    in=find(hFac~=0&isfinite(hFac));
    obc(in)=obc(in)./hFac(in);
    writebin(fout,obc,1,prec,t-1)
end
% }}}

% }}}

% }}}