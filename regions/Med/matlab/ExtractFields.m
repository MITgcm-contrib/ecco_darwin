% Mediteranean Sea for Louisa Giannoudi and Aleka Pavlidou, HCMR
% lats 30.2N to 47.3N, lons -6E to 42E
% (example extraction on face 1)

% This code is best viewed using a "folding" package with the opening
% and closing folds marked by, respectively, "% {{{" and "% }}}".
% For emacs, I use folding.el available here:
% https://github.com/jaalto/project-emacs--folding-mode/blob/master/folding.el

% {{{ define desired region
region_name='Med';
minlat=30.2;
maxlat=47.3;
minlon=-6;
maxlon=42;
NX=270;
prec='real*4';
% }}}

% {{{ extract indices for desired region
pin='/nobackup/dmenemen/llc/llc_270/grid/';
fnm=[pin 'Depth.data'];
[fld fc ix jx] = ...
    quikread_llc(fnm,NX,1,prec,pin,minlat,maxlat,minlon,maxlon);
quikpcolor(fld'); caxis([0 1])
RF=-readbin([pin 'RF.data'],51);
kx=1:min(find(RF(2:end)>mmax(fld)));
[nx ny]=size(fld); nz=length(kx);
suf1=['_' int2str(nx) 'x' int2str(ny)];
suf2=[suf1 'x' int2str(nz)];
% }}}

% {{{ get and save grid information
close all
pout=['/nobackup/dmenemen/ecco_darwin/' region_name '/grid/'];
eval(['mkdir ' pout])
eval(['cd ' pout])
for fnm={'Depth','RAC','XC','YC','hFacC','hFacS','hFacW','AngleCS', ...
         'AngleSN','XG','YG','RAZ','DXC','DYC','DXG','DYG'}
    fin=[pin fnm{1} '.data'];
    switch fnm{1}
      case{'hFacC','hFacS','hFacW'}
        fout=[fnm{1} suf2];
        fld=read_llc_fkij(fin,NX,fc,kx,ix,jx);
      otherwise
        fout=[fnm{1} suf1];
        fld=read_llc_fkij(fin,NX,fc,1,ix,jx);
    end
    writebin(fout,fld);
end
for fnm={'RC','DRF'}
    fin=[pin fnm{1} '.data'];
    fout=[fnm{1} '_' int2str(nz)];
    fld=readbin(fin,nz);
    writebin(fout,fld);
end
% }}}

% {{{ get and save regional fields

pin='/nobackup/dcarrol2/v05_latest/darwin3/run/diags/';
pout=['/nobackup/dmenemen/ecco_darwin/' region_name '/'];
% {{{ get and save scalar 2D fields
fin={'budget/average_2d.'};
fot={'ETAN'};
eval(['mkdir ' pout fot{1}])
eval(['cd ' pout fot{1}])
dnm=dir([pin fin{1} '*.data']);
for t=1:length(dnm)
    fnm=[dnm(t).folder '/' dnm(t).name];
    l=strfind(fnm,'.000');
    ts=str2num(fnm((l+1):(l+10)));
    dy=ts2dte(ts,1200,1992,1,1,30);
    fout=[fot{1} suf1 '.' dy];
    fld=read_llc_fkij(fnm,NX,fc,1,ix,jx);
    writebin(fout,fld);
end
% }}}
% {{{ get and save scalar 3D fields
fin={'monthly/'};
for fot={'THETA', 'DIC', 'NO3', 'NO2', 'NH4', 'PO4', 'FeT', 'SiO2', ...
         'DOC', 'DON', 'DOP', 'DOFe', 'POC', 'PON', 'POP', 'POFe', ...
         'POSi', 'PIC', 'ALK', 'O2', 'c1', 'c2', 'c3', 'c4', 'c5', ...
         'c6', 'c7', 'Chl1', 'Chl2', 'Chl3', 'Chl4', 'Chl5'}
    eval(['mkdir ' pout fot{1}])
    eval(['cd ' pout fot{1}])
    dnm=dir([pin fin{1} fot{1} '.*.data']);
    for t=1:length(dnm)
        fnm=[dnm(t).folder '/' dnm(t).name];
        l=strfind(fnm,'.000');
        ts=str2num(fnm((l+1):(l+10)));
        dy=ts2dte(ts,1200,1992,1,1,30);
        fout=[fot{1} suf2 '.' dy];
        fld=read_llc_fkij(fnm,NX,fc,kx,ix,jx);
        writebin(fout,fld);
    end
end
eval(['mkdir ' pout 'SALT'])
eval(['cd ' pout 'SALT'])
dnm=dir([pin fin{1} 'SALTanom*.data']);
for t=1:length(dnm)
    fnm=[dnm(t).folder '/' dnm(t).name];
    l=strfind(fnm,'.000');
    ts=str2num(fnm((l+1):(l+10)));
    dy=ts2dte(ts,1200,1992,1,1,30);
    fout=['SALT' suf2 '.' dy];
    fld=read_llc_fkij(fnm,NX,fc,kx,ix,jx);
    writebin(fout,fld+35);
end
% }}}
% {{{ get and save vector 3D fields
% Note that zonal velocity is U in faces 1/2 and V in faces 4/5
% and meridional velocity is V in faces 1/2 and -U in faces 4/5.
fin={'budget/average_velmass_3d.'};
eval(['mkdir ' pout 'U'])
eval(['mkdir ' pout 'V'])
eval(['cd ' pout])
dnm=dir([pin fin{1} '*.data']);
kxu=kx;    % kx indices for UVELMASS (1st variable in average_velmass_3d)
kxv=kx+50; % kx indices for VVELMASS (2nd variable in average_velmass_3d)
for t=1:length(dnm)
    fnm=[dnm(t).folder '/' dnm(t).name];
    l=strfind(fnm,'.000');
    ts=str2num(fnm((l+1):(l+10)));
    dy=ts2dte(ts,1200,1992,1,1,30);
    foutu=['U/U' suf2 '.' dy];
    foutv=['V/V' suf2 '.' dy];
    fldu=read_llc_fkij(fnm,NX,fc,kxu,ix,jx);
    fldv=read_llc_fkij(fnm,NX,fc,kxv,ix,jx);
    writebin(foutu,fldu);
    writebin(foutv,fldv);
end
% }}}

pin='/nobackup/hzhang1/forcing/era_xx_it42_v2/';
pout=['/nobackup/dmenemen/ecco_darwin/' region_name '/run_template/'];
eval(['mkdir ' pout])
pout=['/nobackup/dmenemen/ecco_darwin/' region_name '/run_template/era_xx_it42_v2/'];
eval(['mkdir ' pout])
% {{{ get and save scalar surface forcing
for fot={'aqh', 'atemp', 'lwdn', 'preci', 'swdn', 'uwind', 'vwind'}
    dnm=dir([pin 'EXF' fot{1} '*']);
    for d=1:length(dnm)
        fnm=[dnm(d).folder '/' dnm(d).name];
        fout=[pout region_name suf1 '_' dnm(d).name];
        nt=dnm(d).bytes/NX/NX/13/4;
        for t=1:nt
            fld=read_llc_fkij(fnm,NX,fc,t,ix,jx);
            writebin(fout,fld,1,'real*4',t-1);
        end
    end
end
% }}}
% {{{ get and save mixing coefficients
for fot={'diffkr','kapgm','kapredi'}
    fnm=[pin 'llc270_it42_' fot{1} '.data'];
    fout=[pout region_name suf1 '_' fot{1}];
    fld=read_llc_fkij(fnm,NX,fc,kx,ix,jx);
    writebin(fout,fld);
end
% }}}
% {{{ get and save river discharge
fnm=['/nobackup/hzhang1/forcing/era-interim/' ...
     'runoff-2d-Fekete-1deg-mon-V4-SMOOTH.bin'];
fout=[pout region_name suf1 '_Fekete_runoff'];
fld=read_llc_fkij(fnm,NX,fc,1:12,ix,jx);
writebin(fout,fld);
% }}}
% {{{ get and save iron dust
fnm=['/nobackupp19/dmenemen/public/llc_270/ecco_darwin_v5/input/' ...
     'darwin_forcing/llc270_Mahowald_2009_soluble_iron_dust.bin'];
fout=[pout region_name suf1 '_iron_dust'];
fld=read_llc_fkij(fnm,NX,fc,1:12,ix,jx);
writebin(fout,fld);
% }}}

% }}}
