% Gulf of Mexico region extracted for Jessica Zaiss on November 15, 2024
% lats 60S to 40S, lons 30E to 100E
% Based on ecco_darwin/v06/1deg/readme_v4r5_v2.txt
% and ecco_darwin/v06/1deg/readme_darwin_v4r5.txt
% (example extraction on face 5, i.e., rotated UV fields)

% This code is best viewed using a "folding" package with the opening
% and closing folds marked by, respectively, "% {{{" and "% }}}".
% For emacs, I use folding.el available here:
% https://github.com/jaalto/project-emacs--folding-mode/blob/master/folding.el

% {{{ define desired region
region_name='GoM';
minlon=-99;
maxlon=-79;
minlat=18;
maxlat=31;
NX=90;
DeltaT=3600;
prec='real*4';
% }}}

% {{{ extract indices for desired region
pin='/nobackup/dcarrol2/grid/ECCO_V4r5/';
fnm=[pin 'Depth.data'];
[fld fc ix jx] = ...
    quikread_llc(fnm,NX,1,prec,pin,minlat,maxlat,minlon,maxlon);
fnm=[pin 'XC.data'];
[xc fc ix jx] = ...
    quikread_llc(fnm,NX,1,prec,pin,minlat,maxlat,minlon,maxlon);
fnm=[pin 'YC.data'];
[yc fc ix jx] = ...
    quikread_llc(fnm,NX,1,prec,pin,minlat,maxlat,minlon,maxlon);
fld(find(~fld))=nan;
clf
pcolorcen(360+xc',yc',fld');
plotland('k-',12)
colorbar('h')
title('depth (m)')

RF=-readbin([pin 'RF.data'],51);
kx=1:min(find(RF(2:end)>mmax(fld)));
[nx ny]=size(fld); nz=length(kx);
suf1=['_' int2str(nx) 'x' int2str(ny)];
suf2=[suf1 'x' int2str(nz)];
close all
% }}}

% {{{ get and save grid information
pout=['/nobackup/dmenemen/ecco_darwin/' region_name '/grid/'];
eval(['mkdir ' pout])
eval(['cd ' pout])
% {{{ grid cell center
for fnm={'Depth','RAC','XC','YC','hFacC'}
    fin=[pin fnm{1} '.data'];
    switch fnm{1}
      case{'hFacC'}
        fld=read_llc_fkij(fin,NX,fc,kx,ix,jx);
        fout=[fnm{1} suf2];
      otherwise
        fld=read_llc_fkij(fin,NX,fc,1,ix,jx);
        fout=[fnm{1} suf1];
    end
    writebin(fout,fld);
end
% }}}
% {{{ Southwest corner (vorticity) points, no direction
fld=zeros(nx,ny);
for fnm={'XG','YG','RAZ'}
    fin=[pin fnm{1} '.data'];
    fout=[fnm{1} suf1];
    for f=1:length(fc)
        switch fc(f)
          case {1,2}
            fld=read_llc_fkij(fin,NX,fc,1,ix,jx);
          case {4,5}
            fld=read_llc_fkij(fin,NX,fc,1,ix,jx-1); % <<<<<<<<
        end
    end
    writebin(fout,fld);
end
% }}}
% {{{ DXC/DYC: West or South edge points, no direction
fnx='DXC';
fny='DYC';
finx=[pin fnx '.data'];
finy=[pin fny '.data'];
foutx=[fnx suf1];
fouty=[fny suf1];
fldx=zeros(nx,ny);
fldy=zeros(nx,ny);
for f=1:length(fc)
    switch fc(f)
      case {1,2}
        fldx=read_llc_fkij(finx,NX,fc,1,ix,jx);
        fldy=read_llc_fkij(finy,NX,fc,1,ix,jx);
      case {4,5}
        fldx=read_llc_fkij(finy,NX,fc,1,ix,jx);
        fldy=read_llc_fkij(finx,NX,fc,1,ix,jx); % <<<<<<<<
    end
end
writebin(foutx,fldx);
writebin(fouty,fldy);
% }}}
% {{{ DXG/DYG: Southwest edge points, no direction
fnx='DXG';
fny='DYG';
finx=[pin fnx '.data'];
finy=[pin fny '.data'];
foutx=[fnx suf1];
fouty=[fny suf1];
for f=1:length(fc)
    switch fc(f)
      case {1,2}
        fldx=read_llc_fkij(finx,NX,fc,1,ix,jx);
        fldy=read_llc_fkij(finy,NX,fc,1,ix,jx);
      case {4,5}
        fldx=read_llc_fkij(finy,NX,fc,1,ix,jx-1); % <<<<<<<<
        fldy=read_llc_fkij(finx,NX,fc,1,ix,jx);
    end
end
writebin(foutx,fldx);
writebin(fouty,fldy);
% }}}
% {{{ RAW/RAS: West or South edge points, no direction
fnx='RAW';
fny='RAS';
finx=[pin fnx '.data'];
finy=[pin fny '.data'];
foutx=[fnx suf1];
fouty=[fny suf1];
fldx=zeros(nx,ny);
fldy=zeros(ny,ny);
for f=1:length(fc)
    switch fc(f)
      case {1,2}
        fldx=read_llc_fkij(finx,NX,fc,1,ix,jx);
        fldy=read_llc_fkij(finy,NX,fc,1,ix,jx);
      case {4,5}
        fldx=read_llc_fkij(finy,NX,fc,1,ix,jx);
        fldy=read_llc_fkij(finx,NX,fc,1,ix,jx-1); % <<<<<<<<
    end
end
writebin(foutx,fldx);
writebin(fouty,fldy);
% }}}
% {{{ hFacW/S: West or South edge points, no direction
fnx='hFacW';
fny='hFacS';
finx=[pin fnx '.data'];
finy=[pin fny '.data'];
foutx=[fnx suf2];
fouty=[fny suf2];
fldx=zeros(ny,ny,nz);
fldy=zeros(ny,ny,nz);
for f=1:length(fc)
    switch fc(f)
      case {1,2}
        fldx=read_llc_fkij(finx,NX,fc,kx,ix,jx);
        fldy=read_llc_fkij(finy,NX,fc,kx,ix,jx);
      case {4,5}
        fldx=read_llc_fkij(finy,NX,fc,kx,ix,jx-1); % <<<<<<<<
        fldy=read_llc_fkij(finx,NX,fc,kx,ix,jx);
    end
end
writebin(foutx,fldx);
writebin(fouty,fldy);
% }}}
% {{{ vertical
for fnm={'RC','DRF'}
    fin=[pin fnm{1} '.data'];
    fout=[fnm{1} '_' int2str(nz)];
    fld=readbin(fin,nz);
    writebin(fout,fld);
end
% }}}
% }}}

% {{{ get and save regional fields
% 1-deg monthly-mean diagnostics from BaseRun10 (our v06 baseline)
pin='/nobackupnfs1/hzhang1/testing_darwin68y/darwin3/BaseRun10/diags/';
pout=['/nobackup/dmenemen/ecco_darwin/' region_name '/'];
% {{{ get and save scalar 2D fields
fin={'diags_state/state_2d_set1.'};
fot={'ETAN'};
eval(['mkdir ' pout fot{1}])
eval(['cd ' pout fot{1}])
dnm=dir([pin fin{1} '*.data']);
for t=1:length(dnm)
    fnm=[dnm(t).folder '/' dnm(t).name];
    l=strfind(fnm,'.000');
    ts=str2num(fnm((l+1):(l+10)));
    dy=ts2dte(ts,DeltaT,1992,1,1,30);
    fout=[fot{1} suf1 '.' dy];
    fld=read_llc_fkij(fnm,NX,fc,1,ix,jx);
    writebin(fout,fld);
end
fin={'monthly/'};
for fot={'apCO2', 'fugfCO2', 'CO2_flux'}
    eval(['mkdir ' pout fot{1}])
    eval(['cd ' pout fot{1}])
    dnm=dir([pin fin{1} fot{1} '.*.data']);
    for t=1:length(dnm)
        fnm=[dnm(t).folder '/' dnm(t).name];
        l=strfind(fnm,'.000');
        ts=str2num(fnm((l+1):(l+10)));
        dy=ts2dte(ts,DeltaT,1992,1,1,30);
        fout=[fot{1} suf1 '.' dy];
        fld=read_llc_fkij(fnm,NX,fc,1,ix,jx);
        writebin(fout,fld);
    end
end
% }}}
% {{{ get and save scalar 3D fields
fin={'monthly/'};
fld=zeros(nx,ny,nz);
for fot={'THETA', 'DIC', 'NO3', 'NO2', 'NH4', 'PO4', 'FeT', 'SiO2', ...
         'DOC', 'DON', 'DOP', 'DOFe', 'POC', 'PON', 'POP', 'POFe', ...
         'POSi', 'PIC', 'ALK', 'O2', 'c1', 'c2', 'c3', 'c4', 'c5', ...
         'c6', 'c7', 'Chl1', 'Chl2', 'Chl3', 'Chl4', 'Chl5', 'pH'}
    eval(['mkdir ' pout fot{1}])
    eval(['cd ' pout fot{1}])
    dnm=dir([pin fin{1} fot{1} '.*.data']);
    for t=1:length(dnm)
        fnm=[dnm(t).folder '/' dnm(t).name];
        l=strfind(fnm,'.000');
        ts=str2num(fnm((l+1):(l+10)));
        dy=ts2dte(ts,DeltaT,1992,1,1,30);
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
    dy=ts2dte(ts,DeltaT,1992,1,1,30);
    fout=['SALT' suf2 '.' dy];
    fld=read_llc_fkij(fnm,NX,fc,kx,ix,jx);
    writebin(fout,fld+35);
end
% }}}
% {{{ get and save vector 3D fields
% Note that zonal velocity is U in faces 1/2 and V in faces 4/5
% and meridional velocity is V in faces 1/2 and -U in faces 4/5.
fin='diags_state/trsp_3d_set1.';
eval(['mkdir ' pout 'U'])
eval(['mkdir ' pout 'V'])
eval(['cd ' pout])
dnm=dir([pin fin '*.data']);
kxu=kx;    % kx indices for UVELMASS (1st variable in trsp_3d_set1)
kxv=kx+50; % kx indices for VVELMASS (2nd variable in trsp_3d_set1)
for t=1:length(dnm)
    fnm=[dnm(t).folder '/' dnm(t).name];
    l=strfind(fnm,'.000');
    ts=str2num(fnm((l+1):(l+10)));
    dy=ts2dte(ts,DeltaT,1992,1,1,30);
    foutu=['U/U' suf2 '.' dy];
    foutv=['V/V' suf2 '.' dy];
    fldu=read_llc_fkij(fnm,NX,fc,kxv,ix,jx);
    fldv=read_llc_fkij(fnm,NX,fc,kxu,ix,jx-1);
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
fld=zeros(nx,ny);
for fot={'aqh', 'atemp', 'lwdn', 'preci', 'swdn'}
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
fld=zeros(nx,ny,nz);
for fot={'diffkr','kapgm','kapredi'}
    fnm=[pin 'llc270_it42_' fot{1} '.data'];
    fout=[pout region_name suf1 '_' fot{1}];
    fld=read_llc_fkij(fnm,NX,fc,kx,ix,jx);
    writebin(fout,fld);
end
% }}}
% {{{ get and save river discharge
fld=zeros(nx,ny,12);
fnm=['/nobackup/hzhang1/forcing/era-interim/' ...
     'runoff-2d-Fekete-1deg-mon-V4-SMOOTH.bin'];
fout=[pout region_name suf1 '_Fekete_runoff'];
fld=read_llc_fkij(fnm,NX,fc,1:12,ix,jx);
writebin(fout,fld);
% }}}
% {{{ get and save iron dust
fnm=['/nobackup/dcarrol2/forcing/iron_dust/LLC_270/' ...
     'llc270_Mahowald_2009_soluble_iron_dust.bin'];
fout=[pout region_name suf1 '_iron_dust'];
fld=read_llc_fkij(fnm,NX,fc,1:12,ix,jx);
writebin(fout,fld);
% }}}
% {{{ get and save vector surface forcing
% Note that zonal velocity is U in faces 1/2 and V in faces 4/5
% and meridional velocity is V in faces 1/2 and -U in faces 4/5;
% but no index shift is needed for V in faces 4/5 because
% EXFuwind and EXfvwind are provided at tracer points.
fldu=zeros(nx,ny);
fldv=zeros(nx,ny);
dnmu=dir([pin 'EXFuwind*']);
dnmv=dir([pin 'EXFvwind*']);
for d=1:length(dnmu)
    fnmu=[dnmu(d).folder '/' dnmu(d).name];
    fnmv=[dnmv(d).folder '/' dnmv(d).name];
    foutu=[pout region_name suf1 '_' dnmu(d).name];
    foutv=[pout region_name suf1 '_' dnmv(d).name];
    nt=dnmu(d).bytes/NX/NX/13/4;
    for t=1:nt
        switch fc(f)
          case {1,2}
            fldu=read_llc_fkij(fnmu,NX,fc,t,ix,jx);
            fldv=read_llc_fkij(fnmv,NX,fc,t,ix,jx);
          case {4,5}
            fldu=read_llc_fkij(fnmv,NX,fc,t,ix,jx);
            fldv=read_llc_fkij(fnmu,NX,fc,t,ix,jx);
        end
        writebin(foutu,fldu,1,'real*4',t-1);
        writebin(foutv,fldv,1,'real*4',t-1);
    end
end
% }}}

% }}}
