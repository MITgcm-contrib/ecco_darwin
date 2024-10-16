% Gulf of Alaska region extracted for Takamitsu Ito on October 31, 2023
% lats 42N to 62N, lons -160E to -120E
% (example extraction on faces 4+5)

% This code is best viewed using a "folding" package with the opening
% and closing folds marked by, respectively, "% {{{" and "% }}}".
% For emacs, I use folding.el available here:
% https://github.com/jaalto/project-emacs--folding-mode/blob/master/folding.el

% {{{ define desired region
region_name='GoA';
minlat=41.3;
maxlat=60.5;
minlon=-160.3;
maxlon=-122.5;
NX=270;
prec='real*4';
% }}}

% {{{ extract indices for desired region
pin='/nobackup/dmenemen/llc/llc_270/grid/';
fnm=[pin 'Depth.data'];
[tmp fc ix jx] = ...
    quikread_llc(fnm,NX,1,prec,pin,minlat,maxlat,minlon,maxlon);
m(1)=0;
for f=1:length(fc)
    m(f+1)=length(ix{fc(f)});
end
n=length(jx{fc(1)});
fld=zeros(sum(m),n);
for f=1:length(fc)
    fld((sum(m(1:f))+1):sum(m(1:(f+1))),:)=tmp{fc(f)};
end
quikpcolor(fld'); caxis([0 1])
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
        fout=[fnm{1} suf2];
        fld=zeros(nx,ny,nz);
        for f=1:length(fc)
            fld((sum(m(1:f))+1):sum(m(1:(f+1))),:,:) = ...
                read_llc_fkij(fin,NX,fc(f),kx,ix{fc(f)},jx{fc(f)});
        end
      otherwise
        fout=[fnm{1} suf1];
        fld=zeros(nx,ny);
        for f=1:length(fc)
            fld((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
                read_llc_fkij(fin,NX,fc(f),1,ix{fc(f)},jx{fc(f)});
        end
    end
    writebin(fout,fld);
end
% }}}

% {{{ AngleCS and AngleSN at grid cell centers
% need to be rotated because U/V vectors are rotated
fnx='AngleCS';
fny='AngleSN';
finx=[pin fnx '.data'];
finy=[pin fny '.data'];
foutx=[fnx suf1];
fouty=[fny suf1];
fldx=zeros(nx,ny);
fldy=zeros(nx,ny);
for f=1:length(fc)
    switch fc(f)
      case {1,2}
        fldx((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finx,NX,fc(f),1,ix{fc(f)},jx{fc(f)});
        fldy((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finy,NX,fc(f),1,ix{fc(f)},jx{fc(f)});
      case {4,5}
        fldx((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            -read_llc_fkij(finy,nx,fc(f),1,ix{fc(f)},jx{fc(f)}); % <<<<<<<<
        fldy((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finx,nx,fc(f),1,ix{fc(f)},jx{fc(f)});
    end
end
writebin(foutx,fldx);
writebin(fouty,fldy);
% }}}

% {{{ Southwest corner (vorticity) points, no direction
fld=zeros(nx,ny);
for fnm={'XG','YG','RAZ'}
    fin=[pin fnm{1} '.data'];
    fout=[fnm{1} suf1];
    for f=1:length(fc)
        switch fc(f)
          case {1,2}
            fld((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
                read_llc_fkij(fin,NX,fc(f),1,ix{fc(f)},jx{fc(f)});
          case {4,5}
            fld((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
                read_llc_fkij(fin,NX,fc(f),1,ix{fc(f)},jx{fc(f)}-1); % <<<<<<<<
        end
    end
    writebin(fout,fld);
end
% }}}

% {{{ DXF/DYF: Grid cell center, no direction
fnx='DXF';
fny='DYF';
finx=[pin fnx '.data'];
finy=[pin fny '.data'];
foutx=[fnx suf1];
fouty=[fny suf1];
fldx=zeros(nx,ny);
fldy=zeros(nx,ny);
for f=1:length(fc)
    switch fc(f)
      case {1,2}
        fldx((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finx,NX,fc(f),1,ix{fc(f)},jx{fc(f)});
        fldy((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finy,NX,fc(f),1,ix{fc(f)},jx{fc(f)});
      case {4,5}
        fldx((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finy,NX,fc(f),1,ix{fc(f)},jx{fc(f)});
        fldy((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finx,NX,fc(f),1,ix{fc(f)},jx{fc(f)});
    end
end
writebin(foutx,fldx);
writebin(fouty,fldy);
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
        fldx((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finx,NX,fc(f),1,ix{fc(f)},jx{fc(f)});
        fldy((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finy,NX,fc(f),1,ix{fc(f)},jx{fc(f)});
      case {4,5}
        fldx((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finy,NX,fc(f),1,ix{fc(f)},jx{fc(f)});
        fldy((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finx,NX,fc(f),1,ix{fc(f)},jx{fc(f)}-1); % <<<<<<<<
    end
end
writebin(foutx,fldx);
writebin(fouty,fldy);
% }}}

% {{{ DXV/DYU: Southwest corner (vorticity) points, no direction
fnx='DXV';
fny='DYU';
finx=[pin fnx '.data'];
finy=[pin fny '.data'];
foutx=[fnx suf1];
fouty=[fny suf1];
for f=1:length(fc)
    switch fc(f)
      case {1,2}
        fldx((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finx,NX,fc(f),1,ix{fc(f)},jx{fc(f)});
        fldy((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finy,NX,fc(f),1,ix{fc(f)},jx{fc(f)});
      case {4,5}
        fldx((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finy,NX,fc(f),1,ix{fc(f)},jx{fc(f)}-1); % <<<<<<<<
        fldy((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finx,NX,fc(f),1,ix{fc(f)},jx{fc(f)}-1); % <<<<<<<<
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
        fldx((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finx,NX,fc(f),1,ix{fc(f)},jx{fc(f)});
        fldy((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finy,NX,fc(f),1,ix{fc(f)},jx{fc(f)});
      case {4,5}
        fldx((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finy,NX,fc(f),1,ix{fc(f)},jx{fc(f)}-1); % <<<<<<<<
        fldy((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finx,NX,fc(f),1,ix{fc(f)},jx{fc(f)});
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
        fldx((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finx,NX,fc(f),1,ix{fc(f)},jx{fc(f)});
        fldy((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finy,NX,fc(f),1,ix{fc(f)},jx{fc(f)});
      case {4,5}
        fldx((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finy,NX,fc(f),1,ix{fc(f)},jx{fc(f)});
        fldy((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finx,NX,fc(f),1,ix{fc(f)},jx{fc(f)}-1); % <<<<<<<<
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
        fldx((sum(m(1:f))+1):sum(m(1:(f+1))),:,:) = ...
            read_llc_fkij(finx,NX,fc(f),kx,ix{fc(f)},jx{fc(f)});
        fldy((sum(m(1:f))+1):sum(m(1:(f+1))),:,:) = ...
            read_llc_fkij(finy,NX,fc(f),kx,ix{fc(f)},jx{fc(f)});
      case {4,5}
        fldx((sum(m(1:f))+1):sum(m(1:(f+1))),:,:) = ...
            read_llc_fkij(finy,NX,fc(f),kx,ix{fc(f)},jx{fc(f)}-1); % <<<<<<<<
        fldy((sum(m(1:f))+1):sum(m(1:(f+1))),:,:) = ...
            read_llc_fkij(finx,NX,fc(f),kx,ix{fc(f)},jx{fc(f)});
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
pin='/nobackup/dcarrol2/v05_latest/darwin3/run/diags/';
pout=['/nobackup/dmenemen/ecco_darwin/' region_name '/'];
% {{{ get and save scalar 2D fields
fin={'budget/average_2d.'};
fot={'ETAN'};
eval(['mkdir ' pout fot{1}])
eval(['cd ' pout fot{1}])
dnm=dir([pin fin{1} '*.data']);
fld=zeros(ny,ny);
for t=1:length(dnm)
    fnm=[dnm(t).folder '/' dnm(t).name];
    l=strfind(fnm,'.000');
    ts=str2num(fnm((l+1):(l+10)));
    dy=ts2dte(ts,1200,1992,1,1,30);
    fout=[fot{1} suf1 '.' dy];
    for f=1:length(fc)
        fld((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(fnm,NX,fc(f),1,ix{fc(f)},jx{fc(f)});
    end
    writebin(fout,fld);
end
fin={'monthly/'};
for fot={'apCO2', 'fugCO2', 'CO2_flux','SIarea','SIheff'}
    eval(['mkdir ' pout fot{1}])
    eval(['cd ' pout fot{1}])
    dnm=dir([pin fin{1} fot{1} '.*.data']);
    for t=1:length(dnm)
        fnm=[dnm(t).folder '/' dnm(t).name];
        l=strfind(fnm,'.000');
        ts=str2num(fnm((l+1):(l+10)));
        dy=ts2dte(ts,1200,1992,1,1,30);
        fout=[fot{1} suf1 '.' dy];
        for f=1:length(fc)
            fld((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
                read_llc_fkij(fnm,NX,fc(f),1,ix{fc(f)},jx{fc(f)});
        end
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
        dy=ts2dte(ts,1200,1992,1,1,30);
        fout=[fot{1} suf2 '.' dy];
        for f=1:length(fc)
            fld((sum(m(1:f))+1):sum(m(1:(f+1))),:,:) = ...
                read_llc_fkij(fnm,NX,fc(f),kx,ix{fc(f)},jx{fc(f)});
        end
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
    for f=1:length(fc)
        fld((sum(m(1:f))+1):sum(m(1:(f+1))),:,:) = ...
            read_llc_fkij(fnm,NX,fc(f),kx,ix{fc(f)},jx{fc(f)});
    end
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
fldu=zeros(nx,ny,nz);
fldv=zeros(nx,ny,nz);
for t=1:length(dnm)
    fnm=[dnm(t).folder '/' dnm(t).name];
    l=strfind(fnm,'.000');
    ts=str2num(fnm((l+1):(l+10)));
    dy=ts2dte(ts,1200,1992,1,1,30);
    foutu=['U/U' suf2 '.' dy];
    foutv=['V/V' suf2 '.' dy];
    for f=1:length(fc)
        switch fc(f)
          case {1,2}
            fldu((sum(m(1:f))+1):sum(m(1:(f+1))),:,:) = ...
                read_llc_fkij(fnm,NX,fc(f),kxu,ix{fc(f)},jx{fc(f)});
            fldv((sum(m(1:f))+1):sum(m(1:(f+1))),:,:) = ...
                read_llc_fkij(fnm,NX,fc(f),kxv,ix{fc(f)},jx{fc(f)});
          case {4,5}
            fldu((sum(m(1:f))+1):sum(m(1:(f+1))),:,:) = ...
                read_llc_fkij(fnm,NX,fc(f),kxv,ix{fc(f)},jx{fc(f)});
            fldv((sum(m(1:f))+1):sum(m(1:(f+1))),:,:) = - ...
                read_llc_fkij(fnm,NX,fc(f),kxu,ix{fc(f)},jx{fc(f)}-1);
        end
    end
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
            for f=1:length(fc)
                fld((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
                    read_llc_fkij(fnm,NX,fc(f),t,ix{fc(f)},jx{fc(f)});
            end
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
    for f=1:length(fc)
        fld((sum(m(1:f))+1):sum(m(1:(f+1))),:,:) = ...
            read_llc_fkij(fnm,NX,fc(f),kx,ix{fc(f)},jx{fc(f)});
    end
    writebin(fout,fld);
end
% }}}
% {{{ get and save river discharge
fld=zeros(nx,ny,12);
fnm=['/nobackup/hzhang1/forcing/era-interim/' ...
     'runoff-2d-Fekete-1deg-mon-V4-SMOOTH.bin'];
fout=[pout region_name suf1 '_Fekete_runoff'];
for f=1:length(fc)
    fld((sum(m(1:f))+1):sum(m(1:(f+1))),:,:) = ...
        read_llc_fkij(fnm,NX,fc(f),1:12,ix{fc(f)},jx{fc(f)});
end
writebin(fout,fld);
% }}}
% {{{ get and save iron dust
fld=zeros(nx,ny,12);
fnm=['/nobackupp19/dmenemen/public/llc_270/ecco_darwin_v5/input/' ...
     'darwin_forcing/llc270_Mahowald_2009_soluble_iron_dust.bin'];
fout=[pout region_name suf1 '_iron_dust'];
for f=1:length(fc)
    fld((sum(m(1:f))+1):sum(m(1:(f+1))),:,:) = ...
        read_llc_fkij(fnm,NX,fc(f),1:12,ix{fc(f)},jx{fc(f)});
end
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
        for f=1:length(fc)
            switch fc(f)
              case {1,2}
                fldu((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
                    read_llc_fkij(fnmu,NX,fc(f),t,ix{fc(f)},jx{fc(f)});
                fldv((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
                    read_llc_fkij(fnmv,NX,fc(f),t,ix{fc(f)},jx{fc(f)});
              case {4,5}
                fldu((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
                    read_llc_fkij(fnmv,NX,fc(f),t,ix{fc(f)},jx{fc(f)});
                fldv((sum(m(1:f))+1):sum(m(1:(f+1))),:) = - ...
                    read_llc_fkij(fnmu,NX,fc(f),t,ix{fc(f)},jx{fc(f)});
            end
            writebin(foutu,fldu,1,'real*4',t-1);
            writebin(foutv,fldv,1,'real*4',t-1);
        end
    end
end
% }}}

% }}}