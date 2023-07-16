% California Coastal System for Carmen's kelp project
% extract monthly-mean ETAN, Salt, Theta, U, V, and PTRACERS
% extract January 1992 (ts=2232) to February 2023 (ts=819504)
% lats 20N to 45.5N, lons 131.5W to 105W
% (example extraction on face 4 and 5, i.e., rotated UV fields)

% This code is best viewed using a "folding" package with the opening
% and closing folds marked by, respectively, "% {{{" and "% }}}".
% For emacs, I use folding.el available here:
% https://github.com/jaalto/project-emacs--folding-mode/blob/master/folding.el

% {{{ define desired region
region_name='CCS_kelp';
minlat=20;
maxlat=45.5;
minlon=-131.5;
maxlon=-105;
nx=270;
prec='real*4';
% }}}

% {{{ extract indices for desired region
pin='/nobackup/dmenemen/llc/llc_270/grid/';
fnm=[pin 'Depth.data'];
[tmp fc ix jx] = ...
    quikread_llc(fnm,nx,1,prec,pin,minlat,maxlat,minlon,maxlon);
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
suf1=['_' int2str(sum(m)) 'x' int2str(n)];
suf2=[suf1 'x' int2str(length(kx))];
% }}}

% {{{ get and save grid information
close all
pout=['/nobackup/dmenemen/' region_name '/grid/'];
eval(['mkdir ' pout])
eval(['cd ' pout])

% {{{ grid cell center
for fnm={'Depth','RAC','XC','YC','hFacC'}
    fin=[pin fnm{1} '.data'];
    switch fnm{1}
      case{'hFacC'}
        fout=[fnm{1} suf2];
        fld=zeros(sum(m),n,length(kx));
        for f=1:length(fc)
            fld((sum(m(1:f))+1):sum(m(1:(f+1))),:,:) = ...
                read_llc_fkij(fin,nx,fc(f),kx,ix{fc(f)},jx{fc(f)});
        end
      otherwise
        fout=[fnm{1} suf1];
        fld=zeros(sum(m),n);
        for f=1:length(fc)
            fld((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
                read_llc_fkij(fin,nx,fc(f),1,ix{fc(f)},jx{fc(f)});
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
fldx=zeros(sum(m),n);
fldy=zeros(sum(m),n);
for f=1:length(fc)
    switch fc(f)
      case {1,2}
        fldx((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finx,nx,fc(f),1,ix{fc(f)},jx{fc(f)});
        fldy((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finy,nx,fc(f),1,ix{fc(f)},jx{fc(f)});
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
fld=zeros(sum(m),n);
for fnm={'XG','YG','RAZ'}
    fin=[pin fnm{1} '.data'];
    fout=[fnm{1} suf1];
    for f=1:length(fc)
        switch fc(f)
          case {1,2}
            fld((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
                read_llc_fkij(fin,nx,fc(f),1,ix{fc(f)},jx{fc(f)});
          case {4,5}
            fld((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
                read_llc_fkij(fin,nx,fc(f),1,ix{fc(f)},jx{fc(f)}-1); % <<<<<<<<
        end
    end
    writebin(fout,fld);
end
% }}}

% {{{ West edge points, no direction
fnx='DXC';
fny='DYC';
finx=[pin fnx '.data'];
finy=[pin fny '.data'];
foutx=[fnx suf1];
fouty=[fny suf1];
fldx=zeros(sum(m),n);
fldy=zeros(sum(m),n);
for f=1:length(fc)
    switch fc(f)
      case {1,2}
        fldx((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finx,nx,fc(f),1,ix{fc(f)},jx{fc(f)});
        fldy((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finy,nx,fc(f),1,ix{fc(f)},jx{fc(f)});
      case {4,5}
        fldx((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finy,nx,fc(f),1,ix{fc(f)},jx{fc(f)});
        fldy((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finx,nx,fc(f),1,ix{fc(f)},jx{fc(f)}-1); % <<<<<<<<
    end
end
writebin(foutx,fldx);
writebin(fouty,fldy);
% }}}

% {{{ Southwest edge points, no direction
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
            read_llc_fkij(finx,nx,fc(f),1,ix{fc(f)},jx{fc(f)});
        fldy((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finy,nx,fc(f),1,ix{fc(f)},jx{fc(f)});
      case {4,5}
        fldx((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finy,nx,fc(f),1,ix{fc(f)},jx{fc(f)}-1); % <<<<<<<<
        fldy((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(finx,nx,fc(f),1,ix{fc(f)},jx{fc(f)});
    end
end
writebin(foutx,fldx);
writebin(fouty,fldy);
% }}}
% }}}

% {{{ get and save regional fields

pin='/nobackup/dcarrol2/v05_latest/darwin3/run/diags/';
pout=['/nobackup/dmenemen/' region_name '/'];
% {{{ get and save scalar 2D fields
fin={'budget/average_2d.'};
fot={'ETAN'};
eval(['mkdir ' pout fot{1}])
eval(['cd ' pout fot{1}])
dnm=dir([pin fin{1} '*.data']);
fld=zeros(sum(m),n);
for t=1:length(dnm)
    fnm=[dnm(t).folder '/' dnm(t).name];
    l=strfind(fnm,'.000');
    ts=str2num(fnm((l+1):(l+10)));
    dy=ts2dte(ts,1200,1992,1,1,30);
    fout=[fot{1} suf1 '.' dy];
    for f=1:length(fc)
        fld((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
            read_llc_fkij(fnm,nx,fc(f),1,ix{fc(f)},jx{fc(f)});
    end
    writebin(fout,fld);
end
% }}}
% {{{ get and save scalar 3D fields
fin={'monthly/'};
fld=zeros(sum(m),n,length(kx));
for fot={'THETA', 'DIC', 'NO3', 'NO2', 'NH4', 'PO4', 'FeT', 'SiO2', ...
         'DOC', 'DON', 'DOP', 'DOFe', 'POC', 'PON', 'POP', 'POFe', ...
         'POSi', 'PIC', 'ALK', 'O2', 'c1', 'c2', 'c3', 'c4', 'c5', ...
         'c6', 'c7', 'Chl1', 'Chl2', 'Chl3', 'Chl4', 'Chl5'}
    eval(['mkdir ' pout fot{1}])
    eval(['cd ' pout fot{1}])
    dnm=dir([pin fin{1} fot{1} '*.data']);
    for t=1:length(dnm)
        fnm=[dnm(t).folder '/' dnm(t).name];
        l=strfind(fnm,'.000');
        ts=str2num(fnm((l+1):(l+10)));
        dy=ts2dte(ts,1200,1992,1,1,30);
        fout=[fot{1} suf2 '.' dy];
        for f=1:length(fc)
            fld((sum(m(1:f))+1):sum(m(1:(f+1))),:,:) = ...
                read_llc_fkij(fnm,nx,fc(f),kx,ix{fc(f)},jx{fc(f)});
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
            read_llc_fkij(fnm,nx,fc(f),kx,ix{fc(f)},jx{fc(f)});
    end
    writebin(fout,fld+35);
end
% }}}
% {{{ get and save vector 3D fields
% Note that zonal velocity is U in faces 1/2 and V in faces 4/5
% and meridional velocity is V in faces 1/2 and -U in faces 4/5;
% but not index shift is needed for V in faces 4/5 because
% EXFuwind and EXfvwind are provided at tracer points.
fin={'budget/average_velmass_3d.'};
eval(['mkdir ' pout 'U'])
eval(['mkdir ' pout 'V'])
eval(['cd ' pout])
dnm=dir([pin fin{1} '*.data']);
fld=zeros(sum(m),n,length(kx));
kxu=kx;    % kx indices for UVELMASS (1st variable in average_velmass_3d)
kxv=kx+50; % kx indices for VVELMASS (2nd variable in average_velmass_3d)
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
                read_llc_fkij(fnm,nx,fc(f),kxu,ix{fc(f)},jx{fc(f)});
            fldv((sum(m(1:f))+1):sum(m(1:(f+1))),:,:) = ...
                read_llc_fkij(fnm,nx,fc(f),kxv,ix{fc(f)},jx{fc(f)});
          case {4,5}
            fldu((sum(m(1:f))+1):sum(m(1:(f+1))),:,:) = ...
                read_llc_fkij(fnm,nx,fc(f),kxv,ix{fc(f)},jx{fc(f)});
            fldv((sum(m(1:f))+1):sum(m(1:(f+1))),:,:) = - ...
                read_llc_fkij(fnm,nx,fc(f),kxu,ix{fc(f)},jx{fc(f)});
        end
        writebin(foutu,fldu);
        writebin(foutv,fldv);
    end
end
% }}}

pin='/nobackup/hzhang1/forcing/era_xx_it42_v2/';
pout=['/nobackup/dmenemen/' region_name '/run_template/era_xx_it42_v2/'];
% {{{ get and save scalar surface forcing
fld=zeros(sum(m),n);
for fot={'aqh', 'atemp', 'lwdn', 'preci', 'swdn', 'uwind', 'vwind'}
    dnm=dir([pin 'EXF' fot{1} '*']);
    for d=1:length(dnm)
        fnm=[dnm(d).folder '/' dnm(d).name];
        fout=[pout region_name suf1 '_' dnm(d).name];
        nt=dnm(d).bytes/nx/nx/13/4;
        for t=1:nt
            for f=1:length(fc)
                fld((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
                    read_llc_fkij(fnm,nx,fc(f),t,ix{fc(f)},jx{fc(f)});
            end
            writebin(fout,fld,1,'real*4',t-1);
        end
    end
end
% }}}
% {{{ get and save mixing coefficients
fld=zeros(sum(m),n,length(kx));
for fot={'diffkr','kapgm','kapredi'}
    fnm=[pin 'llc270_it42_' fot{1} '.data'];
    fout=[pout region_name suf1 '_' fot{1}];
    for f=1:length(fc)
        fld((sum(m(1:f))+1):sum(m(1:(f+1))),:,:) = ...
            read_llc_fkij(fnm,nx,fc(f),kx,ix{fc(f)},jx{fc(f)});
    end
    writebin(fout,fld);
end
% }}}
% {{{ get and save river discharge
fld=zeros(sum(m),n,12);
fnm=['/nobackup/hzhang1/forcing/era-interim/' ...
     'runoff-2d-Fekete-1deg-mon-V4-SMOOTH.bin'];
fout=[pout region_name suf1 '_Fekete_runoff'];
for f=1:length(fc)
    fld((sum(m(1:f))+1):sum(m(1:(f+1))),:,:) = ...
        read_llc_fkij(fnm,nx,fc(f),1:12,ix{fc(f)},jx{fc(f)});
end
writebin(fout,fld);
% }}}
% {{{ get and save iron dust
fld=zeros(sum(m),n,12);
fnm=['/nobackupp19/dmenemen/public/llc_270/ecco_darwin_v5/input/' ...
     'darwin_forcing/llc270_Mahowald_2009_soluble_iron_dust.bin'];
fout=[pout region_name suf1 '_iron_dust'];
for f=1:length(fc)
    fld((sum(m(1:f))+1):sum(m(1:(f+1))),:,:) = ...
        read_llc_fkij(fnm,nx,fc(f),1:12,ix{fc(f)},jx{fc(f)});
end
writebin(fout,fld);
% }}}
% {{{ get and save vector surface forcing
% note that zonal velocity is U in faces 1/2 and V in faces 4/5
% and meridional velocity is V in faces 1/2 and -U in faces 4/5
fldu=zeros(sum(m),n);
fldv=zeros(sum(m),n);
dnmu=dir([pin 'EXFuwind*']);
dnmv=dir([pin 'EXFvwind*']);
for d=1:length(dnmu)
    fnmu=[dnmu(d).folder '/' dnmu(d).name];
    fnmv=[dnmv(d).folder '/' dnmv(d).name];
    foutu=[pout region_name suf1 '_' dnmu(d).name];
    foutv=[pout region_name suf1 '_' dnmv(d).name];
    nt=dnmu(d).bytes/nx/nx/13/4;
    for t=1:nt
        for f=1:length(fc)
            switch fc(f)
              case {1,2}
                fldu((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
                    read_llc_fkij(fnmu,nx,fc(f),t,ix{fc(f)},jx{fc(f)});
                fldv((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
                    read_llc_fkij(fnmv,nx,fc(f),t,ix{fc(f)},jx{fc(f)});
              case {4,5}
                fldu((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
                    read_llc_fkij(fnmv,nx,fc(f),t,ix{fc(f)},jx{fc(f)});
                fldv((sum(m(1:f))+1):sum(m(1:(f+1))),:) = - ...
                    read_llc_fkij(fnmu,nx,fc(f),t,ix{fc(f)},jx{fc(f)}-1);
            end
            writebin(foutu,fldu,1,'real*4',t-1);
            writebin(foutv,fldv,1,'real*4',t-1);
        end
    end
end
% }}}

% }}}
