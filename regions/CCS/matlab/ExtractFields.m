% extract 01-Jan-1992 to 31-Dec-1992
% lats -1 to 7N, lons -8E to 3.5E
% (example extraction on face 1)

% define desired region
nx=270;
prec='real*4';
region_name='GulfGuinea';
minlat=-13.5;
maxlat=6.75;
minlon=-9;
maxlon=14;
mints=dte2ts('02-Jan-1992',1200,1992,1,1);
maxts=dte2ts('31-Dec-1992',1200,1992,1,1);

% extract indices for desired region
gdir='/nobackup/dmenemen/llc/llc_270/grid/';
fnam=[gdir 'Depth.data'];
[fld fc ix jx] = ...
    quikread_llc(fnam,nx,1,prec,gdir,minlat,maxlat,minlon,maxlon);
quikpcolor(fld')
RF=-readbin([gdir 'RF.data'],51);
kx=1:min(find(RF(2:end)>mmax(fld)));

% Get and save grid information
close all
pout=['/nobackup/dmenemen/' region_name '/grid/'];
eval(['mkdir ' pout])
eval(['cd ' pout])
suf1=['_' int2str(length(ix)) 'x' int2str(length(jx))];
suf2=[suf1 'x' int2str(length(kx))];
for fnm={'Depth','DXC','DXG','DYC','DYG','hFacC','hFacS', ...
         'hFacW','RAC','RAS','RAW','RAZ','rLowC','rLowS', ...
         'rLowW','XC','XG','YC','YG'}
    fin=[gdir fnm{1} '.data'];
    switch fnm{1}
      case{'hFacC','hFacS','hFacW'}
        fld=read_llc_fkij(fin,nx,fc,kx,ix,jx);
        fout=[fnm{1} suf2];
      otherwise
        fld=read_llc_fkij(fin,nx,fc,1,ix,jx);
        fout=[fnm{1} suf1];
    end
    writebin(fout,fld);
end

% get and save regional fields
pin='/nobackup/dmenemen/CMS/ecco_darwin/MITgcm/run/';
pout=['/nobackup/dmenemen/' region_name '/'];

% 2D fields
for fnm={'ETAN'}
    eval(['mkdir ' pout fnm{1}])
    eval(['cd ' pout fnm{1}])
    for ts=mints:72:maxts
        fin=[pin fnm{1} '.' myint2str(ts,10) '.data'];
        dy=ts2dte(ts,1200,1992,1,1);
        fout=[fnm{1} '_' int2str(length(ix)) 'x' int2str(length(jx)) '.' dy];
        fld=read_llc_fkij(fin,nx,fc,1,ix,jx);
        writebin(fout,fld);
    end
end

% 3D fields
for fnm={'SALTanom','THETA','UVELMASS','VVELMASS'}
    eval(['mkdir ' pout fnm{1}])
    eval(['cd ' pout fnm{1}])
    for ts=mints:72:maxts
        fin=[pin fnm{1} '.' myint2str(ts,10) '.data'];
        dy=ts2dte(ts,1200,1992,1,1);
        fout=[fnm{1} '_' int2str(length(ix)) 'x' int2str(length(jx)) ...
              'x' int2str(length(kx)) '.' dy];
        for k=1:length(kx);
            fld=read_llc_fkij(fin,nx,fc,kx(k),ix,jx);
            writebin(fout,fld,1,'real*4',k-1);
        end
    end
end
