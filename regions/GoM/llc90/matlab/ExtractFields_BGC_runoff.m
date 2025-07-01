% Gulf of Mexico region extracted for Jessica Zaiss on November 15, 2024
% lats 60S to 40S, lons 30E to 100E
% Based on ecco_darwin/v06/1deg/readme_v4r5_v2.txt
% and ecco_darwin/v06/1deg/readme_darwin_v4r5.txt
% (example extraction on face 5, i.e., rotated UV fields)

% This code is best viewed using a "folding" package with the opening
% and closing folds marked by, respectively, "% {{{" and "% }}}".
% For emacs, I use folding.el available here:
% https://github.com/jaalto/project-emacs--folding-mode/blob/master/folding.el

addpath(genpath('/nobackup/dcarrol2/MATLAB'));

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
[nx ny]=size(fld);
suf1=['_' int2str(nx) 'x' int2str(ny)];
% }}}

% {{{ get and save regional fields
% 1-deg monthly-mean diagnostics from BaseRun10 (our v06 baseline)
%pout=['/nobackup/dmenemen/ecco_darwin/' region_name '/run_template/'];
pout=['/nobackup/dcarrol2/cutout/GoM/BGC_runoff/'];

eval(['mkdir ' pout])
% {{{ get and save biogeochemical river discharge
pin='/nobackup/rsavelli/LOAC/ECCO_V4r5/bgc_runoff/';

nm={'DOC';'DON';'DOP';'DIN';'DIP';'DSi';'POC';'PP';'PN';'DIC'};
    
for yr=1991:2023
        
	for i = 1:length(nm)

        fnm=[pin nm{i} '_ECCO_V4r5_' int2str(yr)];
        fout=[pout region_name suf1 '_' nm{i} '_' int2str(yr)];
        
        fdir=dir(fnm);
        nlevs=fdir.bytes/4/NX/NX/13;
        fld=read_llc_fkij(fnm,NX,fc,1:nlevs,ix,jx);
        writebin(fout,fld);

	end

end
% }}}

% }}}
