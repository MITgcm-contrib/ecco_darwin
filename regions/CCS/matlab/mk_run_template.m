%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build runtime input files needed to run the Gulf of Guinea
% sub-region of llc_270.  Run after running ExtractFields.m.

% {{{ Define desired region and initialize some variables 
NX=270;
prec='real*4';
region_name='GulfGuinea';
minlat=-13.5;
maxlat=6.75;
minlon=-9;
maxlon=14;
mints=dte2ts('02-Jan-1992',1200,1992,1,1);
maxts=dte2ts('31-Dec-1992',1200,1992,1,1);
[status hostname]=system('echo $HOSTNAME');
if startsWith(hostname,'LMC-051773')
    gdir= '/Users/dmenemen/projects/llc/llc270/grid/';
    pin =['/Users/dmenemen/projects/' region_name '/'];
    pout=['/Users/dmenemen/projects/' region_name '/run_template/'];
else
    gdir= '/nobackup/dmenemen/llc/llc_270/grid/';
    pin =['/nobackup/dmenemen/' region_name '/'];
    pout=['/nobackup/dmenemen/' region_name '/run_template/'];
end
eval(['mkdir ' pout])
% }}}

% {{{ Extract indices for desired region 
fnam=[gdir 'Depth.data'];
[fld fc ix jx]=quikread_llc(fnam,NX,1,prec,gdir,minlat,maxlat,minlon,maxlon);
quikpcolor(fld')
[nx ny]=size(fld);
RF=-readbin([gdir 'RF.data'],51);
kx=1:min(find(RF(2:end)>max(fld(:))));
kx=1:49;
nz=length(kx);
% }}}

% {{{ Make bathymetry file 
close all
suf=['_' int2str(nx) 'x' int2str(ny)];
writebin([pout 'BATHY' suf  '_' region_name],-fld);
% }}}

% {{{ Make delYFile 
nx=69; ny=66;
pout=['~//projects/GulfGuinea/run_template/'];
fnm='~/projects/GulfGuinea/grid/YG_69x66';
tmp=readbin(fnm,[nx ny]);
delY=zeros(ny,1);
delY(1:ny-1)=diff(tmp(1,:));
delY(ny)=2*delY(ny-1)-delY(ny-2);
writebin([pout 'delYFile'],delY);
% }}}

% {{{ Generate initial conditions 
eval(['!cp ' pin 'ETAN/ETAN_69x66.02-Jan-1992 ' pout])
eval(['!cp ' pin 'THETA/THETA_69x66x49.02-Jan-1992 ' pout])
eval(['!cp ' pin 'UVELMASS/UVELMASS_69x66x49.02-Jan-1992 ' pout])
eval(['!cp ' pin 'VVELMASS/VVELMASS_69x66x49.02-Jan-1992 ' pout])
SALT=readbin([pin 'SALTanom/SALTanom_69x66x49.02-Jan-1992'],[nx ny nz])+35;
writebin([pout 'SALT_69x66x49.02-Jan-1992'],SALT)
% }}}

% {{{ Generate U/V/T/S lateral boundary conditions
for fld={'THETA','SALTanom','UVELMASS','VVELMASS'}
    pnm=[pin fld{1} '/'];
    fnm=dir([pnm '*' fld{1} '*']);
    for t=1:length(fnm)
        fin=[pnm fnm(t).name];
        disp(fin)
        if fld{1}(1)=='S'
            tmp=readbin(fin,[nx ny nz])+35;
        else
            tmp=readbin(fin,[nx ny nz]); 
        end
        fout=[pout fld{1}(1) '_West']; % western boundary condition
        if fld{1}(1)=='U'
            writebin(fout,squeeze(tmp(2,:,:)),1,prec,t-1)
        else
            writebin(fout,squeeze(tmp(1,:,:)),1,prec,t-1)
        end
        fout=[pout fld{1}(1) '_South']; % southern boundary condition
        if fld{1}=='V'
            writebin(fout,squeeze(tmp(:,2,:)),1,prec,t-1)
        else
            writebin(fout,squeeze(tmp(:,1,:)),1,prec,t-1)
        end
    end
end
% }}}
