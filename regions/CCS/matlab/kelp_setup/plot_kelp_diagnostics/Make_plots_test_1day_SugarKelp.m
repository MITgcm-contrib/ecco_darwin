

% Dimensions of the file
nx=80;
ny=100;
nz=48;
% Reading info

typp=1;
precc='real*4';
%skipp1=2*50;
mformm='ieee-be';
sk=0;
% Creating a fake array wth zeros to add to 
% the existing data.ptracers 
%ptrc4=ones(llc,nf,nz,trcs);

%fnam_c='/nobackupp18/mmanizza/Kelp/CCS/run_test/diags/hourly/';

inputDir='/nobackupp18/mmanizza/Kelp/CCS/run_test/diags/hourly/';

fXC='/nobackupp18/mmanizza/Kelp/CCS/run_test/XC.data';
fYC='/nobackupp18/mmanizza/Kelp/CCS/run_test/YC.data';
fmask='/nobackupp18/mmanizza/Kelp/CCS/run_test/hFacC.data';

lon = readbin(fXC,[nx ny],1);
lat = readbin(fYC,[nx ny],1);
smask = readbin(fmask,[nx ny],1);

MacFiles = dir([inputDir 'MaC.0*data']);
MaNFiles = dir([inputDir 'MaN.0*data']);

TFiles = dir([inputDir 'THETA.0*data']);
PFiles = dir([inputDir 'PAR.0*data']);
NFiles = dir([inputDir 'NO3.0*data']);

% Reading the grid info
% ptrc31=readbin(fnam,[llc,llc*nf,nz,ptrs]);


NT=24;
ocean=1;

for jt=1:NT
jt;
MacFiles(jt).name;

% surface fiedls only
cc = readbin([inputDir MacFiles(jt).name],[nx ny],1,'real*4');
tt = readbin([inputDir TFiles(jt).name],[nx ny],1,'real*4');
ll= readbin([inputDir PFiles(jt).name],[nx ny],1,'real*4');
nn=readbin([inputDir NFiles(jt).name],[nx ny],1,'real*4');

MaC(:,:,jt)=cc;
temp(:,:,jt)=tt;
par(:,:,jt)=ll;
no3(:,:,jt)=nn;

sumC=0;
sumT=0;
sumP=0;
sumN=0;
NC=0;

for jx=1:nx
for jy=1:ny

if lat(jx,jy)  > 33 &  lat(jx,jy) < 34 
%if lon(jx,jy)  < -119 & lon(jx,jy) > -120 

if  smask(jx,jy) == ocean

sumC = sumC + MaC(jx,jy,jt);
sumT = sumT + temp(jx,jy,jt);
sumP = sumP + par(jx,jy,jt);
sumN = sumN + no3(jx,jy,jt);
NC = NC + 1;

end


%end % if lon
end  % if lat

end
end

MaCave(jt)=sumC/NC;

Tave(jt)=sumT/NC;

PARave(jt)=sumP/NC;

NO3ave(jt)=sumN/NC;

end


% MaC
% MaN
% MaNum
% MaFra






return

%smask=readbin(fnmsk,[llc,llc*nf,1],1);

%ptrc19=ptrc31(:,:,:,1:19);
%size(ptrc31)
%ptrc2031=ptrc31(:,:,:,20:31);
%size(ptrc2031)



% ptrc31=readbin(fnam,[llc,llc*nf,nz,ptrs]);

ddim=4; % Dimesions for concat

% Appending the new additional four
% fake tracers
ptrc35=cat(ddim,ptrc31,ptrc4);
size(ptrc35)

fnamout='pickup_ptracers_35_test_CCS.0000000001.data';
%pprecc='real*8';
writebin(fnamout,ptrc35,typp,precc);








