

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
MaNumFiles = dir([inputDir 'MaNum.0*data']);

TFiles = dir([inputDir 'THETA.0*data']);
PFiles = dir([inputDir 'PAR.0*data']);
NFiles = dir([inputDir 'NO3.0*data']);

% Reading the grid info
% ptrc31=readbin(fnam,[llc,llc*nf,nz,ptrs]);


NT=72;
ocean=1;

for jt=1:NT
jt
MacFiles(jt).name;

% surface fiedls only
cc = readbin([inputDir MacFiles(jt).name],[nx ny],1,'real*4');
cn = readbin([inputDir MaNFiles(jt).name],[nx ny],1,'real*4');

manu = readbin([inputDir MaNumFiles(jt).name],[nx ny],1,'real*4');

tt = readbin([inputDir TFiles(jt).name],[nx ny],1,'real*4');
ll= readbin([inputDir PFiles(jt).name],[nx ny],1,'real*4');
nn=readbin([inputDir NFiles(jt).name],[nx ny],1,'real*4');

MaC(:,:,jt)=cc;
MaN(:,:,jt)=cn;
Manum(:,:,jt)=manu;

temp(:,:,jt)=tt;
par(:,:,jt)=ll;
no3(:,:,jt)=nn;

sumC=0;
sumCN=0;
sumT=0;
sumP=0;
sumN=0;
sumNum=0;
NC=0;

for jx=1:nx
for jy=1:ny

if lat(jx,jy)  > 33 &  lat(jx,jy) < 34 
if lon(jx,jy)  < -119 & lon(jx,jy) > -120 

if  smask(jx,jy) == ocean

sumC = sumC + MaC(jx,jy,jt);
sumCN = sumCN + MaN(jx,jy,jt);
sumT = sumT + temp(jx,jy,jt);
sumP = sumP + par(jx,jy,jt);
sumN = sumN + no3(jx,jy,jt);

sumNum = sumNum +  Manum(jx,jy,jt);

NC = NC + 1;

end


end % if lon
end  % if lat

end
end

MaCave(jt)=sumC/NC;

MaNave(jt)=sumCN/NC;

MaNumave(jt)=sumNum/NC;

Tave(jt)=sumT/NC;

PARave(jt)=sumP/NC;

NO3ave(jt)=sumN/NC;

end


% MaC
% MaN
% MaNum
% MaFra

%timehours

subplot(6,1,1)
plot(PARave,'linewidth',2)
grid on
%xlabel('Time (hours)')
ylabel('(W/m2)')
legend('Surface PAR')
ylim([0 600])
title('ECCO-Darwin CCS (33-34 N,119-120 W)')
%

subplot(6,1,2)
plot(Tave,'linewidth',2)
grid on
%xlabel('Time (hours)')
ylabel('(deg C)')
legend('SST')


subplot(6,1,3)
plot(NO3ave,'linewidth',2)
grid on
xlabel('Time (hours)')
ylabel('(mmol/m3)')
legend('Surface [NO3]')



%
subplot(6,1,4)
plot(MaCave,'linewidth',2)
grid on

legend('Surf. Kelp Carbon')
%
subplot(6,1,5)
plot(MaNave,'linewidth',2)
grid on

legend('Surf. Kelp Nitrogen')

%

subplot(6,1,6)
plot(MaNumave,'linewidth',2)
grid on

xlabel('Time (hours)')
legend('# Surf. Plants/Unit Area')




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








