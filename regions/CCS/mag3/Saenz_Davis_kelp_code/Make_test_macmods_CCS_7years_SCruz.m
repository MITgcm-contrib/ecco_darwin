
% files to read
% to get lat,lon, land mask

fXC = '/home/mmanizza/work/Kelp/ED/out2/XC.data';
fYC = '/home/mmanizza/work/Kelp/ED/out2/YC.data';
fmask = '/home/mmanizza/work/Kelp/ED/out2/hFacC.data';

% Lists  of files
% to read

%file_PAR='File_PAR.txt';
%list_PAR=textread(file_PAR,'%s');

%file_THETA='File_THETA.txt';
%list_THETA=textread(file_THETA,'%s');

file_NO3='File_NO3.txt';
list_NO3=textread(file_NO3,'%s');
%
file_gT='File_gT.txt';
list_gT=textread(file_gT,'%s');

file_gQ='File_gQ.txt';
list_gQ=textread(file_gQ,'%s');

file_gE='File_gE.txt';
list_gE=textread(file_gE,'%s');

file_gH='File_gH.txt';
list_gH=textread(file_gH,'%s');

file_Bmag='File_Bmag.txt';
list_Bmag=textread(file_Bmag,'%s');

file_Nmag='File_Nmag.txt';
list_Nmag=textread(file_Nmag,'%s');

file_wvp='File_wv_prd.txt';
list_wvp=textread(file_wvp,'%s');

file_wvh='File_wv_hght.txt';
list_wvh=textread(file_wvh,'%s');

file_mrt='File_mort.txt';
list_mrt=textread(file_mrt,'%s');

file_mu='File_mu_mag.txt';
list_mu=textread(file_mu,'%s');


% Dimensions
% of CCS regional model
L1=100;
L2=80;


% file to read setting
typp=1;
precc='real*4';
% this parameter controls the depth to select when reading 
%the files
% if sk=0; then you are reading the surface values
% meaning the first vertical model level
sk=0; 
mformm='ieee-be';

% Read longitude
lon = readbin(fXC,[L2 L1],typp,precc,sk,mformm);

% Read latitude
lat = readbin(fYC,[L2 L1],typp,precc,sk,mformm);

% Read mask
mask = readbin(fmask,[L2 L1],typp,precc,sk,mformm);

% Values of land/ocean mask
ocean=1;
land=0;

%return

% Number of files to process
NT=12*24;

%NT=1

for jt=1:NT


% read PAR

%f_PAR = list_PAR(jt);
%f_PAR=char(f_PAR);
%ppar = readbin(f_PAR,[L2 L1],typp,precc,sk,mformm);
%PAR(:,:,jt) = ppar;

% read THETA 

%f_THETA = list_THETA(jt);
%f_THETA=char(f_THETA);
%tt = readbin(f_THETA,[L2 L1],typp,precc,sk,mformm);
%THETA(:,:,jt) =tt; 

% read NO3

f_NO3 = list_NO3(jt);
f_NO3=char(f_NO3);
n3 = readbin(f_NO3,[L2 L1],typp,precc,sk,mformm);
NO3(:,:,jt) =n3;

% read gT

f_gT = list_gT(jt);
f_gT=char(f_gT);
tt = readbin(f_gT,[L2 L1],typp,precc,sk,mformm);
gT(:,:,jt) = tt;

% read gQ

f_gQ = list_gQ(jt);
f_gQ=char(f_gQ);
qq = readbin(f_gQ,[L2 L1],typp,precc,sk,mformm);
gQ(:,:,jt) = qq;

% read gE

f_gE = list_gE(jt);
f_gE=char(f_gE);
ee = readbin(f_gE,[L2 L1],typp,precc,sk,mformm);
gE(:,:,jt) = ee;

% read gH

f_gH = list_gH(jt);
f_gH=char(f_gH);
hh = readbin(f_gH,[L2 L1],typp,precc,sk,mformm);
gH(:,:,jt) = hh;



% read Bmag 

f_B = list_Bmag(jt);
f_B=char(f_B);
BB = readbin(f_B,[L2 L1],typp,precc,sk,mformm);
Bmag(:,:,jt) = BB;

f_N = list_Nmag(jt);
f_N=char(f_N);
NN = readbin(f_N,[L2 L1],typp,precc,sk,mformm);
Nmag(:,:,jt) = NN;

f_wvp = list_wvp(jt);
f_wvp=char(f_wvp);
WVP = readbin(f_wvp,[L2 L1],typp,precc,sk,mformm);
wvp(:,:,jt) = WVP;

f_wvh = list_wvh(jt);
f_wvh=char(f_wvh);
WVH = readbin(f_wvh,[L2 L1],typp,precc,sk,mformm);
wvh(:,:,jt) = WVH;

f_mrt = list_mrt(jt);
f_mrt=char(f_mrt);
MRT = readbin(f_mrt,[L2 L1],typp,precc,sk,mformm);
mmrt(:,:,jt) = MRT;

f_mu = list_mu(jt);
f_mu=char(f_mu);
mm = readbin(f_mu,[L2 L1],typp,precc,sk,mformm);
mu(:,:,jt) = mm;

end

%return

% Now computing an average of the seasonal
% cycle in the SD region

Nlim=36.5;
Slim=35.5;

time = [1:1:NT];



for jt=1:NT
jt

sum_Bmag=0;
sum_Nmag=0;
sum_mort=0;
sum_wvh=0;
sum_wvp=0;


sum_gT=0;
sum_gQ=0;
sum_gE=0;
sum_gH=0;
sum_mu=0;

%sum_par=0;
%sum_theta=0;
sum_no3=0;

N=0;

for ji=1:L2
for jj=1:L1

if lat(ji,jj) < Nlim  & lat(ji,jj) > Slim 

if lon(ji,jj) > -123  & lon(ji,jj) < -122

if mask(ji,jj) == ocean 

sum_Bmag = sum_Bmag + Bmag(ji,jj,jt);

sum_Nmag = sum_Nmag + Nmag(ji,jj,jt);

sum_mort = sum_mort + mmrt(ji,jj,jt);

sum_mort = sum_mort + mmrt(ji,jj,jt);

sum_wvp = sum_wvp + wvp(ji,jj,jt);

sum_wvh = sum_wvh + wvh(ji,jj,jt);

sum_gQ = sum_gQ + gQ(ji,jj,jt);

sum_gE = sum_gE + gE(ji,jj,jt);

sum_gH = sum_gH + gH(ji,jj,jt);

sum_gT = sum_gT + gT(ji,jj,jt);

sum_mu = sum_mu + mu(ji,jj,jt);

%sum_par =sum_par + PAR(ji,jj,jt);

%sum_theta = sum_theta + THETA(ji,jj,jt);


if  NO3(ji,jj,jt) < 0

sum_no3 = sum_no3 + 0;

else

sum_no3 = sum_no3 + NO3(ji,jj,jt);

end

N = N +1;

end % loop mask

end % loop lon limits

end % loop lat limits



end % loop jj

end % loop ji


% Storing averages


Bmag_ave(jt) = sum_Bmag/N;

Nmag_ave(jt) = sum_Nmag/N;

mort_ave(jt) = sum_mort/N;

wvp_ave(jt) = sum_wvp/N;

wvh_ave(jt) = sum_wvh/N;

gT_ave(jt) = sum_gT/N;

gQ_ave(jt) = sum_gQ/N;

gE_ave(jt) = sum_gE/N;

gH_ave(jt) = sum_gH/N;

mu_ave(jt) = sum_mu/N;

%par_ave(jt) = sum_par/N;

%theta_ave(jt) = sum_theta/N;

no3_ave(jt) = sum_no3/N;


end % time loop

return

time=(1992:1/12:1999);
time=time(1:end-1);

% plot setting
fntsz=14;
tckline=2;

% max growth rate Macroc. Pyr.
mumax=0.2;

xrmn=1992;
xrmx=1998;

xrmn=1992;
xrmx=1998;

% PLOT SST
subplot(4,2,1)
plot(time,theta_ave,'b','linewidth',2)
grid on
hold on
xlim([xrmn xrmx])
ylabel('(deg C)','FontSize', fntsz,'FontWeight','bold')
title('Temperature')
%xlabel('Time (Years)')
set(gca,'FontSize',fntsz)
set(gca,'XTickLabel',[]);
set(gca,'LineWidth',tckline)

% PLOT PAR
subplot(4,2,3)
plot(time,par_ave,'b','linewidth',2)
grid on
hold on
%ylim([yrmn yrmx])
xlim([xrmn xrmx])

title('PAR')
ylabel('W/m2','FontSize', fntsz,'FontWeight','bold')
%xlabel('Time (Years)')
set(gca,'FontSize',fntsz)
set(gca,'XTickLabel',[]);
set(gca,'LineWidth',tckline)

% PLOT  NO3
subplot(4,2,5)
plot(time,no3_ave,'b','linewidth',2)
grid on
hold on
ylim([0 1])
xlim([xrmn xrmx])

title('[NO_3]')
ylabel('mmol/m3','FontSize', fntsz,'FontWeight','bold')
%xlabel('Time (Years)')
set(gca,'FontSize',fntsz)
set(gca,'XTickLabel',[]);
set(gca,'LineWidth',tckline)

% PLOT GROWTH RATE
subplot(4,2,7)
plot(time,mu_ave,'b','linewidth',2)
grid on
hold on
ylim([0 mumax])
xlim([xrmn xrmx])
title('Growth Rate (\mu) ')
ylabel('1/day)','FontSize', fntsz,'FontWeight','bold')
xlabel('Time (Years)')
set(gca,'FontSize',fntsz)
%set(gca,'XTickLabel',[]);
set(gca,'LineWidth',tckline)

% PLOT gT
subplot(4,2,2)
plot(time,gT_ave,'r','linewidth',2)
grid on
title('gT')
% 
set(gca,'LineWidth',tckline)
set(gca,'FontSize',fntsz)
set(gca,'XTickLabel',[]);
ylim([0 1])
xlim([xrmn xrmx])


% PLOT gQ
subplot(4,2,4)
plot(time,gQ_ave,'m','linewidth',2)
grid on
title('gQ')
set(gca,'LineWidth',tckline)
set(gca,'FontSize',fntsz)
set(gca,'XTickLabel',[]);
ylim([0 1])
xlim([xrmn xrmx])

% PLOT gE
subplot(4,2,6)
plot(time,gE_ave,'g','linewidth',2)
grid on
title('gE')
set(gca,'LineWidth',tckline)
set(gca,'FontSize',fntsz)
set(gca,'XTickLabel',[]);
ylim([0 1])
xlim([xrmn xrmx])

% PLOT gH
subplot(4,2,8)
plot(time,gH_ave,'b','linewidth',2)
grid on
title('gH')
set(gca,'LineWidth',tckline)
set(gca,'FontSize',fntsz)
ylim([0 1])
xlim([xrmn xrmx])
xlabel('Time (Years)')
















