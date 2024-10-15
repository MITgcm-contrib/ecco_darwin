clear all
close all

%%%% load field data %%%%
years=[1.596 3.415 3.523 3.696 3.713 3.804 4.085 4.3 5.383];
days2=years*365;
days=round(days2);
dates=days-582;
area_data=[10.75 85.77 93.34 106.44 107.75 115.11 139.84 160.85 296.96];%cm2
SGR_data=[0.004691 0.002191 0.002124 .002024 .002015 .001967 .001831 .001740 .001398];

%%%%%% Zooxanthellae %%%%%%%%%%%%%%%

Cells_data=10^6*[1.7 1.2 2.15 1.3 0.9 0.8 0.9 0.85 0.7 0.8 1.3 1.4];
dates2= [1 46 74 105 135 164 195 225 256 286 316 346 ] ;


% a= datenum(2003, 01, 01);
% 
% b=datenum(2003,4 , 15 );  % (year,month,day)
% 
% date=b-a+1


data=load('Host.dat');
n=data(:,1);
H=data(:,2);
Jx=data(:,3);
Jn=data(:,4);
JHG=data(:,5);
JHT=data(:,6);
dH=data(:,7);
dt=data(:,8);
Vol=data(:,9);
Area=data(:,10);
SGR=data(:,11);


a=1./dt;

data2=load('Symbiont.dat');%data2b=load('Symbiont1.dat');
n=data2(:,1);%nb=data2b(:,1);
S=data2(:,2);
JL=data2(:,3);
JCO2=data2(:,4);
JeL=data2(:,5);
JNPQ=data2(:,6);
JCP=data2(:,7);%JCPb=data2b(:,7);
CROS=data2(:,8);
JST=data2(:,9);
JSG=data2(:,10);
dS=data2(:,11);
Cells=data2(:,12);


data3=load('Symbiosis.dat');%data3b=load('Symbiosis1.dat');
n=data3(:,1);
SH=data3(:,2);%SHb=data3b(:,2);
pC=data3(:,3);
pN=data3(:,4);
JeC=data3(:,5);
rNH=data3(:,6);
rCH=data3(:,7);
rNS=data3(:,8);
rCS=data3(:,9);

data4=load('Forcing.dat');
n=data4(:,1);
X=data4(:,2);
DIN=data4(:,3);
L=data4(:,4);
T=data4(:,5);T=T-273.15;
kHT=data4(:,6);
kST=data4(:,7);
% 

JCPm=13.3937;
yCL=0.1;
yC=0.8;
jSGm=0.0210;
nNS=0.13;
nNX=0.2;
jHGm=0.0244;
nNH=0.18;

% for i=1:max(n)
% 
% pl(i)=log(min(JL(i)* yCL, JCPm)/min((JCO2(i) + rCH(i))*H(i)/S(i) + rCS(i),JCPm));
% sl(i)=log(min((pN(i)*H(i)/S(i) + rNS(i))/nNS, jSGm)/min(yC*JCP(i), jSGm));
% hl(i)=log(min((Jn (i) + nNX*Jx(i) + rNH(i)) /nNH, jHGm)/min(yC*pC(i)*S(i)/H(i) + Jx(i), jHGm));
% end


figure(1)
plot(n(1:365/dt(1)),Area(1:365/dt(1)),'linewidth',2)
hold on
%plot(dates,area_data,'s','linewidth',2)
%scatter(dates/dt(1),area_data,100,'o','filled')
xlabel('Time(days)')
ylabel('Colony surface area (cm^2)')
legend('Model','Field data')
set(gca,'fontsize',16)
set(gca,'xtick',[1/dt(1);32/dt(1);60/dt(1);91/dt(1);121/dt(1);152/dt(1);182/dt(1);213/dt(1);244/dt(1);274/dt(1);305/dt(1);335/dt(1);366/dt(1)]);
set(gca,'XTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D','J'});


figure(2)
%plot(n(1:365/dt(1)),Cells(1:365/dt(1)),'b','linewidth',2)
plot(n,Cells,'b','linewidth',2)
hold on
%scatter(dates2/dt(1),Cells_data,100,'ro','filled')
ylabel('Zooxanthellae Density (cells cm^-^2)')
legend('Model','Field data')
%datetick('x','mmm')
set(gca,'fontsize',16)
set(gca,'xtick',[1/dt(1);183/dt(1);367/dt(1);548/dt(1);732/dt(1);913/dt(1);1097/dt(1);1278/dt(1);1462/dt(1);1644/dt(1)]);
set(gca,'XTickLabel',{'Jan08','Jul08','Jan09','Jul09','Jan10','Jul10','Jan11','Jul11','Jan12','Jul12'});
    
%set(gca,'xtick',[1/dt(1);32/dt(1);60/dt(1);91/dt(1);121/dt(1);152/dt(1);182/dt(1);213/dt(1);244/dt(1);274/dt(1);305/dt(1);335/dt(1);366/dt(1)]);
%set(gca,'XTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D','J'});

figure(3)
plot(SH)