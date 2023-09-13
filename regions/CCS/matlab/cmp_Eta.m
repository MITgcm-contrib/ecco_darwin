clear


% Compare regional results to the global simulation
nx=80; ny=100; nz=48;
p1='/nobackup/dmenemen/ecco_darwin/CCS_kelp/';
grid1=[p1 'grid/'];
rac1=readbin([grid1 'RAC_80x100'],[nx ny]);
hfac1=readbin([grid1 'hFacC_80x100x48'],[nx ny]);
tmp2d1=rac1.*hfac1;
area1=sum(tmp2d1(:));


months=360;
fld='ETAN';

eta0=1:months;

t0=datenum(1992,1,16);
for mn=1:months
	mydisp(mn)
	ts=(datenum(1992,mn+1,1)-t0)*72;

    fn=[p1 fld '/' fld '_80x100.' ts2dte(ts,1200,1992,1,16,30)];
    fl1=readbin(fn,[nx ny]);

	tmp2d=fl1.*tmp2d1;
	eta0(mn)=sum(tmp2d(:));
end %mn
eta0=eta0/area1;


var1='dynstat_eta_mean';

%ETAN from STDOUT.0000
fn4='/nobackup/dmenemen/ecco_darwin/darwin3/run_CCS_kelp_2/STDOUT.0000'
tmp=eta0*nan;
val=mitgcmhistory(fn4,'time_tsnumber','dynstat_eta_mean');
for mn=1:months
	ts1=(datenum(1992,mn  ,1)-t0)*72;
	ts2=(datenum(1992,mn+1,1)-t0)*72;
	ix=find(val(:,1)>ts1 & val(:,1)<=ts2);
	tmp(mn)=mean(val(ix,2));
end
eta4=tmp;

%ETAN from STDOUT.0000
fn5='/nobackup/dmenemen/ecco_darwin/darwin3/run/STDOUT.0000'
tmp=eta0*nan;
val=mitgcmhistory(fn5,'time_tsnumber','dynstat_eta_mean');
for mn=1:months
	ts1=(datenum(1992,mn  ,1)-t0)*72;
	ts2=(datenum(1992,mn+1,1)-t0)*72;
	ix=find(val(:,1)>ts1 & val(:,1)<=ts2);
	tmp(mn)=mean(val(ix,2));
end
eta5=tmp;

%plot
figure
tt=(.5:months)/12+1992; 
subplot(211)
    plot(tt,eta0, tt,eta5, tt,eta4, 'linew',2)
        grid on
        title('Eta')
	legend('cutout','fix2','fix1')
%        xlim([1992 2018])
subplot(212)
    h=plot(tt,eta0-eta5, tt,eta0-eta4,'linew',2);
        grid on
	hold on
	plot(tt,0*tt,'k')
	legend('cutout - fix2', 'cutout - fix1')

