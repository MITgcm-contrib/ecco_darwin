% Compare regional results to the global simulation
nx=80; ny=100; nz=48;
p1='/nobackup/dmenemen/ecco_darwin/CCS_kelp/';
p2='/nobackupnfs1/hzhang1/testing_CCS/darwin3/run/diags/';

grid1=[p1 'grid/'];
grid2='/nobackupnfs1/hzhang1/testing_CCS/darwin3/GRID/';


rac1=readbin([grid1 'RAC_80x100'],[nx ny]);
rac2=readbin([grid2 'RAC.data'  ],[nx ny]);
%diff

drf=readbin([grid2 'DRF.data'  ],nz);
%common

hfac1=readbin([grid1 'hFacC_80x100x48'],[nx ny nz]);
hfac2=readbin([grid2 'hFacC.data'     ],[nx ny nz]);
%same, also for Depth

tmp3d1=zeros(nx,ny,nz); tmp3d2=zeros(nx,ny,nz); 
for k=1:nz
tmp3d1(:,:,k)=rac1.*hfac1(:,:,k).*drf(k);
tmp3d2(:,:,k)=rac2.*hfac2(:,:,k).*drf(k);
end
vol1=sum(tmp3d1(:));
vol2=sum(tmp3d2(:));
%almost same


months=240;
flds={'ALK','c1','c2','c3','c4','c5','c6','c7','Chl1','Chl2','Chl3', ...
         'Chl4','Chl5','DIC','DOC','DOFe','DON','DOP','FeT','NH4', ...
         'NO2','NO3','O2','PIC','PO4','POC','POFe','PON','POP','POSi', ...
         'SALT','SiO2','THETA'};

work1=zeros(months, length(flds));
work2=zeros(months, length(flds));

t0=datenum(1992,1,16);
for mn=1:months
	mydisp(mn)
	ts=(datenum(1992,mn+1,1)-t0)*72;

	for f=1:length(flds)
	fld=flds{f};
	
    fn=[p1 fld '/' fld '_80x100x48.' ts2dte(ts,1200,1992,1,16,30)];
    fl1=readbin(fn,[nx ny nz]);
    if strcmp(fld,'SALT')
        fl2=35+rdmds([p2 'monthly/SALTanom'],ts);
    else
        fl2=rdmds([p2 'monthly/' fld],ts);
    end

	tmp3d=fl1.*tmp3d1;
	work1(mn,f)=sum(tmp3d(:));
	tmp3d=fl2.*tmp3d2;
	work2(mn,f)=sum(tmp3d(:));

	end %f
end %mn

	work1=work1./vol1;
	work2=work2./vol2;

%plot
tt=(.5:months)/12+1992; 
for i=1:length(flds)
    [p,f]=ind2sub([9 3],i);
    figure(f)
    subplot(3,3,p)
    plot(tt,work1(:,i), tt,work2(:,i), 'linew',2)
        grid on
        title(flds{i})
	if p==1; legend('cutout','region'); end
%        xlim([1992 2018])
end

