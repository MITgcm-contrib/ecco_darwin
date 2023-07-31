% Compare regional results to the global simulation
nx=10; ny=22; nz=46;
s1='_10x22'; s2='_10x22x46';
p1='~/Links/Box/Public/LR17/GlobalCutout/';
p2='~/mitgcm/darwin3/run/diags/';

grid1=[p1 'grid/'];
grid2=[p2 '../'];

% compare RAC
rac1=readbin([grid1 'RAC' s1],[nx ny]);
rac2=readbin([grid2 'RAC.data'  ],[nx ny]);
disp('Percent difference, min/max/mean for RAC')
disp(minmax(100*(rac2-rac1)./rac1))

% read vertical grid thickness
drf=readbin([grid2 'DRF.data'  ],nz);

% compare hHacC
hfac1=readbin([grid1 'hFacC' s2],[nx ny nz]);
hfac2=readbin([grid2 'hFacC.data'     ],[nx ny nz]);
disp('Difference, min/max/mean for hFacC')
disp(minmax(hfac2-hfac1))

% compare area
tmp2d1=rac1.*hfac1(:,:,1);
tmp2d2=rac2.*hfac2(:,:,1);
area1=sum(tmp2d1(:));
area2=sum(tmp2d2(:));
disp(['Area1: ' num2str(area1) '; Area2: ' num2str(area2)])
disp(['Area1-Area2: ' num2str(area1-area2)])

% compare volume
tmp3d1=zeros(nx,ny,nz); tmp3d2=zeros(nx,ny,nz); 
for k=1:nz
    tmp3d1(:,:,k)=rac1.*hfac1(:,:,k).*drf(k);
    tmp3d2(:,:,k)=rac2.*hfac2(:,:,k).*drf(k);
end
vol1=sum(tmp3d1(:));
vol2=sum(tmp3d2(:));
disp(['Volume1: ' num2str(vol1) '; Volume2: ' num2str(vol2)])
disp(['Volume1-Volume2: ' num2str(vol1-vol2)])

% compute surface-averaged and volume-averaged time series
months=359;
flds={'ALK','c1','c2','c3','c4','c5','c6','c7','Chl1','Chl2','Chl3', ...
      'Chl4','Chl5','DIC','DOC','DOFe','DON','DOP','FeT','NH4', ...
      'NO2','NO3','O2','PIC','PO4','POC','POFe','PON','POP','POSi', ...
      'SALT','SiO2','THETA'};

work1=zeros(months, length(flds));
work2=zeros(months, length(flds));
work3=zeros(months, length(flds));
work4=zeros(months, length(flds));

t0=datenum(1992,1,16);
for mn=1:months
    mydisp(mn)
    ts=(datenum(1992,mn+1,1)-t0)*72;

    for f=1:length(flds)
	fld=flds{f};
	
        fn=[p1 fld '/' fld s2 '.' ts2dte(ts,1200,1992,1,16,30)];
        fl1=readbin(fn,[nx ny nz]);
        if strcmp(fld,'SALT')
            fl2=35+rdmds([p2 'monthly/SALTanom'],ts);
        else
            fl2=rdmds([p2 'monthly/' fld],ts);
        end

	tmp2d=fl1(:,:,1).*tmp2d1;
	work1(mn,f)=sum(tmp2d(:));
	tmp2d=fl2(:,:,1).*tmp2d2;
	work2(mn,f)=sum(tmp2d(:));
	tmp3d=fl1.*tmp3d1;
	work3(mn,f)=sum(tmp3d(:));
	tmp3d=fl2.*tmp3d2;
	work4(mn,f)=sum(tmp3d(:));

    end %f
end %mn

work1=work1./area1;
work2=work2./area2;
work3=work3./vol1;
work4=work4./vol2;

% plot figures
eval(['cd ' p1 '..'])
mkdir Regional
cd Regional
tt=(.5:months)/12+1992;

% plot surface-averaged time series 
for i=1:length(flds)
    [p,f]=ind2sub([9 3],i);
    figure(f)
    subplot(3,3,p)
    plot(tt,work1(:,i), tt,work2(:,i))
    grid on
    title(['surface ' flds{i}])
    if p==1; legend('cutout','region'); end
end
figure(1), print -dpdf SurfTimSer_1
figure(2), print -dpdf SurfTimSer_2
figure(3), print -dpdf SurfTimSer_3
figure(4), print -dpdf SurfTimSer_4

% plot volume-integrated time series 
for i=1:length(flds)
    [p,f]=ind2sub([9 3],i);
    figure(f)
    subplot(3,3,p)
    plot(tt,work3(:,i), tt,work4(:,i))
    grid on
    title(flds{i})
    if p==1; legend('cutout','region'); end
end
figure(1), print -dpdf VolTimSer_1
figure(2), print -dpdf VolTimSer_2
figure(3), print -dpdf VolTimSer_3
figure(4), print -dpdf VolTimSer_4

% plot volume-integrated first/last year
it1=1:12; it2=337:348;
for i=1:length(flds)
    [p,f]=ind2sub([9 3],i);
    figure(f)
    subplot(3,3,p)
    plot(it1,work3(it1,i), it1,work4(it1,i), ...
         it1,work3(it2,i), it1,work4(it2,i))
    grid on
    title(flds{i})
    if p==1
        legend('cutout 1992','region 1992','cutout 2019','region 2019')
    end
end
figure(1), print -dpdf VolTimSer_Yr_1
figure(2), print -dpdf VolTimSer_Yr_2
figure(3), print -dpdf VolTimSer_Yr_3
figure(4), print -dpdf VolTimSer_Yr_4
