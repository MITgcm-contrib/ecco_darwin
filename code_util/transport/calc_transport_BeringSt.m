clear

%model
dirModel='/nobackup/dcarrol2/v05_latest/darwin3/run/';
dirDiagB=[dirModel 'diags/budget/']; 
dirDiagM=[dirModel 'diags/monthly/']; 

%grid
gcmfaces_global; global mygrid
nx=270;ny=nx*13;nz=50;
dirGrid='/nobackup/hzhang1/llc_1080/MITgcm/DM_270/GRID_up/';
nF=5;fileFormat='compact';
grid_load(dirGrid,nF,fileFormat);
if ~isfield(mygrid,'LINES_MASKS');
    [lonPairs,latPairs,names]=gcmfaces_lines_pairs;
    gcmfaces_lines_transp(lonPairs,latPairs,names);
end;
LINES_MASKS=mygrid.LINES_MASKS(1); %BS

%time series
TT=12*29; %1992-2020 month
fld={'NO3', 'PO4', 'SiO2'};
TP=zeros(TT,3);

for mn=1:TT
disp(mn)

	ts=(datenum(1992,mn+1,1)-datenum(1992,1,1) )*72;
	for i=1:3
	fldU1=readbin([dirDiagB 'average_' fld{i} '_3d.' myint2str(ts,10) '.data'],[nx ny nz],1,'real*4',2 -1);
	fldV1=readbin([dirDiagB 'average_' fld{i} '_3d.' myint2str(ts,10) '.data'],[nx ny nz],1,'real*4',3 -1);
	fldU=convert2gcmfaces(fldU1);
	fldV=convert2gcmfaces(fldV1);
	fldU=fldU.*mygrid.mskW; 
	fldV=fldV.*mygrid.mskS;
	tmp=1e-6*calc_transports(fldU,fldV,LINES_MASKS);
	TP(mn,i)=nansum(tmp);
	end

end %mn
	
%plot
figure
yrs=1992:2020;
for i=1:3
amc=TP(:,i);
subplot(3,1,i)
amc_yr=reshape(amc,[12 29]);
amc_yr=squeeze(mean(amc_yr,1));
plot((.5:TT)/12+yrs(1),amc,yrs+.5,amc_yr,'linewidth',2)
title(['Bering St transport [Sv*' fld{i} ']'])
grid on
xlim([yrs(1) yrs(end)+1])
end

