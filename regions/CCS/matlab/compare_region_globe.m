% Compare regional results to the global simulation
nx=80; ny=100;
p1='/nobackup/dmenemen/ecco_darwin/CCS_kelp/';
p2='/nobackup/dmenemen/ecco_darwin/darwin3/run/diags/';

% Compare some monthly-mean values
ts=807552;
for fld={'ALK','c1','c2','c3','c4','c5','c6','c7','Chl1','Chl2','Chl3', ...
         'Chl4','Chl5','DIC','DOC','DOFe','DON','DOP','FeT','NH4', ...
         'NO2','NO3','O2','PIC','PO4','POC','POFe','PON','POP','POSi', ...
         'SALT','SiO2','THETA'}
    fn=[p1 fld{1} '/' fld{1} '_80x100x48.' ts2dte(ts,1200,1992,1,16,30)];
    fl1=readbin(fn,[nx ny]);
    if strcmp(fld{1},'SALT')
        fl2=35+rdmds([p2 'monthly/SALTanom'],ts);
    else
        fl2=rdmds([p2 'monthly/' fld{1}],ts);
    end
    in=find(fl1==0); fl1(in)=nan; fl2(in)=nan;
    clf, subplot(121), pcolorcen(fl1'), colorbar('h')
    title([fld{1} ' ' ts2dte(ts,1200,1992,1,16,30)])
    subplot(122), pcolorcen(fl2(:,:,1)'-fl1'), colorbar('h')
    in=find(~isnan(fl1));
    title(['regional-global: ' num2str((rms(fl2(in)-fl1(in)))/rms(fl1(in)))])
    pause
end
