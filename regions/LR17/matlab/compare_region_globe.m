% Compare regional results to the global simulation
nx=10; ny=22; nz=46;
s1='_10x22.'; s2='_10x22x46.';
p1='~/Links/Box/Public/LR17/';
p2='~/mitgcm/darwin3/run/diags/';

% Compare some monthly-mean values
ts=5472;
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
    clf, subplot(121), mypcolor(fl1'), colorbar('h')
    title([fld{1} ' ' ts2dte(ts,1200,1992,1,16,30)])
    subplot(122), mypcolor(fl2(:,:,1)'-fl1'), colorbar('h')
    in=find(~isnan(fl1));
    title(['regional-global: ' num2str((rms(fl2(in)-fl1(in)))/rms(fl1(in)))])
    pause
end

% Compare some monthly-mean values
fld={'SALT'};
for ts=[1152 3240 5472 7632 9864 12024 14256 16488 18648 20880 23040]
    fn=[p1 fld{1} '/' fld{1} s2 ts2dte(ts,1200,1992,1,16,30)];
    fl1=readbin(fn,[nx ny]);
    if strcmp(fld{1},'SALT')
        fl2=35+rdmds([p2 'monthly/SALTanom'],ts);
    else
        fl2=rdmds([p2 'monthly/' fld{1}],ts);
    end
    in=find(fl1==0); fl1(in)=nan; fl2(in)=nan;
    clf, colormap(jet), subplot(121), mypcolor(fl1'), colorbar('h')
    title([fld{1} ' ' ts2dte(ts,1200,1992,1,16,30)])
    subplot(122), mypcolor(fl2(:,:,1)'-fl1'); caxis([-1 1]*.3); colorbar('h')
    in=find(~isnan(fl1));
    title(['regional-global: ' num2str((rms(fl2(in)-fl1(in)))/rms(fl1(in)))])
    pause
end

% Compare some monthly-mean values
fld={'THETA'};
for ts=[1152 3240 5472 7632 9864 12024 14256 16488 18648 20880 23040]
    fn=[p1 fld{1} '/' fld{1} s2 ts2dte(ts,1200,1992,1,16,30)];
    fl1=readbin(fn,[nx ny]);
    if strcmp(fld{1},'SALT')
        fl2=35+rdmds([p2 'monthly/SALTanom'],ts);
    else
        fl2=rdmds([p2 'monthly/' fld{1}],ts);
    end
    in=find(fl1==0); fl1(in)=nan; fl2(in)=nan;
    clf, colormap('default'), subplot(121), mypcolor(fl1'),caxis([10 14]), colorbar('h')
    title([fld{1} ' ' ts2dte(ts,1200,1992,1,16,30)])
    subplot(122), mypcolor(fl2(:,:,1)'-fl1'); caxis([-1 1]*3); colorbar('h')
    in=find(~isnan(fl1));
    title(['regional-global: ' num2str((rms(fl2(in)-fl1(in)))/rms(fl1(in)))])
    pause
end
