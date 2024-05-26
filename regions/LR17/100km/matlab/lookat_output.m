pin='~/mitgcm/MITgcm/run/';
pot='~/links/Box/Public/LR17/figs/';

fld={'Eta','S','T','U','V'};
ts=11520:1440:15840;
nx=936; ny=875;

for f=1:length(fld)
    clf
    for s=1:length(ts)
        fin=[pin fld{f} '.' myint2str(ts(s),10) '.data'];
        eval(['tmp=readbin(''' fin ''',[nx,ny]);'])
        tmp(find(~tmp))=nan;
        subplot(2,2,s)
        mypcolor(tmp');
        colorbar
        title([fld{f} ' at day ' int2str(ts(s)/1440)])
    end
    eval(['print -djpeg ' pot fld{f}])
end
