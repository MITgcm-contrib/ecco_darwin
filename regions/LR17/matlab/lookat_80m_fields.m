nx=936;
ny=875;
for timestep=0:1440:7200
for timestep=11520
    for fld={'Eta','S','T','U','V'}
        tmp=readbin([fld{1},'.',myint2str(timestep,10),'.data'],[nx ny]);
        tmp(find(tmp==0))=nan;
        pcolorcen(tmp')
        title([fld{1} ' day ' int2str(timestep/1440)])
        colorbar
        eval(['print -djpeg ' fld{1} '_day_' myint2str(timestep/1440,3)])
    end
end
