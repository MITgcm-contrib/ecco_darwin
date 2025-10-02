% compare regional simulation to global cutout
clear
nx=20;
ny=15;
nz=47;

% regional simulation directory
reg_dir='~/mitgcm/darwin3_68y/run/diags/monthly/';

% global cutout directory
cut_dir='~/Box/Public/GoM/llc90/';

% compare ETAN
fld='ETAN';
reg_files=dir([reg_dir fld '*.data']);
cut_files=dir([cut_dir fld '/' fld '*']);
for m=1:length(reg_files)
    fnm=[reg_files(m).folder '/' reg_files(m).name];
    rtmp=readbin([reg_files(m).folder '/' reg_files(m).name],[nx ny]);
    rtmp(find(~rtmp))=nan;
    ctmp=readbin([cut_files(m).folder '/' cut_files(m).name],[nx ny]);
    ctmp(find(~ctmp))=nan;
    clf
    subplot(311)
    pcolorcen(rtmp');
    colormap(cmap)
    caxis([-1 1])
    colorbar
    title(reg_files(m).name)
    subplot(312)
    pcolorcen(ctmp');
    colormap(cmap)
    caxis([-1 1])
    colorbar
    title(cut_files(m).name)
    subplot(313)
    pcolorcen(rtmp'-ctmp');
    colormap(cmap)
    caxis([-1 1])
    colorbar
    title('regional-cutout')
    pause(.1)
end

% compare surface UVELMASS
fld='UVELMASS';
reg_files=dir([reg_dir fld '*.data']);
cut_files=dir([cut_dir 'U/U*']);
for m=1:length(reg_files)
    fnm=[reg_files(m).folder '/' reg_files(m).name];
    rtmp=readbin([reg_files(m).folder '/' reg_files(m).name],[nx ny]);
    rtmp(find(~rtmp))=nan;
    ctmp=readbin([cut_files(m).folder '/' cut_files(m).name],[nx ny]);
    ctmp(find(~ctmp))=nan;
    clf
    subplot(311)
    pcolorcen(rtmp');
    caxis([-1 1])
    colorbar
    title(reg_files(m).name)
    subplot(312)
    pcolorcen(ctmp');
    caxis([-1 1])
    colorbar
    title(cut_files(m).name)
    subplot(313)
    pcolorcen(rtmp'-ctmp');
    caxis([-1 1])
    colorbar
    title('regional-cutout')
    pause
end

% compare surface VVELMASS
fld='VVELMASS';
reg_files=dir([reg_dir fld '*.data']);
cut_files=dir([cut_dir 'V/V*']);
for m=1:length(reg_files)
    fnm=[reg_files(m).folder '/' reg_files(m).name];
    rtmp=readbin([reg_files(m).folder '/' reg_files(m).name],[nx ny]);
    rtmp(find(~rtmp))=nan;
    ctmp=readbin([cut_files(m).folder '/' cut_files(m).name],[nx ny]);
    ctmp(find(~ctmp))=nan;
    clf
    subplot(311)
    pcolorcen(rtmp');
    caxis([-1 1])
    colorbar
    title(reg_files(m).name)
    subplot(312)
    pcolorcen(ctmp');
    caxis([-1 1])
    colorbar
    title(cut_files(m).name)
    subplot(313)
    pcolorcen(rtmp'-ctmp');
    caxis([-1 1])
    colorbar
    title('regional-cutout')
    pause
end
