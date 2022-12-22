clc;clear;

ieee = 'b';           % big-endian format
accuracy = 'float32'; % this is single-precision (='real*4')

fns=dir('edsE*');

for ii=1:length(fns)
    fid = fopen(fns(ii).name, 'r', ieee);
    dat=fread(fid, accuracy); % surface (2m) air temp
    fclose(fid)

    lon=-179.5:1:179.5;
    lat=-89.5:1:89.5;
    datr=reshape(dat,[360,180,12]);

%     pcolor(lon,lat,datr(:,:,5)')

    dati=squeeze(mean(mean(datr(:,lat>70,:))));
    dati2=[dati(end);dati;dati(1)];
    ti=(-15:30:390)';
    datii=interp1(ti,dati2,1:366);
    plot(1:366,datii)
    hold on
    
    fid = fopen(['i_' fns(ii).name], 'w', ieee);
    fwrite(fid,datii,accuracy)
    fclose(fid)
%     pause
end