clc;clear;

for ii=1:31
    
    ieee = 'b';           % big-endian format
    accuracy = 'float32'; % this is single-precision (='real*4')
    fid = fopen(['ptracer' num2str(ii,'%02d') '_ini'], 'r', ieee);
    dat=fread(fid,accuracy);
    plot(dat)
    fclose(fid)
    pause
end