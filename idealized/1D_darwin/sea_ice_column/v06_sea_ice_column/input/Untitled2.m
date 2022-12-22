clc;clear;

ieee = 'b';           % big-endian format
accuracy = 'float32'; % this is single-precision (='real*4')
fns=dir('i_edsE*');
for ii=1:length(fns)
    fid = fopen(fns(ii).name, 'r', ieee);
    dat(:,ii)=fread(fid, accuracy); % surface (2m) air temp
    fclose(fid)
end

plot(sum(dat,2)/0.2)