clc;clear;

ieee = 'b';           % big-endian format
accuracy = 'float32'; % this is single-precision (='real*4')

fid = fopen('atemp_1x1_one_year', 'r', ieee);
dat=fread(fid, accuracy); % surface (2m) air temp
fclose(fid)

plot(1:366,dat)
hold on
plot(1:366,dat*1.05)
xlim([1 366])
ylabel('K')
title('Surface air temperature')

fid = fopen('atemp_1x1_one_year_5p', 'w', ieee);
fwrite(fid,dat*1.05,accuracy)
fclose(fid)

fid = fopen('atemp_1x1_one_year_5p', 'r', ieee);
dat=fread(fid, accuracy); % surface (2m) air temp
fclose(fid)
