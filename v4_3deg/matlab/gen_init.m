% Generate initial condition files
pin='../../../../MITgcm/verification/tutorial_global_oce_biogeo/input/';

% Uvel, Vvel, Theta, Salt, Eta
fin=[pin 'pickup.0005184000.data'];
tmp=readbin(fin,[128,64,123],1,'real*8');
writebin('../data/Uvel.0005184000',tmp(:,:,1:15));
writebin('../data/Vvel.0005184000',tmp(:,:,16:30));
writebin('../data/Theta.0005184000',tmp(:,:,31:45));
writebin('../data/Salt.0005184000',tmp(:,:,46:60));
writebin('../data/Eta.0005184000',tmp(:,:,121));

% Area, Heff
Ice=tmp(:,:,31);
Ice(find(Ice>-1.8))=0;
Ice(find(Ice))=1;
writebin('../data/Area.0005184000',Ice);
writebin('../data/Heff.0005184000',Ice);

% DIC, Alk, PO4, DOP, O2
fin=[pin 'pickup_ptracers.0005184000.data'];
tmp=readbin(fin,[128,64,75],1,'real*8');
writebin('../data/DIC.0005184000',tmp(:,:,1:15));
writebin('../data/Alk.0005184000',tmp(:,:,16:30));
writebin('../data/PO4.0005184000',tmp(:,:,31:45));
writebin('../data/DOP.0005184000',tmp(:,:,46:60));
writebin('../data/O2.0005184000',tmp(:,:,61:75));

% individual ptracers extracted from
% /nobackup/dcarrol2/v4_3deg/pickup/v4_3deg_pickup_ptracers_optimized.0000000001.data
% Type is real*4, dimensions are [128 64 15 40]
% Use xpolate to fill in some suprious zeros
fn1='v4_3deg_pickup_ptracers_optimized.0000000001.data';
for k=1:15
    tmp=readbin(fn1,[128 64],1,'real*4',31*15+k-1);
    in=find(tmp==0);
    for t=1:39
        tmp=readbin(fn1,[128 64],1,'real*4',(t-1)*15+k-1);
        tmp(in)=nan;
        tmp=xpolate(tmp);
        fn2=['../data/ptracers_optimized_' myint2str(i) '.0000000001'];
        writebin(fn2,tmp,1,'real*4',k-1)
    end
end

% convert Dustin's iron file from 3x3 grid to 2.8125x2.8125 grid
lat1=[-90 -88.5:3:90 90];
lon1=-178.5:3:(180+360);
dx=2.8125;
lat2=(-90+dx/2):dx:90;
lon2=(dx/2):dx:360;
fnm1='../data/llc270_Mahowald_2009_soluble_iron_dust.bin';
fnm2='../data/3deg_Mahowald_2009_soluble_iron_dust.bin';
for m=1:12
    tmp=readbin(fnm1,[120 60],1,'real*4',m-1);
    tmp=[tmp(:,end) tmp tmp(:,1)];
    tmp1=[tmp;tmp];
    tmp2=interp2(lat1',lon1,tmp1,lat2',lon2);
    writebin(fnm2,tmp2,1,'real*4',m-1);
end
