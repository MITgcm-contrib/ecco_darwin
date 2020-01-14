% Generate initial condition files
pin='../../../../MITgcm/verification/tutorial_global_oce_biogeo/input/';

% Uvel, Vvel, Theta, Salt, Eta
fin=[pin 'pickup.0005184000.data'];
tmp=readbin(fin,[128,64,123],1,'real*8');
writebin('Uvel.0005184000',tmp(:,:,1:15));
writebin('Vvel.0005184000',tmp(:,:,16:30));
writebin('Theta.0005184000',tmp(:,:,31:45));
writebin('Salt.0005184000',tmp(:,:,46:60));
writebin('Eta.0005184000',tmp(:,:,121));

% DIC, Alk, PO4, DOP, O2
fin=[pin 'pickup_ptracers.0005184000.data'];
tmp=readbin(fin,[128,64,75],1,'real*8');
writebin('DIC.0005184000',tmp(:,:,1:15));
writebin('Alk.0005184000',tmp(:,:,16:30));
writebin('PO4.0005184000',tmp(:,:,31:45));
writebin('DOP.0005184000',tmp(:,:,46:60));
writebin('O2.0005184000',tmp(:,:,61:75));
