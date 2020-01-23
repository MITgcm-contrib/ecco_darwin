% Generate initial condition files
pin='../../../../MITgcm/verification/tutorial_global_oce_biogeo/input/';

% Uvel, Vvel, Theta, Salt, Eta
fin=[pin 'pickup.0005184000.data'];
tmp=readbin(fin,[128,64,123],1,'real*8');
writebin('../input/Uvel.0005184000',tmp(:,:,1:15));
writebin('../input/Vvel.0005184000',tmp(:,:,16:30));
writebin('../input/Theta.0005184000',tmp(:,:,31:45));
writebin('../input/Salt.0005184000',tmp(:,:,46:60));
writebin('../input/Eta.0005184000',tmp(:,:,121));

% Area, Heff
Ice=tmp(:,:,31);
Ice(find(Ice>-1.8))=0;
Ice(find(Ice))=1;
writebin('../input/Area.0005184000',Ice);
writebin('../input/Heff.0005184000',Ice);

% DIC, Alk, PO4, DOP, O2
fin=[pin 'pickup_ptracers.0005184000.data'];
tmp=readbin(fin,[128,64,75],1,'real*8');
writebin('../input/DIC.0005184000',tmp(:,:,1:15));
writebin('../input/Alk.0005184000',tmp(:,:,16:30));
writebin('../input/PO4.0005184000',tmp(:,:,31:45));
writebin('../input/DOP.0005184000',tmp(:,:,46:60));
writebin('../input/O2.0005184000',tmp(:,:,61:75));
