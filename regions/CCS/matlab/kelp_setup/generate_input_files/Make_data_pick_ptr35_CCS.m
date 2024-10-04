

% Dimensions of the file
llc=80;
nf=100;
nz=48;
ptrs=31;
trcs=4;
% Reading info

typp=1;
precc='real*8';
skipp1=2*50;
mformm='ieee-be';
sk=0;
% Creating a fake array wth zeros to add to 
% the existing data.ptracers 
ptrc4=ones(llc,nf,nz,trcs);
ptrc4=ptrc4.*1e-16;

fnam='/nobackupp18/mmanizza/Kelp/CCS/run_test/pickup_ptracers.ckptA.data';
% Reading the file data.ptracers with 31 tracers

ptrc31=readbin(fnam,[llc nf nz ptrs],typp,precc,sk,mformm);

%ptrc19=ptrc31(:,:,:,1:19);
%size(ptrc31)
%ptrc2031=ptrc31(:,:,:,20:31);
%size(ptrc2031)



% ptrc31=readbin(fnam,[llc,llc*nf,nz,ptrs]);

ddim=4; % Dimesions for concat

% Appending the new additional four
% fake tracers
ptrc35=cat(ddim,ptrc31,ptrc4);
size(ptrc35)

fnamout='pickup_ptracers_35_test_CCS.0000000001.data';
%pprecc='real*8';
writebin(fnamout,ptrc35,typp,precc);








