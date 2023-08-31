to fix volume drift, need to convert UVELMASS to UVEL

diagnostics_fill_state.F
uVelMass = uVel * hFacW

update_r_star.F
hFacW = h0FacW * rStarFacW

calc_r_star.F
      rStarAreaWeight = .TRUE.
C-    Area-weighted average consistent with KE (& vert. advection):
      IF ( vectorInvariantMomentum .AND.
     &     (selectKEscheme.EQ.1 .OR. selectKEscheme.EQ.3)
     &   ) rStarAreaWeight =.FALSE.

       IF ( rStarAreaWeight ) THEN
           tmpfldW = rSurfW(i,j,bi,bj) - rLowW(i,j,bi,bj)
           rStarFacW(i,j,bi,bj) =
     &       ( 0.5 _d 0 *( etaFld(i-1,j,bi,bj)*rA(i-1,j,bi,bj)
     &                    +etaFld(i,j,bi,bj)*rA(i,j,bi,bj)
     &                   )*recip_rAw(i,j,bi,bj)
     &        +tmpfldW )/tmpfldW
       ELSE
C-     Simple average
           tmpfldW = rSurfW(i,j,bi,bj) - rLowW(i,j,bi,bj)
           rStarFacW(i,j,bi,bj) =
     &       ( 0.5 _d 0 *( etaFld(i-1,j,bi,bj) + etaFld(i,j,bi,bj) )
     &        +tmpfldW )/tmpfldW

which, for "rStarAreaWeight=.TRUE." becomes:
rStarFacW(i,j) = ( ( eta(i-1,j)*rA(i-1,j) + eta(i,j)*rA(i,j) )
                   / rAw(i,j) / 2 + Depth(i,j) )
                   / Depth(i,j)

% {{{ Horizontal velocity
RAC  =readbin([pin '/grid/RAC'   suf1],[nx ny]);
RAS  =readbin([pin '/grid/RAS'   suf1],[nx ny]);
RAW  =readbin([pin '/grid/RAW'   suf1],[nx ny]);
Depth=readbin([pin '/grid/Depth' suf1],[nx ny]);
hFacS=readbin([pin '/grid/hFacS' suf2],[nx ny nz]);
hFacW=readbin([pin '/grid/hFacW' suf2],[nx ny nz]);
fne=dir([pin 'ETAN/ETAN*01T000000']);

% {{{ U
disp('U')
fnm=dir([pin 'U/U*01T000000']);
for t=1:length(fnm)
    eta=readbin([fne(t).folder '/' fne(t).name],[nx ny]);
    tmp=readbin([fnm(t).folder '/' fnm(t).name],[nx ny nz]);
    
    % western boundary condition
    fout=[region_name suf2 '_U_West'];
    rStarFac=( ( eta(1,:).*RAC(1,:) + eta(2,:).*RAC(2,:) ) ...
                ./ RAW(2,:) / 2 + Depth(2,:) ) ./ Depth(2,:);
    hFac=squeeze(hFacW(2,:,:));
    for k=1:nz
        hFac(:,k)=hFac(:,k).*rStarFac';
    end
    obc=squeeze(tmp(2,:,:));
    in=find(hFac);
    obc(in)=obc(in)./hFac(in);
    writebin(fout,obc,1,prec,t-1)
