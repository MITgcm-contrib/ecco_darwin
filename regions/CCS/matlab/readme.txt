to fix volume drift, need to convert UVELMASS to UVEL

diagnostics_fill_state.F
UVELMASS = uVel * hFacW

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

