#include "SEDIMENT_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: SEDIMENT_BED_COLUMN

C !INTERFACE:
      SUBROUTINE SEDIMENT_BED_COLUMN( i, j, bi, bj, tauB, colD, dtS,
     U                                depN, eroN, myThid )

C !DESCRIPTION:
C     Layered bed of one column (sedNlay > 0), after the ROMS Community
C     Sediment Transport Modeling System (Warner et al. 2008; Sherwood
C     et al. 2018, GMD 11, Sects. 2.3-2.6; Moriarty et al. 2017,
C     Biogeosciences 14, Eq. 2). Layers l = 1..sedNlay, top down, each
C     with a mass per class sedLay [kg/m^2] and a bulk critical stress
C     sedTcb [Pa]; layer 1 is the active layer. Each step:
C     a) the deposit depN(n) enters layer 1 at the fresh-deposit
C        critical stress sedTcbNew (mass-weighted mixing)
C     b) the active layer is set to za = max( sedZaMin, sedZaCoef*
C        (tau_b - tau_ce) ) [m] (tau_ce of the surface, dry density of
C        the bed mixture): mass is drawn up from the layers below when
C        it is thinner, and the excess is moved down when it is thicker
C        (into layer 2, or into a new layer 2 once layer 2 holds more
C        than sedLayMass, the two deepest layers then merging)
C     c) every layer's tau_cb relaxes toward the equilibrium profile
C        tau_cb,eq = exp( (ln z_m - sedTcbOffset)/sedTcbSlope ) at the
C        mass depth z_m of its middle (Sherwood et al. Eq. 4), with the
C        time scale sedTcons from below (consolidation) and sedTswell
C        from above (swelling)
C     d) erosion from the active layer only, class by class:
C        E_n = sedEroM(n)*(1 - sedPorosity)*f_n*(tau_b/tau_ce,n - 1)*dt,
C        f_n the mass fraction of class n in layer 1, limited to its mass
C        there; tau_ce,n = max( Pc*tau_cb + (1-Pc)*sedTauCE(n),
C        sedTauCE(n) ) (Sherwood et al. Eq. 6), Pc going linearly from 0
C        to 1 as the cohesive fraction of layer 1 goes from sedFcLo to
C        sedFcHi. No erosion from columns shallower than sedDryDepth.
C     sedBed(n) is updated to the sum over the layers. Mass is only moved
C     between layers, so the column budget is unchanged.

C !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "SEDIMENT.h"

C !INPUT/OUTPUT PARAMETERS:
C     tauB   :: bed shear stress [Pa]
C     colD   :: water depth of the column [m]
C     dtS    :: time step [s]
C     depN   :: deposit of each class in this step [kg/m^2]
C     eroN   :: (out) erosion of each class in this step [kg/m^2]
      INTEGER i, j, bi, bj
      _RL tauB, colD, dtS
      _RL depN(sedMaxCl), eroN(sedMaxCl)
      INTEGER myThid
CEOP

#ifdef ALLOW_SEDIMENT
C !LOCAL VARIABLES:
      INTEGER l, n, nl, L2
      _RL mb(sedMaxLay,sedMaxCl), tc(sedMaxLay), ml(sedMaxLay)
      _RL m1, ds, sm, vol, rhoD, fc, Pc, tcp, tce, za, Ma
      _RL need, xm, xn, fr, zm, zmid, teq, mDeep1, mDeep2, ero

      L2 = sedNlay
      DO l=1,L2
       tc(l) = sedTcb(i,j,l,bi,bj)
       DO n=1,sedNcl
        mb(l,n) = sedLay(i,j,l,n,bi,bj)
       ENDDO
      ENDDO

C--   a) deposition into the active layer
      m1 = 0. _d 0
      ds = 0. _d 0
      DO n=1,sedNcl
       m1 = m1 + mb(1,n)
       ds = ds + depN(n)
      ENDDO
      IF ( m1+ds.GT.zeroRL )
     &  tc(1) = ( m1*tc(1) + ds*sedTcbNew )/( m1 + ds )
      DO n=1,sedNcl
       mb(1,n) = mb(1,n) + depN(n)
      ENDDO
      m1 = m1 + ds

C--   b) active layer thickness (surface bulk critical stress, dry
C     density of the whole bed)
      sm  = 0. _d 0
      vol = 0. _d 0
      DO l=1,L2
       DO n=1,sedNcl
        sm  = sm  + mb(l,n)
        vol = vol + mb(l,n)/sedBulkRho(n)
       ENDDO
      ENDDO
      rhoD = sedBulkRho(1)
      IF ( vol.GT.zeroRL ) rhoD = sm/vol
      tce = tc(1)
      IF ( m1.GT.zeroRL ) THEN
       fc  = 0. _d 0
       tcp = 0. _d 0
       DO n=1,sedNcl
        IF ( sedCohesive(n) ) fc = fc + mb(1,n)/m1
        tcp = tcp + sedTauCE(n)*mb(1,n)/m1
       ENDDO
       Pc  = MIN( oneRL, MAX( zeroRL,
     &            ( fc - sedFcLo )/( sedFcHi - sedFcLo ) ) )
       tce = MAX( Pc*tc(1) + ( oneRL - Pc )*tcp, tcp )
      ENDIF
      za = MAX( sedZaMin, sedZaCoef*( tauB - tce ) )
      sedZa(i,j,bi,bj) = za
      Ma = za*rhoD

      IF ( L2.GT.1 .AND. m1.LT.Ma ) THEN
C-    thin active layer: draw mass up from the layers below
       need = Ma - m1
       DO l=2,L2
        ml(l) = 0. _d 0
        DO n=1,sedNcl
         ml(l) = ml(l) + mb(l,n)
        ENDDO
        IF ( need.GT.zeroRL .AND. ml(l).GT.zeroRL ) THEN
         xm = MIN( need, ml(l) )
         tc(1) = ( m1*tc(1) + xm*tc(l) )/( m1 + xm )
         IF ( xm.GE.ml(l) ) THEN
          DO n=1,sedNcl
           mb(1,n) = mb(1,n) + mb(l,n)
           mb(l,n) = 0. _d 0
          ENDDO
         ELSE
          fr = xm/ml(l)
          DO n=1,sedNcl
           xn = mb(l,n)*fr
           mb(1,n) = mb(1,n) + xn
           mb(l,n) = mb(l,n) - xn
          ENDDO
         ENDIF
         m1   = m1 + xm
         need = need - xm
        ENDIF
       ENDDO
C-    close the gaps left by emptied layers
       nl = 1
       DO l=2,L2
        ml(l) = 0. _d 0
        DO n=1,sedNcl
         ml(l) = ml(l) + mb(l,n)
        ENDDO
        IF ( ml(l).GT.zeroRL ) THEN
         nl = nl + 1
         IF ( nl.LT.l ) THEN
          tc(nl) = tc(l)
          DO n=1,sedNcl
           mb(nl,n) = mb(l,n)
           mb(l,n)  = 0. _d 0
          ENDDO
         ENDIF
        ENDIF
       ENDDO

      ELSEIF ( L2.GT.1 .AND. m1.GT.Ma ) THEN
C-    thick active layer: move the excess down into layer 2
       xm = m1 - Ma
       fr = xm/m1
       ml(2) = 0. _d 0
       DO n=1,sedNcl
        ml(2) = ml(2) + mb(2,n)
       ENDDO
       IF ( ml(2).GT.zeroRL .AND. ml(2)+xm.GT.sedLayMass
     &      .AND. L2.GT.2 ) THEN
C-    layer 2 is full: merge the two deepest layers, shift the others
C     down and start a new layer 2
        mDeep1 = 0. _d 0
        mDeep2 = 0. _d 0
        DO n=1,sedNcl
         mDeep1 = mDeep1 + mb(L2-1,n)
         mDeep2 = mDeep2 + mb(L2,n)
        ENDDO
        IF ( mDeep1+mDeep2.GT.zeroRL )
     &    tc(L2) = ( mDeep1*tc(L2-1) + mDeep2*tc(L2) )
     &            /( mDeep1 + mDeep2 )
        DO n=1,sedNcl
         mb(L2,n) = mb(L2,n) + mb(L2-1,n)
        ENDDO
        DO l=L2-1,3,-1
         tc(l) = tc(l-1)
         DO n=1,sedNcl
          mb(l,n) = mb(l-1,n)
         ENDDO
        ENDDO
        DO n=1,sedNcl
         mb(2,n) = 0. _d 0
        ENDDO
        ml(2) = 0. _d 0
       ENDIF
       tc(2) = ( ml(2)*tc(2) + xm*tc(1) )/( ml(2) + xm )
       DO n=1,sedNcl
        xn = mb(1,n)*fr
        mb(2,n) = mb(2,n) + xn
        mb(1,n) = mb(1,n) - xn
       ENDDO
       m1 = m1 - xm
      ENDIF

C--   c) consolidation and swelling toward the equilibrium profile
      zm = 0. _d 0
      DO l=1,L2
       ml(l) = 0. _d 0
       DO n=1,sedNcl
        ml(l) = ml(l) + mb(l,n)
       ENDDO
       IF ( ml(l).GT.zeroRL ) THEN
        zmid = zm + halfRL*ml(l)
        teq  = EXP( ( LOG(zmid) - sedTcbOffset )/sedTcbSlope )
        teq  = MIN( sedTcbMax, MAX( sedTcbMin, teq ) )
        IF ( tc(l).LT.teq ) THEN
         tc(l) = teq + ( tc(l) - teq )*EXP( -dtS/sedTcons )
        ELSE
         tc(l) = teq + ( tc(l) - teq )*EXP( -dtS/sedTswell )
        ENDIF
        zm = zm + ml(l)
       ENDIF
      ENDDO

C--   d) erosion from the active layer (mixed bed)
      DO n=1,sedNcl
       eroN(n) = 0. _d 0
      ENDDO
      m1 = ml(1)
      IF ( m1.GT.zeroRL .AND. colD.GT.sedDryDepth ) THEN
       fc = 0. _d 0
       DO n=1,sedNcl
        IF ( sedCohesive(n) ) fc = fc + mb(1,n)/m1
       ENDDO
       Pc = MIN( oneRL, MAX( zeroRL,
     &           ( fc - sedFcLo )/( sedFcHi - sedFcLo ) ) )
       DO n=1,sedNcl
        IF ( mb(1,n).GT.zeroRL ) THEN
         tce = MAX( Pc*tc(1) + ( oneRL - Pc )*sedTauCE(n), sedTauCE(n) )
         IF ( tauB.GT.tce ) THEN
          ero = sedEroM(n)*( oneRL - sedPorosity )*( mb(1,n)/m1 )
     &         *( tauB/tce - oneRL )*dtS
          ero = MIN( ero, mb(1,n) )
          mb(1,n) = mb(1,n) - ero
          eroN(n) = ero
         ENDIF
        ENDIF
       ENDDO
      ENDIF

C--   store; the bed total of each class is the sum over the layers
      DO n=1,sedNcl
       sedBed(i,j,n,bi,bj) = 0. _d 0
      ENDDO
      DO l=1,L2
       sedTcb(i,j,l,bi,bj) = tc(l)
       DO n=1,sedNcl
        sedLay(i,j,l,n,bi,bj) = mb(l,n)
        sedBed(i,j,n,bi,bj)   = sedBed(i,j,n,bi,bj) + mb(l,n)
       ENDDO
      ENDDO
#endif /* ALLOW_SEDIMENT */

      RETURN
      END
