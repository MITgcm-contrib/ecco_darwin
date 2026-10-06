#include "WAD_OPTIONS.h"
#ifdef ALLOW_SEAICE
# include "SEAICE_OPTIONS.h"
#endif

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: WAD_INIT_VARIA

C !INTERFACE:
      SUBROUTINE WAD_INIT_VARIA( myThid )

C     !DESCRIPTION:
C     Initialize pkg/wad variables and apply the film floor to the
C     initial free surface. The initial face mask is computed later from
C     etaH (S/R WAD_UPDATE_MASKS, called from INITIALISE_VARIA just
C     before UPDATE_R_STAR).

C     !USES:
      IMPLICIT NONE
C     == Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SURFACE.h"
#include "DYNVARS.h"
#include "WAD.h"
#ifdef ALLOW_SEAICE
# include "SEAICE_SIZE.h"
# include "SEAICE_GRID.h"
#endif

C     !INPUT/OUTPUT PARAMETERS:
C     myThid ::  my Thread Id number
      INTEGER myThid
CEOP

#ifdef ALLOW_WAD
C     !LOCAL VARIABLES:
      INTEGER i, j, bi, bj, ks

      DO bj=myByLo(myThid), myByHi(myThid)
       DO bi=myBxLo(myThid), myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          wadMaskW(i,j,bi,bj) = oneRS
          wadMaskS(i,j,bi,bj) = oneRS
          wadDryC (i,j,bi,bj) = zeroRS
C-    no Qsw seen yet: the first WAD_DRY_FORCING call takes Qsw as is
          wadQswIn (i,j,bi,bj) = zeroRS
          wadQswOut(i,j,bi,bj) = -1. _d 30
          wadRelaxW(i,j,bi,bj) = oneRL
          wadRelaxS(i,j,bi,bj) = oneRL
          wadAdvFacW(i,j,bi,bj) = oneRL
          wadAdvFacS(i,j,bi,bj) = oneRL
          wadIceFacW(i,j,bi,bj) = oneRL
          wadIceFacS(i,j,bi,bj) = oneRL
          wadSIsolU0(i,j,bi,bj) = oneRS
          wadSIsolV0(i,j,bi,bj) = oneRS
#if (defined ALLOW_SEAICE) && (defined SEAICE_CGRID)
C-    static ice-solver masks (set in SEAICE_INIT_VARIA, called before)
          IF ( useSEAICE ) THEN
           wadSIsolU0(i,j,bi,bj) = seaiceMaskU(i,j,bi,bj)
           wadSIsolV0(i,j,bi,bj) = seaiceMaskV(i,j,bi,bj)
          ENDIF
#endif
          wadOutFac(i,j,bi,bj) = oneRL
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      _BEGIN_MASTER(myThid)
      wadVolAdj     = 0. _d 0
      wadVolAdjStep = 0. _d 0
      wadCumIn      = 0. _d 0
      wadNflips     = 0
      wadNlimit     = 0
      _END_MASTER(myThid)
      _BARRIER

#ifdef NONLIN_FRSURF
C--   Surface-level NLFS: INI_NLFS_VARS keeps the surface above the bed
C     of shallower neighbours (+hFacInf*drF), which would stop a channel
C     from draining below the level of the adjacent flats. With WAD the
C     faces take care of that: keep only the column's own limit.
      IF ( select_rStar.EQ.0 .AND. nonlinFreeSurf.GT.0 ) THEN
       DO bj=myByLo(myThid), myByHi(myThid)
        DO bi=myBxLo(myThid), myBxHi(myThid)
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           ks = kSurfC(i,j,bi,bj)
           IF ( ks.LE.Nr ) Rmin_surf(i,j,bi,bj) =
     &       MAX( rF(ks+1), R_low(i,j,bi,bj) ) + hFacInf*drF(ks)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDIF
#endif /* NONLIN_FRSURF */

C--   Initial eta can be below the bed (e.g. uniform eta over a beach):
C     raise it to the film
      CALL WAD_APPLY_FLOOR( .FALSE., startTime, nIter0, myThid )

C--   Upwind side of each face (wadUpwindFace) for the initial mask
C     call: from the initial velocity; replaced by the pickup on restart
      DO bj=myByLo(myThid), myByHi(myThid)
       DO bi=myBxLo(myThid), myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          ks = 1
          IF ( select_rStar.EQ.0 ) ks = MIN( kSurfW(i,j,bi,bj), Nr )
          wadUpW(i,j,bi,bj) = 0. _d 0
          IF ( uVel(i,j,ks,bi,bj).GT.zeroRL ) wadUpW(i,j,bi,bj) = 1.
          IF ( uVel(i,j,ks,bi,bj).LT.zeroRL ) wadUpW(i,j,bi,bj) = -1.
          ks = 1
          IF ( select_rStar.EQ.0 ) ks = MIN( kSurfS(i,j,bi,bj), Nr )
          wadUpS(i,j,bi,bj) = 0. _d 0
          IF ( vVel(i,j,ks,bi,bj).GT.zeroRL ) wadUpS(i,j,bi,bj) = 1.
          IF ( vVel(i,j,ks,bi,bj).LT.zeroRL ) wadUpS(i,j,bi,bj) = -1.
         ENDDO
        ENDDO
       ENDDO
      ENDDO

C--   Restart: restore the volume-budget scalars and the upwind sides
C     (masks are recomputed)
      IF ( .NOT.( startTime.EQ.baseTime .AND. nIter0.EQ.0
     &            .AND. pickupSuff.EQ.' ' ) ) THEN
        CALL WAD_READ_PICKUP( nIter0, myThid )
      ENDIF

#endif /* ALLOW_WAD */

      RETURN
      END
