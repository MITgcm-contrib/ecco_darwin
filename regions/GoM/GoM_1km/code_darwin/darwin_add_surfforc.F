#include "DARWIN_OPTIONS.h"
#include "EXF_OPTIONS.h"

CBOP
C !ROUTINE: DARWIN_ADD_SURFFORC
C !INTERFACE: ==========================================================
      SUBROUTINE DARWIN_ADD_SURFFORC(
     U                        gPtr,
     I                        surfForc, fac,
     I                        bi, bj, iMin, iMax, jMin, jMax,
     I                        myIter, myTime, myThid )

C !DESCRIPTION:
C     Add a per-area surface flux to a tracer tendency.
C
C     By default the whole flux goes into level 1, which is what this
C     routine has always done. Set darwinRunoffMixD in data.darwin to a
C     depth in metres and it is instead spread over the top of the water
C     column, weighted by cell thickness, with the column integral left
C     unchanged.
C
C     Why the option exists: with a 1 m top cell, the depth over which
C     a river's dissolved load is deposited is a numerical choice, not
C     a physical one. rStar spreads the river's WATER over the whole
C     column while the load lands in level 1 alone, and in a column of
C     only a few metres the two are badly mismatched. Spreading over
C     about 10 m is closer to a real buoyant plume.
C
C     gPtr is declared here over Nr levels. The callers pass
C     gPtr(1-OLx,1-OLy,1,iTrc) of a contiguous (x,y,Nr,nDarwin)
C     array, so
C     levels 2..Nr follow in memory and need no change at the call sites
C     -- the same sequence association DARWIN_ADD_BOTFORC already relies
C     on.

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "GRID.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#ifdef ALLOW_EXF
#include "EXF_INTERP_SIZE.h"
#endif
#include "DARWIN_EXF_PARAMS.h"

C!INPUT PARAMETERS: ====================================================
C  surfForc  :: surface forcing field to apply
C  fac       :: factor to multiply surface forcing by
C  myThid    :: thread number
      _RL surfForc(1-OLx:sNx+OLx, 1-OLy:sNy+OLy)
      _RL fac
      _RL myTime
      INTEGER iMin, iMax, jMin, jMax, bi, bj, myIter, myThid

C!INPUT/OUTPUT PARAMETERS: =============================================
C  gPtr      :: tracer tendency array to update
      _RL gPtr(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, Nr)
CEOP

#ifdef ALLOW_DARWIN

C!LOCAL VARIABLES: =====================================================
C  i,j,k         :: loop indices
C  darwinRunoffMixD is read from data.darwin (namelist
C                   darwin_forcing_params); 0 reproduces the
C                   level-1-only behaviour EXACTLY.
C  kSprd         :: deepest level receiving any of the flux
C  hSum          :: thickness actually used, m (>= darwinRunoffMixD,
C                   or the whole column where that is shallower)
      INTEGER i,j,k,kTop
      INTEGER kSprd(1-OLx:sNx+OLx, 1-OLy:sNy+OLy)
      _RL tmp2d(1-OLx:sNx+OLx, 1-OLy:sNy+OLy)
C  tmp3d         :: spreading-branch work array; FREESURF_RESCALE_G
C                   takes (x,y,Nr) and rescales level k, so it must not
C                   be given the 2-D tmp2d for k > 1
      _RL tmp3d(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, Nr)
      _RL hSum(1-OLx:sNx+OLx, 1-OLy:sNy+OLy)
      _RL dzk

      IF ( darwinRunoffMixD .LE. 0. _d 0 ) THEN

C===  Level 1 only: the original code, unchanged.
      DO j=jMin,jMax
       DO i=iMin,iMax
        tmp2d(i,j) = surfForc(i,j)*fac
     &             * recip_drF(1)*recip_hFacC(i,j,1,bi,bj)
       ENDDO
      ENDDO
#ifdef NONLIN_FRSURF
C-    Account for change in level thickness
      IF (nonlinFreeSurf.GT.0) THEN
        CALL FREESURF_RESCALE_G(
     I                           bi, bj, 1,
     U                           tmp2d,
     I                           myThid )
      ENDIF
#endif /* NONLIN_FRSURF */
      DO j=jMin,jMax
       DO i=iMin,iMax
        gPtr(i,j,1) = gPtr(i,j,1) + tmp2d(i,j)
       ENDDO
      ENDDO

      ELSE

C===  Spread over the top darwinRunoffMixD metres.
C     Whole cells are taken until the accumulated thickness reaches the
C     target, so hSum is the thickness really used and dividing by it
C     makes the column integral identical to the level-1-only case:
C     SUM_k ( surfForc*fac/hSum )*dz(k) = surfForc*fac.
      kTop = 1
      DO j=jMin,jMax
       DO i=iMin,iMax
        hSum(i,j)  = 0. _d 0
        kSprd(i,j) = 0
        DO k=1,Nr
         IF ( hFacC(i,j,k,bi,bj) .GT. 0. _d 0 .AND.
     &        hSum(i,j) .LT. darwinRunoffMixD ) THEN
          dzk = drF(k)*hFacC(i,j,k,bi,bj)
          hSum(i,j)  = hSum(i,j) + dzk
          kSprd(i,j) = k
         ENDIF
        ENDDO
        kTop = MAX( kTop, kSprd(i,j) )
       ENDDO
      ENDDO

      DO k=1,kTop
       DO j=jMin,jMax
        DO i=iMin,iMax
         IF ( k .LE. kSprd(i,j) .AND. hSum(i,j) .GT. 0. _d 0 ) THEN
          tmp3d(i,j,k) = surfForc(i,j)*fac/hSum(i,j)
         ELSE
          tmp3d(i,j,k) = 0. _d 0
         ENDIF
        ENDDO
       ENDDO
#ifdef NONLIN_FRSURF
C-    Account for change in level thickness
       IF (nonlinFreeSurf.GT.0) THEN
         CALL FREESURF_RESCALE_G(
     I                            bi, bj, k,
     U                            tmp3d,
     I                            myThid )
       ENDIF
#endif /* NONLIN_FRSURF */
       DO j=jMin,jMax
        DO i=iMin,iMax
         gPtr(i,j,k) = gPtr(i,j,k) + tmp3d(i,j,k)
        ENDDO
       ENDDO
      ENDDO

      ENDIF

#endif /* ALLOW_DARWIN */

      RETURN
      END
