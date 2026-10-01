#include "WAD_OPTIONS.h"
#ifdef ALLOW_PTRACERS
# include "PTRACERS_OPTIONS.h"
#endif

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: WAD_DRY_FORCING

C !INTERFACE:
      SUBROUTINE WAD_DRY_FORCING( bi, bj, myThid )

C !DESCRIPTION:
C     No surface forcing on dry columns: the surface heat, salt and
C     fresh-water forcing (surfaceForcingT/S), the penetrating short-wave
C     (Qsw) and the passive-tracer surface forcing (surfaceForcingPTr)
C     are scaled by a column factor
C       fD = MIN(1, MAX(0, (D - wadCritDepth)/wadCritDepth))
C     (0 on dry columns, 1 once D >= 2*wadCritDepth). On a dry flat these
C     fluxes would otherwise act on the few cm of the film (salinity,
C     heat or brine spikes); physically they go into the exposed bed.
C     EmPmR is not changed (it moves volume with useRealFreshWaterFlux).
C     Called at the end of S/R EXTERNAL_FORCING_SURF.

C !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "FFIELDS.h"
#include "WAD.h"
#ifdef ALLOW_PTRACERS
# include "PTRACERS_SIZE.h"
# include "PTRACERS_PARAMS.h"
# include "PTRACERS_FIELDS.h"
#endif

C !INPUT PARAMETERS:
      INTEGER bi, bj
      INTEGER myThid
CEOP

#ifdef ALLOW_WAD
C !LOCAL VARIABLES:
      INTEGER i, j
#ifdef ALLOW_PTRACERS
      INTEGER iTr
#endif
      _RL fD(1-OLx:sNx+OLx,1-OLy:sNy+OLy)

      IF ( .NOT.( useWAD .AND. wadDryForcing ) ) RETURN

      DO j=1-OLy,sNy+OLy
       DO i=1-OLx,sNx+OLx
        fD(i,j) = MIN( oneRL, MAX( zeroRL,
     &     ( etaH(i,j,bi,bj) + Ro_surf(i,j,bi,bj) - R_low(i,j,bi,bj)
     &       - wadCritDepth )/wadCritDepth ) )
        surfaceForcingT(i,j,bi,bj) = surfaceForcingT(i,j,bi,bj)*fD(i,j)
        surfaceForcingS(i,j,bi,bj) = surfaceForcingS(i,j,bi,bj)*fD(i,j)
C-    Qsw is not reloaded every step (only with EXF or periodic forcing):
C     if it still holds our last output, scale the last input again
        IF ( Qsw(i,j,bi,bj).NE.wadQswOut(i,j,bi,bj) )
     &    wadQswIn(i,j,bi,bj) = Qsw(i,j,bi,bj)
        Qsw(i,j,bi,bj) = wadQswIn(i,j,bi,bj)*fD(i,j)
        wadQswOut(i,j,bi,bj) = Qsw(i,j,bi,bj)
       ENDDO
      ENDDO
#ifdef ALLOW_PTRACERS
      IF ( usePTRACERS ) THEN
       DO iTr=1,PTRACERS_numInUse
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          surfaceForcingPTr(i,j,bi,bj,iTr) =
     &       surfaceForcingPTr(i,j,bi,bj,iTr)*fD(i,j)
         ENDDO
        ENDDO
       ENDDO
      ENDIF
#endif /* ALLOW_PTRACERS */
#endif /* ALLOW_WAD */

      RETURN
      END
