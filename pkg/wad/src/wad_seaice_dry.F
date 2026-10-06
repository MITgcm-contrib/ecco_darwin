#include "WAD_OPTIONS.h"
#ifdef ALLOW_SEAICE
# include "SEAICE_OPTIONS.h"
#endif

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: WAD_SEAICE_DRY

C !INTERFACE:
      SUBROUTINE WAD_SEAICE_DRY( afterGrowth, myThid )

C !DESCRIPTION:
C     No sea-ice thermodynamics on dry columns, consistent with
C     S/R WAD_DRY_FORCING. SEAICE_GROWTH exchanges heat, salt and fresh
C     water between ice and ocean; WAD_DRY_FORCING (called later, in
C     EXTERNAL_FORCING_SURF) scales the ocean side by the column factor
C       fD = MIN(1, MAX(0, (D - wadCritDepth)/wadCritDepth)),
C     so the ice side must be scaled by the same factor. Without this a
C     supercooled film on a dry column makes ice at every step while
C     the film never warms (runaway ice thickness on exposed flats).
C     Called with afterGrowth=.FALSE. just before SEAICE_GROWTH (saves
C     HEFF, AREA, HSNOW, HSALT and EmPmR) and with afterGrowth=.TRUE.
C     just after it, which sets X = X0 + fD*(X - X0). Ice on a dry
C     column neither grows nor melts; EmPmR there keeps its value from
C     before SEAICE_GROWTH (atmospheric fresh water only).
C     With wadIceTopMelt, ice that thinned during the step keeps the whole
C     change (grounded ice melting from the atmosphere, its meltwater in
C     EmPmR); only growth is ramped off.

C !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SURFACE.h"
#include "DYNVARS.h"
#include "FFIELDS.h"
#include "WAD.h"
#ifdef ALLOW_SEAICE
# include "SEAICE_SIZE.h"
# include "SEAICE_PARAMS.h"
# include "SEAICE.h"
#endif

C !INPUT PARAMETERS:
C     afterGrowth :: .FALSE.: save the state ; .TRUE.: apply the ramp
      LOGICAL afterGrowth
      INTEGER myThid
CEOP

#if ( defined ALLOW_WAD && defined ALLOW_SEAICE )
C !LOCAL VARIABLES:
      INTEGER i, j, bi, bj
      _RL fD

      IF ( .NOT.( useWAD .AND. wadDryForcing ) ) RETURN

      DO bj=myByLo(myThid), myByHi(myThid)
       DO bi=myBxLo(myThid), myBxHi(myThid)
        IF ( .NOT.afterGrowth ) THEN
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           wadSIheff(i,j,bi,bj) = HEFF(i,j,bi,bj)
           wadSIarea(i,j,bi,bj) = AREA(i,j,bi,bj)
           wadSIsnow(i,j,bi,bj) = HSNOW(i,j,bi,bj)
           wadSIempr(i,j,bi,bj) = EmPmR(i,j,bi,bj)
#ifdef SEAICE_VARIABLE_SALINITY
           wadSIsalt(i,j,bi,bj) = HSALT(i,j,bi,bj)
#endif
          ENDDO
         ENDDO
        ELSE
         DO j=1,sNy
          DO i=1,sNx
           fD = MIN( oneRL, MAX( zeroRL,
     &        ( etaH(i,j,bi,bj) + Ro_surf(i,j,bi,bj) - R_low(i,j,bi,bj)
     &          - wadCritDepth )/wadCritDepth ) )
           IF ( wadIceTopMelt .AND.
     &          HEFF(i,j,bi,bj).LT.wadSIheff(i,j,bi,bj) ) fD = oneRL
           IF ( fD.LT.oneRL ) THEN
            HEFF(i,j,bi,bj)  = wadSIheff(i,j,bi,bj)
     &           + fD*( HEFF(i,j,bi,bj)  - wadSIheff(i,j,bi,bj) )
            AREA(i,j,bi,bj)  = wadSIarea(i,j,bi,bj)
     &           + fD*( AREA(i,j,bi,bj)  - wadSIarea(i,j,bi,bj) )
            HSNOW(i,j,bi,bj) = wadSIsnow(i,j,bi,bj)
     &           + fD*( HSNOW(i,j,bi,bj) - wadSIsnow(i,j,bi,bj) )
            EmPmR(i,j,bi,bj) = wadSIempr(i,j,bi,bj)
     &           + fD*( EmPmR(i,j,bi,bj) - wadSIempr(i,j,bi,bj) )
#ifdef SEAICE_VARIABLE_SALINITY
            HSALT(i,j,bi,bj) = wadSIsalt(i,j,bi,bj)
     &           + fD*( HSALT(i,j,bi,bj) - wadSIsalt(i,j,bi,bj) )
#endif
            IF ( useRealFreshWaterFlux ) sIceLoad(i,j,bi,bj) =
     &           HEFF(i,j,bi,bj)*SEAICE_rhoIce
     &         + HSNOW(i,j,bi,bj)*SEAICE_rhoSnow
           ENDIF
          ENDDO
         ENDDO
        ENDIF
       ENDDO
      ENDDO
#endif /* ALLOW_WAD and ALLOW_SEAICE */

      RETURN
      END
