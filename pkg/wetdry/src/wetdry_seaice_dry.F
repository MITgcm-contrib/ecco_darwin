#include "WETDRY_OPTIONS.h"
#ifdef ALLOW_SEAICE
# include "SEAICE_OPTIONS.h"
#endif

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: WETDRY_SEAICE_DRY

C !INTERFACE:
      SUBROUTINE WETDRY_SEAICE_DRY( afterGrowth, myThid )

C !DESCRIPTION:
C     No sea-ice thermodynamics on dry columns, consistent with
C     S/R WETDRY_DRY_FORCING. SEAICE_GROWTH exchanges heat, salt and fresh
C     water between ice and ocean; WETDRY_DRY_FORCING (called later, in
C     EXTERNAL_FORCING_SURF) scales the ocean side by the column factor
C       fD = MIN(1, MAX(0, (D - wdCritDepth)/wdCritDepth)),
C     so the ice side must be scaled by the same factor. Without this a
C     supercooled film on a dry column makes ice at every step while
C     the film never warms (runaway ice thickness on exposed flats).
C     Called with afterGrowth=.FALSE. just before SEAICE_GROWTH (saves
C     HEFF, AREA, HSNOW, HSALT and EmPmR) and with afterGrowth=.TRUE.
C     just after it, which sets X = X0 + fD*(X - X0). Ice on a dry
C     column neither grows nor melts; EmPmR there keeps its value from
C     before SEAICE_GROWTH (atmospheric fresh water only).
C     With wdIceTopMelt, ice that thinned during the step keeps the whole
C     change (grounded ice melting from the atmosphere, its meltwater in
C     EmPmR); only growth is ramped off.
C     With wdTopMeltSheet, the melt of ice frozen to the bed (wdBond = 1)
C     stays on the ice: EmPmR keeps its value from before SEAICE_GROWTH
C     and the melted ice + snow mass goes to wdMeltSh, which the next
C     OVERFLOOD_STEP adds to the sheet (ice salinity SEAICE_salt0 must
C     be 0, checked in OVERFLOOD_CHECK). With wdTopMeltFloat, floating
C     ice keeps its top melt the same way (melt ponds on all ice).

C !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SURFACE.h"
#include "DYNVARS.h"
#include "FFIELDS.h"
#include "WETDRY.h"
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

#if ( defined ALLOW_WETDRY && defined ALLOW_SEAICE )
C !LOCAL VARIABLES:
      INTEGER i, j, bi, bj
      _RL fD, dm

      IF ( .NOT.( useWETDRY .AND. wdDryForcing ) ) RETURN

      DO bj=myByLo(myThid), myByHi(myThid)
       DO bi=myBxLo(myThid), myBxHi(myThid)
        IF ( .NOT.afterGrowth ) THEN
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           wdSIheff(i,j,bi,bj) = HEFF(i,j,bi,bj)
           wdSIarea(i,j,bi,bj) = AREA(i,j,bi,bj)
           wdSIsnow(i,j,bi,bj) = HSNOW(i,j,bi,bj)
           wdSIempr(i,j,bi,bj) = EmPmR(i,j,bi,bj)
#ifdef SEAICE_VARIABLE_SALINITY
           wdSIsalt(i,j,bi,bj) = HSALT(i,j,bi,bj)
#endif
          ENDDO
         ENDDO
        ELSE
         DO j=1,sNy
          DO i=1,sNx
           IF ( wdTopMeltSheet ) THEN
            wdMeltSh(i,j,bi,bj) = zeroRL
            IF ( ( wdBond(i,j,bi,bj).GE.oneRL .OR. wdTopMeltFloat )
     &     .AND. HEFF(i,j,bi,bj).LT.wdSIheff(i,j,bi,bj) ) THEN
             dm = ( wdSIheff(i,j,bi,bj) - HEFF(i,j,bi,bj) )
     &          *SEAICE_rhoIce
     &          + MAX( zeroRL, wdSIsnow(i,j,bi,bj) - HSNOW(i,j,bi,bj) )
     &           *SEAICE_rhoSnow
             dm = MIN( dm, wdMeltAtm(i,j,bi,bj) )
             wdMeltSh(i,j,bi,bj) = dm/deltaTClock
C-    only the atmospheric melt leaves the ocean side
             EmPmR(i,j,bi,bj) = EmPmR(i,j,bi,bj) + wdMeltSh(i,j,bi,bj)
            ENDIF
           ENDIF
           fD = MIN( oneRL, MAX( zeroRL,
     &        ( etaH(i,j,bi,bj) + Ro_surf(i,j,bi,bj) - R_low(i,j,bi,bj)
     &          - wdCritDepth )/wdCritDepth ) )
           IF ( wdIceTopMelt .AND.
     &          HEFF(i,j,bi,bj).LT.wdSIheff(i,j,bi,bj) ) fD = oneRL
           IF ( fD.LT.oneRL ) THEN
            HEFF(i,j,bi,bj)  = wdSIheff(i,j,bi,bj)
     &           + fD*( HEFF(i,j,bi,bj)  - wdSIheff(i,j,bi,bj) )
            AREA(i,j,bi,bj)  = wdSIarea(i,j,bi,bj)
     &           + fD*( AREA(i,j,bi,bj)  - wdSIarea(i,j,bi,bj) )
            HSNOW(i,j,bi,bj) = wdSIsnow(i,j,bi,bj)
     &           + fD*( HSNOW(i,j,bi,bj) - wdSIsnow(i,j,bi,bj) )
            EmPmR(i,j,bi,bj) = wdSIempr(i,j,bi,bj)
     &           + fD*( EmPmR(i,j,bi,bj) - wdSIempr(i,j,bi,bj) )
#ifdef SEAICE_VARIABLE_SALINITY
            HSALT(i,j,bi,bj) = wdSIsalt(i,j,bi,bj)
     &           + fD*( HSALT(i,j,bi,bj) - wdSIsalt(i,j,bi,bj) )
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
#ifdef ALLOW_DIAGNOSTICS
      IF ( afterGrowth .AND. wdTopMeltSheet .AND. useDiagnostics )
     &  CALL DIAGNOSTICS_FILL( wdMeltSh, 'WDmeltS ', 0, 1, 0, 1, 1,
     &                         myThid )
#endif
#endif /* ALLOW_WETDRY and ALLOW_SEAICE */

      RETURN
      END
