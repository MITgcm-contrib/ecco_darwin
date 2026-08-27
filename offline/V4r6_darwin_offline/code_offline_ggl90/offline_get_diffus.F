#include "OFFLINE_OPTIONS.h"
#ifdef ALLOW_DIC
# include "DIC_OPTIONS.h"
#endif
#ifdef ALLOW_DARWIN
# include "DARWIN_OPTIONS.h"
#endif
#ifdef ALLOW_GMREDI
# include "GMREDI_OPTIONS.h"
#endif

CBOP
C     !ROUTINE: OFFLINE_GET_SURFFORCING
C     !INTERFACE:
      SUBROUTINE OFFLINE_GET_DIFFUS( myTime, myIter, myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE OFFLINE_GET_DIFFUS
C     | o Interpolate in time diffusivity fields that have
C     |   been loaded from file
C     *==========================================================*
C     *==========================================================*
C     \ev
C
C     Patched (local override) to add a GGL90diffKr interpolation
C     block, mirroring the existing KPPdiffKzS block. This routine is
C     called from do_oceanic_phys.F strictly after GGL90_CALC (and
C     KPP_CALC) run within the same timestep, so the file-read value
C     written here overwrites whatever GGL90_CALC's own prognostic
C     TKE-based computation produced from the (nonexistent, in
C     offline mode) local surface forcing -- the same mechanism that
C     already makes the KPP block below work. GGL90_CALC_DIFF (called
C     later, from calc_diffusivity.F during thermodynamics) then reads
C     this GGL90diffKr value and adds it into KappaRx exactly as it
C     would for an online run; no changes to GGL90_CALC_DIFF,
C     calc_diffusivity.F or do_oceanic_phys.F are needed. GGL90 has no
C     non-local "ghat" flux analog, so there is no GGL90 counterpart
C     to KPPghat here.

C     !USES:
      IMPLICIT NONE
C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DYNVARS.h"
#include "FFIELDS.h"
#include "OFFLINE.h"
#include "OFFLINE_SWITCH.h"
#ifdef ALLOW_GMREDI
# include "GMREDI.h"
#endif
#ifdef ALLOW_KPP
# include "KPP.h"
#endif
#ifdef ALLOW_GGL90
# include "GGL90.h"
#endif

C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     myTime  :: current time in simulation
C     myIter  :: current iteration number in simulation
C     myThid  :: my Thread Id number
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

C     !LOCAL VARIABLES:
      INTEGER i,j,k
      INTEGER bi,bj
      _RL aWght, bWght
#ifdef ALLOW_AUTODIFF
      _RL locTime
      INTEGER intimeP, intime0, intime1
#endif /* ALLOW_AUTODIFF */

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

#ifdef ALLOW_AUTODIFF
C--   Re-compute the weights (bWght, aWght) to simplify dependencies
C     (since they are not stored on tapes)
      locTime = myTime - offlineTimeOffset
      CALL GET_PERIODIC_INTERVAL(
     O                  intimeP, intime0, intime1, bWght, aWght,
     I                  offlineForcingCycle, offlineForcingPeriod,
     I                  deltaToffline, locTime, myThid )
#endif /* ALLOW_AUTODIFF */

C--   Interpolate Diffusivity Components:
      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
#ifndef ALLOW_AUTODIFF
        bWght = offline_Wght(1,bi,bj)
        aWght = offline_Wght(2,bi,bj)
#endif /* ndef ALLOW_AUTODIFF */

        IF ( Wvelfile .NE. ' '  .AND. myIter.EQ.nIter0 ) THEN
         DO k=1,Nr
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            wVel(i,j,k,bi,bj) = bWght*wvel0(i,j,k,bi,bj)
     &                        + aWght*wvel1(i,j,k,bi,bj)
           ENDDO
          ENDDO
         ENDDO
        ENDIF
        IF ( offlineLoadConvec ) THEN
         DO k=1,Nr
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            IVDConvCount(i,j,k,bi,bj) = bWght*conv0(i,j,k,bi,bj)
     &                                + aWght*conv1(i,j,k,bi,bj)
           ENDDO
          ENDDO
         ENDDO
        ENDIF
#ifdef ALLOW_GMREDI
        IF ( offlineLoadGMRedi ) THEN
         DO k=1,Nr
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            Kwx(i,j,k,bi,bj)  = bWght*gmkx0(i,j,k,bi,bj)
     &                        + aWght*gmkx1(i,j,k,bi,bj)
            Kwy(i,j,k,bi,bj)  = bWght*gmky0(i,j,k,bi,bj)
     &                        + aWght*gmky1(i,j,k,bi,bj)
            Kwz(i,j,k,bi,bj)  = bWght*gmkz0(i,j,k,bi,bj)
     &                        + aWght*gmkz1(i,j,k,bi,bj)
#ifndef GM_NON_UNITY_DIAGONAL
            Kux(i,j,k,bi,bj) = GM_isopycK
            Kvy(i,j,k,bi,bj) = GM_isopycK
#endif
           ENDDO
          ENDDO
         ENDDO
        ENDIF
#endif /* ALLOW_GMREDI */
#ifdef ALLOW_KPP
        IF ( offlineLoadKPP ) THEN
         DO k=1,Nr
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            KPPdiffKzS(i,j,k,bi,bj) = bWght*kdfs0(i,j,k,bi,bj)
     &                              + aWght*kdfs1(i,j,k,bi,bj)
C-- Note: for convenience, the array KPPghat will contain
C         the product ghat*diffKzS (and not ghat alone).
            KPPghat(i,j,k,bi,bj) = bWght*kght0(i,j,k,bi,bj)
     &                           + aWght*kght1(i,j,k,bi,bj)
           ENDDO
          ENDDO
         ENDDO
        ENDIF
#endif /* ALLOW_KPP */
#ifdef ALLOW_GGL90
        IF ( offlineLoadGGL90 ) THEN
         DO k=1,Nr
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            GGL90diffKr(i,j,k,bi,bj) = bWght*ggl90dfk0(i,j,k,bi,bj)
     &                               + aWght*ggl90dfk1(i,j,k,bi,bj)
           ENDDO
          ENDDO
         ENDDO
        ENDIF
#endif /* ALLOW_GGL90 */

C--   Interpolate surface forcing
#if ( (defined ALLOW_DIC) || (defined ALLOW_DARWIN) )
#ifdef ALLOW_OLD_VIRTUALFLUX
        IF ( SFluxFile.NE.' ' ) THEN
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            surfaceForcingS(i,j,bi,bj) = bWght*sflx0(i,j,bi,bj)
     &                                 + aWght*sflx1(i,j,bi,bj)
            surfaceForcingS(i,j,bi,bj) = surfaceForcingS(i,j,bi,bj)
     &                                  *mass2rUnit
           ENDDO
          ENDDO
        ENDIF
#endif /* ALLOW_OLD_VIRTUALFLUX */
#endif /* ALLOW_DIC or ALLOW_DARWIN */

C--    kept from older version:
c           surfaceForcingT(i,j,bi,bj) = bWght*hflx0(i,j,bi,bj)
c    &                                 + aWght*hflx1(i,j,bi,bj)
c           surfaceForcingT(i,j,bi,bj) = surfaceForcingT(i,j,bi,bj)
c    &                         *recip_Cp*mass2rUnit
c           ICEM(i,j,bi,bj) = bWght*icem0(i,j,bi,bj)
c    &                      + aWght*icem1(i,j,bi,bj)

C--   end bi,bj loops
       ENDDO
      ENDDO

      RETURN
      END
