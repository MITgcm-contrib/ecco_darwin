#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

CBOP
C     !ROUTINE: SWFRAC2D
C     !INTERFACE:
      SUBROUTINE SWFRAC2D(
     I                  imax, fact,
     U                  swdk,
     I                  bi, bj, myTime, myIter, myThid )
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE SWFRAC2D
C     | o Compute solar short-wave flux penetration.
C     *==========================================================*
C     | Compute fraction of solar short-wave flux penetrating to
C     | specified depth, swdk, due to exponential decay in
C     | Jerlov water type jwtype.
C     | Reference : Two band solar absorption model of Paulson
C     |             and Simpson (1977, JPO, 7, 952-956)
C     | Notes
C     | =====
C     | Below 200m the solar penetration gets set to zero,
C     | otherwise the limit for the exponent (+/- 5678) needs to
C     | be taken care of.
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE
C     == Global data ==
#include "SIZE.h"
#include "FFIELDS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     input arguments
C     imax    :: number of vertical grid points
C     fact    :: scale  factor to apply to depth array
C     bi,bj     :: Current tile indices
C     myTime  :: Current time in simulation
C     myIter  :: Current iteration number in simulation
C     myThid  :: My Thread Id. number
      INTEGER imax
      _RL     fact
      INTEGER bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
C     input/output arguments
C     swdk    :: on input: vertical depth for desired sw fraction
C               (fact*swdk) is negative distance (m) from surface
C     swdk    :: on output: short wave (radiation) fractional decay
      _RL     swdk(0:sNx+1,0:sNy+1,imax)

C     !LOCAL VARIABLES:
C     === Local variables ===
C     max number of different water types
      INTEGER   nwtype
      PARAMETER(nwtype=5)
      _RL facz
      _RL rfac(nwtype),a1(nwtype),a2(nwtype)
      INTEGER i, j, m
CEOP

C     Jerlov water type :
C                  I          IA         IB         II         III
C     jwtype :     1          2          3          4          5
      DATA rfac / 0.58 _d 0, 0.62 _d 0, 0.67 _d 0, 0.77 _d 0, 0.78 _d 0/
      DATA a1   / 0.35 _d 0, 0.6  _d 0, 1.0  _d 0, 1.5  _d 0, 1.4  _d 0/
      DATA a2   / 23.0 _d 0, 20.0 _d 0, 17.0 _d 0, 14.0 _d 0, 7.9  _d 0/

      DO j=0,sNy+1
       DO i=0,sNx+1
        DO m = 1,imax
          facz = fact*swdk(i,j,m)
          IF ( facz .LT. -200. _d 0 ) THEN
            swdk(i,j,m) = 0. _d 0
          ELSE
            swdk(i,j,m) = rfac(NINT(waterTyp2d(i,j,bi,bj)))
     &         * exp( facz/a1(NINT(waterTyp2d(i,j,bi,bj))) )
     &         + (1. _d 0 - rfac(NINT(waterTyp2d(i,j,bi,bj))))
     &         * exp( facz/a2(NINT(waterTyp2d(i,j,bi,bj))) )
          ENDIF
        ENDDO
       ENDDO
      ENDDO

      RETURN
      END
