#include "WETDRY_OPTIONS.h"
#ifdef ALLOW_SEAICE
# include "SEAICE_OPTIONS.h"
#endif

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: WETDRY_SEAICE_FREESLIP

C !INTERFACE:
      SUBROUTINE WETDRY_SEAICE_FREESLIP( etaZ, zetaZ, myThid )

C !DESCRIPTION:
C     Free slip at the ice faces pkg/wetdry closes (dry faces, wdIceSeal,
C     wdIceGround, ice thinner than wdIceMinHeff), as at land.
C     With SEAICE_no_slip=.FALSE. S/R SEAICE_CALC_VISCOSITIES zeroes the
C     shear viscosity at Z points with the static mask HEFFM only. At a
C     face pkg/wetdry closes all four HEFFM are 1, so etaZ stays the average
C     of the four cells, including the stiff grounded fast ice: the partly
C     covered cell next to it feels a no-slip wall with the fast ice's
C     viscosity, which its own weak ice and area-scaled drag cannot
C     balance (Colville freshet: 2-10 m/s ice next to the fast-ice edge).
C     Here etaZ and zetaZ are also zeroed at Z points next to a closed
C     face (seaiceMaskU/V, set in S/R WETDRY_UPDATE_MASKS).
C     Called from S/R SEAICE_CALC_VISCOSITIES (free-slip branch).

C !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#ifdef ALLOW_SEAICE
# include "SEAICE_SIZE.h"
# include "SEAICE_GRID.h"
#endif

C !INPUT/OUTPUT PARAMETERS:
C     etaZ, zetaZ :: shear and bulk viscosity at Z points (masked here)
      _RL etaZ (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL zetaZ(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      INTEGER myThid
CEOP

#if ( defined ALLOW_WETDRY && defined ALLOW_SEAICE && defined SEAICE_CGRID )
C !LOCAL VARIABLES:
      INTEGER i, j, bi, bj
      _RL maskZ

      IF ( .NOT.useWETDRY ) RETURN

      DO bj=myByLo(myThid), myByHi(myThid)
       DO bi=myBxLo(myThid), myBxHi(myThid)
        DO j=1-OLy+1,sNy+OLy-1
         DO i=1-OLx+1,sNx+OLx-1
          maskZ = seaiceMaskU(i,j,bi,bj)*seaiceMaskU(i,j-1,bi,bj)
     &          * seaiceMaskV(i,j,bi,bj)*seaiceMaskV(i-1,j,bi,bj)
          etaZ (i,j,bi,bj) = etaZ (i,j,bi,bj)*maskZ
          zetaZ(i,j,bi,bj) = zetaZ(i,j,bi,bj)*maskZ
         ENDDO
        ENDDO
       ENDDO
      ENDDO
#endif /* ALLOW_WETDRY and ALLOW_SEAICE and SEAICE_CGRID */

      RETURN
      END
