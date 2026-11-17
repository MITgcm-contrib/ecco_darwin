#include "SEAICE_OPTIONS.h"

CBOP
C     !ROUTINE: SEAICE_OCNVEL_AVG
C     !INTERFACE:
      SUBROUTINE SEAICE_OCNVEL_AVG( myThid )
C     !DESCRIPTION:
C     Nares 2026-11-15 (Dustin: fix (a)). Thickness-weighted mean of the ocean
C     velocity over the top SIdragDepth metres (dz x hFacW/hFacS), stored in
C     uOcnSI/vOcnSI and used for the ice-ocean drag instead of level 1.
C     Layers partly inside the depth range count with their overlap only.
C     Cells with no wet water in the range (e.g. under the ice shelf) keep
C     the level-1 value, as in the original code. r* stretching is ignored.
C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "SEAICE_OCNVEL.h"
C     !INPUT PARAMETERS:
      INTEGER myThid
CEOP
      INTEGER i, j, k, bi, bj
      _RL zTop, zBot, dz, w
      _RL sU(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sV(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL hU(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL hV(1-OLx:sNx+OLx,1-OLy:sNy+OLy)

      IF ( usingPCoords ) THEN
       STOP 'SEAICE_OCNVEL_AVG: z-coordinates only'
      ENDIF
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          sU(i,j) = 0. _d 0
          sV(i,j) = 0. _d 0
          hU(i,j) = 0. _d 0
          hV(i,j) = 0. _d 0
         ENDDO
        ENDDO
        DO k=1,Nr
         zTop = rF(1) - rF(k)
         zBot = rF(1) - rF(k+1)
         IF ( zTop .LT. SIdragDepth ) THEN
          dz = MIN( zBot, SIdragDepth ) - zTop
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            w = dz*hFacW(i,j,k,bi,bj)
            sU(i,j) = sU(i,j) + w*uVel(i,j,k,bi,bj)
            hU(i,j) = hU(i,j) + w
            w = dz*hFacS(i,j,k,bi,bj)
            sV(i,j) = sV(i,j) + w*vVel(i,j,k,bi,bj)
            hV(i,j) = hV(i,j) + w
           ENDDO
          ENDDO
         ENDIF
        ENDDO
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          IF ( hU(i,j) .GT. 0. _d 0 ) THEN
           uOcnSI(i,j,bi,bj) = sU(i,j)/hU(i,j)
          ELSE
           uOcnSI(i,j,bi,bj) = uVel(i,j,1,bi,bj)
          ENDIF
          IF ( hV(i,j) .GT. 0. _d 0 ) THEN
           vOcnSI(i,j,bi,bj) = sV(i,j)/hV(i,j)
          ELSE
           vOcnSI(i,j,bi,bj) = vVel(i,j,1,bi,bj)
          ENDIF
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      RETURN
      END
