#include "SEDIMENT_OPTIONS.h"

CBOP
C !ROUTINE: SEDIMENT_INIT_FIXED
C !INTERFACE:
      SUBROUTINE SEDIMENT_INIT_FIXED( myThid )
C !DESCRIPTION:
C     Fixed initialisation of pkg sediment: concentration of the added
C     fluid in each cell (sedRiverC or the sedRiverCFile map), and
C     register the diagnostics
C !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "SEDIMENT.h"
C !INPUT PARAMETERS:
      INTEGER myThid
CEOP

#ifdef ALLOW_SEDIMENT
C !LOCAL VARIABLES:
      INTEGER i, j, bi, bj, n
      _RL tmp2d(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      DO n=1,sedNcl
       IF ( sedRiverCFile(n).NE.' ' ) THEN
        CALL READ_FLD_XY_RL( sedRiverCFile(n), ' ', tmp2d, 0, myThid )
        CALL EXCH_XY_RL( tmp2d, myThid )
       ENDIF
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           sedRivC(i,j,bi,bj,n) = sedRiverC(n)
           IF ( sedRiverCFile(n).NE.' ' )
     &       sedRivC(i,j,bi,bj,n) = tmp2d(i,j,bi,bj)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) CALL SEDIMENT_DIAGNOSTICS_INIT( myThid )
#endif
#endif /* ALLOW_SEDIMENT */

      RETURN
      END
