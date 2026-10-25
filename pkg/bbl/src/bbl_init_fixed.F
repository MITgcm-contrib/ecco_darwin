#include "BBL_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: BBL_INIT_FIXED

C !INTERFACE:
      SUBROUTINE BBL_INIT_FIXED( myThid )

C     !DESCRIPTION:
C     Initialize BBL variables that are kept fixed during the run.

C     !USES:
      IMPLICIT NONE
#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#include "GRID.h"
#include "BBL.h"

C     !INPUT/OUTPUT PARAMETERS:
C     myThid ::  my Thread Id number
      INTEGER myThid
CEOP

C     !LOCAL VARIABLES:
C     === Local variables ===
      INTEGER bi, bj, i, j

C     Bottom level of the BBL: kLowC, or 0 (no BBL) in ice-shelf cavity
C     columns (top wet level kSurfC > 1) if bbl_skipIceCavity
      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          bbl_kBot(i,j,bi,bj) = kLowC(i,j,bi,bj)
          IF ( bbl_skipIceCavity .AND. kLowC(i,j,bi,bj).GT.0
     &         .AND. kSurfC(i,j,bi,bj).GT.1 ) bbl_kBot(i,j,bi,bj) = 0
         ENDDO
        ENDDO
       ENDDO
      ENDDO

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN
        CALL BBL_DIAGNOSTICS_INIT( myThid )
      ENDIF
#endif

      RETURN
      END
