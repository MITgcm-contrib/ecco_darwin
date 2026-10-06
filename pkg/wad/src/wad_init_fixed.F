#include "WAD_OPTIONS.h"
#ifdef ALLOW_SEAICE
# include "SEAICE_OPTIONS.h"
#endif

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: WAD_INIT_FIXED

C !INTERFACE:
      SUBROUTINE WAD_INIT_FIXED( myThid )

C     !DESCRIPTION:
C     Save the static (topography-only) W and S masks, before any WAD
C     face closing, and register WAD diagnostics.

C     !USES:
      IMPLICIT NONE
C     == Global variables ===
#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#include "GRID.h"
#include "WAD.h"
#ifdef ALLOW_SEAICE
# include "SEAICE_SIZE.h"
# include "SEAICE_GRID.h"
#endif

C     !INPUT/OUTPUT PARAMETERS:
C     myThid ::  my Thread Id number
      INTEGER myThid
CEOP

#ifdef ALLOW_WAD
C     !LOCAL VARIABLES:
      INTEGER i, j, k, bi, bj

      DO bj=myByLo(myThid), myByHi(myThid)
       DO bi=myBxLo(myThid), myBxHi(myThid)
        DO k=1,Nr
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           wadMaskW0(i,j,k,bi,bj) = maskW(i,j,k,bi,bj)
           wadMaskS0(i,j,k,bi,bj) = maskS(i,j,k,bi,bj)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

C-    extra Manning roughness (river ice, ice jams)
      DO bj=myByLo(myThid), myByHi(myThid)
       DO bi=myBxLo(myThid), myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          wadManEx(i,j,bi,bj) = 0. _d 0
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      IF ( wadManningFile.NE.' ' ) THEN
        CALL READ_FLD_XY_RS( wadManningFile, ' ', wadManEx, 0, myThid )
        CALL EXCH_XY_RS( wadManEx, myThid )
      ENDIF

#ifdef ALLOW_SEAICE
C-    static sea-ice velocity masks (SEAICE_INIT_FIXED is called before)
      DO bj=myByLo(myThid), myByHi(myThid)
       DO bi=myBxLo(myThid), myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          wadSIMaskU0(i,j,bi,bj) = oneRS
          wadSIMaskV0(i,j,bi,bj) = oneRS
          IF ( useSEAICE ) THEN
           wadSIMaskU0(i,j,bi,bj) = SIMaskU(i,j,bi,bj)
           wadSIMaskV0(i,j,bi,bj) = SIMaskV(i,j,bi,bj)
          ENDIF
         ENDDO
        ENDDO
       ENDDO
      ENDDO
#endif /* ALLOW_SEAICE */

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN
        CALL WAD_DIAGNOSTICS_INIT( myThid )
      ENDIF
#endif

#endif /* ALLOW_WAD */

      RETURN
      END
