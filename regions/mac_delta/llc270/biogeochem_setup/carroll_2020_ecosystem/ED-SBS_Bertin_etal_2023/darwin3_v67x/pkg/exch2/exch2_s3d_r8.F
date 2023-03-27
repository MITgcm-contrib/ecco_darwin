#include "CPP_EEOPTIONS.h"
#include "W2_OPTIONS.h"

CBOP
C     !ROUTINE: EXCH2_S3D_R8

C     !INTERFACE:
      SUBROUTINE EXCH2_S3D_R8(
     U                       phi,
     I                       myNz, myThid )

C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE EXCH2_S3D_R8
C     | o Handle Simple exchanges (= that ignore corners)
C     |   for _R8, 3-dim scalar arrays with overlap size = 1
C     *==========================================================*

C     !USES:
      IMPLICIT NONE
C     === Global data ===
#include "SIZE.h"
#include "EEPARAMS.h"
c#include "W2_EXCH2_SIZE.h"
c#include "W2_EXCH2_TOPOLOGY.h"

C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     phi    :: Array with overlap regions are to be exchanged
C     myNz   :: 3rd dimension of array to exchange
C     myThid :: My thread id.
      INTEGER myNz
      _R8 phi(0:sNx+1,0:sNy+1,myNz,nSx,nSy)
      INTEGER myThid

C     !LOCAL VARIABLES:
C     == Local variables ==
C     OL[wens]       :: Overlap extents in west, east, north, south.
C     exchWidth[XY]  :: Extent of regions that will be exchanged.
      INTEGER OLw, OLe, OLn, OLs, exchWidthX, exchWidthY

CEOP

      OLw        = 1
      OLe        = 1
      OLn        = 1
      OLs        = 1
      exchWidthX = 1
      exchWidthY = 1

      CALL EXCH2_R81_CUBE( phi, .FALSE., 'T ',
     I            OLw, OLe, OLs, OLn, myNz,
     I            exchWidthX, exchWidthY,
     I            EXCH_IGNORE_CORNERS, myThid )

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
