#include "WAD_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: WAD_WRITE_PICKUP

C !INTERFACE:
      SUBROUTINE WAD_WRITE_PICKUP( permPickup,
     &                    suff, myTime, myIter, myThid )

C !DESCRIPTION:
C     Write the pkg/wad pickup file "pickup_wad.{suff}": the global
C     volume-budget scalars wadVolAdj, wadCumIn and wadVol0. The face
C     masks are recomputed from etaH at restart (S/R WAD_UPDATE_MASKS).

C !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "WAD.h"

C !INPUT PARAMETERS:
C     permPickup :: write a permanent pickup
C     suff       :: suffix for pickup file (eg. ckptA or 0000000010)
C     myTime     :: model time
C     myIter     :: time-step number
C     myThid     :: thread number
      LOGICAL permPickup
      CHARACTER*(*) suff
      _RL myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef ALLOW_WAD
C !LOCAL VARIABLES:
      _RL budget(3)
      CHARACTER*(MAX_LEN_FNAM) fn

      budget(1) = wadVolAdj
      budget(2) = wadCumIn
      budget(3) = wadVol0
      CALL WRITE_GLVEC_RL( 'pickup_wad.', suff, budget, 3,
     &                     myIter, myThid )
C-    upwind side of each face (wadUpwindFace), for an exact restart
      WRITE(fn,'(A,A)') 'pickup_wadface.', suff
      CALL WRITE_REC_3D_RL( fn, precFloat64, 1, wadUpW, 1,
     &                      myIter, myThid )
      CALL WRITE_REC_3D_RL( fn, precFloat64, 1, wadUpS, 2,
     &                      myIter, myThid )
#endif /* ALLOW_WAD */

      RETURN
      END
