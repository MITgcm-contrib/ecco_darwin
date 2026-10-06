#include "WAD_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: WAD_READ_PICKUP

C !INTERFACE:
      SUBROUTINE WAD_READ_PICKUP( myIter, myThid )

C !DESCRIPTION:
C     Read the pkg/wad pickup file "pickup_wad.{suff}" (see
C     S/R WAD_WRITE_PICKUP) to restore the volume-budget scalars.

C !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "WAD.h"

C !INPUT PARAMETERS:
C     myIter :: time-step number of the pickup
C     myThid :: thread number
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef ALLOW_WAD
C !LOCAL VARIABLES:
      _RL budget(3)
      CHARACTER*(10) suff
      CHARACTER*(MAX_LEN_FNAM) fn, filNam
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      LOGICAL fileExist

      IF ( pickupSuff.EQ.' ' ) THEN
        IF ( rwSuffixType.EQ.0 ) THEN
          WRITE(suff,'(I10.10)') myIter
        ELSE
          CALL RW_GET_SUFFIX( suff, startTime, myIter, myThid )
        ENDIF
      ELSE
        suff = pickupSuff
      ENDIF
      CALL READ_GLVEC_RL( 'pickup_wad.', suff, budget, 3, 1, myThid )

      _BEGIN_MASTER(myThid)
      wadVolAdj = budget(1)
      wadCumIn  = budget(2)
      wadVol0   = budget(3)
      _END_MASTER(myThid)
      _BARRIER

C-    upwind side of each face; pickups written before this file existed
C     keep the side from the pickup velocity (restart not bit-exact)
      WRITE(fn,'(A,A)') 'pickup_wadface.', suff
      CALL MDS_CHECK4FILE( fn, '.data', 'WAD_READ_PICKUP',
     &                     filNam, fileExist, .TRUE., myThid )
      IF ( fileExist ) THEN
        CALL READ_REC_3D_RL( fn, precFloat64, 1, wadUpW, 1,
     &                       myIter, myThid )
        CALL READ_REC_3D_RL( fn, precFloat64, 1, wadUpS, 2,
     &                       myIter, myThid )
        CALL EXCH_UV_XY_RL( wadUpW, wadUpS, .FALSE., myThid )
      ELSE
        _BEGIN_MASTER(myThid)
        WRITE(msgBuf,'(3A)') '** WARNING ** WAD_READ_PICKUP: no ',
     &   fn(1:25), ': upwind face sides from the pickup velocity'
        CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
        _END_MASTER(myThid)
      ENDIF
C-    relaxed sea-ice load factor (wadLoadDepth); without the file it is
C     reset to its depth ramp (restart not bit-exact)
      IF ( wadLoadDepth.GT.wadCritDepth ) THEN
       WRITE(fn,'(A,A)') 'pickup_wadload.', suff
       CALL MDS_CHECK4FILE( fn, '.data', 'WAD_READ_PICKUP',
     &                      filNam, fileExist, .TRUE., myThid )
       IF ( fileExist ) THEN
        CALL READ_REC_3D_RL( fn, precFloat64, 1, wadLoadFac, 1,
     &                       myIter, myThid )
        CALL EXCH_XY_RL( wadLoadFac, myThid )
        _BEGIN_MASTER(myThid)
        wadLoadRead = .TRUE.
        _END_MASTER(myThid)
        _BARRIER
       ENDIF
      ENDIF
#endif /* ALLOW_WAD */

      RETURN
      END
