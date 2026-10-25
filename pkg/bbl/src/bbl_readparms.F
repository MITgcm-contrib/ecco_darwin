#include "BBL_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: BBL_READPARMS

C !INTERFACE:
      SUBROUTINE BBL_READPARMS( myThid )

C     !DESCRIPTION:
C     Initialize BBL variables and constants.

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "BBL.h"

C     !INPUT PARAMETERS:
      INTEGER myThid
CEOP

#ifdef ALLOW_BBL

C     !LOCAL VARIABLES:
C     msgBuf     :: Informational/error message buffer
C     iUnit      :: Work variable for IO unit number
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER iUnit

      NAMELIST /BBL_PARM01/
     &     bbl_wvel, bbl_hvel, bbl_initEta,
     &     bbl_useNofSpeed, bbl_nofCoeff, bbl_fMin,
     &     bbl_selectEntrain, bbl_selectDetrain, bbl_uDet,
     &     bbl_gpTaper, bbl_tauMix,
     &     bbl_coldStart, bbl_doPtracers, bbl_skipIceCavity,
     &     bbl_thetaFile, bbl_saltFile, bbl_etaFile

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      IF ( .NOT.useBBL ) THEN
C-    pkg BBL is not used
        _BEGIN_MASTER(myThid)
C-    Track pkg activation status:
C     print a (weak) warning if data.bbl is found
         CALL PACKAGES_UNUSED_MSG( 'useBBL', ' ', ' ' )
        _END_MASTER(myThid)
        RETURN
      ENDIF

      _BEGIN_MASTER(myThid)

C--   Default values for BBL
      bbl_wvel      = 1. _d -5
      bbl_hvel      = 1. _d 0
      bbl_initEta   = 20. _d 0
      bbl_useNofSpeed = .FALSE.
      bbl_nofCoeff  = 1. _d 0
      bbl_fMin      = 1. _d -5
      bbl_selectEntrain = 0
      bbl_selectDetrain = 0
      bbl_uDet      = 1. _d -2
      bbl_gpTaper   = 0. _d 0
      bbl_tauMix    = 0. _d 0
      bbl_coldStart = .FALSE.
      bbl_doPtracers = usePTRACERS
      bbl_skipIceCavity = .FALSE.

C-    file names for initial conditions:
      bbl_thetaFile = ' '
      bbl_saltFile  = ' '
      bbl_etaFile   = ' '

      WRITE(msgBuf,'(A)') 'BBL_READPARMS: opening data.bbl'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT , 1)
      CALL OPEN_COPY_DATA_FILE(
     I                     'data.bbl', 'BBL_READPARMS',
     O                     iUnit,
     I                     myThid )

C     Read parameters from open data file
      READ(UNIT=iUnit,NML=BBL_PARM01)
      WRITE(msgBuf,'(A)')
     &    'BBL_READPARMS: finished reading data.bbl'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &    SQUEEZE_RIGHT , 1)
C     Close the open data file
#ifdef SINGLE_DISK_IO
      CLOSE(iUnit)
#else
      CLOSE(iUnit,STATUS='DELETE')
#endif /* SINGLE_DISK_IO */

      _END_MASTER(myThid)

C--   Everyone else must wait for the parameters to be loaded
      _BARRIER

#endif /* ALLOW_BBL */

      RETURN
      END
