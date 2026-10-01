#include "WAD_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: WAD_READPARMS

C !INTERFACE:
      SUBROUTINE WAD_READPARMS( myThid )

C     !DESCRIPTION:
C     Read pkg/wad (wetting and drying) parameters from data.wad

C     !USES:
      IMPLICIT NONE
C     == Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "WAD.h"

C     !INPUT PARAMETERS:
      INTEGER myThid
CEOP

#ifdef ALLOW_WAD
C     !LOCAL VARIABLES:
C     msgBuf     :: Informational/error message buffer
C     iUnit      :: Work variable for IO unit number
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER iUnit

      NAMELIST /WAD_PARM01/
     &       wadMinDepth, wadCritDepth, wadDragDepth, wadAdvDepth,
     &       wadGMDepth,
     &       wadKPPDepth, wadKPPCapDepth,
     &       wadSmoothWidth, wadMonFreq, wadMaxSpeed,
     &       wadManningN, wadDragMax, wadMaxFroude,
     &       wadManningFile, wadManningT1, wadManningT2,
     &       wadConserveVol, wadCarryVel,
     &       wadUpwindFace, wadDryForcing

      IF ( .NOT.useWAD ) THEN
C-    pkg WAD is not used
        _BEGIN_MASTER(myThid)
C-    Track pkg activation status:
C     print a (weak) warning if data.wad is found
         CALL PACKAGES_UNUSED_MSG( 'useWAD', ' ', ' ' )
        _END_MASTER(myThid)
        RETURN
      ENDIF

      _BEGIN_MASTER(myThid)

C--   Default values for WAD
      wadMinDepth    = 0.05 _d 0
      wadCritDepth   = 0.10 _d 0
C-    off by default: the relaxation strands water on receding beaches
      wadDragDepth   = 0. _d 0
C-    momentum advection is unstable in thin water on 2-D flats
      wadAdvDepth    = 0.5 _d 0
C-    GM/Redi makes large tracer extrema in thin r* layers
      wadGMDepth     = 2. _d 0
C-    optional: the KPP non-local flux fraction is capped at 1 in any
C     case (WAD_KPP_TAPER), which removes the T/S extrema it made in
C     thin r* columns
      wadKPPDepth    = 0. _d 0
      wadKPPCapDepth = 20. _d 0
      wadSmoothWidth = 0. _d 0
C-    depth-dependent (Manning) bottom drag: off
      wadManningN    = 0. _d 0
      wadDragMax     = 0.05 _d 0
      wadMaxFroude   = 0. _d 0
      wadManningFile = ' '
      wadManningT1   = 0. _d 0
      wadManningT2   = 0. _d 0
      wadMonFreq     = monitorFreq
C-    coastal flows never reach this (transient fronts with wadCarryVel
C     reach ~10 m/s); a thin-water instability passes it within steps
      wadMaxSpeed    = 30. _d 0
      wadConserveVol = .TRUE.
      wadCarryVel    = .FALSE.
      wadDryForcing  = .TRUE.
C-    the minimum (NLFS) or mean (r*) surface throttles the shoreline flux
      wadUpwindFace  = .TRUE.

      WRITE(msgBuf,'(A)') 'WAD_READPARMS: opening data.wad'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      CALL OPEN_COPY_DATA_FILE(
     I                     'data.wad', 'WAD_READPARMS',
     O                     iUnit,
     I                     myThid )

C     Read parameters from open data file
      READ(UNIT=iUnit,NML=WAD_PARM01)
      WRITE(msgBuf,'(A)')
     &    'WAD_READPARMS: finished reading data.wad'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
C     Close the open data file
#ifdef SINGLE_DISK_IO
      CLOSE(iUnit)
#else
      CLOSE(iUnit,STATUS='DELETE')
#endif /* SINGLE_DISK_IO */

      _END_MASTER(myThid)

C--   Everyone else must wait for the parameters to be loaded
      _BARRIER

#endif /* ALLOW_WAD */

      RETURN
      END
