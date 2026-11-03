#include "DEBUG_OPTIONS.h"

CBOP
C     !ROUTINE: DEBUG_MSG
C     !INTERFACE:
      SUBROUTINE DEBUG_MSG(
     I                text,
     I                myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE DEBUG_MSG
C     | o Prints to STDOUT the text argument
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE

C     === Global data ===
#include "SIZE.h"
#include "EEPARAMS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     text      :: test to write
C     myThid    :: my Thread Id number
      CHARACTER*(*) text
      INTEGER myThid

C     !LOCAL VARIABLES:
C     msgBuf    :: Informational/error message buffer
      CHARACTER*(MAX_LEN_MBUF) msgBuf
CEOP

      _BEGIN_MASTER( myThid )
      WRITE(msgBuf,'(A11,A60)') 'DEBUG_MSG: ',text
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
C-    Kerguelen STEP4 smoke-test diagnostic (2026-08-07): local override
C     of the upstream pkg/debug/debug_msg.F, enabling the FLUSH the
C     original leaves commented out ("not standard f77"). Needed because
C     FORT_BUFFERED=1 on this system otherwise delays STDOUT writes well
C     behind actual execution, hiding exactly where the model is when a
C     job dies mid-run. ifort supports FLUSH as an extension. Not a
C     permanent behavior change -- only meaningful when debugMode=.TRUE.
C     in eedata, which is off by default.
      CALL FLUSH(standardMessageUnit)
      _END_MASTER( myThid )

      RETURN
      END
