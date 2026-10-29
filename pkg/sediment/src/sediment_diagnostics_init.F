#include "SEDIMENT_OPTIONS.h"

CBOP
C !ROUTINE: SEDIMENT_DIAGNOSTICS_INIT
C !INTERFACE:
      SUBROUTINE SEDIMENT_DIAGNOSTICS_INIT( myThid )
C !DESCRIPTION:
C     Register the sediment diagnostics: bed mass, deposition and erosion
C     fluxes of each class (SEDBEDnn, SEDDEPnn, SEDEROnn) and the bed
C     shear stress (SEDTAUB)
C !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "SEDIMENT.h"
C !INPUT PARAMETERS:
      INTEGER myThid
CEOP

#if (defined ALLOW_SEDIMENT) && (defined ALLOW_DIAGNOSTICS)
C !LOCAL VARIABLES:
      INTEGER       diagNum, n
      CHARACTER*8   diagName
      CHARACTER*16  diagCode
      CHARACTER*16  diagUnits
      CHARACTER*(80) diagTitle

      diagCode  = 'SM P    L1      '
      DO n=1,sedNcl
        WRITE(diagName,'(A,I2.2)') 'SEDBED', n
        WRITE(diagTitle,'(A,I2)') 'Bed sediment mass, class', n
        diagUnits = 'kg/m^2          '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
        WRITE(diagName,'(A,I2.2)') 'SEDFLF', n
        WRITE(diagTitle,'(A,I2)') 'Fluff (fresh deposit) mass, class', n
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
        WRITE(diagName,'(A,I2.2)') 'SEDDEP', n
        WRITE(diagTitle,'(A,I2)') 'Deposition flux to the bed, class', n
        diagUnits = 'kg/m^2/s        '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
        WRITE(diagName,'(A,I2.2)') 'SEDERO', n
        WRITE(diagTitle,'(A,I2)') 'Erosion flux from the bed, class', n
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      ENDDO
      diagName  = 'SEDTAUB '
      diagTitle = 'Bed shear stress (current, or max. wave+current)'
      diagUnits = 'Pa              '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'SEDTAUW '
      diagTitle = 'Wave bed shear stress'
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'SEDTCB  '
      diagTitle = 'Bulk critical stress of the active bed layer'
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'SEDUBW  '
      diagTitle = 'Near-bed wave orbital velocity'
      diagUnits = 'm/s             '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'SEDHSW  '
      diagTitle = 'Significant wave height (parametric waves)'
      diagUnits = 'm               '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'SEDZA   '
      diagTitle = 'Active bed layer thickness'
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'SEDACT  '
      diagTitle = 'Mass of the active bed layer (all classes)'
      diagUnits = 'kg/m^2          '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
#endif /* ALLOW_SEDIMENT and ALLOW_DIAGNOSTICS */

      RETURN
      END
