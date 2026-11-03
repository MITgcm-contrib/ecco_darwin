#include "SEDIMENT_OPTIONS.h"

CBOP
C !ROUTINE: SEDIMENT_READ_PICKUP
C !INTERFACE:
      SUBROUTINE SEDIMENT_READ_PICKUP( myIter, myThid )
C !DESCRIPTION:
C     Read the bed and cumulative budgets written by SEDIMENT_WRITE_PICKUP
C !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "SEDIMENT.h"
C !INPUT PARAMETERS:
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef ALLOW_SEDIMENT
C !LOCAL VARIABLES:
      CHARACTER*(MAX_LEN_FNAM) fn, filNam
      CHARACTER*(10) suff
      LOGICAL fileExist
      INTEGER n, l, i, j, bi, bj
      _RL tmp2d(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL cum(3*sedMaxCl)

      IF ( pickupSuff.EQ.' ' ) THEN
        IF ( rwSuffixType.EQ.0 ) THEN
          WRITE(suff,'(I10.10)') myIter
        ELSE
          CALL RW_GET_SUFFIX( suff, startTime, myIter, myThid )
        ENDIF
      ELSE
        suff = pickupSuff
      ENDIF
      WRITE(fn,'(2A)') 'pickup_sediment.', suff
      CALL READ_REC_3D_RL( fn, precFloat64, sedMaxCl,
     &                     sedBed, 1, myIter, myThid )
C-    fluff part of the bed (absent in pickups written before it existed:
C     then the whole bed is consolidated)
      WRITE(fn,'(2A)') 'pickup_sedfluff.', suff
      CALL MDS_CHECK4FILE( fn, '.data', 'SEDIMENT_READ_PICKUP',
     &                     filNam, fileExist, .TRUE., myThid )
      IF ( fileExist ) CALL READ_REC_3D_RL( fn, precFloat64, sedMaxCl,
     &                     sedFluff, 1, myIter, myThid )
C-    layered bed (absent: built from the bed totals)
      IF ( sedNlay.GT.0 ) THEN
       WRITE(fn,'(2A)') 'pickup_sedbed.', suff
       CALL MDS_CHECK4FILE( fn, '.data', 'SEDIMENT_READ_PICKUP',
     &                      filNam, fileExist, .TRUE., myThid )
       IF ( fileExist ) THEN
        DO n=1,sedMaxCl
         DO l=1,sedMaxLay
          CALL READ_REC_XY_RL( fn, tmp2d, (n-1)*sedMaxLay+l,
     &                         myIter, myThid )
          DO bj=myByLo(myThid), myByHi(myThid)
           DO bi=myBxLo(myThid), myBxHi(myThid)
            DO j=1-OLy,sNy+OLy
             DO i=1-OLx,sNx+OLx
              sedLay(i,j,l,n,bi,bj) = tmp2d(i,j,bi,bj)
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        WRITE(fn,'(2A)') 'pickup_sedtcb.', suff
        DO l=1,sedMaxLay
         CALL READ_REC_XY_RL( fn, tmp2d, l, myIter, myThid )
         DO bj=myByLo(myThid), myByHi(myThid)
          DO bi=myBxLo(myThid), myBxHi(myThid)
           DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
             sedTcb(i,j,l,bi,bj) = tmp2d(i,j,bi,bj)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ELSE
        CALL SEDIMENT_BED_FROM_TOTAL( myThid )
       ENDIF
      ENDIF
      CALL READ_GLVEC_RL( 'pickup_sedimentcum.', suff, cum, 3*sedMaxCl,
     &                    1, myThid )
      _BEGIN_MASTER(myThid)
      DO n=1,sedMaxCl
        sedCumRiv(n) = cum(n)
        sedCumDep(n) = cum(sedMaxCl+n)
        sedCumEro(n) = cum(2*sedMaxCl+n)
      ENDDO
      _END_MASTER(myThid)
      _BARRIER
#endif /* ALLOW_SEDIMENT */

      RETURN
      END
