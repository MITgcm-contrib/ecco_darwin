#include "SEDIMENT_OPTIONS.h"

CBOP
C !ROUTINE: SEDIMENT_WRITE_PICKUP
C !INTERFACE:
      SUBROUTINE SEDIMENT_WRITE_PICKUP( permPickup,
     &                    suff, myTime, myIter, myThid )
C !DESCRIPTION:
C     Write the bed (pickup_sediment.*) and the cumulative budgets
C     (pickup_sedimentcum.*)
C !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "SEDIMENT.h"
C !INPUT PARAMETERS:
      LOGICAL permPickup
      CHARACTER*(*) suff
      _RL myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef ALLOW_SEDIMENT
C !LOCAL VARIABLES:
      CHARACTER*(MAX_LEN_FNAM) fn
      INTEGER n, ilnb, l, i, j, bi, bj
      _RL tmp2d(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      EXTERNAL ILNBLNK
      INTEGER ILNBLNK
      _RL cum(3*sedMaxCl)

      ilnb = MAX( 1, ILNBLNK(suff) )
      WRITE(fn,'(2A)') 'pickup_sediment.', suff(1:ilnb)
      CALL WRITE_REC_3D_RL( fn, precFloat64, sedMaxCl,
     &                      sedBed, 1, myIter, myThid )
C-    fluff part of the bed (separate file: older pickups stay readable)
      WRITE(fn,'(2A)') 'pickup_sedfluff.', suff(1:ilnb)
      CALL WRITE_REC_3D_RL( fn, precFloat64, sedMaxCl,
     &                      sedFluff, 1, myIter, myThid )
C-    layered bed: masses and bulk critical stresses (separate files,
C     one 2-D record per layer and class: the mdsio buffer holds 2*Nr
C     levels at most)
      IF ( sedNlay.GT.0 ) THEN
       WRITE(fn,'(2A)') 'pickup_sedbed.', suff(1:ilnb)
       DO n=1,sedMaxCl
        DO l=1,sedMaxLay
         DO bj=myByLo(myThid), myByHi(myThid)
          DO bi=myBxLo(myThid), myBxHi(myThid)
           DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
             tmp2d(i,j,bi,bj) = sedLay(i,j,l,n,bi,bj)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
         CALL WRITE_REC_XY_RL( fn, tmp2d, (n-1)*sedMaxLay+l,
     &                         myIter, myThid )
        ENDDO
       ENDDO
       WRITE(fn,'(2A)') 'pickup_sedtcb.', suff(1:ilnb)
       DO l=1,sedMaxLay
        DO bj=myByLo(myThid), myByHi(myThid)
         DO bi=myBxLo(myThid), myBxHi(myThid)
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            tmp2d(i,j,bi,bj) = sedTcb(i,j,l,bi,bj)
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        CALL WRITE_REC_XY_RL( fn, tmp2d, l, myIter, myThid )
       ENDDO
      ENDIF
      DO n=1,sedMaxCl
        cum(n)              = sedCumRiv(n)
        cum(sedMaxCl+n)     = sedCumDep(n)
        cum(2*sedMaxCl+n)   = sedCumEro(n)
      ENDDO
      CALL WRITE_GLVEC_RL( 'pickup_sedimentcum.', suff, cum, 3*sedMaxCl,
     &                     myIter, myThid )
#endif /* ALLOW_SEDIMENT */

      RETURN
      END
