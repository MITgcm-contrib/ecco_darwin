#include "SEDIMENT_OPTIONS.h"

CBOP
C !ROUTINE: SEDIMENT_INIT_VARIA
C !INTERFACE:
      SUBROUTINE SEDIMENT_INIT_VARIA( myThid )
C !DESCRIPTION:
C     Initial bed (sedBedFile or clean), zero budgets; on a restart
C     read them from the pickup
C !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SEDIMENT.h"
C !INPUT PARAMETERS:
      INTEGER myThid
CEOP

#ifdef ALLOW_SEDIMENT
C !LOCAL VARIABLES:
      INTEGER i, j, n, l, bi, bj
      _RL tmp2d(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      DO bj=myByLo(myThid), myByHi(myThid)
       DO bi=myBxLo(myThid), myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          sedTauB(i,j,bi,bj) = 0. _d 0
          sedTauW(i,j,bi,bj) = 0. _d 0
          sedZa  (i,j,bi,bj) = 0. _d 0
          sedUbw (i,j,bi,bj) = 0. _d 0
          sedPerW(i,j,bi,bj) = 0. _d 0
          sedDirW(i,j,bi,bj) = 0. _d 0
          sedHsW (i,j,bi,bj) = 0. _d 0
          sedWub0(i,j,bi,bj) = 0. _d 0
          sedWub1(i,j,bi,bj) = 0. _d 0
          sedWpe0(i,j,bi,bj) = 0. _d 0
          sedWpe1(i,j,bi,bj) = 0. _d 0
          sedWdi0(i,j,bi,bj) = 0. _d 0
          sedWdi1(i,j,bi,bj) = 0. _d 0
         ENDDO
        ENDDO
        DO n=1,sedMaxCl
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           sedBed(i,j,n,bi,bj) = 0. _d 0
           sedFluff(i,j,n,bi,bj) = 0. _d 0
          ENDDO
         ENDDO
         DO l=1,sedMaxLay
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            sedLay(i,j,l,n,bi,bj) = 0. _d 0
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        DO l=1,sedMaxLay
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           sedTcb(i,j,l,bi,bj) = sedTcbNew
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      _BEGIN_MASTER(myThid)
      DO n=1,sedMaxCl
        sedCumRiv(n) = 0. _d 0
        sedCumDep(n) = 0. _d 0
        sedCumEro(n) = 0. _d 0
      ENDDO
      sedErrMax = 0. _d 0
      sedWrec   = 0
      _END_MASTER(myThid)

      IF ( startTime.EQ.baseTime .AND. nIter0.EQ.0
     &     .AND. pickupSuff.EQ.' ' ) THEN
C-    initial bed from files (kg/m^2), masked to the ocean
       DO n=1,sedNcl
        IF ( sedBedFile(n).NE.' ' ) THEN
         CALL READ_FLD_XY_RL( sedBedFile(n), ' ', tmp2d, 0, myThid )
         DO bj=myByLo(myThid), myByHi(myThid)
          DO bi=myBxLo(myThid), myBxHi(myThid)
           DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
             sedBed(i,j,n,bi,bj) = MAX( zeroRL, tmp2d(i,j,bi,bj) )
     &                            *maskInC(i,j,bi,bj)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDIF
       ENDDO
       IF ( sedNlay.GT.0 ) CALL SEDIMENT_BED_FROM_TOTAL( myThid )
      ELSE
        CALL SEDIMENT_READ_PICKUP( nIter0, myThid )
      ENDIF
      _BARRIER
#endif /* ALLOW_SEDIMENT */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: SEDIMENT_BED_FROM_TOTAL
C !INTERFACE:
      SUBROUTINE SEDIMENT_BED_FROM_TOTAL( myThid )
C !DESCRIPTION:
C     Layered bed from the bed totals sedBed (an initial bed file, or a
C     pickup written without layers): all of it in the deepest layer,
C     at the equilibrium bulk critical stress of its mid mass depth; the
C     active layer is drawn from it in the first step
C !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "SEDIMENT.h"
C !INPUT PARAMETERS:
      INTEGER myThid
CEOP

#ifdef ALLOW_SEDIMENT
C !LOCAL VARIABLES:
      INTEGER i, j, n, l, bi, bj
      _RL sm, teq

      DO bj=myByLo(myThid), myByHi(myThid)
       DO bi=myBxLo(myThid), myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          sm = 0. _d 0
          DO n=1,sedNcl
           DO l=1,sedNlay
            sedLay(i,j,l,n,bi,bj) = 0. _d 0
           ENDDO
           sedLay(i,j,sedNlay,n,bi,bj) = sedBed(i,j,n,bi,bj)
           sm = sm + sedBed(i,j,n,bi,bj)
          ENDDO
          IF ( sm.GT.zeroRL ) THEN
           teq = EXP( ( LOG(halfRL*sm) - sedTcbOffset )/sedTcbSlope )
           sedTcb(i,j,sedNlay,bi,bj) = MIN( sedTcbMax,
     &                                     MAX( sedTcbMin, teq ) )
          ENDIF
         ENDDO
        ENDDO
       ENDDO
      ENDDO
#endif /* ALLOW_SEDIMENT */

      RETURN
      END
