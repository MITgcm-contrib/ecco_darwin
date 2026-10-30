#include "CPP_OPTIONS.h"

CBOP
C     !ROUTINE: LR17_RIVERS
C     !INTERFACE:
      SUBROUTINE LR17_RIVERS( myTime, myIter, myThid )
C     !DESCRIPTION:
C     Time-varying river sources for LR17: adds each river's discharge
C     (data.rivers, read on the first call) to addMass, shared equally
C     among the wet cells of a (2*rivHalf+1)^2 block centered on its
C     source cell (rivHalf = 0, the default: the source cell alone).
C     Spreading lowers the divergence of the single source cell, where
C     salinity undershot below zero (2026-10-24). Called from
C     LOAD_FIELDS_DRIVER right after addMass is reset, so other
C     contributions still add on top.
C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "FFIELDS.h"
#include "LR17_RIVERS.h"
C     !INPUT PARAMETERS:
      _RL     myTime
      INTEGER myIter, myThid
CEOP
#ifdef ALLOW_ADDFLUID
      INTEGER iUnit, r, n, i, j, k, bi, bj, iG, jG
      _RL     tRec, w, q, cnt
      LOGICAL done, counted
      SAVE done, counted
      DATA done /.FALSE./, counted /.FALSE./
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      NAMELIST /RIVERS_PARM/ rivN, rivNrec, rivT0, rivDt, rivI, rivJ,
     &                       rivQ, rivHalf

      IF ( .NOT.done ) THEN
        _BARRIER
        _BEGIN_MASTER(myThid)
        rivN    = 0
        rivNrec = 0
        rivT0   = 0. _d 0
        rivDt   = 86400. _d 0
        rivHalf = 0
        CALL OPEN_COPY_DATA_FILE( 'data.rivers', 'LR17_RIVERS',
     O                            iUnit, myThid )
        READ(UNIT=iUnit,NML=RIVERS_PARM)
        CLOSE(iUnit,STATUS='DELETE')
        IF ( rivN.GT.rivMax .OR. rivNrec.GT.rivRecMax
     &       .OR. rivNrec.LT.1 ) THEN
          WRITE(msgBuf,'(A,2I6)')
     &      'LR17_RIVERS: bad rivN, rivNrec =', rivN, rivNrec
          CALL PRINT_ERROR( msgBuf, myThid )
          STOP 'ABNORMAL END: S/R LR17_RIVERS'
        ENDIF
        WRITE(msgBuf,'(A,I3,A,I5,A)') 'LR17_RIVERS: ', rivN,
     &    ' rivers, ', rivNrec, ' records'
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
        _END_MASTER(myThid)
        _BARRIER
        done = .TRUE.
      ENDIF
      IF ( selectAddFluid.EQ.0 .OR. rivN.EQ.0 ) RETURN

C--   wet cells (top wet level) in each river's block, over all tiles
      IF ( .NOT.counted ) THEN
       DO r=1,rivN
        cnt = 0. _d 0
        DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)
          DO j=1,sNy
           jG = myYGlobalLo - 1 + (bj-1)*sNy + j
           DO i=1,sNx
            iG = myXGlobalLo - 1 + (bi-1)*sNx + i
            k = kSurfC(i,j,bi,bj)
            IF ( ABS(iG-rivI(r)).LE.rivHalf
     &           .AND. ABS(jG-rivJ(r)).LE.rivHalf
     &           .AND. k.LE.Nr ) cnt = cnt + maskC(i,j,k,bi,bj)
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        _GLOBAL_SUM_RL( cnt, myThid )
        IF ( cnt.LT.0.5 _d 0 ) THEN
          WRITE(msgBuf,'(A,I3,A)') 'LR17_RIVERS: river ', r,
     &      ' has no wet cell in its block'
          CALL PRINT_ERROR( msgBuf, myThid )
          STOP 'ABNORMAL END: S/R LR17_RIVERS'
        ENDIF
        _BEGIN_MASTER(myThid)
        rivCells(r) = cnt
        WRITE(msgBuf,'(A,I3,A,I3,A,I3)') 'LR17_RIVERS: river ', r,
     &    ' shared by', NINT(cnt), ' wet cells, rivHalf =', rivHalf
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
        _END_MASTER(myThid)
       ENDDO
       _BARRIER
       counted = .TRUE.
      ENDIF

C--   record weights (linear between records, held at the ends)
      tRec = ( myTime - rivT0 ) / rivDt
      n = INT( tRec ) + 1
      w = tRec - (n-1)
      IF ( tRec.LT.0. _d 0 ) THEN
        n = 1
        w = 0. _d 0
      ELSEIF ( n.GE.rivNrec ) THEN
        n = rivNrec - 1
        w = 1. _d 0
      ENDIF
      IF ( rivNrec.EQ.1 ) THEN
        n = 1
        w = 0. _d 0
      ENDIF

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1,sNy
         jG = myYGlobalLo - 1 + (bj-1)*sNy + j
         DO i=1,sNx
          iG = myXGlobalLo - 1 + (bi-1)*sNx + i
          DO r=1,rivN
           IF ( ABS(iG-rivI(r)).LE.rivHalf
     &          .AND. ABS(jG-rivJ(r)).LE.rivHalf ) THEN
            q = rivQ(r,n)
            IF ( rivNrec.GT.1 ) q = (1.-w)*rivQ(r,n) + w*rivQ(r,n+1)
            k = kSurfC(i,j,bi,bj)
            IF ( k.LE.Nr ) addMass(i,j,k,bi,bj) =
     &           addMass(i,j,k,bi,bj) + rhoConst*q*maskC(i,j,k,bi,bj)
     &                                  /rivCells(r)
           ENDIF
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
#endif /* ALLOW_ADDFLUID */
      RETURN
      END
