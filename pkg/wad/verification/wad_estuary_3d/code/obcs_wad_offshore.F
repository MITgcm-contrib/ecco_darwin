#include "OBCS_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: OBCS_WAD_OFFSHORE

C !INTERFACE:
      SUBROUTINE OBCS_WAD_OFFSHORE( futureIter, myThid )

C !DESCRIPTION:
C     Offshore water mass at the western open boundary of
C     the wad_* experiments (called at the end of OBCS_CALC, after the
C     Flather tide has set OBWu and the zero-gradient OBWt/s/ptr).
C     - T, S and pTracers are upwinded: where the boundary flow enters
C       the domain (OBWu > 0) they take the offshore value, elsewhere
C       they keep the zero-gradient (interior) value. The offshore
C       profile of each boundary row is its state at the first call
C       (initial conditions on a cold start; the restart state after a
C       pickup, since OBCS boundary values are not in the pickup).
C     Parameters in the optional file data.wadob, namelist WADOB_PARM01:
C       obUpwindTS  :: .TRUE. (default): upwind T/S/pTracers as above;
C                      .FALSE.: zero-gradient (the previous behaviour)
C       obTempOff, obSaltOff :: uniform offshore temperature / salinity,
C                      replacing the captured profile (default: unset)

C !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "OBCS_PARAMS.h"
#include "OBCS_GRID.h"
#include "OBCS_FIELDS.h"
#ifdef ALLOW_PTRACERS
# include "PTRACERS_SIZE.h"
# include "PTRACERS_PARAMS.h"
# include "OBCS_PTRACERS.h"
#endif

C !INPUT PARAMETERS:
      INTEGER futureIter
      INTEGER myThid
CEOP

#if (defined ALLOW_OBCS) && (defined ALLOW_OBCS_WEST)
C !LOCAL VARIABLES:
      INTEGER bi, bj, j, k, iUnit
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      LOGICAL exst
#ifdef ALLOW_PTRACERS
      INTEGER iTr
#endif
C     offshore state (saved between calls)
      LOGICAL obInit, obUpwindTS
      _RL obTempOff, obSaltOff
      _RL obOffT(1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL obOffS(1-OLy:sNy+OLy,Nr,nSx,nSy)
      SAVE obInit
      COMMON /OBCS_WAD_OFFSHORE_L/ obUpwindTS
      COMMON /OBCS_WAD_OFFSHORE_R/ obTempOff, obSaltOff, obOffT, obOffS
#ifdef ALLOW_PTRACERS
      _RL obOffP(1-OLy:sNy+OLy,Nr,nSx,nSy,PTRACERS_num)
      COMMON /OBCS_WAD_OFFSHORE_P/ obOffP
#endif
      DATA obInit /.FALSE./
      NAMELIST /WADOB_PARM01/ obUpwindTS, obTempOff, obSaltOff

C--   First call: parameters, and the offshore profile of each row
      IF ( .NOT.obInit ) THEN
       _BARRIER
       _BEGIN_MASTER(myThid)
       obUpwindTS = .TRUE.
       obTempOff  = UNSET_RL
       obSaltOff  = UNSET_RL
       INQUIRE( FILE='data.wadob', EXIST=exst )
       IF ( exst ) THEN
        CALL OPEN_COPY_DATA_FILE( 'data.wadob', 'OBCS_WAD_OFFSHORE',
     O                            iUnit, myThid )
        READ(UNIT=iUnit, NML=WADOB_PARM01)
#ifdef SINGLE_DISK_IO
        CLOSE(iUnit)
#else
        CLOSE(iUnit,STATUS='DELETE')
#endif
       ENDIF
       WRITE(msgBuf,'(A,L2)')
     &  'OBCS_WAD_OFFSHORE: obUpwindTS=', obUpwindTS
       CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
       _END_MASTER(myThid)
       _BARRIER
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         DO k=1,Nr
          DO j=1-OLy,sNy+OLy
           obOffT(j,k,bi,bj) = OBWt(j,k,bi,bj)
           obOffS(j,k,bi,bj) = OBWs(j,k,bi,bj)
           IF ( obTempOff.NE.UNSET_RL ) obOffT(j,k,bi,bj) = obTempOff
           IF ( obSaltOff.NE.UNSET_RL ) obOffS(j,k,bi,bj) = obSaltOff
#ifdef ALLOW_PTRACERS
           DO iTr=1,PTRACERS_numInUse
            obOffP(j,k,bi,bj,iTr) = OBWptr(j,k,bi,bj,iTr)
           ENDDO
#endif
          ENDDO
         ENDDO
        ENDDO
       ENDDO
       _BARRIER
       _BEGIN_MASTER(myThid)
       obInit = .TRUE.
       _END_MASTER(myThid)
       _BARRIER
      ENDIF

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         IF ( OB_Iw(j,bi,bj).NE.OB_indexNone ) THEN
C--   T, S, pTracers: offshore water where the flow enters
          IF ( obUpwindTS .AND. OBWu(j,1,bi,bj).GT.zeroRL ) THEN
           DO k=1,Nr
            OBWt(j,k,bi,bj) = obOffT(j,k,bi,bj)
            OBWs(j,k,bi,bj) = obOffS(j,k,bi,bj)
#ifdef ALLOW_PTRACERS
            DO iTr=1,PTRACERS_numInUse
             OBWptr(j,k,bi,bj,iTr) = obOffP(j,k,bi,bj,iTr)
            ENDDO
#endif
           ENDDO
          ENDIF
         ENDIF
        ENDDO
       ENDDO
      ENDDO
#endif /* ALLOW_OBCS and ALLOW_OBCS_WEST */

      RETURN
      END
