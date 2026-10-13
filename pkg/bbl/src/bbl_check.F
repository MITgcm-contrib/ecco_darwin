#include "BBL_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: BBL_CHECK

C !INTERFACE:
      SUBROUTINE BBL_CHECK( myThid )

C     !DESCRIPTION:
C     Check dependances with other packages

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "BBL.h"

C     !INPUT PARAMETERS:
C     myThid   :: my Thread Id number
      INTEGER myThid
CEOP

#ifdef ALLOW_BBL
C     !LOCAL VARIABLES:
C     msgBuf   :: Informational/error message buffer
      CHARACTER*(MAX_LEN_MBUF) msgBuf

      _BEGIN_MASTER(myThid)

       WRITE(msgBuf,'(A)') 'BBL_CHECK: #define BBL'
       CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT , 1)

C--  Print out some key parameters :
       CALL WRITE_0D_RL( bbl_wvel, INDEX_NONE,
     &  'bbl_wvel =', ' /* vertical entrainment (m/s) */')
       CALL WRITE_0D_RL( bbl_hvel, INDEX_NONE,
     &  'bbl_hvel =', ' /* horizontal velocity (m/s) */')
       CALL WRITE_0D_RL( bbl_initEta, INDEX_NONE,
     &  'bbl_initEta =', ' /* default thickness (m) */')
       CALL WRITE_0D_L ( bbl_useNofSpeed, INDEX_NONE,
     &  'bbl_useNofSpeed =', ' /* speed = nofCoeff*gp*slope/f */')
       CALL WRITE_0D_RL( bbl_nofCoeff, INDEX_NONE,
     &  'bbl_nofCoeff =', ' /* downslope fraction of Nof speed */')
       CALL WRITE_0D_RL( bbl_fMin, INDEX_NONE,
     &  'bbl_fMin =', ' /* lower bound of |f| (1/s) */')
       CALL WRITE_0D_I ( bbl_selectEntrain, INDEX_NONE,
     &  'bbl_selectEntrain =', ' /* 0: none, 1: Turner 1986 */')
       CALL WRITE_0D_I ( bbl_selectDetrain, INDEX_NONE,
     &  'bbl_selectDetrain =', ' /* 0: bbl_wvel, 1: plume-like */')
       CALL WRITE_0D_RL( bbl_uDet, INDEX_NONE,
     &  'bbl_uDet =', ' /* no detrainment above this outflow (m/s) */')
       CALL WRITE_0D_RL( bbl_gpTaper, INDEX_NONE,
     &  'bbl_gpTaper =', ' /* g-prime taper scale (m/s^2) */')
       CALL WRITE_0D_RL( bbl_tauMix, INDEX_NONE,
     &  'bbl_tauMix =', ' /* dissolution timescale (s) */')
       CALL WRITE_0D_L ( bbl_coldStart, INDEX_NONE,
     &  'bbl_coldStart =', ' /* restart without pickup_bbl */')
       CALL WRITE_0D_L ( bbl_doPtracers, INDEX_NONE,
     &  'bbl_doPtracers =', ' /* BBL carries pkg/ptracers */')

C--  Check parameters:
C    flux-form exchange in BBL_CALC_RHS needs 2 halo points
      IF ( OLx.LT.2 .OR. OLy.LT.2 ) THEN
        WRITE(msgBuf,'(A)') 'BBL_CHECK: pkg/bbl needs OLx,OLy >= 2'
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R BBL_CHECK'
      ENDIF
      IF ( bbl_useNofSpeed .AND.
     &     ( bbl_nofCoeff.LE.0. .OR. bbl_fMin.LE.0. ) ) THEN
        WRITE(msgBuf,'(2A)') 'BBL_CHECK: bbl_useNofSpeed needs ',
     &       'bbl_nofCoeff > 0 and bbl_fMin > 0'
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R BBL_CHECK'
      ENDIF

      IF ( bbl_selectDetrain.LT.0 .OR. bbl_selectDetrain.GT.1 .OR.
     &     ( bbl_selectDetrain.EQ.1 .AND. bbl_uDet.LE.0. ) ) THEN
        WRITE(msgBuf,'(2A)') 'BBL_CHECK: bbl_selectDetrain must be ',
     &       '0 or 1, with bbl_uDet > 0'
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R BBL_CHECK'
      ENDIF
      IF ( bbl_gpTaper.LT.0. .OR. bbl_tauMix.LT.0. ) THEN
        WRITE(msgBuf,'(2A)') 'BBL_CHECK: bbl_gpTaper and ',
     &       'bbl_tauMix must be >= 0'
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R BBL_CHECK'
      ENDIF

C    entrainment scales with the BBL speed (w_E = E(Ri)*U); with a
C    prescribed speed of O(1 m/s) it thickens the BBL without bound
C    (in cs32x15 tests: bbl_eta ~150 m, summed transport ~500 Sv)
      IF ( bbl_selectEntrain.GE.1 .AND. .NOT.bbl_useNofSpeed ) THEN
        WRITE(msgBuf,'(2A)') '** WARNING ** BBL_CHECK: entrainment ',
     &       'with a constant BBL speed (bbl_hvel); consider'
        CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(2A)') '** WARNING ** BBL_CHECK: ',
     &       'bbl_useNofSpeed=.TRUE. (see pkg/bbl/README.md)'
        CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
      ENDIF
      IF ( bbl_doPtracers .AND. .NOT.usePTRACERS ) THEN
        WRITE(msgBuf,'(2A)') 'BBL_CHECK: bbl_doPtracers=.TRUE. ',
     &       'requires usePTRACERS=.TRUE.'
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R BBL_CHECK'
      ENDIF
#ifndef ALLOW_PTRACERS
      IF ( bbl_doPtracers ) THEN
        WRITE(msgBuf,'(2A)') 'BBL_CHECK: bbl_doPtracers=.TRUE. ',
     &       'requires pkg/ptracers to be compiled'
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R BBL_CHECK'
      ENDIF
#endif
      IF ( bbl_selectEntrain.LT.0 .OR. bbl_selectEntrain.GT.1 ) THEN
        WRITE(msgBuf,'(A,I4)')
     &       'BBL_CHECK: invalid bbl_selectEntrain =', bbl_selectEntrain
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R BBL_CHECK'
      ENDIF

      _END_MASTER(myThid)

#endif /* ALLOW_BBL */

      RETURN
      END
