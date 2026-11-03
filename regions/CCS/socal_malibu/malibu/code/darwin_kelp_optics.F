#include "DARWIN_OPTIONS.h"

C malibu_v01 [K2]: read the kelp blade optics (data.kelpoptics, namelist
C DARWIN_KELP_OPTICS). Without the file, mode 0 (non-spectral) is used.

CBOP
C !ROUTINE: DARWIN_KELP_OPTICS_INIT
C !INTERFACE: ==========================================================
      SUBROUTINE DARWIN_KELP_OPTICS_INIT( myThid )
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#ifdef ALLOW_RADTRANS
#include "RADTRANS_SIZE.h"
#include "DARWIN_KELP_OPTICS.h"
#endif
      INTEGER myThid
CEOP
#if defined(ALLOW_RADTRANS) && defined(DARWIN_ALLOW_MACROALGAE)
      INTEGER iUnit, l
      LOGICAL exst
      _RL tot
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER ios
      NAMELIST /DARWIN_KELP_OPTICS/ kelp_optics_mode, kelp_proj,
     &   kelp_sfa, kelp_absBlade, kelp_refBlade, kelp_trnBlade
      NAMELIST /DARWIN_TSS_OPTICS/ tss_optics_mode, tss_aStar,
     &   tss_bStar, tss_bbRatio
      NAMELIST /DARWIN_BBP_OPTICS/ bbp_eta, bbp_lambda0
      NAMELIST /DARWIN_CDOM_RIVER/ cdom_riv_frac

      _BEGIN_MASTER(myThid)
      kelp_optics_mode = 0
      kelp_proj = 0.5 _d 0
      kelp_sfa = 1. _d 0 / 58. _d 0
      tss_optics_mode = 0
      tss_bbRatio = 0. _d 0
      bbp_eta = 0. _d 0
      bbp_lambda0 = 550. _d 0
      cdom_riv_frac = 0. _d 0
      DO l = 1, nlam
       tss_aStar(l) = 0. _d 0
       tss_bStar(l) = 0. _d 0
       kelp_absBlade(l) = 0. _d 0
       kelp_refBlade(l) = 0. _d 0
       kelp_trnBlade(l) = 0. _d 0
      ENDDO
      INQUIRE( FILE='data.kelpoptics', EXIST=exst )
      IF ( exst ) THEN
       CALL OPEN_COPY_DATA_FILE(
     I                   'data.kelpoptics', 'DARWIN_KELP_OPTICS_INIT',
     O                   iUnit,
     I                   myThid )
       READ(UNIT=iUnit, NML=DARWIN_KELP_OPTICS)
C      [K4] TSS optics group is optional
       READ(UNIT=iUnit, NML=DARWIN_TSS_OPTICS, IOSTAT=ios)
       IF ( ios.NE.0 ) tss_optics_mode = 0
C      [K6], [K7] optional groups (rewind: a missing group leaves EOF)
       REWIND(iUnit)
       READ(UNIT=iUnit, NML=DARWIN_BBP_OPTICS, IOSTAT=ios)
       IF ( ios.NE.0 ) bbp_eta = 0. _d 0
       REWIND(iUnit)
       READ(UNIT=iUnit, NML=DARWIN_CDOM_RIVER, IOSTAT=ios)
       IF ( ios.NE.0 ) cdom_riv_frac = 0. _d 0
#ifdef SINGLE_DISK_IO
       CLOSE(iUnit)
#else
       CLOSE(iUnit,STATUS='DELETE')
#endif /* SINGLE_DISK_IO */
      ENDIF
      IF ( kelp_optics_mode.EQ.1 ) THEN
       DO l = 1, nlam
        tot = kelp_absBlade(l) + kelp_refBlade(l) + kelp_trnBlade(l)
        IF ( ABS(tot - 1. _d 0).GT.1. _d -3 ) THEN
         WRITE(msgBuf,'(A,I3,A,F8.4)')
     &    'DARWIN_KELP_OPTICS_INIT: band', l, ' A+R+T =', tot
         CALL PRINT_ERROR( msgBuf, myThid )
         STOP 'ABNORMAL END: S/R DARWIN_KELP_OPTICS_INIT'
        ENDIF
       ENDDO
      ENDIF
      WRITE(msgBuf,'(A,I2,A,F6.3,A,F8.5)')
     & 'DARWIN_KELP_OPTICS_INIT: mode', kelp_optics_mode,
     & ' proj', kelp_proj, ' sfa [m2/gDW]', kelp_sfa
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,I2,A,F7.4)')
     & 'DARWIN_KELP_OPTICS_INIT: TSS optics mode', tss_optics_mode,
     & ' bb ratio', tss_bbRatio
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,F6.3,A,F6.1,A,1PE10.3)')
     & 'DARWIN_KELP_OPTICS_INIT: bbp eta', bbp_eta,
     & ' lambda0', bbp_lambda0, ' river CDOM/DOP', cdom_riv_frac
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      _END_MASTER(myThid)
      _BARRIER
#endif
      RETURN
      END
