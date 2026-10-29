#include "SEDIMENT_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: SEDIMENT_READPARMS

C !INTERFACE:
      SUBROUTINE SEDIMENT_READPARMS( myThid )

C !DESCRIPTION:
C     Read the sediment parameters from data.sediment (namelist SEDIMENT_PARM01)

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
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER iUnit, n, ios

      NAMELIST /SEDIMENT_PARM01/
     &       sedNcl, sedPtrIdx, sedWs, sedTauCE, sedEroM, sedTauCD,
     &       sedRiverC, sedRiverCFile, sedBulkRho, sedBedFile,
     &       sedDryDepth,
     &       sedDryDeposit,
     &       sedRefDepth, sedMonFreq, sedConsTime, sedTauCEf, sedEroMf
      NAMELIST /SEDIMENT_PARM02/
     &       sedNlay, sedCohesive, sedPorosity, sedFcLo, sedFcHi,
     &       sedTcbNew, sedTcbOffset, sedTcbSlope, sedTcbMin, sedTcbMax,
     &       sedTcons, sedTswell, sedZaCoef, sedZaMin, sedLayMass
      NAMELIST /SEDIMENT_PARM03/
     &       sedWaveMode, sedWaveUbFile, sedWavePerFile, sedWaveDirFile,
     &       sedWaveHsFile,
     &       sedWaveNrec, sedWaveDt, sedWindSpeed, sedWindDir, sedFetch,
     &       sedWaveGamma, sedZ0, sedWaveIce

      IF ( .NOT.useSEDIMENT ) THEN
        _BEGIN_MASTER(myThid)
         CALL PACKAGES_UNUSED_MSG( 'useSEDIMENT', ' ', ' ' )
        _END_MASTER(myThid)
        RETURN
      ENDIF

      _BEGIN_MASTER(myThid)

C--   Defaults: one cohesive mud class in pTracer 1
      sedNcl = 1
      DO n=1,sedMaxCl
        sedPtrIdx(n)  = n
        sedWs(n)      = 0.5 _d -3
        sedTauCE(n)   = 0.2 _d 0
        sedEroM(n)    = 5. _d -5
        sedTauCD(n)   = 0. _d 0
        sedRiverC(n)  = 0. _d 0
        sedBulkRho(n) = 500. _d 0
        sedBedFile(n) = ' '
        sedRiverCFile(n) = ' '
        sedTauCEf(n)  = UNSET_RL
        sedEroMf(n)   = UNSET_RL
      ENDDO
      sedConsTime = 0. _d 0
      sedDryDepth = 0.10 _d 0
      sedDryDeposit = .FALSE.
      sedRefDepth = 0.5 _d 0
      sedMonFreq  = monitorFreq
C-    layered bed (off by default)
      sedNlay     = 0
      DO n=1,sedMaxCl
        sedCohesive(n) = .FALSE.
      ENDDO
      sedPorosity  = 0. _d 0
      sedFcLo      = 0.03 _d 0
      sedFcHi      = 0.2  _d 0
      sedTcbNew    = 0.05 _d 0
      sedTcbOffset = 0. _d 0
      sedTcbSlope  = 2. _d 0
      sedTcbMin    = 0.05 _d 0
      sedTcbMax    = 5. _d 0
      sedTcons     = 86400. _d 0
      sedTswell    = 864000. _d 0
      sedZaCoef    = 0.007 _d 0
      sedZaMin     = 0.5 _d -3
      sedLayMass   = 1. _d 0
C-    waves (off by default)
      sedWaveMode    = 0
      sedWaveUbFile  = ' '
      sedWavePerFile = ' '
      sedWaveDirFile = ' '
      sedWaveHsFile  = ' '
      sedWaveNrec    = 1
      sedWaveDt      = 0. _d 0
      sedWindSpeed   = 0. _d 0
      sedWindDir     = 0. _d 0
      sedFetch       = 1. _d 4
      sedWaveGamma   = 0.5 _d 0
      sedZ0          = 2. _d -4
      sedWaveIce     = .TRUE.

      WRITE(msgBuf,'(A)') 'SEDIMENT_READPARMS: opening data.sediment'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      CALL OPEN_COPY_DATA_FILE(
     I                     'data.sediment', 'SEDIMENT_READPARMS',
     O                     iUnit,
     I                     myThid )
      READ(UNIT=iUnit,NML=SEDIMENT_PARM01)
C-    optional namelists (absent in older data.sediment files: end of
C     file, ios < 0)
      REWIND(iUnit)
      READ(UNIT=iUnit,NML=SEDIMENT_PARM02,IOSTAT=ios)
      IF ( ios.GT.0 ) THEN
        WRITE(msgBuf,'(A)')
     &   'SEDIMENT_READPARMS: error reading SEDIMENT_PARM02'
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R SEDIMENT_READPARMS'
      ENDIF
      REWIND(iUnit)
      READ(UNIT=iUnit,NML=SEDIMENT_PARM03,IOSTAT=ios)
      IF ( ios.GT.0 ) THEN
        WRITE(msgBuf,'(A)')
     &   'SEDIMENT_READPARMS: error reading SEDIMENT_PARM03'
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R SEDIMENT_READPARMS'
      ENDIF
C-    fluff layer: unset fluff parameters take the bed values
      DO n=1,sedMaxCl
        IF ( sedTauCEf(n).EQ.UNSET_RL ) sedTauCEf(n) = sedTauCE(n)
        IF ( sedEroMf(n) .EQ.UNSET_RL ) sedEroMf(n)  = sedEroM(n)
      ENDDO
      WRITE(msgBuf,'(A)')
     &    'SEDIMENT_READPARMS: finished reading data.sediment'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
#ifdef SINGLE_DISK_IO
      CLOSE(iUnit)
#else
      CLOSE(iUnit,STATUS='DELETE')
#endif /* SINGLE_DISK_IO */

      _END_MASTER(myThid)
      _BARRIER

#endif /* ALLOW_SEDIMENT */

      RETURN
      END
