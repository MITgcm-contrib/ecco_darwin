#include "SEDIMENT_OPTIONS.h"
#ifdef ALLOW_PTRACERS
# include "PTRACERS_OPTIONS.h"
#endif

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: SEDIMENT_CHECK

C !INTERFACE:
      SUBROUTINE SEDIMENT_CHECK( myThid )

C !DESCRIPTION:
C     Check sediment parameters and print them

C !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "SEDIMENT.h"
#ifdef ALLOW_PTRACERS
# include "PTRACERS_SIZE.h"
# include "PTRACERS_PARAMS.h"
#endif

C !INPUT PARAMETERS:
      INTEGER myThid
CEOP

#ifdef ALLOW_SEDIMENT
C !LOCAL VARIABLES:
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER n, m, errCount

      _BEGIN_MASTER(myThid)
      errCount = 0

      WRITE(msgBuf,'(A)') 'SEDIMENT_CHECK: #define ALLOW_SEDIMENT'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )

#ifndef ALLOW_PTRACERS
      WRITE(msgBuf,'(A)') 'SEDIMENT_CHECK: needs pkg/ptracers'
      CALL PRINT_ERROR( msgBuf, myThid )
      errCount = errCount + 1
#else
      IF ( .NOT.usePTRACERS ) THEN
        WRITE(msgBuf,'(A)') 'SEDIMENT_CHECK: needs usePTRACERS=.TRUE.'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF
      IF ( sedNcl.LT.1 .OR. sedNcl.GT.sedMaxCl ) THEN
        WRITE(msgBuf,'(A,I4,A,I4)') 'SEDIMENT_CHECK: sedNcl =', sedNcl,
     &   ' must be in 1 ..', sedMaxCl
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ELSE
       DO n=1,sedNcl
        IF ( sedPtrIdx(n).LT.1 .OR.
     &       sedPtrIdx(n).GT.PTRACERS_numInUse ) THEN
          WRITE(msgBuf,'(A,I2,A,I4,A)') 'SEDIMENT_CHECK: class', n,
     &     ': sedPtrIdx =', sedPtrIdx(n), ' is not a pTracer in use'
          CALL PRINT_ERROR( msgBuf, myThid )
          errCount = errCount + 1
        ENDIF
        IF ( sedWs(n).LT.zeroRL .OR. sedTauCE(n).LE.zeroRL .OR.
     &       sedEroM(n).LT.zeroRL .OR. sedRiverC(n).LT.zeroRL .OR.
     &       sedBulkRho(n).LE.zeroRL ) THEN
          WRITE(msgBuf,'(A,I2,A)') 'SEDIMENT_CHECK: class', n,
     &     ': need sedWs, sedEroM, sedRiverC >= 0 and'
     &     //' sedTauCE, sedBulkRho > 0'
          CALL PRINT_ERROR( msgBuf, myThid )
          errCount = errCount + 1
        ENDIF
       ENDDO
       DO n=2,sedNcl
        DO m=1,n-1
         IF ( sedPtrIdx(n).EQ.sedPtrIdx(m) ) THEN
          WRITE(msgBuf,'(A,I2,A,I2,A)') 'SEDIMENT_CHECK: classes', m,
     &     ' and', n, ' use the same pTracer'
          CALL PRINT_ERROR( msgBuf, myThid )
          errCount = errCount + 1
         ENDIF
        ENDDO
       ENDDO
      ENDIF
#endif /* ALLOW_PTRACERS */
C-    layered bed
      IF ( sedNlay.LT.0 .OR. sedNlay.GT.sedMaxLay ) THEN
        WRITE(msgBuf,'(A,I4,A,I4)') 'SEDIMENT_CHECK: sedNlay =',
     &   sedNlay, ' must be in 0 ..', sedMaxLay
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF
      IF ( sedNlay.GT.0 ) THEN
       IF ( sedConsTime.GT.zeroRL ) THEN
        WRITE(msgBuf,'(2A)') 'SEDIMENT_CHECK: the fluff layer',
     &   ' (sedConsTime > 0) and the layered bed (sedNlay > 0) exclude'
        CALL PRINT_ERROR( msgBuf, myThid )
        WRITE(msgBuf,'(2A)') 'SEDIMENT_CHECK: each other: fresh',
     &   ' deposit takes sedTcbNew and consolidates with sedTcons'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
       ENDIF
       IF ( sedPorosity.LT.zeroRL .OR. sedPorosity.GE.oneRL .OR.
     &      sedFcHi.LE.sedFcLo .OR. sedTcbSlope.EQ.zeroRL .OR.
     &      sedTcons.LE.zeroRL .OR. sedTswell.LE.zeroRL .OR.
     &      sedZaMin.LE.zeroRL .OR. sedLayMass.LE.zeroRL .OR.
     &      sedTcbMin.GT.sedTcbMax ) THEN
        WRITE(msgBuf,'(2A)') 'SEDIMENT_CHECK: need 0 <= sedPorosity',
     &   ' < 1, sedFcHi > sedFcLo, sedTcbSlope /= 0, sedTcbMin <='
        CALL PRINT_ERROR( msgBuf, myThid )
        WRITE(msgBuf,'(2A)') 'SEDIMENT_CHECK: sedTcbMax and sedTcons,',
     &   ' sedTswell, sedZaMin, sedLayMass > 0'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
       ENDIF
      ENDIF
C-    waves
      IF ( sedWaveMode.LT.0 .OR. sedWaveMode.GT.3 ) THEN
        WRITE(msgBuf,'(A,I4,A)') 'SEDIMENT_CHECK: sedWaveMode =',
     &   sedWaveMode, ' must be 0, 1, 2 or 3'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF
      IF ( sedWaveMode.EQ.1 .AND. ( sedWaveUbFile.EQ.' ' .OR.
     &     sedWavePerFile.EQ.' ' .OR. sedWaveNrec.LT.1 .OR.
     &     ( sedWaveNrec.GT.1 .AND. sedWaveDt.LE.zeroRL ) ) ) THEN
        WRITE(msgBuf,'(2A)') 'SEDIMENT_CHECK: sedWaveMode=1 needs',
     &   ' sedWaveUbFile, sedWavePerFile, sedWaveNrec >= 1'
        CALL PRINT_ERROR( msgBuf, myThid )
        WRITE(msgBuf,'(A)')
     &   'SEDIMENT_CHECK: and sedWaveDt > 0 when sedWaveNrec > 1'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF
      IF ( sedWaveMode.EQ.3 .AND. ( sedWaveHsFile.EQ.' ' .OR.
     &     sedWavePerFile.EQ.' ' .OR. sedWaveNrec.LT.1 .OR.
     &     ( sedWaveNrec.GT.1 .AND. sedWaveDt.LE.zeroRL ) ) ) THEN
        WRITE(msgBuf,'(2A)') 'SEDIMENT_CHECK: sedWaveMode=3 needs',
     &   ' sedWaveHsFile, sedWavePerFile, sedWaveNrec >= 1'
        CALL PRINT_ERROR( msgBuf, myThid )
        WRITE(msgBuf,'(A)')
     &   'SEDIMENT_CHECK: and sedWaveDt > 0 when sedWaveNrec > 1'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF
      IF ( sedWaveMode.EQ.2 .AND. ( sedFetch.LE.zeroRL .OR.
     &     sedWindSpeed.LT.zeroRL ) ) THEN
        WRITE(msgBuf,'(A)')
     &   'SEDIMENT_CHECK: sedWaveMode=2 needs sedFetch > 0, U10 >= 0'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF
      IF ( sedWaveMode.GT.0 .AND. sedZ0.LE.zeroRL ) THEN
        WRITE(msgBuf,'(A)') 'SEDIMENT_CHECK: need sedZ0 > 0'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF

C-    surface-level NLFS: the surface cell refers to hFac_surfC, which
C     SEDIMENT_STEP does not use (the bed/suspended budget would drift)
      IF ( nonlinFreeSurf.GT.0 .AND. select_rStar.EQ.0 ) THEN
        WRITE(msgBuf,'(A)')
     &   'SEDIMENT_CHECK: select_rStar=0 is not supported (use r*)'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF

      IF ( errCount.GT.0 ) THEN
        WRITE(msgBuf,'(A,I3,A)')
     &   'SEDIMENT_CHECK: detected', errCount, ' fatal error(s)'
        CALL PRINT_ERROR( msgBuf, myThid )
        CALL ALL_PROC_DIE( 0 )
        STOP 'ABNORMAL END: S/R SEDIMENT_CHECK'
      ENDIF

C--   Print the parameters
      CALL WRITE_0D_I( sedNcl, INDEX_NONE, 'sedNcl =',
     &  ' /* number of sediment classes */')
      CALL WRITE_1D_I( sedPtrIdx, sedNcl, INDEX_K, 'sedPtrIdx =',
     &  ' /* pTracer index of each class */')
      CALL WRITE_1D_RL( sedWs, sedNcl, INDEX_K, 'sedWs =',
     &  ' /* settling velocity (m/s) */')
      CALL WRITE_1D_RL( sedTauCE, sedNcl, INDEX_K, 'sedTauCE =',
     &  ' /* critical stress for erosion (Pa) */')
      CALL WRITE_1D_RL( sedEroM, sedNcl, INDEX_K, 'sedEroM =',
     &  ' /* erosion rate constant (kg/m2/s) */')
      CALL WRITE_1D_RL( sedTauCD, sedNcl, INDEX_K, 'sedTauCD =',
     &  ' /* critical stress for deposition (Pa, <=0: none) */')
      CALL WRITE_1D_RL( sedRiverC, sedNcl, INDEX_K, 'sedRiverC =',
     &  ' /* concentration of added fluid (kg/m3) */')
      DO n=1,sedNcl
       IF ( sedRiverCFile(n).NE.' ' ) THEN
        CALL WRITE_0D_C( sedRiverCFile(n), -1, INDEX_NONE,
     &   'sedRiverCFile =', ' /* map of the added-fluid conc. */')
       ENDIF
      ENDDO
      CALL WRITE_1D_RL( sedBulkRho, sedNcl, INDEX_K, 'sedBulkRho =',
     &  ' /* dry bulk density of deposit (kg/m3) */')
      CALL WRITE_0D_L( sedDryDeposit, INDEX_NONE, 'sedDryDeposit =',
     &  ' /* deposit from drying columns (films) too */')
      CALL WRITE_0D_RL( sedDryDepth, INDEX_NONE, 'sedDryDepth =',
     &  ' /* no erosion from shallower columns (m) */')
      CALL WRITE_0D_RL( sedConsTime, INDEX_NONE, 'sedConsTime =',
     &  ' /* fluff consolidation time (s); <= 0: no fluff layer */')
      CALL WRITE_1D_RL( sedTauCEf, sedNcl, INDEX_K, 'sedTauCEf =',
     &  ' /* critical stress for erosion of the fluff (Pa) */')
      CALL WRITE_1D_RL( sedEroMf, sedNcl, INDEX_K, 'sedEroMf =',
     &  ' /* erosion rate constant of the fluff (kg/m2/s) */')
      CALL WRITE_0D_RL( sedRefDepth, INDEX_NONE, 'sedRefDepth =',
     &  ' /* bed stress from the lowest sedRefDepth (m) */')
      CALL WRITE_0D_RL( sedMonFreq, INDEX_NONE, 'sedMonFreq =',
     &  ' /* %SEDIMENT_MON frequency (s) */')
      CALL WRITE_0D_I( sedNlay, INDEX_NONE, 'sedNlay =',
     &  ' /* bed layers (0: single store) */')
      IF ( sedNlay.GT.0 ) THEN
      CALL WRITE_1D_L( sedCohesive, sedNcl, INDEX_K, 'sedCohesive =',
     &  ' /* cohesive classes */')
      CALL WRITE_0D_RL( sedPorosity, INDEX_NONE, 'sedPorosity =',
     &  ' /* bed porosity */')
      CALL WRITE_0D_RL( sedFcLo, INDEX_NONE, 'sedFcLo =',
     &  ' /* cohesive fraction where Pc starts > 0 */')
      CALL WRITE_0D_RL( sedFcHi, INDEX_NONE, 'sedFcHi =',
     &  ' /* cohesive fraction where Pc reaches 1 */')
      CALL WRITE_0D_RL( sedTcbNew, INDEX_NONE, 'sedTcbNew =',
     &  ' /* bulk critical stress of fresh deposit (Pa) */')
      CALL WRITE_0D_RL( sedTcbOffset, INDEX_NONE, 'sedTcbOffset =',
     &  ' /* equilibrium tau_cb profile offset (ln kg/m2) */')
      CALL WRITE_0D_RL( sedTcbSlope, INDEX_NONE, 'sedTcbSlope =',
     &  ' /* equilibrium tau_cb profile slope (ln kg/m2) */')
      CALL WRITE_0D_RL( sedTcbMin, INDEX_NONE, 'sedTcbMin =',
     &  ' /* minimum equilibrium tau_cb (Pa) */')
      CALL WRITE_0D_RL( sedTcbMax, INDEX_NONE, 'sedTcbMax =',
     &  ' /* maximum equilibrium tau_cb (Pa) */')
      CALL WRITE_0D_RL( sedTcons, INDEX_NONE, 'sedTcons =',
     &  ' /* consolidation time scale (s) */')
      CALL WRITE_0D_RL( sedTswell, INDEX_NONE, 'sedTswell =',
     &  ' /* swelling time scale (s) */')
      CALL WRITE_0D_RL( sedZaCoef, INDEX_NONE, 'sedZaCoef =',
     &  ' /* active layer growth with excess stress (m/Pa) */')
      CALL WRITE_0D_RL( sedZaMin, INDEX_NONE, 'sedZaMin =',
     &  ' /* minimum active layer thickness (m) */')
      CALL WRITE_0D_RL( sedLayMass, INDEX_NONE, 'sedLayMass =',
     &  ' /* mass that starts a new bed layer (kg/m2) */')
      ENDIF
      CALL WRITE_0D_I( sedWaveMode, INDEX_NONE, 'sedWaveMode =',
     &  ' /* waves: 0 none, 1 files, 2 parametric */')
      IF ( sedWaveMode.EQ.2 ) THEN
      CALL WRITE_0D_RL( sedWindSpeed, INDEX_NONE, 'sedWindSpeed =',
     &  ' /* 10 m wind speed (m/s) */')
      CALL WRITE_0D_RL( sedWindDir, INDEX_NONE, 'sedWindDir =',
     &  ' /* wind blows toward (deg ccw from +x) */')
      CALL WRITE_0D_RL( sedFetch, INDEX_NONE, 'sedFetch =',
     &  ' /* fetch (m) */')
      ENDIF
      IF ( sedWaveMode.GT.0 ) THEN
      CALL WRITE_0D_RL( sedWaveGamma, INDEX_NONE, 'sedWaveGamma =',
     &  ' /* breaking: Hs <= gamma*depth */')
      CALL WRITE_0D_RL( sedZ0, INDEX_NONE, 'sedZ0 =',
     &  ' /* bed roughness length for waves (m) */')
      CALL WRITE_0D_L( sedWaveIce, INDEX_NONE, 'sedWaveIce =',
     &  ' /* waves damped by sea-ice area */')
      ENDIF
      _END_MASTER(myThid)
#endif /* ALLOW_SEDIMENT */

      RETURN
      END
