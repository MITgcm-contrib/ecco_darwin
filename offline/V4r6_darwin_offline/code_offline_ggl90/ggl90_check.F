#include "GGL90_OPTIONS.h"

C     Patched (local override) to relax the implicitViscosity requirement
C     when running under pkg/offline -- see the header comment below at
C     the actual change for why. Everything else is verbatim from the
C     stock darwin3 pkg/ggl90/ggl90_check.F (backport_ckpt68y).

      SUBROUTINE GGL90_CHECK( myThid )
C     *==========================================================*
C     | SUBROUTINE GGL90_CHECK                                   |
C     | o Validate basic package setup and inter-package         |
C     | dependencies.                                            |
C     *==========================================================*
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GGL90.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_PARAMS.h"
#endif

C     === Routine arguments ===
C     myThid   :: Number of this instance of GGL90_CHECK
      INTEGER myThid

#ifdef ALLOW_GGL90

C     === Local variables ===
C     errCount :: error counter
C     msgBuf   :: Informational/error message buffer
      INTEGER errCount
      CHARACTER*(MAX_LEN_MBUF) msgBuf

      _BEGIN_MASTER(myThid)
      errCount = 0

      WRITE(msgBuf,'(A)') 'GGL90_CHECK: #define ALLOW_GGL90'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )

C     Convective adjustment might be needed even with ggl90
      IF ( cAdjFreq.NE.0. .OR. ivdc_kappa.NE.0. ) THEN
         WRITE(msgBuf,'(A)')
     &  'GGL90_CHECK: Some form of convection has been enabled'
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                       SQUEEZE_RIGHT, myThid )
      ENDIF

#ifdef ALLOW_GGL90_HORIZDIFF
#ifdef ALLOW_GGL90_SMOOTH
      IF ( GGL90diffTKEh .NE. 0. _d 0 ) THEN
         WRITE(msgBuf,'(A)')
     &      'Horizontal diffusion of TKE and smoothing are both on'
         CALL PRINT_ERROR( msgBuf, myThid )
         errCount = errCount + 1
      ENDIF
#endif
#else
      IF ( GGL90diffTKEh .NE. 0. _d 0 ) THEN
         WRITE(msgBuf,'(2A)')
     &      'Define ALLOW_GGL90_HORIZDIFF in GGL90_OPTIONS.h',
     &      ' to get diffusion of TKE'
         CALL PRINT_ERROR( msgBuf, myThid )
         errCount = errCount + 1
      ENDIF
#endif
#ifndef ALLOW_GGL90_LANGMUIR
      IF (useLANGMUIR ) THEN
         WRITE(msgBuf,'(2A)')
     &      'Define ALLOW_GGL90_LANGMUIR in GGL90_OPTIONS.h',
     &      ' to use Langmuir Parameterization'
         CALL PRINT_ERROR( msgBuf, myThid )
         errCount = errCount + 1
      ENDIF
#endif

C     GGL90 needs implicit vertical diffusion and viscosity
      IF (.NOT.implicitDiffusion ) THEN
         WRITE(msgBuf,'(A)')
     &      'GGL90 needs implicitDiffusion to be enabled'
         CALL PRINT_ERROR( msgBuf, myThid )
         errCount = errCount + 1
      ENDIF
C     PATCHED: stock code requires implicitViscosity=.TRUE.
C     unconditionally, but set_parms.F forces implicitViscosity =
C     momViscosity .AND. implicitViscosity, and momViscosity = momStepping
C     .AND. momViscosity -- so with pkg/offline's OFFLINE_RESET_PARMS
C     setting momStepping=.FALSE. (momentum is read from file, never
C     stepped), implicitViscosity is unconditionally forced back to
C     .FALSE. regardless of what the namelist sets it to. There is no way
C     to satisfy this check as written while using pkg/offline. This is
C     not a real numerical requirement in our case: GGL90's own computed
C     vertical viscosity is for the MOMENTUM equation specifically, which
C     is never solved here (velocity comes from the converted archive);
C     the vertical DIFFUSIVITY it computes for tracers is itself entirely
C     overwritten every step by GGL90diffKrFile (data.off, see
C     offline_get_diffus.F) regardless of how GGL90 internally applies
C     implicit vs explicit treatment. Discovered on the first-ever run
C     attempt combining GGL90 with pkg/offline -- apparently a novel
C     combination (ecco_darwin uses GGL90 online, where momStepping is
C     genuinely .TRUE.; Steph's original offline config used KPP instead).
      IF (.NOT.implicitViscosity .AND. .NOT.useOFFLINE) THEN
         WRITE(msgBuf,'(A)')
     &      'GGL90 needs implicitViscosity to be enabled'
         CALL PRINT_ERROR( msgBuf, myThid )
         errCount = errCount + 1
      ENDIF

C     Limitation of Langmuir Parameterization implementation:
      IF ( useLANGMUIR .AND. useAbsVorticity ) THEN
         WRITE(msgBuf,'(2A)') 'GGL90_CHECK: ',
     &      'with useAbsVorticity, missing Coriolis'
         CALL PRINT_ERROR( msgBuf, myThid )
         WRITE(msgBuf,'(2A)') 'GGL90_CHECK: ',
     &      ' contrib from Langmuir Parameterization'
         CALL PRINT_ERROR( msgBuf, myThid )
         errCount = errCount + 1
      ENDIF
      IF ( useLANGMUIR .AND. useCDscheme ) THEN
         WRITE(msgBuf,'(2A)') 'GGL90_CHECK: ',
     &      'with useCDscheme, missing Coriolis'
         CALL PRINT_ERROR( msgBuf, myThid )
         WRITE(msgBuf,'(2A)') 'GGL90_CHECK: ',
     &      ' contrib from Langmuir Parameterization'
         CALL PRINT_ERROR( msgBuf, myThid )
         errCount = errCount + 1
      ENDIF

C     Compatibility with the use of other pkgs:
      IF ( useKPP ) THEN
         WRITE(msgBuf,'(A)')
     &      'GGL90 and KPP cannot be turned on at the same time'
         CALL PRINT_ERROR( msgBuf, myThid )
         errCount = errCount + 1
      ENDIF

      IF ( usePP81 ) THEN
         WRITE(msgBuf,'(A)')
     &      'GGL90 and PP81 cannot be turned on at the same time'
         CALL PRINT_ERROR( msgBuf, myThid )
         errCount = errCount + 1
      ENDIF

      IF ( useMY82 ) THEN
         WRITE(msgBuf,'(A)')
     &      'GGL90 and MY82 cannot be turned on at the same time'
         CALL PRINT_ERROR( msgBuf, myThid )
         errCount = errCount + 1
      ENDIF

#ifdef ALLOW_AUTODIFF
C     Not everything works stably with AD code
      IF ( useGGL90inAdMode .AND. adMxlMaxFlag .GT. 1 ) THEN
       WRITE(msgBuf,'(2A,L1,A,I2,A)') '** WARNING ** GGL90_CHECK: ',
     &      'useGGL90inAdMode = ', useGGL90inAdMode,
     &      ' with adMxlMaxFlag =', adMxlMaxFlag, ' (> 1)'
       CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
       WRITE(msgBuf,'(2A)') '** WARNING ** GGL90_CHECK: ',
     &      'tends to be unstable, use (ad)mxlMaxFlag = 0 or 1.'
       CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      ENDIF
      IF ( useGGL90inAdMode .AND. adMxlMaxFlag .NE. mxlMaxFlag ) THEN
       WRITE(msgBuf,'(2A,I2,A)') '** WARNING ** GGL90_CHECK: ',
     &      'In AD-mode only, using adMxlMaxFlag =', adMxlMaxFlag,
     &      ' (instead of mxlMaxFlag).'
       CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      ENDIF
#endif

C-    Stop cleanly when error was detected
      IF ( errCount.GE.1 ) THEN
         WRITE(msgBuf,'(A,I3,A)')
     &       'GGL90_CHECK: detected', errCount,' fatal error(s)'
         CALL PRINT_ERROR( msgBuf, myThid )
         CALL ALL_PROC_DIE( 0 )
         STOP 'ABNORMAL END: S/R GGL90_CHECK'
      ENDIF
      _END_MASTER(myThid)

#endif /* ALLOW_GGL90 */

      RETURN
      END
