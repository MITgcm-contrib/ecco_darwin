#include "WAD_OPTIONS.h"
#include "MOM_COMMON_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: WAD_CHECK

C !INTERFACE:
      SUBROUTINE WAD_CHECK( myThid )

C     !DESCRIPTION:
C     Check pkg/wad parameters and dependencies on the main model
C     configuration (non-linear free surface with r* coordinate).

C     !USES:
      IMPLICIT NONE
C     == Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#ifdef ALLOW_GENERIC_ADVDIFF
# include "GAD.h"
#endif
#ifdef ALLOW_PTRACERS
# include "PTRACERS_SIZE.h"
# include "PTRACERS_PARAMS.h"
#endif
#include "WAD.h"

C     !INPUT PARAMETERS:
C     myThid   :: my Thread Id number
      INTEGER myThid
CEOP

#ifdef ALLOW_WAD
C     !LOCAL VARIABLES:
C     msgBuf   :: Informational/error message buffer
C     errCount :: number of fatal parameter errors
C     maxRcol  :: largest resting column thickness in the domain [m]
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER errCount
#ifdef ALLOW_PTRACERS
      INTEGER iTr, nTrAB
#endif
      INTEGER i, j, bi, bj
      _RL maxRcol

C--   Largest resting column thickness (all threads, all procs)
      maxRcol = 0. _d 0
      DO bj=myByLo(myThid), myByHi(myThid)
       DO bi=myBxLo(myThid), myBxHi(myThid)
        DO j=1,sNy
         DO i=1,sNx
          IF ( maskInC(i,j,bi,bj).NE.zeroRS ) maxRcol = MAX( maxRcol,
     &         Ro_surf(i,j,bi,bj) - R_low(i,j,bi,bj) )
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      _GLOBAL_MAX_RL( maxRcol, myThid )

      _BEGIN_MASTER(myThid)

      WRITE(msgBuf,'(A)') 'WAD_CHECK: #define ALLOW_WAD'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )

C--  Print out some key parameters :
      CALL WRITE_0D_RL( wadMinDepth, INDEX_NONE, 'wadMinDepth =',
     &  ' /* thickness of the film kept in dry columns (m) */')
      CALL WRITE_0D_RL( wadCritDepth, INDEX_NONE, 'wadCritDepth =',
     &  ' /* face opens when surface > crest + this (m) */')
      CALL WRITE_0D_RL( wadDragDepth, INDEX_NONE, 'wadDragDepth =',
     &  ' /* velocity relaxed to 0 below this depth (m) */')
      CALL WRITE_0D_RL( wadAdvDepth, INDEX_NONE, 'wadAdvDepth =',
     &  ' /* momentum advection ramped off below this depth (m) */')
      CALL WRITE_0D_RL( wadGMDepth, INDEX_NONE, 'wadGMDepth =',
     &  ' /* GM/Redi tapered off below this depth (m) */')
      CALL WRITE_0D_RL( wadKPPDepth, INDEX_NONE, 'wadKPPDepth =',
     &  ' /* KPP non-local flux ramped off below this depth (m) */')
      CALL WRITE_0D_RL( wadKPPCapDepth, INDEX_NONE, 'wadKPPCapDepth =',
     &  ' /* KPP non-local fraction capped at 1 below this (m) */')
      CALL WRITE_0D_RL( wadManningN, INDEX_NONE, 'wadManningN =',
     &  ' /* Manning n for depth-dependent drag (<=0: off) */')
      CALL WRITE_0D_RL( wadDragMax, INDEX_NONE, 'wadDragMax =',
     &  ' /* upper limit of that drag coefficient */')
      CALL WRITE_0D_RL( wadMaxFroude, INDEX_NONE, 'wadMaxFroude =',
     &  ' /* cap |u| at this times sqrt(g D) (0: off) */')
      CALL WRITE_0D_RL( wadSmoothWidth, INDEX_NONE, 'wadSmoothWidth =',
     &  ' /* width of smooth face mask (m) */')
      CALL WRITE_0D_RL( wadMonFreq, INDEX_NONE, 'wadMonFreq =',
     &  ' /* WAD monitor frequency (s) */')
      CALL WRITE_0D_RL( wadMaxSpeed, INDEX_NONE, 'wadMaxSpeed =',
     &  ' /* stop if an open-face speed exceeds this (m/s) */')
      CALL WRITE_0D_L( wadConserveVol, INDEX_NONE, 'wadConserveVol =',
     &  ' /* take film-floor volume from wet neighbours */')
      CALL WRITE_0D_L( wadCarryVel, INDEX_NONE, 'wadCarryVel =',
     &  ' /* opened face takes upstream u* */')
      CALL WRITE_0D_L( wadUpwindFace, INDEX_NONE, 'wadUpwindFace =',
     &  ' /* NLFS face thickness from upwind surface */')
      CALL WRITE_0D_L( wadDryForcing, INDEX_NONE, 'wadDryForcing =',
     &  ' /* no surface forcing on dry columns */')
      CALL WRITE_0D_C( wadManningFile, -1, INDEX_NONE,
     &  'wadManningFile =', ' /* extra Manning roughness map */')
      CALL WRITE_0D_RL( wadManningT1, INDEX_NONE, 'wadManningT1 =',
     &  ' /* extra roughness full until this time [s] */')
      CALL WRITE_0D_RL( wadManningT2, INDEX_NONE, 'wadManningT2 =',
     &  ' /* extra roughness gone at this time [s] */')

C--  Check parameters:
      errCount = 0
#ifndef ALLOW_BOTTOMDRAG_ROUGHNESS
      IF ( wadManningN.GT.zeroRL .OR. wadManningFile.NE.' ' ) THEN
        WRITE(msgBuf,'(2A)') 'WAD_CHECK: wadManningN/File needs',
     &   ' #define ALLOW_BOTTOMDRAG_ROUGHNESS (MOM_COMMON_OPTIONS.h)'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF
#else
      IF ( (wadManningN.GT.zeroRL .OR. wadManningFile.NE.' ')
     &     .AND. zRoughBot.GT.zeroRL ) THEN
        WRITE(msgBuf,'(2A)') 'WAD_CHECK: wadManningN and zRoughBot',
     &   ' both set the bottom drag coefficient; use one'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF
#endif

#ifndef NONLIN_FRSURF
      WRITE(msgBuf,'(A)')
     &  'WAD_CHECK: pkg/wad needs #define NONLIN_FRSURF'
      CALL PRINT_ERROR( msgBuf, myThid )
      errCount = errCount + 1
#endif
      IF ( .NOT.usingZCoords ) THEN
        WRITE(msgBuf,'(A)')
     &  'WAD_CHECK: pkg/wad only works with height coordinates'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF
      IF ( nonlinFreeSurf.LT.3 ) THEN
C-    the CG2D matrix must be rebuilt every step to follow the face mask
        WRITE(msgBuf,'(A,I3,A)') 'WAD_CHECK: nonlinFreeSurf=',
     &   nonlinFreeSurf, ' but pkg/wad needs nonlinFreeSurf >= 3'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF
      IF ( select_rStar.LT.0 .OR. selectSigmaCoord.NE.0 ) THEN
        WRITE(msgBuf,'(A)')
     &  'WAD_CHECK: pkg/wad needs select_rStar >= 0 and no sigma coord.'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF
      IF ( select_rStar.EQ.0 .AND. hFacInf*drF(1).GT.wadMinDepth ) THEN
C-    surface-level NLFS clips the surface at Rmin_surf (hFacInf) without
C     telling WAD: the film must be thicker than that clip
        WRITE(msgBuf,'(A,1PE10.3,A)') 'WAD_CHECK: select_rStar=0 needs'
     &   //' hFacInf <=', wadMinDepth/drF(1), ' (wadMinDepth/drF(1))'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF
      IF ( .NOT.exactConserv ) THEN
        WRITE(msgBuf,'(A)') 'WAD_CHECK: pkg/wad needs exactConserv=T'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF
      IF ( .NOT.implicitDiffusion ) THEN
        WRITE(msgBuf,'(A)')
     &  'WAD_CHECK: pkg/wad needs implicitDiffusion=T (thin films)'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF
      IF ( wadMinDepth.LE.zeroRL ) THEN
        WRITE(msgBuf,'(A)') 'WAD_CHECK: wadMinDepth must be > 0'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF
      IF ( wadCritDepth.LE.wadMinDepth ) THEN
        WRITE(msgBuf,'(A)')
     &  'WAD_CHECK: wadCritDepth must be > wadMinDepth'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF
      IF ( wadCarryVel .AND. Nr.GT.1 .AND. .NOT.momImplVertAdv ) THEN
C-    a face opening at full speed into a thin r* column drives a large
C     vertical velocity through cm-thin layers: explicit vertical
C     advection of momentum blows up
        WRITE(msgBuf,'(A)')
     &  'WAD_CHECK: wadCarryVel with Nr > 1 needs momImplVertAdv=.TRUE.'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF
C-    Combinations that pkg/wad does not handle
C-    useRealFreshWaterFlux: supported (evaporation in
C     the outflow limiter and the dry-column ramp, EmPmR in the budget)
      IF ( useRealFreshWaterFlux ) THEN
        WRITE(msgBuf,'(2A)') 'WAD_CHECK: useRealFreshWaterFlux:',
     &   ' evaporation limited on thin and dry columns (WADEVLIM)'
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
      ENDIF
      IF ( useAbsVorticity .AND. momAdvection .AND.
     &     vectorInvariantMomentum ) THEN
C-    f is then inside the advection term that WAD ramps off in thin water
        WRITE(msgBuf,'(A)') 'WAD_CHECK: useAbsVorticity=.TRUE. would'
     &   //' ramp off Coriolis in thin water; not supported'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF
      IF ( implicDiv2DFlow.LT.oneRL ) THEN
C-    etaH halos are then not updated before WAD_APPLY_FLOOR
        WRITE(msgBuf,'(A)') 'WAD_CHECK: pkg/wad needs implicDiv2DFlow=1'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF
      IF ( applyExchUV_early ) THEN
C-    uVel halos would miss the outflow limiter
        WRITE(msgBuf,'(A)')
     &   'WAD_CHECK: applyExchUV_early=.TRUE. is not supported'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF
      IF ( useThSIce ) THEN
        WRITE(msgBuf,'(A)') 'WAD_CHECK: pkg/thsice is not supported'
     &   //' (ice would grow on dry columns)'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF
      IF ( useSALT_PLUME ) THEN
        WRITE(msgBuf,'(A)') 'WAD_CHECK: pkg/salt_plume is not supported'
     &   //' (brine bypasses the dry-column ramp)'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF
      IF ( useSEAICE ) THEN
C-    sea ice on drying flats (masks, grounding, dry-column growth) is
C     not part of this version of pkg/wad
        WRITE(msgBuf,'(A)') 'WAD_CHECK: pkg/seaice is not supported'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF
      IF ( useGCHEM ) THEN
        WRITE(msgBuf,'(2A)') '** WARNING ** WAD_CHECK: pkg/gchem ',
     &   'air-sea fluxes are not ramped off on dry columns'
        CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
      ENDIF
      IF ( wadSmoothWidth.GT.zeroRL ) THEN
        WRITE(msgBuf,'(A)')
     &  'WAD_CHECK: wadSmoothWidth > 0 is not implemented yet'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF
#ifdef ALLOW_SHELFICE_REMESHING
      IF ( useShelfIce ) THEN
        WRITE(msgBuf,'(A)')
     &  'WAD_CHECK: pkg/wad is not compatible with shelfice remeshing'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF
#endif

C-   The r* "too small" STOP in CALC_R_STAR is bypassed with WAD, so a
C    hFacInf above the film fraction is harmless, but worth noting
      IF ( maxRcol.GT.zeroRL .AND.
     &     hFacInf*maxRcol.GT.wadMinDepth ) THEN
        WRITE(msgBuf,'(2A,1PE12.4,A)') '** WARNING ** WAD_CHECK: ',
     &   'hFacInf >', wadMinDepth/maxRcol,
     &   ' = wadMinDepth/max(Rcol) ; r* STOP is bypassed by WAD'
        CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
      ENDIF

C-   Momentum advection in thin columns (see MITgcm_WAD_design.md):
C    flux form divides by hFacW and blows up at the drying front ;
C    vector-invariant needs upwindShear when Nr > 1
      IF ( momStepping .AND. momAdvection ) THEN
       IF ( .NOT.vectorInvariantMomentum ) THEN
        WRITE(msgBuf,'(2A)') '** WARNING ** WAD_CHECK: flux-form ',
     &   'momentum advection is unstable at drying fronts ;'
        CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(2A)') '** WARNING ** WAD_CHECK: ',
     &   'use vectorInvariantMomentum=.TRUE.'
        CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
       ELSEIF ( Nr.GT.1 .AND. .NOT.upwindShear ) THEN
        WRITE(msgBuf,'(2A)') '** WARNING ** WAD_CHECK: centred ',
     &   'vertical shear is unstable in thin columns ;'
        CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(2A)') '** WARNING ** WAD_CHECK: ',
     &   'use upwindShear=.TRUE.'
        CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
       ENDIF
      ENDIF

C-    explicit quadratic bottom drag on a layer of thickness h is unstable
C     when Cd |u| dt / h > 2 ; thin partial cells in drying columns reach
C     this (e.g. Cd=2.5e-3, |u|=1 m/s, dt=6 s: h < 7.5 mm)
      IF ( momStepping .AND. Nr.GT.1 .AND. selectImplicitDrag.EQ.0
     &     .AND. bottomDragQuadratic.GT.zeroRL ) THEN
        WRITE(msgBuf,'(2A,1PE10.3,A)') '** WARNING ** WAD_CHECK: ',
     &   'explicit quadratic bottom drag is unstable in layers < ',
     &   bottomDragQuadratic*deltaTMom/2. _d 0, ' m (|u| = 1 m/s) ;'
        CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(2A)') '** WARNING ** WAD_CHECK: ',
     &   'with thin partial cells use selectImplicitDrag=2'
        CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(2A)') '** WARNING ** WAD_CHECK: ',
     &   '(#define ALLOW_SOLVE4_PS_AND_DRAG) or a larger hFacMin'
        CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
      ENDIF

#ifdef ALLOW_GENERIC_ADVDIFF
C-   Adams-Bashforth extrapolation of the tracer tendency is inconsistent
C    with the fast volume changes of thin columns: a uniform tracer does
C    not stay uniform (and can blow up with films of a few mm).
      IF ( ( tempStepping .AND. (AdamsBashforthGt.OR.AdamsBashforth_T) )
     & .OR.( saltStepping .AND. (AdamsBashforthGs.OR.AdamsBashforth_S) )
     &   ) THEN
        WRITE(msgBuf,'(2A)') '** WARNING ** WAD_CHECK: ',
     &   'T or S uses an Adams-Bashforth advection scheme ;'
        CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(2A)') '** WARNING ** WAD_CHECK: ',
     &   'prefer a 2-level scheme (e.g. 33, 77) in drying areas'
        CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
      ENDIF
#endif
#ifdef ALLOW_PTRACERS
      IF ( usePTRACERS ) THEN
        nTrAB = 0
        DO iTr=1,PTRACERS_numInUse
          IF ( PTRACERS_AdamsBashGtr(iTr) .OR.
     &         PTRACERS_AdamsBash_Tr(iTr) ) nTrAB = nTrAB + 1
        ENDDO
        IF ( nTrAB.GT.0 ) THEN
          WRITE(msgBuf,'(2A,I4,A)') '** WARNING ** WAD_CHECK: ',
     &     'pTracers:', nTrAB, ' use an Adams-Bashforth scheme ;'
          CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                        SQUEEZE_RIGHT, myThid )
          WRITE(msgBuf,'(2A)') '** WARNING ** WAD_CHECK: ',
     &     'prefer a 2-level scheme (e.g. 33, 77) in drying areas'
          CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                        SQUEEZE_RIGHT, myThid )
        ENDIF
      ENDIF
#endif

      IF ( errCount.GT.0 ) THEN
        WRITE(msgBuf,'(A,I3,A)')
     &       'WAD_CHECK: detected', errCount,' fatal error(s)'
        CALL PRINT_ERROR( msgBuf, myThid )
        CALL ALL_PROC_DIE( 0 )
        STOP 'ABNORMAL END: S/R WAD_CHECK'
      ENDIF

      _END_MASTER(myThid)

#endif /* ALLOW_WAD */

      RETURN
      END
