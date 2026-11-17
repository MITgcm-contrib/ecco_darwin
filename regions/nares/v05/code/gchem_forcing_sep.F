C pan_AO LOCAL MOD of pkg/gchem/gchem_forcing_sep.F
C (darwin3_ckpt69k, copied verbatim 2026-10-03; the ONLY changes are the four
C  TIMER_START/TIMER_STOP pairs marked "pan_AO timer" below. No logic, no
C  arguments, no ordering is altered. Diff against
C  $DARWIN3/pkg/gchem/gchem_forcing_sep.F before trusting that claim.)
C
C WHY THIS MOD EXISTS
C -------------------
C Turning diagnostics on costs this configuration 4.6x wall clock -- 2795 s
C -> 12753 s for the same 288-step model day on the same binary -- which is
C 12 -> 54 days wall per model year and makes production infeasible. MITgcm's
C own section timers localise it to here and nowhere else:
C
C     section              no-diags   with-diags
C     GCHEM_FORCING_SEP      1129.0      9542.2   <- 75% of the run
C       of which DARWIN_FORCING 294.3     1125.5
C     THERMODYNAMICS         1030.2      2472.7
C     DYNAMICS                157.3       131.2
C     DO_THE_MODEL_IO           4.6        13.6   <- writing is free
C
C So 8417 s -- 66% of the entire run -- is inside GCHEM_FORCING_SEP but
C OUTSIDE the one block that was already timed. Everything else in this
C routine is untimed, so the existing instrumentation can only say "not
C darwin_forcing". These four timers split the remainder into the blocks that
C can plausibly carry it.
C
C Note the enabling mechanism, because it is counter-intuitive: MITgcm sets
C useDiagnostics=.FALSE. when no diagnostics list is active. The "fast"
C baseline was therefore not merely discarding output -- it was skipping
C every IF (useDiagnostics) block in the model, including the 33-ptracer
C 3-D DIAGNOSTICS_FILL loop below, which runs every timestep whether or not
C any stream has asked for a GC_Tr* field.
C
C HYPOTHESES, in the order this instrumentation tests them:
C   GC_DARWIN_CONS1/2  global reductions, twice per step. A collective is a
C                      synchronisation point, so rank 0's timer absorbs the
C                      load imbalance of every other rank -- which would also
C                      explain the standing ~22% rank occupancy. Leading
C                      suspect.
C   GC_DIAGFILL        the 33-ptracer x Nr fill loop that useDiagnostics
C                      gates.
C   GC_OBCS_APPLY_PTR  33 OBCS applications per tile per step.
C
C TIMER BUDGET: eesupp/src/timers.F has maxTimers = 50 and this run registers
C 44. These four take it to 48. Adding a fifth timer anywhere is fine; a
C seventh overflows and hits the GOTO 904 error path. Count before adding.
C
C FIXED FORM: every TIMER call below is split across a continuation line, not
C because it must be but because it only just must. Written on one line,
C 'GC_DARWIN_CONS1    [GCHEM_FORCING_SEP]' lands on column 73 at this
C indentation and build 25226744 died on it, while the identically-shaped
C GC_DIAGFILL call one block up landed on exactly column 72 and compiled.
C A padded timer name is one character from the cliff; continue them all.
C
C This mod is instrumentation only and is safe to keep. If it is ever
C removed, remove the whole file -- do not hand-edit the timers out of a copy
C that has drifted from upstream.

#include "GCHEM_OPTIONS.h"
#ifdef ALLOW_DIC
# include "DIC_OPTIONS.h"
#endif
#ifdef ALLOW_BLING
# include "BLING_OPTIONS.h"
#endif
#ifdef ALLOW_DARWIN
# include "DARWIN_OPTIONS.h"
#endif

CBOP
C !ROUTINE: GCHEM_FORCING_SEP
C !INTERFACE: ==========================================================
      SUBROUTINE GCHEM_FORCING_SEP( myTime, myIter, myThid )

C !DESCRIPTION:
C     calls subroutine that will update passive tracers values
C     with a separate timestep. Since GCHEM_FORCING_SEP is now
C     called before DO_FIELDS_BLOCKING_EXCHANGES, the passive
C     tracer values in the halo regions are not up to date and
C     must not be used.

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "PTRACERS_FIELDS.h"
#include "GCHEM.h"
#ifdef ALLOW_DIC
# include "DIC_VARS.h"
#endif /* ALLOW_DIC */
#ifdef ALLOW_BLING
# include "BLING_VARS.h"
#endif /* ALLOW_BLING */
#ifdef ALLOW_DARWIN
# include "DARWIN_SIZE.h"
# include "DARWIN_INDICES.h"
# include "DARWIN_FLUX.h"
#endif

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      _RL myTime
      INTEGER myIter, myThid
CEOP

#ifdef ALLOW_GCHEM
#ifdef GCHEM_SEPARATE_FORCING

C!LOCAL VARIABLES: ====================================================
C  i,j                  :: loop indices
C  bi,bj                :: tile indices
C  k                    :: vertical level
      INTEGER bi,bj,iMin,iMax,jMin,jMax
c     INTEGER i,j
      PARAMETER( iMin = 1 , iMax = sNx )
      PARAMETER( jMin = 1 , jMax = sNy )
#if (defined ALLOW_OBCS) || (defined ALLOW_DIAGNOSTICS)
      INTEGER iTr
#endif
#ifdef ALLOW_DIAGNOSTICS
      CHARACTER*8 diagName
#endif /* ALLOW_DIAGNOSTICS */

#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_ENTER('GCHEM_FORCING_SEP',myThid)
#endif

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN
C--   pan_AO timer 1/4
       CALL TIMER_START('GC_DIAGFILL        [GCHEM_FORCING_SEP]',
     &                  myThid)
C--   fill-in tracer diagnostics before any GChem udate
       DO iTr = 1,gchem_sepFTr_num
        diagName = '        '
        WRITE(diagName,'(A5,A2)') 'GC_Tr', PTRACERS_ioLabel(iTr)
        CALL DIAGNOSTICS_FILL( pTracer(1-OLx,1-OLy,1,1,1,iTr), diagName,
     &                         0, Nr, 0, 1, 1, myThid )
       ENDDO
       CALL TIMER_STOP ('GC_DIAGFILL        [GCHEM_FORCING_SEP]',
     &                  myThid)
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

ccccccccccccccccccccccccc
c global calculations   c
ccccccccccccccccccccccccc
#ifdef ALLOW_DARWIN
      IF ( useDARWIN ) THEN
C--   pan_AO timer 2/4: global reductions, before the tile loop
        CALL TIMER_START('GC_DARWIN_CONS1    [GCHEM_FORCING_SEP]',
     &                   myThid)
        CALL DARWIN_CONSERVE_SURFACE( myTime, myIter, myThid )
        CALL DARWIN_CONS( 1, myTime, myIter, myThid )
        CALL DARWIN_CONS_RESET( myThid )
        CALL TIMER_STOP ('GC_DARWIN_CONS1    [GCHEM_FORCING_SEP]',
     &                   myThid)
      ENDIF
#endif

#ifdef ALLOW_OLD_VIRTUALFLUX
#ifdef ALLOW_DIC
# ifdef ALLOW_AUTODIFF
      IF ( .NOT.useDIC ) STOP 'ABNORMAL END: S/R GCHEM_FORCING_SEP (1)'
# else /* ALLOW_AUTODIFF */
      IF ( useDIC ) THEN
# endif /* ALLOW_AUTODIFF */
c find global surface averages
       gsm_s = 0. _d 0
       gsm_dic = 0. _d 0
       gsm_alk = 0. _d 0
       CALL GCHEM_SURFMEAN(salt,gsm_s,myThid)
       CALL GCHEM_SURFMEAN(
     &             pTracer(1-OLx,1-OLy,1,1,1,1), gsm_dic, myThid )
       CALL GCHEM_SURFMEAN(
     &             pTracer(1-OLx,1-OLy,1,1,1,2), gsm_alk, myThid )
# ifndef ALLOW_AUTODIFF
      ENDIF
# endif /* ALLOW_AUTODIFF */
#endif /* ALLOW_DIC */
#ifdef ALLOW_DARWIN
c     IF ( useDARWIN ) THEN
c find global surface averages
       gsm_s = 0. _d 0
       gsm_dic = 0. _d 0
       gsm_alk = 0. _d 0
       CALL GCHEM_SURFMEAN(salt,gsm_s,myThid)
       CALL GCHEM_SURFMEAN(
     &             pTracer(1-OLx,1-OLy,1,1,1,iDIC), gsm_dic, myThid )
       CALL GCHEM_SURFMEAN(
     &             pTracer(1-OLx,1-OLy,1,1,1,iALK), gsm_alk, myThid )
c     ENDIF
#endif
ccccccccccccccccccccccccccccccccccccccccccc
#endif /* ALLOW_OLD_VIRTUALFLUX */

ccccccccccccccccccccccccc
c chemical forcing      c
ccccccccccccccccccccccccc
C$taf loop = parallel
      DO bj=myByLo(myThid),myByHi(myThid)
C$taf loop = parallel
       DO bi=myBxLo(myThid),myBxHi(myThid)

ccccccccccccccccccccccccccc DIC cccccccccccccccccccccccccccccccc
#ifdef ALLOW_DIC
# ifdef ALLOW_AUTODIFF
        IF (.NOT.useDIC) STOP 'ABNORMAL END: S/R GCHEM_FORCING_SEP (2)'
# else /* ALLOW_AUTODIFF */
        IF ( useDIC ) THEN
# endif /* ALLOW_AUTODIFF */
#ifdef ALLOW_DEBUG
          IF (debugMode) CALL DEBUG_CALL('DIC_BIOTIC_FORCING',myThid)
#endif
#ifdef ALLOW_FE
          CALL DIC_BIOTIC_FORCING( pTracer(1-OLx,1-OLy,1,bi,bj,1),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,2),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,3),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,4),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,5),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,6),
     &                          bi, bj, iMin, iMax, jMin, jMax,
     &                          myTime, myIter, myThid )
#else
#ifdef ALLOW_O2
          CALL DIC_BIOTIC_FORCING( pTracer(1-OLx,1-OLy,1,bi,bj,1),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,2),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,3),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,4),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,5),
     &                          bi, bj, iMin, iMax, jMin, jMax,
     &                          myTime, myIter, myThid )
#else
          CALL DIC_BIOTIC_FORCING( pTracer(1-OLx,1-OLy,1,bi,bj,1),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,2),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,3),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,4),
     &                          bi, bj, iMin, iMax, jMin, jMax,
     &                          myTime, myIter, myThid )
#endif
#endif
# ifndef ALLOW_AUTODIFF
        ENDIF
# endif /* ALLOW_AUTODIFF */
#endif /* ALLOW_DIC */
cccccccccccccccccccccccccc END DIC cccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccc BLING cccccccccccccccccccccccccccccccc
#ifdef ALLOW_BLING
        IF ( useBLING ) THEN
# ifndef USE_BLING_V1
#  ifdef USE_SIBLING
          CALL BLING_MAIN( pTracer(1-OLx,1-OLy,1,bi,bj,1),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,2),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,3),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,4),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,5),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,6),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,7),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,8),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,9),
#   ifdef ADVECT_PHYTO
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,10),
#   endif
     &                          bi, bj, iMin, iMax, jMin, jMax,
     &                          myTime, myIter, myThid )
#  else /* USE_SIBLING */
          CALL BLING_MAIN( pTracer(1-OLx,1-OLy,1,bi,bj,1),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,2),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,3),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,4),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,5),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,6),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,7),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,8),
#   ifdef ADVECT_PHYTO
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,9),
#   endif
     &                          bi, bj, iMin, iMax, jMin, jMax,
     &                          myTime, myIter, myThid )
#  endif /* USE_SIBLING */
# else /* USE_BLING_v1 */
          CALL BLING_MAIN( pTracer(1-OLx,1-OLy,1,bi,bj,1),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,2),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,3),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,4),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,5),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,6),
#  ifdef ADVECT_PHYTO
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,7),
#  endif
     &                          bi, bj, iMin, iMax, jMin, jMax,
     &                          myTime, myIter, myThid )
# endif /* USE_BLING_v1 */
        ENDIF
#endif /* ALLOW_BLING */
cccccccccccccccccccccccccc END BLING cccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccc DARWIN ccccccccccccccccccccccccccccccccccccc
#ifdef ALLOW_DARWIN
        IF ( useDARWIN ) THEN
#ifdef NUT_SUPPLY
c articficial supply of nutrients
#ifdef ALLOW_DEBUG
          IF (debugMode) CALL DEBUG_CALL('DARWIN_NUT_SUPPLY',myThid)
#endif
          CALL DARWIN_NUT_SUPPLY( pTracer(1-OLx,1-OLy,1,bi,bj,iPO4),
     &                          bi, bj, iMin, iMax, jMin, jMax,
     &                          myIter, myTime, myThid )
          CALL DARWIN_NUT_SUPPLY( pTracer(1-OLx,1-OLy,1,bi,bj,iNO3),
     &                          bi, bj, iMin, iMax, jMin, jMax,
     &                          myIter, myTime, myThid )
          CALL DARWIN_NUT_SUPPLY( pTracer(1-OLx,1-OLy,1,bi,bj,iFeT),
     &                          bi, bj, iMin, iMax, jMin, jMax,
     &                          myIter, myTime, myThid )
          CALL DARWIN_NUT_SUPPLY( pTracer(1-OLx,1-OLy,1,bi,bj,iSiO2),
     &                          bi, bj, iMin, iMax, jMin, jMax,
     &                          myIter, myTime, myThid )
#endif
ccccccccccccccc
C darwin_forcing operates on bi,bj part only, but needs to get full
C array because of last (iPtr) index
#ifdef ALLOW_DEBUG
          IF (debugMode) CALL DEBUG_CALL('DARWIN_FORCING',myThid)
#endif
          CALL TIMER_START('DARWIN_FORCING [GCHEM_FORCING_SEP]',myThid)
          CALL DARWIN_FORCING(  pTracer(1-OLx,1-OLy,1,1,1,1),
     &                          bi, bj, iMin, iMax, jMin, jMax,
     &                          myIter, myTime, myThid )
          CALL TIMER_STOP ('DARWIN_FORCING [GCHEM_FORCING_SEP]',myThid)
        ENDIF
#endif /* ALLOW_DARWIN */
cccccccccccccccccccccccccc END DARWIN ccccccccccccccccccccccccccccccccc

#ifdef ALLOW_OBCS
C--   Apply (again) open boundary conditions for each passive tracer
C Note: could skip this 2nd call to OBCS_APPLY if all DIC/DARWIN
C       updates of ptracers were only done in the interior (i.e. with
C       tendency multiplied by maskInC)
        IF ( useOBCS .AND. .NOT.useDIC ) THEN
#ifdef ALLOW_DEBUG
          IF (debugMode) CALL DEBUG_CALL('OBCS_APPLY_PTRACER',myThid)
#endif
C--   pan_AO timer 3/4: 33 tracers x every tile x every step
          CALL TIMER_START('GC_OBCS_APPLY_PTR  [GCHEM_FORCING_SEP]',
     &                     myThid)
          DO iTr = 1,gchem_sepFTr_num
            CALL OBCS_APPLY_PTRACER(
     I                bi, bj, 0, iTr,
     U                pTracer(1-OLx,1-OLy,1,bi,bj,iTr),
     I                myThid )
          ENDDO
          CALL TIMER_STOP ('GC_OBCS_APPLY_PTR  [GCHEM_FORCING_SEP]',
     &                     myThid)
        ENDIF
#endif /* ALLOW_OBCS */

C--   end bi,bj loops.
       ENDDO
      ENDDO

#ifdef ALLOW_DARWIN
      IF ( useDARWIN ) THEN
C--   pan_AO timer 4/4: global reductions, after the tile loop. A collective
C--   here absorbs every other rank's load imbalance into rank 0's number --
C--   read this one together with the rank occupancy, not on its own.
         CALL TIMER_START('GC_DARWIN_CONS2    [GCHEM_FORCING_SEP]',
     &                    myThid)
         CALL DARWIN_CONS( 2, myTime, myIter, myThid )
         CALL DARWIN_ATMOS( myTime, myIter, myThid )
         CALL TIMER_STOP ('GC_DARWIN_CONS2    [GCHEM_FORCING_SEP]',
     &                    myThid)
      ENDIF
#endif /* ALLOW_DARWIN */

#ifdef ALLOW_DIC
# ifdef ALLOW_AUTODIFF
      IF ( .NOT.useDIC ) STOP 'ABNORMAL END: S/R GCHEM_FORCING_SEP (3)'
# else /* ALLOW_AUTODIFF */
      IF ( useDIC ) THEN
# endif /* ALLOW_AUTODIFF */
#ifdef ALLOW_DEBUG
       IF (debugMode) CALL DEBUG_CALL('DIC_ATMOS',myThid)
#endif
       CALL DIC_ATMOS( myTime, myIter, myThid )
# ifdef COMPONENT_MODULE
       CALL DIC_STORE_FLUXCO2( myTime, myIter, myThid )
# endif
# ifdef ALLOW_COST
       CALL DIC_COST( myTime, myIter, myThid )
# endif
# ifndef ALLOW_AUTODIFF
      ENDIF
# endif /* ALLOW_AUTODIFF */
#endif /* ALLOW_DIC */

#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_LEAVE('GCHEM_FORCING_SEP',myThid)
#endif

#endif /* GCHEM_SEPARATE_FORCING */
#endif /* ALLOW_GCHEM */

      RETURN
      END
