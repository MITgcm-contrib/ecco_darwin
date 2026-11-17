
#include "ICEPLUME_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: ICEPLUME_TENDENCY_APPLY_T

C !INTERFACE:
      SUBROUTINE ICEPLUME_TENDENCY_APPLY_T(
     U                    gT_arr,
     I                    iMin,iMax,jMin,jMax,
     I                    k, bi, bj, myTime, myIter, myThid )

C     !DESCRIPTION:
C     As for S/R ICEPLUME_TENDENCY_APPLY_T
C     Add iceplume tendency terms to T (theta) tendency.
C     Routine works for one level at a time.

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "GRID.h"
#include "EEPARAMS.h"
#include "DYNVARS.h"
#include "ICEPLUME.h"


C     !INPUT/OUTPUT PARAMETERS:
C     gT_arr    :: the tendency array
C     k         :: Current vertical level index
C     bi,bj     :: Current tile indices
C     myTime    :: Current time in simulation
C     myIter    :: Current iteration number
C     myThid    :: my Thread Id number
      _RL     gT_arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER iMin, iMax, jMin, jMax
      INTEGER k, bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid

C     !LOCAL VARIABLES:
      INTEGER i, j

#ifdef ALLOW_ICEPLUME

catn: this is the background melt tendency
      IF ( applyIcePlumeBGTendT ) THEN
       DO j=jMin,jMax
        DO i=iMin,iMax
          gT_arr(i,j) = gT_arr(i,j)
     &         + iceplumeBG_TendT(i,j,k,bi,bj)
        ENDDO
       ENDDO
      ENDIF

#endif /* ALLOW_ICEPLUME */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: ICEPLUME_TENDENCY_APPLY_S

C !INTERFACE:
      SUBROUTINE ICEPLUME_TENDENCY_APPLY_S(
     U                    gS_arr,
     I                    iMin,iMax,jMin,jMax,
     I                    k, bi, bj, myTime, myIter, myThid )

C     !DESCRIPTION:
C     As for S/R ICEPLUME_TENDENCY_APPLY_S
C     Add iceplume tendency terms to S tendency.
C     Routine works for one level at a time.

C     !INPUT PARAMETERS:
      IMPLICIT NONE
#include "SIZE.h"
#include "GRID.h"
#include "EEPARAMS.h"
#include "DYNVARS.h"
#include "ICEPLUME.h"

C     !INPUT/OUTPUT PARAMETERS:
C     gS_arr    :: the tendency array
C     k         :: Current vertical level index
C     bi,bj     :: Current tile indices
C     myTime    :: Current time in simulation
C     myIter    :: Current iteration number
C     myThid    :: my Thread Id number
      _RL     gS_arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER iMin, iMax, jMin, jMax
      INTEGER k, bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

C     !LOCAL VARIABLES:
      INTEGER i, j

#ifdef ALLOW_ICEPLUME

catn: this is the background melt tendency
      IF ( applyIcePlumeBGTendS ) THEN
       DO j=jMin,jMax
        DO i=iMin,iMax
          gS_arr(i,j) = gS_arr(i,j)
     &         + iceplumeBG_TendS(i,j,k,bi,bj)
        ENDDO
       ENDDO
      ENDIF

#endif /* ALLOW_ICEPLUME */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: ICEPLUME_TENDENCY_APPLY_PTR
C     !INTERFACE:
      SUBROUTINE ICEPLUME_TENDENCY_APPLY_PTR(
     U                    gPtracer,
     I                    iMin,iMax,jMin,jMax,
     I                    k, bi, bj, iTracer, myTime, myIter, myThid )
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | S/R ICEPLUME_TENDENCY_APPLY_PTR
C     | o Passive-tracer tendency of the plume addMass source,
C     |   the analogue of the T/S plume terms in APPLY_FORCING_T/S:
C     |   water removed (entrainment) leaves at the local
C     |   concentration (no concentration change); water added
C     |   (plume outflow) brings the plume concentration
C     |   ptr_AddMass3D, which includes the discharge's own
C     |   concentration (ptracerMask when useInputPtracers, else 0).
C     |   Background ice-front meltwater (addMass3Dbg) carries no
C     |   tracer (concentration 0), as T/S get their own BG terms.
C     | Called from PTRACERS_APPLY_FORCING (hook added by
C     | utils/copy_iceplume_pkg_to_mitgcm.py).
C     *==========================================================*
C     \ev
C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "ICEPLUME.h"
#ifdef ALLOW_PTRACERS
C     (PTRACERS_SIZE.h is already included by ICEPLUME.h)
#include "PTRACERS_PARAMS.h"
#include "PTRACERS_FIELDS.h"
#endif

C     !INPUT/OUTPUT PARAMETERS:
      _RL     gPtracer(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER iMin, iMax, jMin, jMax
      INTEGER k, bi, bj, iTracer
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP
#if (defined ALLOW_ICEPLUME) && (defined ALLOW_PTRACERS)
#ifdef ALLOW_ADDFLUID
C     !LOCAL VARIABLES:
      INTEGER i, j

      IF ( selectAddFluid.NE.0 ) THEN
       IF ( ( selectAddFluid.GE.1 .AND. nonlinFreeSurf.GT.0 )
     &      .OR. convertFW2Salt.EQ.-1. _d 0 ) THEN
        DO j=0,sNy+1
         DO i=0,sNx+1
          gPtracer(i,j) = gPtracer(i,j)
     &      + ( addMass3Dplume(i,j,k,bi,bj)
     &         *( ptr_AddMass3D(i,j,k,bi,bj,iTracer)
     &           -pTracer(i,j,k,bi,bj,iTracer) )
     &         -addMass3Dbg(i,j,k,bi,bj)*pTracer(i,j,k,bi,bj,iTracer)
     &        )*mass2rUnit
     &       *recip_rA(i,j,bi,bj)
     &       *recip_drF(k)*_recip_hFacC(i,j,k,bi,bj)
         ENDDO
        ENDDO
       ELSE
        DO j=0,sNy+1
         DO i=0,sNx+1
          gPtracer(i,j) = gPtracer(i,j)
     &      + ( addMass3Dplume(i,j,k,bi,bj)
     &         *( ptr_AddMass3D(i,j,k,bi,bj,iTracer)
     &           -PTRACERS_ref(k,iTracer) )
     &         -addMass3Dbg(i,j,k,bi,bj)*PTRACERS_ref(k,iTracer)
     &        )*mass2rUnit
     &       *recip_rA(i,j,bi,bj)
     &       *recip_drF(k)*_recip_hFacC(i,j,k,bi,bj)
         ENDDO
        ENDDO
       ENDIF
      ENDIF
#endif /* ALLOW_ADDFLUID */
#endif /* ALLOW_ICEPLUME and ALLOW_PTRACERS */

      RETURN
      END
