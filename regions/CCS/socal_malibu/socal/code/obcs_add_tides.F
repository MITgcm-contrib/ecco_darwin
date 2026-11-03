#include "OBCS_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: OBCS_ADD_TIDES

C     !INTERFACE:
      SUBROUTINE OBCS_ADD_TIDES( myTime, myIter, myThid )

C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE OBCS_ADD_TIDES
C     | o Modify OB normal flow to add tidal forcing
C     | socal 2-km nest (2026-10-10): tides ramped up over the
C     | first tideRampTime = 2 d, tideRamp = (1-cos(pi t/T))/2,
C     | t = myTime - startTime. Switched on at full amplitude, the
C     | start-up transport of each constituent left a +1.2 m
C     | domain-mean ETAN offset (tides are added after
C     | OBCS_BALANCE_FLOW) and a w-CFL of 11 at the first step.
C     | 2026-10-15: ramp on cold start only, t = myTime - baseTime
C     | (was startTime: every pickup restart re-ramped the tides
C     | and left a lasting domain-mean ETAN step, malibu_v01).
C     |   NOTE that at the moment tidal forcing is applied
C     |   only to "normal" flow.  Code below should eventually
C     |   be augmented to also specify flow parallel to boundary.
C     *==========================================================*

C     !USES:
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "OBCS_PARAMS.h"
#include "OBCS_GRID.h"
#include "OBCS_FIELDS.h"

C     !INPUT/OUTPUT PARAMETERS:
      _RL myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef ALLOW_OBCS
#ifdef ALLOW_OBCS_TIDES

C     !FUNCTIONS:

C     !LOCAL VARIABLES:
C     bi, bj       :: tile indices
C     i,j,k        :: loop indices
C     iB, jB       :: local index of open boundary
C     msgBuf       :: Informational/error message buffer
      INTEGER bi, bj
      INTEGER i, j, k, iB, jB
      INTEGER td
      _RL tideRamp, tideRampTime, tElapsed

#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_ENTER('OBCS_ADD_TIDES',myThid)
#endif

      tideRampTime = 2. _d 0 * 86400. _d 0
      tElapsed = MAX( myTime - baseTime, 0. _d 0 )
      tideRamp = 1. _d 0
      IF ( tElapsed .LT. tideRampTime ) tideRamp =
     &     0.5 _d 0 * ( 1. _d 0 - COS( PI * tElapsed / tideRampTime ) )

C--   Add tidal currents:
      DO td=1,tidalComponents
       IF ( tidalPeriod(td) .GT. 0.D0 ) THEN
        DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)

#ifdef ALLOW_OBCS_EAST
          IF ( tileHasOBE(bi,bj) ) THEN
           DO k=1,Nr
            DO j=1-OLy,sNy+OLy
             iB = OB_Ie(j,bi,bj)
             IF ( iB.NE.OB_indexNone ) THEN
              OBEu(j,k,bi,bj) = OBEu(j,k,bi,bj) + tideRamp *
     &              maskW(iB,j,k,bi,bj) * OBEam(j,td,bi,bj) *
     &              COS( 2.D0 * PI * (myTime-OBEph(j,td,bi,bj)) /
     &                   tidalPeriod(td) )
             ENDIF
            ENDDO
           ENDDO
          ENDIF
#endif /* ALLOW_OBCS_EAST */

#ifdef ALLOW_OBCS_WEST
          IF ( tileHasOBW(bi,bj) ) THEN
           DO k=1,Nr
            DO j=1-OLy,sNy+OLy
             iB = OB_Iw(j,bi,bj)
             IF ( iB.NE.OB_indexNone ) THEN
              OBWu(j,k,bi,bj) = OBWu(j,k,bi,bj) + tideRamp *
     &              maskW(1+iB,j,k,bi,bj) * OBWam(j,td,bi,bj) *
     &              COS( 2.D0 * PI * (myTime-OBWph(j,td,bi,bj)) /
     &                   tidalPeriod(td) )
             ENDIF
            ENDDO
           ENDDO
          ENDIF
#endif /* ALLOW_OBCS_WEST */

#ifdef ALLOW_OBCS_NORTH
          IF ( tileHasOBN(bi,bj) ) THEN
           DO k=1,Nr
            DO i=1-OLx,sNx+OLx
             jB = OB_Jn(i,bi,bj)
             IF ( jB.NE.OB_indexNone ) THEN
              OBNv(i,k,bi,bj) = OBNv(i,k,bi,bj) + tideRamp *
     &              maskS(i,jB,k,bi,bj) * OBNam(i,td,bi,bj) *
     &              COS( 2.D0 * PI * (myTime-OBNph(i,td,bi,bj)) /
     &                   tidalPeriod(td) )
             ENDIF
            ENDDO
           ENDDO
          ENDIF
#endif /* ALLOW_OBCS_NORTH */

#ifdef ALLOW_OBCS_SOUTH
          IF ( tileHasOBS(bi,bj) ) THEN
           DO k=1,Nr
            DO i=1-OLx,sNx+OLx
             jB = OB_Js(i,bi,bj)
             IF ( jB.NE.OB_indexNone ) THEN
              OBSv(i,k,bi,bj) = OBSv(i,k,bi,bj) + tideRamp *
     &              maskS(i,1+jB,k,bi,bj)* OBSam(i,td,bi,bj) *
     &              COS( 2.D0 * PI * (myTime-OBSph(i,td,bi,bj)) /
     &                   tidalPeriod(td) )
             ENDIF
            ENDDO
           ENDDO
          ENDIF
#endif /* ALLOW_OBCS_SOUTH */

         ENDDO
        ENDDO
       ENDIF
      ENDDO

#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_LEAVE('OBCS_ADD_TIDES',myThid)
#endif

#endif /* ALLOW_OBCS_TIDES */
#endif /* ALLOW_OBCS */

      RETURN
      END
