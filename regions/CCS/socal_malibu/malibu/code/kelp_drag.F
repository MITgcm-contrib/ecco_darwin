#include "DARWIN_OPTIONS.h"

C malibu_v01 [KD1]: kelp canopy drag on horizontal momentum.
C
C  Body force per unit mass, after MACMODS (Frieder et al. 2022, Front. Mar.
C  Sci. 9:752951, Eq. 5):  F_D = -1/2 C_D a P |u| u,  with
C    a   = kd_sfa * Bmag       frond surface area per volume [m2 m-3 = m-1]
C          (Bmag = kelp biomass tracer iMagB [gDW m-3]; kd_sfa = 1/58 m2/gDW,
C          as kelp_sfa in data.kelpoptics, G-MACMODS dry_sa)
C    C_D = kd_Cd   (0.0148 in Frieder et al. 2022, citing Utter & Denny 1996)
C    P   = kd_proj (1/2: Frieder et al. 2022 use P = 1/2 I)
C  Bmag is averaged from the two tracer cells next to each U (V) point and
C  |u| includes the 4-point mean of the other velocity component.
C  The drag rate r = 1/2 C_D P a |u| [s-1] enters as r/(1 + r deltaTMom), equal
C  to the explicit term for r deltaTMom << 1 and bounded by 1/deltaTMom for
C  dense canopies, so the term inside Adams-Bashforth cannot overshoot.
C  Runtime switch and parameters in data.kelpdrag (namelist KELP_DRAG_PARM);
C  no file or kd_on = .FALSE. leaves the model unchanged.
C  dir = 1: U point (gArr = gU tendency), dir = 2: V point.

CBOP
C     !ROUTINE: KELP_DRAG_TENDENCY
C     !INTERFACE:
      SUBROUTINE KELP_DRAG_TENDENCY(
     U                 gArr,
     I                 dir, k, bi, bj,
     I                 myTime, myIter, myThid )

      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_FIELDS.h"
#ifdef ALLOW_RADTRANS
#include "RADTRANS_SIZE.h"
#endif
#include "DARWIN_SIZE.h"
#include "DARWIN_INDICES.h"

C     !INPUT/OUTPUT PARAMETERS:
      _RL     gArr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER dir, k, bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef DARWIN_ALLOW_MACROALGAE
C     !LOCAL VARIABLES:
      INTEGER i, j, iUnit, ios
      LOGICAL exst
      _RL a, spd, vbar, r
      CHARACTER*(MAX_LEN_MBUF) msgBuf
C     parameters, read once per process
      LOGICAL kd_init, kd_on
      _RL kd_Cd, kd_proj, kd_sfa
      SAVE kd_init, kd_on, kd_Cd, kd_proj, kd_sfa
      DATA kd_init /.FALSE./
      NAMELIST /KELP_DRAG_PARM/ kd_on, kd_Cd, kd_proj, kd_sfa

      IF ( .NOT. kd_init ) THEN
       kd_on   = .FALSE.
       kd_Cd   = 0.0148 _d 0
       kd_proj = 0.5 _d 0
       kd_sfa  = 1. _d 0 / 58. _d 0
       INQUIRE( FILE='data.kelpdrag', EXIST=exst )
       IF ( exst ) THEN
        CALL OPEN_COPY_DATA_FILE(
     I                   'data.kelpdrag', 'KELP_DRAG_TENDENCY',
     O                   iUnit,
     I                   myThid )
        READ(UNIT=iUnit, NML=KELP_DRAG_PARM, IOSTAT=ios)
        CLOSE(iUnit)
        IF ( ios.NE.0 ) THEN
         WRITE(msgBuf,'(A)')
     &    'KELP_DRAG_TENDENCY: error reading data.kelpdrag'
         CALL PRINT_ERROR( msgBuf, myThid )
         STOP 'ABNORMAL END: S/R KELP_DRAG_TENDENCY'
        ENDIF
       ENDIF
       WRITE(msgBuf,'(A,L2,3(A,1PE12.5))')
     &  'KELP_DRAG: kd_on=', kd_on, ' kd_Cd=', kd_Cd,
     &  ' kd_proj=', kd_proj, ' kd_sfa=', kd_sfa
       CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
       kd_init = .TRUE.
      ENDIF

      IF ( .NOT. kd_on ) RETURN

      IF ( dir.EQ.1 ) THEN
C--    U points: face between tracer cells (i-1,j) and (i,j)
       DO j=0,sNy+1
        DO i=1,sNx+1
         a = kd_sfa * 0.5 _d 0
     &     * ( MAX(pTracer(i-1,j,k,bi,bj,iMagB),0. _d 0)
     &        *maskC(i-1,j,k,bi,bj)
     &       + MAX(pTracer(i  ,j,k,bi,bj,iMagB),0. _d 0)
     &        *maskC(i  ,j,k,bi,bj) )
         IF ( a.GT.0. _d 0 ) THEN
          vbar = 0.25 _d 0 * ( vVel(i-1,j  ,k,bi,bj)
     &                       + vVel(i  ,j  ,k,bi,bj)
     &                       + vVel(i-1,j+1,k,bi,bj)
     &                       + vVel(i  ,j+1,k,bi,bj) )
          spd = SQRT( uVel(i,j,k,bi,bj)**2 + vbar**2 )
          r = 0.5 _d 0 * kd_Cd * kd_proj * a * spd
          gArr(i,j) = gArr(i,j)
     &      - r / ( 1. _d 0 + r*deltaTMom )
     &      * uVel(i,j,k,bi,bj) * maskW(i,j,k,bi,bj)
         ENDIF
        ENDDO
       ENDDO
      ELSE
C--    V points: face between tracer cells (i,j-1) and (i,j)
       DO j=1,sNy+1
        DO i=0,sNx+1
         a = kd_sfa * 0.5 _d 0
     &     * ( MAX(pTracer(i,j-1,k,bi,bj,iMagB),0. _d 0)
     &        *maskC(i,j-1,k,bi,bj)
     &       + MAX(pTracer(i,j  ,k,bi,bj,iMagB),0. _d 0)
     &        *maskC(i,j  ,k,bi,bj) )
         IF ( a.GT.0. _d 0 ) THEN
          vbar = 0.25 _d 0 * ( uVel(i  ,j-1,k,bi,bj)
     &                       + uVel(i+1,j-1,k,bi,bj)
     &                       + uVel(i  ,j  ,k,bi,bj)
     &                       + uVel(i+1,j  ,k,bi,bj) )
          spd = SQRT( vVel(i,j,k,bi,bj)**2 + vbar**2 )
          r = 0.5 _d 0 * kd_Cd * kd_proj * a * spd
          gArr(i,j) = gArr(i,j)
     &      - r / ( 1. _d 0 + r*deltaTMom )
     &      * vVel(i,j,k,bi,bj) * maskS(i,j,k,bi,bj)
         ENDIF
        ENDDO
       ENDDO
      ENDIF
#endif /* DARWIN_ALLOW_MACROALGAE */

      RETURN
      END
