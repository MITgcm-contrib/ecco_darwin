#include "SEDIMENT_OPTIONS.h"
#ifdef ALLOW_SEAICE
# include "SEAICE_OPTIONS.h"
#endif

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: SEDIMENT_WAVES

C !INTERFACE:
      SUBROUTINE SEDIMENT_WAVES( myTime, myIter, myThid )

C !DESCRIPTION:
C     Near-bed wave orbital velocity sedUbw [m/s], period sedPerW [s],
C     direction sedDirW [rad, toward, ccw from +x] and (mode 2)
C     significant height sedHsW [m] of each column, for the wave part of
C     the bed stress in SEDIMENT_STEP.
C     sedWaveMode = 1: from files (e.g. WAVEWATCH III output: UBR, the
C       rms bottom orbital velocity, and a mean or peak period; the
C       direction file is optional, without it waves and current are
C       taken as collinear). Record r holds time (r-1)*sedWaveDt,
C       cyclic over sedWaveNrec records, linear in time between them.
C     sedWaveMode = 2: parametric wind sea, fetch- and depth-limited
C       (Young and Verhagen 1996, Coastal Eng. 29), from sedWindSpeed,
C       sedWindDir and sedFetch with the local water depth h:
C        Hs = 0.241 U^2/g [tanh(A1) tanh(B1/tanh(A1))]^0.87,
C        A1 = 0.493 d^0.75, B1 = 3.13e-3 c^0.57,
C        fp = 0.133 g/U [tanh(A2) tanh(B2/tanh(A2))]^-0.37,
C        A2 = 0.331 d^1.01, B2 = 5.215e-4 c^0.73,
C       d = g h/U^2, c = g F/U^2; Hs <= sedWaveGamma*h (breaking);
C       u_b = pi Hrms/(T sinh(kh)), Hrms = Hs/sqrt(2), T = 1/fp,
C       kh from Guo (2002): kh = x (1 - exp(-x^1.25))^-0.4, x = w^2 h/g.
C     sedWaveMode = 3: Hs, period and (optional) direction from files
C       (e.g. WAVEWATCH III HTSGW, PERPW and DIRPW interpolated to the
C       model grid; direction converted to 'toward, ccw from +x'), with
C       the same records and time interpolation as mode 1; u_b from the
C       local depth as in mode 2 (breaking limit, linear theory).
C     With pkg/seaice and sedWaveIce, u_b is multiplied by (1 - AREA).

C !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SEDIMENT.h"
#ifdef ALLOW_SEAICE
# include "SEAICE_SIZE.h"
# include "SEAICE_PARAMS.h"
# include "SEAICE.h"
#endif

C !INPUT PARAMETERS:
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef ALLOW_SEDIMENT
C !LOCAL VARIABLES:
      INTEGER i, j, k, bi, bj, r0, r1
      _RL colD, U, dlt, chi, A1, B1, A2, B2, hs, fp, per, om, x, kh
      _RL wgt, tt, cx, cy

C--   modes 1, 3: bring the two records around myTime into memory
C     (sedWub0/1 hold the orbital velocity in mode 1, Hs in mode 3)
      wgt = 0. _d 0
      IF ( sedWaveMode.EQ.1 .OR. sedWaveMode.EQ.3 ) THEN
       r0 = 1
       r1 = 1
       IF ( sedWaveNrec.GT.1 .AND. sedWaveDt.GT.zeroRL ) THEN
        tt  = MAX( zeroRL, myTime )/sedWaveDt
        r0  = INT( tt )
        wgt = tt - r0
        r0  = MOD( r0, sedWaveNrec ) + 1
        r1  = MOD( r0, sedWaveNrec ) + 1
       ENDIF
       IF ( r0.NE.sedWrec ) THEN
        IF ( sedWaveMode.EQ.1 ) THEN
         CALL READ_REC_XY_RL( sedWaveUbFile, sedWub0, r0,
     &                        myIter, myThid )
         CALL READ_REC_XY_RL( sedWaveUbFile, sedWub1, r1,
     &                        myIter, myThid )
        ELSE
         CALL READ_REC_XY_RL( sedWaveHsFile, sedWub0, r0,
     &                        myIter, myThid )
         CALL READ_REC_XY_RL( sedWaveHsFile, sedWub1, r1,
     &                        myIter, myThid )
        ENDIF
        CALL READ_REC_XY_RL( sedWavePerFile, sedWpe0, r0,
     &                       myIter, myThid )
        CALL READ_REC_XY_RL( sedWavePerFile, sedWpe1, r1,
     &                       myIter, myThid )
        IF ( sedWaveDirFile.NE.' ' ) THEN
         CALL READ_REC_XY_RL( sedWaveDirFile, sedWdi0, r0,
     &                        myIter, myThid )
         CALL READ_REC_XY_RL( sedWaveDirFile, sedWdi1, r1,
     &                        myIter, myThid )
        ENDIF
        _BARRIER
        _BEGIN_MASTER(myThid)
        sedWrec = r0
        _END_MASTER(myThid)
        _BARRIER
       ENDIF
      ENDIF

      DO bj=myByLo(myThid), myByHi(myThid)
       DO bi=myBxLo(myThid), myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          sedUbw (i,j,bi,bj) = 0. _d 0
          sedPerW(i,j,bi,bj) = 0. _d 0
          sedDirW(i,j,bi,bj) = 0. _d 0
          sedHsW (i,j,bi,bj) = 0. _d 0
         ENDDO
        ENDDO
        DO j=1,sNy
         DO i=1,sNx
          IF ( maskInC(i,j,bi,bj).NE.zeroRS ) THEN

           hs  = 0. _d 0
           per = 0. _d 0
           IF ( sedWaveMode.EQ.1 .OR. sedWaveMode.EQ.3 ) THEN
C-    from files, linear in time
            hs  = MAX( zeroRL, (oneRL-wgt)*sedWub0(i,j,bi,bj)
     &                        + wgt*sedWub1(i,j,bi,bj) )
            per = MAX( zeroRL, (oneRL-wgt)*sedWpe0(i,j,bi,bj)
     &                        + wgt*sedWpe1(i,j,bi,bj) )
            IF ( sedWaveDirFile.NE.' ' ) THEN
C-    interpolate the direction as a unit vector (no 360 deg jump)
             cx = (oneRL-wgt)*COS( deg2rad*sedWdi0(i,j,bi,bj) )
     &          + wgt*COS( deg2rad*sedWdi1(i,j,bi,bj) )
             cy = (oneRL-wgt)*SIN( deg2rad*sedWdi0(i,j,bi,bj) )
     &          + wgt*SIN( deg2rad*sedWdi1(i,j,bi,bj) )
             sedDirW(i,j,bi,bj) = ATAN2( cy, cx )
            ELSE
C-    no direction: flag for collinear waves and current
             sedDirW(i,j,bi,bj) = UNSET_RL
            ENDIF
            sedPerW(i,j,bi,bj) = per
            IF ( sedWaveMode.EQ.1 ) THEN
             sedUbw(i,j,bi,bj) = hs
             hs = 0. _d 0
            ENDIF
           ENDIF

           IF ( sedWaveMode.GE.2 ) THEN
            colD = 0. _d 0
            DO k=1,Nr
             colD = colD + drF(k)*_hFacC(i,j,k,bi,bj)
            ENDDO
            IF ( sedWaveMode.EQ.2 .AND. sedWindSpeed.GT.zeroRL
     &           .AND. colD.GT.zeroRL ) THEN
C-    parametric wind sea
             U   = sedWindSpeed
             dlt = gravity*colD/(U*U)
             chi = gravity*sedFetch/(U*U)
             A1  = 0.493 _d 0*dlt**0.75 _d 0
             B1  = 3.13 _d -3*chi**0.57 _d 0
             hs  = 0.241 _d 0*U*U/gravity
     &           *( TANH(A1)*TANH(B1/TANH(A1)) )**0.87 _d 0
             A2  = 0.331 _d 0*dlt**1.01 _d 0
             B2  = 5.215 _d -4*chi**0.73 _d 0
             fp  = 0.133 _d 0*gravity/U
     &           *( TANH(A2)*TANH(B2/TANH(A2)) )**(-0.37 _d 0)
             per = oneRL/fp
             sedPerW(i,j,bi,bj) = per
             sedDirW(i,j,bi,bj) = deg2rad*sedWindDir
            ENDIF
C-    local depth: breaking limit, then linear theory
            IF ( colD.GT.zeroRL .AND. hs.GT.zeroRL
     &           .AND. per.GT.zeroRL ) THEN
             hs  = MIN( hs, sedWaveGamma*colD )
             om  = 2. _d 0*PI/per
             x   = om*om*colD/gravity
             kh  = x*( oneRL - EXP(-x**1.25 _d 0) )**(-0.4 _d 0)
#ifdef ALLOW_MANGROVE
C-    decay through a mangrove forest (Mendez and Losada 2004)
             IF ( useMANGROVE )
     &         CALL MANGROVE_WAVE_ATT( i, j, bi, bj, colD, kh, hs,
     &                                 myThid )
#endif
             sedHsW(i,j,bi,bj) = hs
             IF ( kh.LT.50. _d 0 ) sedUbw(i,j,bi,bj) =
     &          PI*hs/SQRT(2. _d 0)/( per*SINH(kh) )
            ENDIF
           ENDIF

#ifdef ALLOW_SEAICE
C-    waves damped under sea ice
           IF ( sedWaveIce .AND. useSEAICE ) sedUbw(i,j,bi,bj) =
     &        sedUbw(i,j,bi,bj)*MAX( zeroRL, oneRL - AREA(i,j,bi,bj) )
#endif
          ENDIF
         ENDDO
        ENDDO
       ENDDO
      ENDDO
#endif /* ALLOW_SEDIMENT */

      RETURN
      END
