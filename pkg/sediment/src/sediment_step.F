#include "SEDIMENT_OPTIONS.h"
#ifdef ALLOW_PTRACERS
# include "PTRACERS_OPTIONS.h"
#endif

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: SEDIMENT_STEP

C !INTERFACE:
      SUBROUTINE SEDIMENT_STEP( myTime, myIter, myThid )

C !DESCRIPTION:
C     Sediment source, settling, deposition and erosion, applied to the
C     new pTracer concentrations (kg/m^3) after the tracer time step
C     (operator split; called in FORWARD_STEP before the blocking
C     exchanges). For each class n and water column, with the layer
C     thickness the new tracers refer to, h_k = drF(k)*h0FacC*rStarFacC
C     under r* with staggerTimeStep (hFacC is only updated to it at the
C     next step), else drF(k)*hFacC, and mass m_k = C_k*h_k [kg/m^2]
C     (select_rStar=0 is stopped in SEDIMENT_CHECK):
C     1) the added fluid (addMass, e.g. a river) carries sedRiverC(n):
C        applied as a tracer tendency by SEDIMENT_ADD_TENDENCY (as for T
C        and S); here only the nominal input addMass/rhoConst*sedRiverC
C        *dt is accumulated for the budget monitor
C     2) settling: F_k = min(sedWs*dt*C_k, m_k) leaves each layer
C        downward (explicit, limited so no layer goes negative, which
C        keeps it stable in thin r* layers); from the bottom wet layer
C        the fraction p = 1 (or 1 - tau_b/sedTauCD if sedTauCD > 0)
C        deposits on the bed, the rest stays in the bottom layer
C        (deposition ramped to 0 in drying columns, between
C        2*sedDryDepth and sedDryDepth, unless sedDryDeposit)
C     3) erosion: E = sedEroM*(tau_b/sedTauCE - 1)*dt when
C        tau_b > sedTauCE and the column is deeper than sedDryDepth,
C        limited to the bed mass, into the lowest sedRefDepth of the
C        column (thickness-weighted; bottom wet layer if sedRefDepth<=0)
C     tau_b = rhoConst*(Cd_q*|u_b|^2 + Cd_l*|u_b|), u_b the thickness-
C     weighted mean velocity of the lowest sedRefDepth of the column
C     (cell-centred: mean of the two faces, closed faces as 0). The
C     bottom wet cell alone can be a thin partial cell whose velocity is
C     damped by the drag: its stress
C     is biased low where the bed crosses a layer interface, and those
C     cells become artificial deposition traps.
C     This routine conserves suspended + bed mass in each column; the
C     largest departure (round-off) is reported by the %SEDIMENT_MON
C     monitor, which also reports the nominal river input.

C !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "FFIELDS.h"
#include "SURFACE.h"
#include "SEDIMENT.h"
#ifdef ALLOW_MANGROVE
# include "MANGROVE.h"
#endif
#ifdef ALLOW_PTRACERS
# include "PTRACERS_SIZE.h"
# include "PTRACERS_PARAMS.h"
# include "PTRACERS_FIELDS.h"
#endif

C !INPUT PARAMETERS:
C     myTime :: time at the beginning of the time step
C     myIter :: iteration number at the beginning of the time step
C     myThid :: my Thread Id number
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#if (defined ALLOW_SEDIMENT) && (defined ALLOW_PTRACERS)
C !LOCAL VARIABLES:
      INTEGER i, j, k, n, bi, bj, kb, ip
      _RL dtS, ub, vb, spd2, tauB, colD, p, fl, ero, dep
      _RL acc, wk, eroFl
      _RL tot1, tauW, tauM, Aw, fw, cph, sph, spd
      _RL hk(Nr), m0(Nr), wr(Nr)
      _RL mmC(Nr,sedMaxCl), tot0N(sedMaxCl), rivN(sedMaxCl)
      _RL depN(sedMaxCl), eroN(sedMaxCl)
      _RL rivT(nSx,nSy,sedMaxCl), depT(nSx,nSy,sedMaxCl)
      _RL eroT(nSx,nSy,sedMaxCl), errT
      _RL gSum, errMax
      _RL depF(1-OLx:sNx+OLx,1-OLy:sNy+OLy,sedMaxCl)
      _RL eroF(1-OLx:sNx+OLx,1-OLy:sNy+OLy,sedMaxCl)
      _RL tmp2d(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      LOGICAL doRiver
#ifdef ALLOW_DIAGNOSTICS
      CHARACTER*8 diagName
      LOGICAL  DIAGNOSTICS_IS_ON
      EXTERNAL DIAGNOSTICS_IS_ON
#endif
      LOGICAL  DIFFERENT_MULTIPLE
      EXTERNAL DIFFERENT_MULTIPLE

      dtS = dTtracerLev(1)
C--   wave orbital velocity, period and direction of each column
      IF ( sedWaveMode.GT.0 )
     &  CALL SEDIMENT_WAVES( myTime, myIter, myThid )
      doRiver = selectAddFluid.GE.1
      errT = 0. _d 0

      DO bj=myByLo(myThid), myByHi(myThid)
       DO bi=myBxLo(myThid), myBxHi(myThid)
        DO n=1,sedMaxCl
         rivT(bi,bj,n) = 0. _d 0
         depT(bi,bj,n) = 0. _d 0
         eroT(bi,bj,n) = 0. _d 0
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           depF(i,j,n) = 0. _d 0
           eroF(i,j,n) = 0. _d 0
          ENDDO
         ENDDO
        ENDDO

        DO j=1,sNy
         DO i=1,sNx
          kb = kLowC(i,j,bi,bj)
          IF ( kb.GE.1 .AND. maskInC(i,j,bi,bj).NE.zeroRS ) THEN
           colD = 0. _d 0
           DO k=1,kb
            hk(k) = drF(k)*_hFacC(i,j,k,bi,bj)
#ifdef NONLIN_FRSURF
C-    r* with staggerTimeStep: the new tracers refer to the thickness of
C     the latest rStarFacC (CALC_R_STAR), not yet copied to hFacC
C     (UPDATE_R_STAR, next step). Without staggerTimeStep THERMODYNAMICS
C     runs before CALC_R_STAR and hFacC is the right thickness.
            IF ( nonlinFreeSurf.GT.0 .AND. select_rStar.GT.0
     &           .AND. staggerTimeStep )
     &        hk(k) = drF(k)*h0FacC(i,j,k,bi,bj)*rStarFacC(i,j,bi,bj)
#endif
            colD  = colD + hk(k)
           ENDDO
C--   bed shear stress from the mean velocity of the lowest sedRefDepth
C     of the column (open faces only), or of the bottom wet level
           ub = 0. _d 0
           vb = 0. _d 0
           acc = 0. _d 0
           DO k=kb,1,-1
            wr(k) = 0. _d 0
            wk = hk(k)
            IF ( sedRefDepth.GT.zeroRL )
     &        wk = MIN( hk(k), MAX( zeroRL, sedRefDepth - acc ) )
            IF ( k.LT.kb .AND. sedRefDepth.LE.zeroRL ) wk = 0. _d 0
            IF ( wk.GT.zeroRL ) THEN
C-    cell-centre velocity = mean of the two faces, a closed face
C     counting as u = 0 (dividing by the number of open faces instead
C     doubles the speed in dead-end cells)
             ub = ub + wk*halfRL*
     &          ( uVel(i,  j,k,bi,bj)*_maskW(i,  j,k,bi,bj)
     &          + uVel(i+1,j,k,bi,bj)*_maskW(i+1,j,k,bi,bj) )
             vb = vb + wk*halfRL*
     &          ( vVel(i,j,  k,bi,bj)*_maskS(i,j,  k,bi,bj)
     &          + vVel(i,j+1,k,bi,bj)*_maskS(i,j+1,k,bi,bj) )
             acc = acc + wk
             wr(k) = wk
            ENDIF
           ENDDO
           IF ( acc.GT.zeroRL ) THEN
            ub = ub/acc
            vb = vb/acc
            DO k=1,kb
             wr(k) = wr(k)/acc
            ENDDO
           ELSE
            wr(kb) = oneRL
           ENDIF
           spd2 = ub*ub + vb*vb
           tauB = rhoConst*( bottomDragQuadratic*spd2
     &                     + bottomDragLinear*SQRT(spd2) )
C-    waves: maximum wave + current stress (Soulsby 1995), the wave
C     friction factor from the orbital excursion (Soulsby 1997)
           tauW = 0. _d 0
           IF ( sedWaveMode.GT.0 .AND. sedUbw(i,j,bi,bj).GT.zeroRL
     &          .AND. sedPerW(i,j,bi,bj).GT.zeroRL ) THEN
            Aw   = sedUbw(i,j,bi,bj)*sedPerW(i,j,bi,bj)/(2. _d 0*PI)
            fw   = MIN( 0.3 _d 0, 1.39 _d 0*(Aw/sedZ0)**(-0.52 _d 0) )
            tauW = halfRL*rhoConst*fw*sedUbw(i,j,bi,bj)**2
            tauM = tauB*( oneRL + 1.2 _d 0
     &             *( tauW/(tauB+tauW) )**3.2 _d 0 )
C-    angle between waves and current (collinear without a direction)
            cph = oneRL
            sph = 0. _d 0
            spd = SQRT(spd2)
            IF ( sedDirW(i,j,bi,bj).NE.UNSET_RL
     &           .AND. spd.GT.zeroRL ) THEN
             cph = ( ub*COS(sedDirW(i,j,bi,bj))
     &             + vb*SIN(sedDirW(i,j,bi,bj)) )/spd
             cph = MAX( -oneRL, MIN( oneRL, cph ) )
             sph = SQRT( oneRL - cph*cph )
            ENDIF
            tauB = SQRT( ( tauM + tauW*ABS(cph) )**2 + (tauW*sph)**2 )
           ENDIF
#ifdef ALLOW_MANGROVE
C-    mangrove roots shelter part of the bed
           IF ( useMANGROVE ) tauB = tauB
     &        *MAX( zeroRL, oneRL - mgBedShelter*mgCovTot(i,j,bi,bj) )
#endif
           sedTauB(i,j,bi,bj) = tauB
           sedTauW(i,j,bi,bj) = tauW

C--   A) each class: river input (budget), settling and deposition
           DO n=1,sedNcl
            ip = sedPtrIdx(n)
            tot0N(n) = sedBed(i,j,n,bi,bj)
            DO k=1,kb
             m0(k) = pTracer(i,j,k,bi,bj,ip)*hk(k)
             mmC(k,n) = m0(k)
             tot0N(n) = tot0N(n) + m0(k)
            ENDDO
C-    1) nominal input carried by the added fluid (for the budget)
            rivN(n) = 0. _d 0
#ifdef ALLOW_ADDFLUID
            IF ( doRiver ) THEN
             DO k=1,kb
              IF ( addMass(i,j,k,bi,bj).GT.zeroRL ) rivN(n) = rivN(n)
     &          + addMass(i,j,k,bi,bj)*recip_rhoConst
     &           *recip_rA(i,j,bi,bj)*sedRivC(i,j,bi,bj,n)*dtS
             ENDDO
            ENDIF
#endif
C-    2) settling (from the concentrations before this step), and
C        deposition from the bottom wet layer
            dep = 0. _d 0
            DO k=1,kb
             fl = 0. _d 0
             IF ( hk(k).GT.zeroRL .AND. m0(k).GT.zeroRL )
     &         fl = MIN( sedWs(n)*dtS*m0(k)/hk(k), m0(k) )
             mmC(k,n) = mmC(k,n) - fl
             IF ( k.LT.kb ) THEN
              mmC(k+1,n) = mmC(k+1,n) + fl
             ELSE
              p = oneRL
              IF ( sedTauCD(n).GT.zeroRL )
     &          p = MAX( zeroRL, oneRL - tauB/sedTauCD(n) )
C-    drying columns: no deposition from the film (ramp from 1 at
C     2*sedDryDepth to 0 at sedDryDepth). The film is a numerical
C     minimum depth; a real drained flat keeps mm of water, so its load
C     stays suspended and leaves with the flood (sedDryDeposit=.TRUE.:
C     deposit from the film as from any water)
              IF ( .NOT.sedDryDeposit .AND. sedDryDepth.GT.zeroRL )
     &          p = p*MIN( oneRL, MAX( zeroRL,
     &                     ( colD - sedDryDepth )/sedDryDepth ) )
              dep = p*fl
              mmC(k,n) = mmC(k,n) + fl - dep
             ENDIF
            ENDDO
            depN(n) = dep
            eroN(n) = 0. _d 0
           ENDDO

C--   B) the bed: deposition in, erosion out
           IF ( sedNlay.LE.0 ) THEN
C-    single bed store per class (with the optional fluff layer)
            DO n=1,sedNcl
             sedBed(i,j,n,bi,bj) = sedBed(i,j,n,bi,bj) + depN(n)
C-    fluff layer: deposits enter the fluff, which consolidates into
C     the bed (sedBed stays the total)
             IF ( sedConsTime.GT.zeroRL ) THEN
              sedFluff(i,j,n,bi,bj) = ( sedFluff(i,j,n,bi,bj)
     &            + depN(n) )*EXP( -dtS/sedConsTime )
             ENDIF
C-    fluff first, at its own critical stress and rate; the consolidated
C     bed only once the fluff of the column is gone
             ero = 0. _d 0
             eroFl = 0. _d 0
             IF ( sedConsTime.GT.zeroRL .AND. colD.GT.sedDryDepth
     &            .AND. tauB.GT.sedTauCEf(n) ) THEN
              eroFl = sedEroMf(n)*( tauB/sedTauCEf(n) - oneRL )*dtS
              eroFl = MIN( eroFl, sedFluff(i,j,n,bi,bj) )
              sedFluff(i,j,n,bi,bj) = sedFluff(i,j,n,bi,bj) - eroFl
             ENDIF
             IF ( colD.GT.sedDryDepth .AND. tauB.GT.sedTauCE(n)
     &            .AND. sedFluff(i,j,n,bi,bj).LE.zeroRL ) THEN
              ero = sedEroM(n)*( tauB/sedTauCE(n) - oneRL )*dtS
              ero = MIN( ero, sedBed(i,j,n,bi,bj) - eroFl )
             ENDIF
             eroN(n) = ero + eroFl
             sedBed(i,j,n,bi,bj) = sedBed(i,j,n,bi,bj) - eroN(n)
            ENDDO
           ELSE
C-    layered bed (Sherwood et al. 2018; Warner et al. 2008)
            CALL SEDIMENT_BED_COLUMN( i, j, bi, bj, tauB, colD, dtS,
     U                                depN, eroN, myThid )
           ENDIF

C--   C) erosion into the lowest sedRefDepth of the column, with the
C        same thickness weights as the bed stress: a thin bottom cell
C        (a pocket closed on one side where the bed crosses a layer
C        interface) would otherwise receive all of it, and fast-settling
C        classes drop straight back onto the bed; new concentrations,
C        and the conservation check
           DO n=1,sedNcl
            ip  = sedPtrIdx(n)
            ero = eroN(n)
            IF ( ero.GT.zeroRL ) THEN
             DO k=1,kb
              mmC(k,n) = mmC(k,n) + wr(k)*ero
             ENDDO
            ENDIF
            tot1 = sedBed(i,j,n,bi,bj)
            DO k=1,kb
             IF ( hk(k).GT.zeroRL )
     &         pTracer(i,j,k,bi,bj,ip) = mmC(k,n)/hk(k)
             tot1 = tot1 + mmC(k,n)
            ENDDO
            errT = MAX( errT,
     &                  ABS( tot1 - tot0N(n) )*rA(i,j,bi,bj) )
            rivT(bi,bj,n) = rivT(bi,bj,n) + rivN(n)*rA(i,j,bi,bj)
            depT(bi,bj,n) = depT(bi,bj,n) + depN(n)*rA(i,j,bi,bj)
            eroT(bi,bj,n) = eroT(bi,bj,n) + ero*rA(i,j,bi,bj)
            depF(i,j,n) = depN(n)/dtS
            eroF(i,j,n) = ero/dtS
           ENDDO
          ENDIF
         ENDDO
        ENDDO

#ifdef ALLOW_DIAGNOSTICS
        IF ( useDiagnostics ) THEN
         DO n=1,sedNcl
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            tmp2d(i,j) = sedBed(i,j,n,bi,bj)
           ENDDO
          ENDDO
          WRITE(diagName,'(A,I2.2)') 'SEDBED', n
          CALL DIAGNOSTICS_FILL( tmp2d, diagName,
     &                           0, 1, 2, bi, bj, myThid )
          WRITE(diagName,'(A,I2.2)') 'SEDFLF', n
          CALL DIAGNOSTICS_FILL( sedFluff(1-OLx,1-OLy,n,bi,bj),
     &                           diagName, 0, 1, 1, bi, bj, myThid )
          WRITE(diagName,'(A,I2.2)') 'SEDDEP', n
          CALL DIAGNOSTICS_FILL( depF(1-OLx,1-OLy,n), diagName,
     &                           0, 1, 2, bi, bj, myThid )
          WRITE(diagName,'(A,I2.2)') 'SEDERO', n
          CALL DIAGNOSTICS_FILL( eroF(1-OLx,1-OLy,n), diagName,
     &                           0, 1, 2, bi, bj, myThid )
         ENDDO
         CALL DIAGNOSTICS_FILL( sedTauB, 'SEDTAUB ',
     &                          0, 1, 1, bi, bj, myThid )
         CALL DIAGNOSTICS_FILL( sedTauW, 'SEDTAUW ',
     &                          0, 1, 1, bi, bj, myThid )
         CALL DIAGNOSTICS_FILL( sedUbw,  'SEDUBW  ',
     &                          0, 1, 1, bi, bj, myThid )
         CALL DIAGNOSTICS_FILL( sedHsW,  'SEDHSW  ',
     &                          0, 1, 1, bi, bj, myThid )
         IF ( sedNlay.GT.0 ) THEN
          CALL DIAGNOSTICS_FILL( sedZa, 'SEDZA   ',
     &                           0, 1, 1, bi, bj, myThid )
          CALL DIAGNOSTICS_FILL( sedTcb(1-OLx,1-OLy,1,bi,bj),
     &                   'SEDTCB  ', 0, 1, 2, bi, bj, myThid )
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            tmp2d(i,j) = 0. _d 0
           ENDDO
          ENDDO
          DO n=1,sedNcl
           DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
             tmp2d(i,j) = tmp2d(i,j) + sedLay(i,j,1,n,bi,bj)
            ENDDO
           ENDDO
          ENDDO
          CALL DIAGNOSTICS_FILL( tmp2d, 'SEDACT  ',
     &                           0, 1, 2, bi, bj, myThid )
         ENDIF
        ENDIF
#endif /* ALLOW_DIAGNOSTICS */

       ENDDO
      ENDDO

C--   Global budgets
      DO n=1,sedNcl
        CALL GLOBAL_SUM_TILE_RL( rivT(1,1,n), gSum, myThid )
        _BEGIN_MASTER(myThid)
        sedCumRiv(n) = sedCumRiv(n) + gSum
        _END_MASTER(myThid)
        CALL GLOBAL_SUM_TILE_RL( depT(1,1,n), gSum, myThid )
        _BEGIN_MASTER(myThid)
        sedCumDep(n) = sedCumDep(n) + gSum
        _END_MASTER(myThid)
        CALL GLOBAL_SUM_TILE_RL( eroT(1,1,n), gSum, myThid )
        _BEGIN_MASTER(myThid)
        sedCumEro(n) = sedCumEro(n) + gSum
        _END_MASTER(myThid)
      ENDDO
      errMax = errT
      _GLOBAL_MAX_RL( errMax, myThid )
      _BEGIN_MASTER(myThid)
      sedErrMax = MAX( sedErrMax, errMax )
      _END_MASTER(myThid)
      _BARRIER

C--   Monitor
      IF ( DIFFERENT_MULTIPLE( sedMonFreq, myTime+deltaTClock,
     &                         deltaTClock ) ) THEN
        CALL SEDIMENT_MONITOR( myTime+deltaTClock, myIter+1, myThid )
      ENDIF
#endif /* ALLOW_SEDIMENT and ALLOW_PTRACERS */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: SEDIMENT_MONITOR

C !INTERFACE:
      SUBROUTINE SEDIMENT_MONITOR( myTime, myIter, myThid )

C !DESCRIPTION:
C     Print, for each class, the suspended and bed mass in the interior
C     (maskInC) domain and the cumulative river input, deposition and
C     erosion [kg]. The open-boundary export is
C     cumRiver - (susp + bed - initial total).

C !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SURFACE.h"
#include "SEDIMENT.h"
#ifdef ALLOW_PTRACERS
# include "PTRACERS_SIZE.h"
# include "PTRACERS_PARAMS.h"
# include "PTRACERS_FIELDS.h"
#endif

C !INPUT PARAMETERS:
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#if (defined ALLOW_SEDIMENT) && (defined ALLOW_PTRACERS)
C !LOCAL VARIABLES:
      INTEGER i, j, k, n, bi, bj, ip
      _RL susT(nSx,nSy), bedT(nSx,nSy), gSus, gBed, hk
      CHARACTER*(MAX_LEN_MBUF) msgBuf

      DO n=1,sedNcl
       ip = sedPtrIdx(n)
       DO bj=myByLo(myThid), myByHi(myThid)
        DO bi=myBxLo(myThid), myBxHi(myThid)
         susT(bi,bj) = 0. _d 0
         bedT(bi,bj) = 0. _d 0
         DO j=1,sNy
          DO i=1,sNx
           bedT(bi,bj) = bedT(bi,bj) + sedBed(i,j,n,bi,bj)
     &                  *rA(i,j,bi,bj)*maskInC(i,j,bi,bj)
           DO k=1,Nr
            hk = drF(k)*_hFacC(i,j,k,bi,bj)
#ifdef NONLIN_FRSURF
            IF ( nonlinFreeSurf.GT.0 .AND. select_rStar.GT.0
     &           .AND. staggerTimeStep )
     &        hk = drF(k)*h0FacC(i,j,k,bi,bj)*rStarFacC(i,j,bi,bj)
#endif
            susT(bi,bj) = susT(bi,bj) + pTracer(i,j,k,bi,bj,ip)
     &                   *hk*rA(i,j,bi,bj)*maskInC(i,j,bi,bj)
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
       CALL GLOBAL_SUM_TILE_RL( susT, gSus, myThid )
       CALL GLOBAL_SUM_TILE_RL( bedT, gBed, myThid )
       _BEGIN_MASTER(myThid)
       WRITE(msgBuf,'(A,I10,A,I2,A,1P2E14.6)')
     &  '%SEDIMENT_MON: iter', myIter, ' class', n,
     &  ' susp,bed [kg]=', gSus, gBed
       CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
       WRITE(msgBuf,'(A,I10,A,I2,A,1P3E14.6)')
     &  '%SEDIMENT_MON: iter', myIter, ' class', n,
     &  ' cumRiver,cumDep,cumEro [kg]=',
     &  sedCumRiv(n), sedCumDep(n), sedCumEro(n)
       CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
       _END_MASTER(myThid)
      ENDDO
      _BEGIN_MASTER(myThid)
      WRITE(msgBuf,'(A,I10,A,1PE12.4)')
     &  '%SEDIMENT_MON: iter', myIter, ' max step budget error [kg]=',
     &  sedErrMax
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      _END_MASTER(myThid)
#endif /* ALLOW_SEDIMENT and ALLOW_PTRACERS */

      RETURN
      END
