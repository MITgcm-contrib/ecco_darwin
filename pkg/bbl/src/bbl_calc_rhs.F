#include "BBL_OPTIONS.h"

CBOP
C     !ROUTINE: BBL_CALC_RHS

C     !INTERFACE:
      SUBROUTINE BBL_CALC_RHS(
     I        myTime, myIter, myThid )

C     !DESCRIPTION:
C     Calculate tendency of tracers due to bottom boundary layer.
C
C     The lateral exchange is computed in flux form so that the result
C     does not depend on the order in which faces are visited, and
C     therefore not on the domain decomposition (tiling):
C      1) vertical detrainment, pointwise, over the whole tile + halo;
C      2) decide, for every face, whether it is active and which side
C         is the donor, from the state after step 1; replenish the
C         thickness of interior donor cells, then exchange bbl_eta;
C      3) compute all face volume transports from this frozen state,
C         limit them so that no cell gives away more than half of its
C         bbl volume, nor receives more than half of its residual
C         volume, then update all interior cells at once.
C     Faces between an interior and a halo cell are computed by both
C     neighbouring tiles from identical halo data, hence identically.
C     This requires OLx,OLy >= 2 (see BBL_CHECK).

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "BBL.h"
#ifdef ALLOW_PTRACERS
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "PTRACERS_FIELDS.h"
#include "BBL_PTR.h"
#endif /* ALLOW_PTRACERS */
#ifdef ALLOW_AUTODIFF_TAMC
# include "tamc.h"
#endif

C     !INPUT PARAMETERS:
C     myTime    :: Current time in simulation
C     myIter    :: Current time-step number
C     myThid    :: my Thread Id number
      _RL     myTime
      INTEGER myIter, myThid

C     !OUTPUT PARAMETERS:

C     !LOCAL VARIABLES:
C     bi,bj     :: Tile indices
C     i,j       :: Loop indices
C     id,jd     :: Donor cell indices of a face
C     ir,jr     :: Receiver cell indices of a face
C     kBot      :: k index of bottommost wet grid box
C     kLowC1    :: k index of bottommost (i-1,j) or (i,j-1) cell
C     kLowC2    :: k index of bottommost (i,j) cell
C     kl        :: k index at which to compare 2 cells
C     eta0      :: bbl_eta after vertical detrainment, before replenish
C     dirX      :: zonal face (i-1|i) donor: +1 = (i-1), -1 = (i), 0 = off
C     dirY      :: merid. face (j-1|j) donor: +1 = (j-1), -1 = (j), 0 = off
C     thk       :: thickness of bottommost wet grid cell
C     tLoc,sLoc :: T/S of bottommost wet grid cell
C     rhoLoc    :: in situ density of bottommost wet grid cell
C     rhoBBL    :: in situ density of bottom boundary layer
C     resThk    :: thk - bbl_eta
C     resTheta  :: Theta of residual volume resThk
C     resSalt   :: Salt of residual volume resThk
C     volX,volY :: volume transport through zonal/meridional faces
C     sumOut    :: total volume a cell gives away
C     sumIn     :: total volume a cell receives
C     inTheta   :: received theta content
C     inSalt    :: received salt content
C     bbl_rho1  :: density of the (i-1,j) or (i,j-1) cell
C     bbl_rho2  :: density of the (i,j) cell
C     deltaRho  :: density change
C     deltaDpt  :: depth change
C     isDonor   :: cell is a donor across at least one face
C     etaNew    :: updated bbl_eta
C     tmpVol    :: temporary volume
C     velX,velY :: BBL speed through zonal/meridional faces (m/s)
C     fLoc      :: |f| at a face, bounded below by bbl_fMin
C     uEnt      :: BBL speed used for entrainment (fastest active face)
C     gPrime    :: reduced gravity of the BBL over the cell (m/s^2)
C     riB       :: bulk Richardson number g'*bbl_eta/uEnt^2
C     wEnt      :: entrainment velocity (m/s)
C     dEta      :: entrained thickness
C     ramp      :: 0 to 1 as gPrime goes from 0 to bbl_gpTaper
C     fMix      :: fraction of the BBL dissolved into its cell
C     resetBBL  :: set the BBL to the bottom-cell properties
C     wDet      :: detrainment velocity used in pass 1 (m/s)
C     uOut      :: fastest outflow speed of a cell (m/s)
C     cellT/S   :: cell T/S after this step's BBL tendency
C     kBotDt    :: tracer time step of the bottommost wet level
C     diagFld   :: signed face field for diagnostics
      INTEGER bi, bj
      INTEGER i, j, id, jd, ir, jr, kBot, kLowC1, kLowC2, kl
      _RL     eta0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      INTEGER dirX (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      INTEGER dirY (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL     thk     (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     tLoc    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     sLoc    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     rhoLoc  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     rhoBBL  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     resThk  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     resTheta(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     resSalt (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     volX    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     volY    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     sumOut  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     sumIn   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     inTheta (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
#ifdef ALLOW_PTRACERS
C     iTr       :: passive tracer index
C     inPtr     :: received passive tracer content
C     cellP     :: passive tracer of the cell after this step's tendency
      INTEGER iTr
      _RL     inPtr (1-OLx:sNx+OLx,1-OLy:sNy+OLy,PTRACERS_num)
      _RL     cellP
#endif /* ALLOW_PTRACERS */
      _RL     inSalt  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     bbl_rho1, bbl_rho2, deltaRho, deltaDpt
      _RL     etaNew, tmpVol
      _RL     velX (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL     velY (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL     fLoc
      _RL     uEnt, gPrime, riB, dEta, cellT, cellS, kBotDt
      _RL     wEnt (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     ramp, fMix, wDet, uOut
      LOGICAL isDonor, resetBBL
#ifdef ALLOW_DIAGNOSTICS
      _RL     diagFld (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      LOGICAL doDiagTrp, doDiagVel
      LOGICAL  DIAGNOSTICS_IS_ON
      EXTERNAL DIAGNOSTICS_IS_ON
#endif /* ALLOW_DIAGNOSTICS */
#ifdef BBL_CHECK_BUDGET
C     tileSumT/S :: per-tile sum of BBL T/S tendency x volume x dt
C     tileAbsT/S :: per-tile sum of its absolute value (for scaling)
C     msgBuf     :: Informational/error message buffer
      INTEGER k
      _RL     vol, sumT, sumS, absT, absS
      _RL     tileSumT(nSx,nSy), tileSumS(nSx,nSy)
      _RL     tileAbsT(nSx,nSy), tileAbsS(nSx,nSy)
      CHARACTER*(MAX_LEN_MBUF) msgBuf
#endif /* BBL_CHECK_BUDGET */
CEOP

#ifdef ALLOW_AUTODIFF_TAMC
C     BBL state is overwritten below: store it for the adjoint
CADJ STORE bbl_eta   = comlev1, key = ikey_dynamics, kind = isbyte
CADJ STORE bbl_theta = comlev1, key = ikey_dynamics, kind = isbyte
CADJ STORE bbl_salt  = comlev1, key = ikey_dynamics, kind = isbyte
# ifdef ALLOW_PTRACERS
CADJ STORE bbl_ptr   = comlev1, key = ikey_dynamics, kind = isbyte
# endif
#endif /* ALLOW_AUTODIFF_TAMC */

C     detrainment velocity of pass 1: with bbl_selectDetrain = 1 a dense
C     BBL detrains only where it does not flow out (after the face
C     directions are known, in the replenish loop below)
      wDet = bbl_wvel
      IF ( bbl_selectDetrain .EQ. 1 ) wDet = 0. _d 0

C==== Pass 1: vertical detrainment, face directions, replenishment.
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)

#ifdef ALLOW_PTRACERS
        IF ( bbl_doPtracers ) THEN
         DO iTr=1,PTRACERS_numInUse
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            bbl_TendPtr(i,j,bi,bj,iTr) = 0. _d 0
           ENDDO
          ENDDO
         ENDDO
        ENDIF
#endif /* ALLOW_PTRACERS */

C     Initialize tendency terms, make local copy of
C     bottomost temperature, salinity, in-situ density
C     and in-situ BBL density.
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          bbl_TendTheta(i,j,bi,bj) = 0. _d 0
          bbl_TendSalt (i,j,bi,bj) = 0. _d 0
          kBot        = max(1,kLowC(i,j,bi,bj))
          tLoc(i,j)   = theta(i,j,kBot,bi,bj)
          sLoc(i,j)   = salt (i,j,kBot,bi,bj)
          thk(i,j)    = drF(kBot) * hFacC(i,j,kBot,bi,bj)
          rhoLoc(i,j) = rhoInSitu(i,j,kBot,bi,bj)
          IF ( kBot .EQ. Nr ) THEN
           rhoBBL(i,j) = bbl_rho_nr(i,j,bi,bj)
          ELSE
           rhoBBL(i,j) = rhoInSitu(i,j,kBot+1,bi,bj)
          ENDIF
         ENDDO
        ENDDO

C--   Vertical exchange between BBL and residual volume of bottommost
C     wet grid box. This operation does not change total tracer
C     quantity in bottommost wet grid box. Pointwise, so it is done
C     identically in the halo of this tile and the interior of its
C     neighbour.
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          kBot = kLowC(i,j,bi,bj)
          IF ( kBot .GT. 0 ) THEN
           IF ( bbl_gpTaper.GT.0. .OR. bbl_tauMix.GT.0. ) THEN
C     Continuous variant: detrain at bbl_wvel and dissolve a fraction
C     fMix of the BBL into its cell; fMix goes from 0 for
C     gPrime >= bbl_gpTaper to min(1,dt/bbl_tauMix) for gPrime <= 0.
C     BBL T/S are unchanged, so the cell content is conserved.
            kBotDt = dTtracerLev(kBot)
            gPrime = gravity*recip_rhoConst
     &             * ( rhoBBL(i,j) - rhoLoc(i,j) )
            IF ( bbl_gpTaper.GT.0. ) THEN
             ramp = MIN( 1. _d 0, MAX( 0. _d 0, gPrime/bbl_gpTaper ) )
            ELSEIF ( gPrime.GT.0. ) THEN
             ramp = 1. _d 0
            ELSE
             ramp = 0. _d 0
            ENDIF
            fMix = 1. _d 0
            IF ( bbl_tauMix.GT.0. )
     &           fMix = MIN( 1. _d 0, kBotDt/bbl_tauMix )
            fMix = ( 1. _d 0 - ramp )*fMix
            bbl_eta(i,j,bi,bj) = MAX( 0. _d 0,
     &           bbl_eta(i,j,bi,bj)*( 1. _d 0 - fMix )
     &           - wDet*kBotDt )
            resetBBL = bbl_eta(i,j,bi,bj) .EQ. 0. _d 0
           ELSE
C     If model density is lower than BBL, slowly diffuse upward.
            IF ( rhoLoc(i,j) .LT. rhoBBL(i,j) )
     &            bbl_eta(i,j,bi,bj) = max ( 0. _d 0 ,
     &            bbl_eta(i,j,bi,bj) - wDet * dTtracerLev(kBot) )
C     If model density is higher than BBL then mix instantly.
            resetBBL = rhoLoc(i,j) .GE. rhoBBL(i,j) .OR.
     &                 bbl_eta(i,j,bi,bj) .EQ. 0. _d 0
           ENDIF
           IF ( resetBBL ) THEN
            bbl_theta(i,j,bi,bj) = tLoc(i,j)
            bbl_salt (i,j,bi,bj) = sLoc(i,j)
            bbl_eta  (i,j,bi,bj) = 0. _d 0
#ifdef ALLOW_PTRACERS
            IF ( bbl_doPtracers ) THEN
             DO iTr=1,PTRACERS_numInUse
              bbl_ptr(i,j,bi,bj,iTr) = pTracer(i,j,kBot,bi,bj,iTr)
             ENDDO
            ENDIF
#endif /* ALLOW_PTRACERS */
           ENDIF
          ENDIF
          eta0(i,j,bi,bj) = bbl_eta(i,j,bi,bj)
          dirX(i,j,bi,bj) = 0
          dirY(i,j,bi,bj) = 0
          velX(i,j,bi,bj) = 0. _d 0
          velY(i,j,bi,bj) = 0. _d 0
         ENDDO
        ENDDO

C--   Zonal faces (i-1|i). Heavy BBL water that sits higher than light
C     BBL water flows towards it. Compare the bbl densities at the
C     higher pressure (highest possible density of given t,s);
C     bbl in situ density is stored in k > kLowC indices.
        DO j=1-OLy,sNy+OLy
         DO i=2-OLx,sNx+OLx
          kLowC1 = kLowC(i-1,j,bi,bj)
          kLowC2 = kLowC(i  ,j,bi,bj)
C     maskInW is 0 on faces beyond an open boundary (pkg/obcs), where
C     the halo neighbor is not a real cell; else it is the wet mask
          IF ((kLowC1.GT.0).AND.(kLowC2.GT.0)
     &        .AND.(maskInW(i,j,bi,bj).GT.0.)) THEN
           kl = MAX(kLowC1,kLowC2) + 1
C     parentheses make deltaDpt exactly antisymmetric (and exactly 0
C     for equal BBL tops), so mirror-image faces get opposite signs
           deltaDpt = ( R_low(i-1,j,bi,bj) + eta0(i-1,j,bi,bj) )
     &              - ( R_low(i  ,j,bi,bj) + eta0(i  ,j,bi,bj) )
           IF ( deltaDpt .GT. 0. ) THEN
            IF ( kl .GT. Nr ) THEN
             bbl_rho1 = bbl_rho_nr(i-1,j,bi,bj)
            ELSE
             bbl_rho1 = rhoInSitu(i-1,j,kl,bi,bj)
            ENDIF
            bbl_rho2 = rhoInSitu(i,j,kLowC2,bi,bj)
           ELSE
            bbl_rho1 = rhoInSitu(i-1,j,kLowC1,bi,bj)
            IF ( kl .GT. Nr ) THEN
             bbl_rho2 = bbl_rho_nr(i,j,bi,bj)
            ELSE
             bbl_rho2 = rhoInSitu(i,j,kl,bi,bj)
            ENDIF
           ENDIF
           deltaRho = bbl_rho2 - bbl_rho1
           IF ( (deltaRho*deltaDpt) .LT. 0. ) THEN
            IF ( bbl_rho1 .GT. bbl_rho2 ) THEN
             dirX(i,j,bi,bj) =  1
            ELSE
             dirX(i,j,bi,bj) = -1
            ENDIF
            velX(i,j,bi,bj) = bbl_hvel
            IF ( bbl_useNofSpeed ) THEN
C     downslope fraction of the Nof speed g'*slope/f, g' from the
C     density contrast and slope from the BBL tops of the two cells
             fLoc = MAX( bbl_fMin,
     &            ABS( fCori(i-1,j,bi,bj) + fCori(i,j,bi,bj) )*halfRL )
             velX(i,j,bi,bj) = MIN( bbl_hvel, bbl_nofCoeff
     &            * gravity*recip_rhoConst*ABS(deltaRho)
     &            * ABS(deltaDpt)*recip_dxC(i,j,bi,bj) / fLoc )
            ENDIF
           ENDIF
          ENDIF
         ENDDO
        ENDDO

C--   Meridional faces (j-1|j), same as above.
        DO j=2-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          kLowC1 = kLowC(i,j-1,bi,bj)
          kLowC2 = kLowC(i,j  ,bi,bj)
          IF ((kLowC1.GT.0).AND.(kLowC2.GT.0)
     &        .AND.(maskInS(i,j,bi,bj).GT.0.)) THEN
           kl = MAX(kLowC1,kLowC2) + 1
           deltaDpt = ( R_low(i,j-1,bi,bj) + eta0(i,j-1,bi,bj) )
     &              - ( R_low(i,j  ,bi,bj) + eta0(i,j  ,bi,bj) )
           IF ( deltaDpt .GT. 0. ) THEN
            IF ( kl .GT. Nr ) THEN
             bbl_rho1 = bbl_rho_nr(i,j-1,bi,bj)
            ELSE
             bbl_rho1 = rhoInSitu(i,j-1,kl,bi,bj)
            ENDIF
            bbl_rho2 = rhoInSitu(i,j,kLowC2,bi,bj)
           ELSE
            bbl_rho1 = rhoInSitu(i,j-1,kLowC1,bi,bj)
            IF ( kl .GT. Nr ) THEN
             bbl_rho2 = bbl_rho_nr(i,j,bi,bj)
            ELSE
             bbl_rho2 = rhoInSitu(i,j,kl,bi,bj)
            ENDIF
           ENDIF
           deltaRho = bbl_rho2 - bbl_rho1
           IF ( (deltaRho*deltaDpt) .LT. 0. ) THEN
            IF ( bbl_rho1 .GT. bbl_rho2 ) THEN
             dirY(i,j,bi,bj) =  1
            ELSE
             dirY(i,j,bi,bj) = -1
            ENDIF
            velY(i,j,bi,bj) = bbl_hvel
            IF ( bbl_useNofSpeed ) THEN
             fLoc = MAX( bbl_fMin,
     &            ABS( fCori(i,j-1,bi,bj) + fCori(i,j,bi,bj) )*halfRL )
             velY(i,j,bi,bj) = MIN( bbl_hvel, bbl_nofCoeff
     &            * gravity*recip_rhoConst*ABS(deltaRho)
     &            * ABS(deltaDpt)*recip_dyC(i,j,bi,bj) / fLoc )
            ENDIF
           ENDIF
          ENDIF
         ENDDO
        ENDDO
#if (defined ALLOW_AUTODIFF_TAMC) && (defined BBL_WITH_STABLE_ADJOINT)
C     stable adjoint (cf. GMREDI_WITH_STABLE_ADJOINT): no sensitivity
C     through the BBL speeds, which on steep slopes makes the adjoint
C     grow without bound; the adjoint keeps the transport of BBL water
C     by the speeds of the forward run
        CALL ZERO_ADJ_LOC( 1, velX(1-OLx,1-OLy,bi,bj), myThid )
        CALL ZERO_ADJ_LOC( 1, velY(1-OLx,1-OLy,bi,bj), myThid )
#endif

C--   Replenish thickness of interior donor cells, if needed.
        DO j=1,sNy
         DO i=1,sNx
          isDonor = dirX(i  ,j,bi,bj) .EQ. -1 .OR.
     &              dirX(i+1,j,bi,bj) .EQ.  1 .OR.
     &              dirY(i,j  ,bi,bj) .EQ. -1 .OR.
     &              dirY(i,j+1,bi,bj) .EQ.  1
          IF ( isDonor .AND.
     &         ( bbl_theta(i,j,bi,bj) .EQ. tLoc(i,j) ) .AND.
     &         ( bbl_salt (i,j,bi,bj) .EQ. sLoc(i,j) ) .AND.
     &         ( bbl_eta  (i,j,bi,bj) .LT. bbl_initEta ) )
     &         bbl_eta(i,j,bi,bj) = min ( bbl_initEta, thk(i,j) )
C     plume-like detrainment: a dense BBL detrains at
C     bbl_wvel*max(0, 1 - uOut/bbl_uDet), uOut the fastest outflow
C     speed of the cell, so a stagnant BBL detrains fully and a plume
C     flowing out faster than bbl_uDet keeps its water
          kBot = kLowC(i,j,bi,bj)
          IF ( bbl_selectDetrain.EQ.1 .AND. kBot.GT.0 .AND.
     &         rhoLoc(i,j).LT.rhoBBL(i,j) ) THEN
           uOut = 0. _d 0
           IF ( dirX(i  ,j,bi,bj).EQ.-1 )
     &          uOut = MAX( uOut, velX(i  ,j,bi,bj) )
           IF ( dirX(i+1,j,bi,bj).EQ. 1 )
     &          uOut = MAX( uOut, velX(i+1,j,bi,bj) )
           IF ( dirY(i,j  ,bi,bj).EQ.-1 )
     &          uOut = MAX( uOut, velY(i,j  ,bi,bj) )
           IF ( dirY(i,j+1,bi,bj).EQ. 1 )
     &          uOut = MAX( uOut, velY(i,j+1,bi,bj) )
           bbl_eta(i,j,bi,bj) = MAX( 0. _d 0, bbl_eta(i,j,bi,bj)
     &          - bbl_wvel*dTtracerLev(kBot)
     &          * MAX( 0. _d 0, 1. _d 0 - uOut/bbl_uDet ) )
          ENDIF
         ENDDO
        ENDDO

C--   end bi,bj loops.
       ENDDO
      ENDDO

C     Halo cells now need the replenished thickness of their own tile.
      CALL EXCH_XY_RL( bbl_eta, myThid )

C==== Pass 2a: face volume transports from the frozen state and the
C     limiter factors of interior cells. A cell may give away at most
C     half of its bbl volume and receive at most half of its residual
C     volume, summed over its faces. The factors of halo cells are taken
C     from the tile that owns them (EXCH below): computing them in the
C     halo would need the neighbours of halo cells, which are not valid
C     at cubed-sphere corners, and both tiles must use the same factors
C     on a shared face for the exchange to be conservative.
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)

        CALL BBL_BOTTOM_STATE( bi, bj,
     O       tLoc, sLoc, thk, resThk, resTheta, resSalt, rhoLoc, rhoBBL,
     I       myThid )
        CALL BBL_RAW_VOLUMES( bi, bj,
     I       dirX(1-OLx,1-OLy,bi,bj), dirY(1-OLx,1-OLy,bi,bj),
     I       velX(1-OLx,1-OLy,bi,bj), velY(1-OLx,1-OLy,bi,bj), resThk,
     O       volX, volY, sumOut, sumIn,
     I       myThid )

        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          bbl_sclOut(i,j,bi,bj) = 1. _d 0
          bbl_sclIn (i,j,bi,bj) = 1. _d 0
         ENDDO
        ENDDO
        DO j=1,sNy
         DO i=1,sNx
          tmpVol = bbl_eta(i,j,bi,bj) * rA(i,j,bi,bj) * 0.5 _d 0
          IF ( sumOut(i,j) .GT. tmpVol )
     &         bbl_sclOut(i,j,bi,bj) = tmpVol / sumOut(i,j)
          tmpVol = resThk(i,j) * rA(i,j,bi,bj) * 0.5 _d 0
          IF ( sumIn(i,j) .GT. tmpVol )
     &         bbl_sclIn(i,j,bi,bj) = tmpVol / sumIn(i,j)
         ENDDO
        ENDDO

C--   end bi,bj loops.
       ENDDO
      ENDDO

      CALL EXCH_XY_RL( bbl_sclOut, myThid )
      CALL EXCH_XY_RL( bbl_sclIn , myThid )

C==== Pass 2b: limited face volumes and simultaneous update of all
C     interior cells.
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)

        CALL BBL_BOTTOM_STATE( bi, bj,
     O       tLoc, sLoc, thk, resThk, resTheta, resSalt, rhoLoc, rhoBBL,
     I       myThid )
        CALL BBL_RAW_VOLUMES( bi, bj,
     I       dirX(1-OLx,1-OLy,bi,bj), dirY(1-OLx,1-OLy,bi,bj),
     I       velX(1-OLx,1-OLy,bi,bj), velY(1-OLx,1-OLy,bi,bj), resThk,
     O       volX, volY, sumOut, sumIn,
     I       myThid )
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          wEnt(i,j) = 0. _d 0
         ENDDO
        ENDDO


C--   Apply the limiter to the faces of interior cells.
        DO j=1,sNy
         DO i=1,sNx+1
          IF ( dirX(i,j,bi,bj) .EQ. 1 ) THEN
           volX(i,j) = volX(i,j) * min( bbl_sclOut(i-1,j,bi,bj),
     &                                      bbl_sclIn(i,j,bi,bj) )
          ELSEIF ( dirX(i,j,bi,bj) .EQ. -1 ) THEN
           volX(i,j) = volX(i,j) * min( bbl_sclOut(i,j,bi,bj),
     &                                      bbl_sclIn(i-1,j,bi,bj) )
          ENDIF
         ENDDO
        ENDDO
        DO j=1,sNy+1
         DO i=1,sNx
          IF ( dirY(i,j,bi,bj) .EQ. 1 ) THEN
           volY(i,j) = volY(i,j) * min( bbl_sclOut(i,j-1,bi,bj),
     &                                      bbl_sclIn(i,j,bi,bj) )
          ELSEIF ( dirY(i,j,bi,bj) .EQ. -1 ) THEN
           volY(i,j) = volY(i,j) * min( bbl_sclOut(i,j,bi,bj),
     &                                      bbl_sclIn(i,j-1,bi,bj) )
          ENDIF
         ENDDO
        ENDDO

#ifdef ALLOW_DIAGNOSTICS
C--   Signed face transports (m^3/s) and speeds (m/s), positive towards
C     increasing i or j.
        IF ( useDiagnostics ) THEN
         doDiagTrp = DIAGNOSTICS_IS_ON( 'BBLuTrp ', myThid )
     &          .OR. DIAGNOSTICS_IS_ON( 'BBLvTrp ', myThid )
         doDiagVel = DIAGNOSTICS_IS_ON( 'BBLuVel ', myThid )
     &          .OR. DIAGNOSTICS_IS_ON( 'BBLvVel ', myThid )
         IF ( doDiagTrp ) THEN
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            diagFld(i,j) = 0. _d 0
            IF ( dirX(i,j,bi,bj).EQ.1 ) diagFld(i,j) = volX(i,j)
     &           / dTtracerLev(kLowC(i-1,j,bi,bj))
            IF ( dirX(i,j,bi,bj).EQ.-1 ) diagFld(i,j) = -volX(i,j)
     &           / dTtracerLev(kLowC(i,j,bi,bj))
           ENDDO
          ENDDO
          CALL DIAGNOSTICS_FILL( diagFld, 'BBLuTrp ',
     &                           0, 1, 2, bi, bj, myThid )
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            diagFld(i,j) = 0. _d 0
            IF ( dirY(i,j,bi,bj).EQ.1 ) diagFld(i,j) = volY(i,j)
     &           / dTtracerLev(kLowC(i,j-1,bi,bj))
            IF ( dirY(i,j,bi,bj).EQ.-1 ) diagFld(i,j) = -volY(i,j)
     &           / dTtracerLev(kLowC(i,j,bi,bj))
           ENDDO
          ENDDO
          CALL DIAGNOSTICS_FILL( diagFld, 'BBLvTrp ',
     &                           0, 1, 2, bi, bj, myThid )
         ENDIF
         IF ( doDiagVel ) THEN
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            diagFld(i,j) = dirX(i,j,bi,bj)*velX(i,j,bi,bj)
           ENDDO
          ENDDO
          CALL DIAGNOSTICS_FILL( diagFld, 'BBLuVel ',
     &                           0, 1, 2, bi, bj, myThid )
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            diagFld(i,j) = dirY(i,j,bi,bj)*velY(i,j,bi,bj)
           ENDDO
          ENDDO
          CALL DIAGNOSTICS_FILL( diagFld, 'BBLvVel ',
     &                           0, 1, 2, bi, bj, myThid )
         ENDIF
        ENDIF
#endif /* ALLOW_DIAGNOSTICS */

C--   Tracer tendencies of the bottommost wet cell, and volume and
C     content that each interior cell gives away and receives.
C     Exchanging dVol of donor bbl water against dVol of receiver
C     residual water moves dVol*(bblTracer_d - resTracer_r).
        DO j=1,sNy
         DO i=1,sNx
          sumOut (i,j) = 0. _d 0
          sumIn  (i,j) = 0. _d 0
          inTheta(i,j) = 0. _d 0
          inSalt (i,j) = 0. _d 0
         ENDDO
        ENDDO
#ifdef ALLOW_PTRACERS
        IF ( bbl_doPtracers ) THEN
         DO iTr=1,PTRACERS_numInUse
          DO j=1,sNy
           DO i=1,sNx
            inPtr(i,j,iTr) = 0. _d 0
           ENDDO
          ENDDO
         ENDDO
        ENDIF
#endif /* ALLOW_PTRACERS */
        DO j=1,sNy
         DO i=1,sNx+1
          IF ( volX(i,j) .GT. 0. ) THEN
           IF ( dirX(i,j,bi,bj) .EQ. 1 ) THEN
            id = i-1
            ir = i
           ELSE
            id = i
            ir = i-1
           ENDIF
           CALL BBL_FACE_APPLY(
     I          volX(i,j), id, j, ir, j, bi, bj,
     I          thk, resTheta, resSalt,
     U          sumOut, sumIn, inTheta, inSalt,
     I          myThid )
#ifdef ALLOW_PTRACERS
           IF ( bbl_doPtracers ) CALL BBL_FACE_APPLY_PTR(
     I          volX(i,j), id, j, ir, j, bi, bj, thk, resThk,
     U          inPtr,
     I          myThid )
#endif /* ALLOW_PTRACERS */
          ENDIF
         ENDDO
        ENDDO
        DO j=1,sNy+1
         DO i=1,sNx
          IF ( volY(i,j) .GT. 0. ) THEN
           IF ( dirY(i,j,bi,bj) .EQ. 1 ) THEN
            jd = j-1
            jr = j
           ELSE
            jd = j
            jr = j-1
           ENDIF
           CALL BBL_FACE_APPLY(
     I          volY(i,j), i, jd, i, jr, bi, bj,
     I          thk, resTheta, resSalt,
     U          sumOut, sumIn, inTheta, inSalt,
     I          myThid )
#ifdef ALLOW_PTRACERS
           IF ( bbl_doPtracers ) CALL BBL_FACE_APPLY_PTR(
     I          volY(i,j), i, jd, i, jr, bi, bj, thk, resThk,
     U          inPtr,
     I          myThid )
#endif /* ALLOW_PTRACERS */
          ENDIF
         ENDDO
        ENDDO

C--   Adjust bbl thickness and tracer properties of interior cells.
        DO j=1,sNy
         DO i=1,sNx
          IF ( sumOut(i,j).GT.0. .OR. sumIn(i,j).GT.0. ) THEN
           tmpVol = bbl_eta(i,j,bi,bj)*rA(i,j,bi,bj) - sumOut(i,j)
           etaNew = ( tmpVol + sumIn(i,j) ) / rA(i,j,bi,bj)
           IF ( etaNew .LT. 0.0001 ) THEN
            bbl_eta  (i,j,bi,bj) = 0. _d 0
            bbl_theta(i,j,bi,bj) = tLoc(i,j)
            bbl_salt (i,j,bi,bj) = sLoc(i,j)
#ifdef ALLOW_PTRACERS
            IF ( bbl_doPtracers ) THEN
             kBot = max(1,kLowC(i,j,bi,bj))
             DO iTr=1,PTRACERS_numInUse
              bbl_ptr(i,j,bi,bj,iTr) = pTracer(i,j,kBot,bi,bj,iTr)
             ENDDO
            ENDIF
#endif /* ALLOW_PTRACERS */
           ELSE
            IF ( sumIn(i,j) .GT. 0. ) THEN
             bbl_theta(i,j,bi,bj) =
     &            ( tmpVol*bbl_theta(i,j,bi,bj) + inTheta(i,j) )
     &            / ( etaNew*rA(i,j,bi,bj) )
             bbl_salt (i,j,bi,bj) =
     &            ( tmpVol*bbl_salt (i,j,bi,bj) + inSalt (i,j) )
     &            / ( etaNew*rA(i,j,bi,bj) )
#ifdef ALLOW_PTRACERS
             IF ( bbl_doPtracers ) THEN
              DO iTr=1,PTRACERS_numInUse
               bbl_ptr(i,j,bi,bj,iTr) =
     &            ( tmpVol*bbl_ptr(i,j,bi,bj,iTr) + inPtr(i,j,iTr) )
     &            / ( etaNew*rA(i,j,bi,bj) )
              ENDDO
             ENDIF
#endif /* ALLOW_PTRACERS */
            ENDIF
            bbl_eta(i,j,bi,bj) = etaNew
           ENDIF
          ENDIF
         ENDDO
        ENDDO

C--   Entrainment of residual water into a moving BBL that is denser
C     than its cell: E(Ri) = (0.08-0.1*Ri)/(1+5*Ri) for Ri < 0.8
C     (Turner 1986), w_E = E*U. Mixes within the bottommost cell only,
C     so it does not change the cell's tracer content.
        IF ( bbl_selectEntrain .EQ. 1 ) THEN
         DO j=1,sNy
          DO i=1,sNx
           kBot = kLowC(i,j,bi,bj)
           IF ( kBot.GT.0 .AND. bbl_eta(i,j,bi,bj).GT.0. ) THEN
            uEnt = 0. _d 0
            IF ( dirX(i  ,j,bi,bj).NE.0 )
     &           uEnt = MAX( uEnt, velX(i  ,j,bi,bj) )
            IF ( dirX(i+1,j,bi,bj).NE.0 )
     &           uEnt = MAX( uEnt, velX(i+1,j,bi,bj) )
            IF ( dirY(i,j  ,bi,bj).NE.0 )
     &           uEnt = MAX( uEnt, velY(i,j  ,bi,bj) )
            IF ( dirY(i,j+1,bi,bj).NE.0 )
     &           uEnt = MAX( uEnt, velY(i,j+1,bi,bj) )
            gPrime = gravity*recip_rhoConst
     &             * ( rhoBBL(i,j) - rhoLoc(i,j) )
            IF ( uEnt.GT.0. .AND. gPrime.GT.0. ) THEN
             riB = gPrime*bbl_eta(i,j,bi,bj) / ( uEnt*uEnt )
             IF ( riB .LT. 0.8 _d 0 ) wEnt(i,j) = uEnt
     &            * ( 0.08 _d 0 - 0.1 _d 0*riB )
     &            / ( 1. _d 0 + 5. _d 0*riB )
C     fade out near neutral density (E is largest as gPrime -> 0)
             IF ( bbl_gpTaper.GT.0. ) wEnt(i,j) = wEnt(i,j)
     &            * MIN( 1. _d 0, gPrime/bbl_gpTaper )
            ENDIF
           ENDIF
          ENDDO
         ENDDO
#if (defined ALLOW_AUTODIFF_TAMC) && (defined BBL_WITH_STABLE_ADJOINT)
C     stable adjoint: no sensitivity through the entrainment speed
         CALL ZERO_ADJ_LOC( 1, wEnt, myThid )
#endif
         DO j=1,sNy
          DO i=1,sNx
           kBot = kLowC(i,j,bi,bj)
           IF ( kBot.GT.0 .AND. bbl_eta(i,j,bi,bj).GT.0. ) THEN
            kBotDt = dTtracerLev(kBot)
            dEta = MIN( wEnt(i,j)*kBotDt,
     &                  thk(i,j) - bbl_eta(i,j,bi,bj) )
            IF ( dEta .GT. 0. ) THEN
C     residual water after this step's BBL tendency of the cell
             cellT = tLoc(i,j) + bbl_TendTheta(i,j,bi,bj)*kBotDt
             cellS = sLoc(i,j) + bbl_TendSalt (i,j,bi,bj)*kBotDt
             tmpVol = thk(i,j) - bbl_eta(i,j,bi,bj)
             bbl_theta(i,j,bi,bj) = ( bbl_eta(i,j,bi,bj)
     &            * bbl_theta(i,j,bi,bj) + dEta
     &            * ( cellT*thk(i,j)
     &              - bbl_theta(i,j,bi,bj)*bbl_eta(i,j,bi,bj) )/tmpVol )
     &            / ( bbl_eta(i,j,bi,bj) + dEta )
             bbl_salt (i,j,bi,bj) = ( bbl_eta(i,j,bi,bj)
     &            * bbl_salt (i,j,bi,bj) + dEta
     &            * ( cellS*thk(i,j)
     &              - bbl_salt (i,j,bi,bj)*bbl_eta(i,j,bi,bj) )/tmpVol )
     &            / ( bbl_eta(i,j,bi,bj) + dEta )
#ifdef ALLOW_PTRACERS
             IF ( bbl_doPtracers ) THEN
              DO iTr=1,PTRACERS_numInUse
               cellP = pTracer(i,j,kBot,bi,bj,iTr)
     &               + bbl_TendPtr(i,j,bi,bj,iTr)*PTRACERS_dTLev(kBot)
               bbl_ptr(i,j,bi,bj,iTr) = ( bbl_eta(i,j,bi,bj)
     &            * bbl_ptr(i,j,bi,bj,iTr) + dEta
     &            * ( cellP*thk(i,j)
     &            - bbl_ptr(i,j,bi,bj,iTr)*bbl_eta(i,j,bi,bj) )/tmpVol )
     &            / ( bbl_eta(i,j,bi,bj) + dEta )
              ENDDO
             ENDIF
#endif /* ALLOW_PTRACERS */
             bbl_eta(i,j,bi,bj) = bbl_eta(i,j,bi,bj) + dEta
            ENDIF
           ENDIF
          ENDDO
         ENDDO
        ENDIF
#ifdef ALLOW_DIAGNOSTICS
        IF ( useDiagnostics ) CALL DIAGNOSTICS_FILL( wEnt, 'BBLwEnt ',
     &                                   0, 1, 2, bi, bj, myThid )
#endif /* ALLOW_DIAGNOSTICS */

C--   end bi,bj loops.
       ENDDO
      ENDDO

      CALL EXCH_XY_RL( bbl_eta      , myThid )
      CALL EXCH_XY_RL( bbl_theta    , myThid )
      CALL EXCH_XY_RL( bbl_salt     , myThid )
      CALL EXCH_XY_RL( bbl_TendTheta, myThid )
      CALL EXCH_XY_RL( bbl_TendSalt , myThid )
#ifdef ALLOW_PTRACERS
      IF ( bbl_doPtracers ) THEN
       DO iTr=1,PTRACERS_numInUse
        CALL EXCH_XY_RL( bbl_ptr    (1-OLx,1-OLy,1,1,iTr), myThid )
        CALL EXCH_XY_RL( bbl_TendPtr(1-OLx,1-OLy,1,1,iTr), myThid )
       ENDDO
      ENDIF
#endif /* ALLOW_PTRACERS */

#ifdef BBL_CHECK_BUDGET
C--   BBL lateral exchange only moves T/S between bottom cells, so the
C     volume-weighted sum of its tendencies over the domain should be 0.
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        tileSumT(bi,bj) = 0. _d 0
        tileSumS(bi,bj) = 0. _d 0
        tileAbsT(bi,bj) = 0. _d 0
        tileAbsS(bi,bj) = 0. _d 0
        DO j=1,sNy
         DO i=1,sNx
          k = kLowC(i,j,bi,bj)
          IF ( k .GT. 0 ) THEN
           vol = rA(i,j,bi,bj)*drF(k)*hFacC(i,j,k,bi,bj)
     &         * dTtracerLev(k)
           tileSumT(bi,bj) = tileSumT(bi,bj)
     &                     + bbl_TendTheta(i,j,bi,bj)*vol
           tileSumS(bi,bj) = tileSumS(bi,bj)
     &                     + bbl_TendSalt(i,j,bi,bj)*vol
           tileAbsT(bi,bj) = tileAbsT(bi,bj)
     &                     + ABS(bbl_TendTheta(i,j,bi,bj))*vol
           tileAbsS(bi,bj) = tileAbsS(bi,bj)
     &                     + ABS(bbl_TendSalt(i,j,bi,bj))*vol
          ENDIF
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      CALL GLOBAL_SUM_TILE_RL( tileSumT, sumT, myThid )
      CALL GLOBAL_SUM_TILE_RL( tileSumS, sumS, myThid )
      CALL GLOBAL_SUM_TILE_RL( tileAbsT, absT, myThid )
      CALL GLOBAL_SUM_TILE_RL( tileAbsS, absS, myThid )
      _BEGIN_MASTER( myThid )
      WRITE(msgBuf,'(A,I10,1P4E14.6)')
     &     'BBL_BUDGET it,sumT,absT,sumS,absS',
     &     myIter, sumT, absT, sumS, absS
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      _END_MASTER( myThid )
#ifdef ALLOW_PTRACERS
      IF ( bbl_doPtracers ) THEN
       DO iTr=1,PTRACERS_numInUse
        DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)
          tileSumT(bi,bj) = 0. _d 0
          tileAbsT(bi,bj) = 0. _d 0
          DO j=1,sNy
           DO i=1,sNx
            k = kLowC(i,j,bi,bj)
            IF ( k .GT. 0 ) THEN
             vol = rA(i,j,bi,bj)*drF(k)*hFacC(i,j,k,bi,bj)
     &           * PTRACERS_dTLev(k)
             tileSumT(bi,bj) = tileSumT(bi,bj)
     &                       + bbl_TendPtr(i,j,bi,bj,iTr)*vol
             tileAbsT(bi,bj) = tileAbsT(bi,bj)
     &                       + ABS(bbl_TendPtr(i,j,bi,bj,iTr))*vol
            ENDIF
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        CALL GLOBAL_SUM_TILE_RL( tileSumT, sumT, myThid )
        CALL GLOBAL_SUM_TILE_RL( tileAbsT, absT, myThid )
        _BEGIN_MASTER( myThid )
        WRITE(msgBuf,'(A,I10,I4,1P2E14.6)')
     &       'BBL_BUDGET_PTR it,iTr,sum,abs', myIter, iTr, sumT, absT
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
        _END_MASTER( myThid )
       ENDDO
      ENDIF
#endif /* ALLOW_PTRACERS */
#endif /* BBL_CHECK_BUDGET */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: BBL_FACE_APPLY

C     !INTERFACE:
      SUBROUTINE BBL_FACE_APPLY(
     I        dVol, id, jd, ir, jr, bi, bj,
     I        thk, resTheta, resSalt,
     U        sumOut, sumIn, inTheta, inSalt,
     I        myThid )

C     !DESCRIPTION:
C     Add the effect of a volume dVol of donor (id,jd) bbl water
C     exchanged against receiver (ir,jr) residual water: T/S tendency
C     of both bottommost wet cells, and the volume and T/S content
C     given away / received, for whichever of the two cells is in the
C     tile interior.

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "BBL.h"

C     !INPUT/OUTPUT PARAMETERS:
C     dVol      :: exchanged volume (m^3)
C     id,jd     :: donor cell indices
C     ir,jr     :: receiver cell indices
C     bi,bj     :: tile indices
C     thk       :: thickness of bottommost wet grid cell
C     resTheta  :: Theta of residual volume
C     resSalt   :: Salt of residual volume
C     sumOut    :: total volume a cell gives away
C     sumIn     :: total volume a cell receives
C     inTheta   :: received theta content
C     inSalt    :: received salt content
C     myThid    :: my Thread Id number
      _RL     dVol
      INTEGER id, jd, ir, jr, bi, bj
      _RL     thk     (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     resTheta(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     resSalt (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     sumOut  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     sumIn   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     inTheta (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     inSalt  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER myThid
CEOP

C     !LOCAL VARIABLES:
C     tendT,tendS :: T/S content moved per unit time
C     inD, inR    :: donor / receiver cell is in the tile interior
C     dtD, dtR    :: tracer time step at donor / receiver bottom level
      _RL     tendT, tendS, dtD, dtR
      LOGICAL inD, inR

      inD = id.GE.1 .AND. id.LE.sNx .AND. jd.GE.1 .AND. jd.LE.sNy
      inR = ir.GE.1 .AND. ir.LE.sNx .AND. jr.GE.1 .AND. jr.LE.sNy

C     The tendency of each cell is applied over the tracer time step of
C     its own level, so both get exactly dVol*(bblTracer_d-resTracer_r)
C     (conservative, and consistent with the bbl_eta update).
      dtD = dTtracerLev(kLowC(id,jd,bi,bj))
      dtR = dTtracerLev(kLowC(ir,jr,bi,bj))

      IF ( inD ) THEN
       tendT = dVol * (bbl_theta(id,jd,bi,bj) - resTheta(ir,jr)) / dtD
       tendS = dVol * (bbl_salt (id,jd,bi,bj) - resSalt (ir,jr)) / dtD
       bbl_TendTheta(id,jd,bi,bj) = bbl_TendTheta(id,jd,bi,bj)
     &      - tendT / rA(id,jd,bi,bj) / thk(id,jd)
       bbl_TendSalt (id,jd,bi,bj) = bbl_TendSalt (id,jd,bi,bj)
     &      - tendS / rA(id,jd,bi,bj) / thk(id,jd)
       sumOut(id,jd) = sumOut(id,jd) + dVol
      ENDIF
      IF ( inR ) THEN
       tendT = dVol * (bbl_theta(id,jd,bi,bj) - resTheta(ir,jr)) / dtR
       tendS = dVol * (bbl_salt (id,jd,bi,bj) - resSalt (ir,jr)) / dtR
       bbl_TendTheta(ir,jr,bi,bj) = bbl_TendTheta(ir,jr,bi,bj)
     &      + tendT / rA(ir,jr,bi,bj) / thk(ir,jr)
       bbl_TendSalt (ir,jr,bi,bj) = bbl_TendSalt (ir,jr,bi,bj)
     &      + tendS / rA(ir,jr,bi,bj) / thk(ir,jr)
       sumIn  (ir,jr) = sumIn  (ir,jr) + dVol
       inTheta(ir,jr) = inTheta(ir,jr) + dVol*bbl_theta(id,jd,bi,bj)
       inSalt (ir,jr) = inSalt (ir,jr) + dVol*bbl_salt (id,jd,bi,bj)
      ENDIF

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: BBL_BOTTOM_STATE

C     !INTERFACE:
      SUBROUTINE BBL_BOTTOM_STATE( bi, bj,
     O       tLoc, sLoc, thk, resThk, resTheta, resSalt, rhoLoc, rhoBBL,
     I       myThid )

C     !DESCRIPTION:
C     T/S, thickness and density of the bottommost wet cell, and
C     thickness and T/S of its residual (non-BBL) volume.

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "BBL.h"

C     !INPUT/OUTPUT PARAMETERS:
      INTEGER bi, bj, myThid
      _RL     tLoc    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     sLoc    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     thk     (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     resThk  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     resTheta(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     resSalt (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     rhoLoc  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     rhoBBL  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
CEOP

C     !LOCAL VARIABLES:
      INTEGER i, j, kBot

        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          kBot        = max(1,kLowC(i,j,bi,bj))
          tLoc(i,j)   = theta(i,j,kBot,bi,bj)
          sLoc(i,j)   = salt (i,j,kBot,bi,bj)
          thk(i,j)    = drF(kBot) * hFacC(i,j,kBot,bi,bj)
          resThk(i,j) = thk(i,j) - bbl_eta(i,j,bi,bj)
          IF ( kLowC(i,j,bi,bj).GT.0 .AND. resThk(i,j).GT.0. ) THEN
           resTheta(i,j) = ( tLoc(i,j)*thk(i,j)
     &          - bbl_theta(i,j,bi,bj)*bbl_eta(i,j,bi,bj) )/resThk(i,j)
           resSalt (i,j) = ( sLoc(i,j)*thk(i,j)
     &          - bbl_salt (i,j,bi,bj)*bbl_eta(i,j,bi,bj) )/resThk(i,j)
          ELSE
           resThk  (i,j) = 0. _d 0
           resTheta(i,j) = tLoc(i,j)
           resSalt (i,j) = sLoc(i,j)
          ENDIF
          rhoLoc(i,j) = rhoInSitu(i,j,kBot,bi,bj)
          IF ( kBot .EQ. Nr ) THEN
           rhoBBL(i,j) = bbl_rho_nr(i,j,bi,bj)
          ELSE
           rhoBBL(i,j) = rhoInSitu(i,j,kBot+1,bi,bj)
          ENDIF
         ENDDO
        ENDDO


      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: BBL_RAW_VOLUMES

C     !INTERFACE:
      SUBROUTINE BBL_RAW_VOLUMES( bi, bj,
     I       dirX, dirY, velX, velY, resThk,
     O       volX, volY, sumOut, sumIn,
     I       myThid )

C     !DESCRIPTION:
C     Volume exchanged through each active face before the limiter, and
C     the total volume each cell gives away / receives.

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "BBL.h"

C     !INPUT/OUTPUT PARAMETERS:
      INTEGER bi, bj, myThid
      INTEGER dirX    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER dirY    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     velX    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     velY    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     resThk  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     volX    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     volY    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     sumOut  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     sumIn   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
CEOP

C     !LOCAL VARIABLES:
      INTEGER i, j, id, jd, ir, jr

      DO j=1-OLy,sNy+OLy
       DO i=1-OLx,sNx+OLx
        volX  (i,j) = 0. _d 0
        volY  (i,j) = 0. _d 0
        sumOut(i,j) = 0. _d 0
        sumIn (i,j) = 0. _d 0
       ENDDO
      ENDDO

C--   Unlimited volume transport through each active face: each face
C     on its own may move at most half the donor bbl volume and half
C     the receiver residual volume (as in the original scheme).
        DO j=1-OLy,sNy+OLy
         DO i=2-OLx,sNx+OLx
          IF ( dirX(i,j) .NE. 0 ) THEN
           IF ( dirX(i,j) .EQ. 1 ) THEN
            id = i-1
            ir = i
           ELSE
            id = i
            ir = i-1
           ENDIF
           volX(i,j) = max( 0. _d 0, min(
     &          bbl_eta(id,j,bi,bj) * rA(id,j,bi,bj) * 0.5 _d 0,
     &          resThk(ir,j) * rA(ir,j,bi,bj) * 0.5 _d 0,
     &          dyG(i,j,bi,bj)*bbl_eta(id,j,bi,bj)*velX(i,j)
     &          *dTtracerLev(kLowC(id,j,bi,bj)) ))
           sumOut(id,j) = sumOut(id,j) + volX(i,j)
           sumIn (ir,j) = sumIn (ir,j) + volX(i,j)
          ENDIF
         ENDDO
        ENDDO
        DO j=2-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          IF ( dirY(i,j) .NE. 0 ) THEN
           IF ( dirY(i,j) .EQ. 1 ) THEN
            jd = j-1
            jr = j
           ELSE
            jd = j
            jr = j-1
           ENDIF
           volY(i,j) = max( 0. _d 0, min(
     &          bbl_eta(i,jd,bi,bj) * rA(i,jd,bi,bj) * 0.5 _d 0,
     &          resThk(i,jr) * rA(i,jr,bi,bj) * 0.5 _d 0,
     &          dxG(i,j,bi,bj)*bbl_eta(i,jd,bi,bj)*velY(i,j)
     &          *dTtracerLev(kLowC(i,jd,bi,bj)) ))
           sumOut(i,jd) = sumOut(i,jd) + volY(i,j)
           sumIn (i,jr) = sumIn (i,jr) + volY(i,j)
          ENDIF
         ENDDO
        ENDDO

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: BBL_FACE_APPLY_PTR

C     !INTERFACE:
      SUBROUTINE BBL_FACE_APPLY_PTR(
     I        dVol, id, jd, ir, jr, bi, bj, thk, resThk,
     U        inPtr,
     I        myThid )

C     !DESCRIPTION:
C     Same as BBL_FACE_APPLY, for the passive tracers: a volume dVol of
C     donor (id,jd) bbl water is exchanged against receiver (ir,jr)
C     residual water; tendency of both bottommost wet cells (over the
C     passive tracer time step of their level) and tracer content
C     received by the receiver, for cells in the tile interior.

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "BBL.h"
#ifdef ALLOW_PTRACERS
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "PTRACERS_FIELDS.h"
#include "BBL_PTR.h"
#endif /* ALLOW_PTRACERS */

C     !INPUT/OUTPUT PARAMETERS:
C     dVol      :: exchanged volume (m^3)
C     id,jd     :: donor cell indices
C     ir,jr     :: receiver cell indices
C     bi,bj     :: tile indices
C     thk       :: thickness of bottommost wet grid cell
C     resThk    :: thickness of residual (non-BBL) volume
C     inPtr     :: received passive tracer content
C     myThid    :: my Thread Id number
      _RL     dVol
      INTEGER id, jd, ir, jr, bi, bj
      _RL     thk     (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     resThk  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
#ifdef ALLOW_PTRACERS
      _RL     inPtr   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,PTRACERS_num)
#else
      _RL     inPtr(1)
#endif
      INTEGER myThid
CEOP

#ifdef ALLOW_PTRACERS
C     !LOCAL VARIABLES:
C     iTr         :: passive tracer index
C     kD, kR      :: bottommost wet level of donor / receiver
C     resP        :: passive tracer of the receiver residual volume
C     tendP       :: tracer content moved per unit time
C     inD, inR    :: donor / receiver cell is in the tile interior
      INTEGER iTr, kD, kR
      _RL     resP, tendP
      LOGICAL inD, inR

      inD = id.GE.1 .AND. id.LE.sNx .AND. jd.GE.1 .AND. jd.LE.sNy
      inR = ir.GE.1 .AND. ir.LE.sNx .AND. jr.GE.1 .AND. jr.LE.sNy
      kD  = max(1,kLowC(id,jd,bi,bj))
      kR  = max(1,kLowC(ir,jr,bi,bj))

      DO iTr=1,PTRACERS_numInUse
C     same expressions as resTheta in BBL_BOTTOM_STATE
       IF ( resThk(ir,jr) .GT. 0. ) THEN
        resP = ( pTracer(ir,jr,kR,bi,bj,iTr)*thk(ir,jr)
     &       - bbl_ptr(ir,jr,bi,bj,iTr)*bbl_eta(ir,jr,bi,bj) )
     &       / resThk(ir,jr)
       ELSE
        resP = pTracer(ir,jr,kR,bi,bj,iTr)
       ENDIF
       IF ( inD ) THEN
        tendP = dVol * (bbl_ptr(id,jd,bi,bj,iTr) - resP)
     &        / PTRACERS_dTLev(kD)
        bbl_TendPtr(id,jd,bi,bj,iTr) = bbl_TendPtr(id,jd,bi,bj,iTr)
     &       - tendP / rA(id,jd,bi,bj) / thk(id,jd)
       ENDIF
       IF ( inR ) THEN
        tendP = dVol * (bbl_ptr(id,jd,bi,bj,iTr) - resP)
     &        / PTRACERS_dTLev(kR)
        bbl_TendPtr(ir,jr,bi,bj,iTr) = bbl_TendPtr(ir,jr,bi,bj,iTr)
     &       + tendP / rA(ir,jr,bi,bj) / thk(ir,jr)
        inPtr(ir,jr,iTr) = inPtr(ir,jr,iTr)
     &       + dVol*bbl_ptr(id,jd,bi,bj,iTr)
       ENDIF
      ENDDO
#endif /* ALLOW_PTRACERS */

      RETURN
      END
