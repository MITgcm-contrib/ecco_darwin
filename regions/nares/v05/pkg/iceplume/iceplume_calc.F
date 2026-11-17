#include "ICEPLUME_OPTIONS.h"

CBOP
C     !ROUTINE: ICEPLUME_CALC
C     !INTERFACE:
      SUBROUTINE ICEPLUME_CALC(
     I     myTime, myIter,
     I     myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE ICEPLUME_CALC
C     | o Send ambient conditions to plume model
C     | o Calculate source/sink terms to parameterise movement of
C     |   water and tracers between vertical layers by resulting plume
C     | o Calculate melt rates and tendencies due to melting of
C     |   ice front in none-plume locations
C     |   
C     *==========================================================*

C     | atn 09.Oct.2021: added notes
C     | Because the thermodynamic tendencies code was copied from
C     | pkg/icefront , the sign conventions are as follows:
C     | Outward (pos) = leaving the ocean, inward (neg) = entering ocean.
C     | Freshwater (kg/m2/s): inward (negative) freshwater flux implies 
C     | glacier melting due to outward (positive) heat flux (W/m^2)

C     | kiki 27. October 2021: started to implement new plume option 6&7
C     | frictionless plume following Wells and Worster, 2008 to simulate
C     | purely melt water driven plumes (see Jackson et al., 2020, GRL)
C     | and truncated line plume by Jackson et al., 2017
catn I flipped so that truncated = 6, buoy = 7, in case Kiki removes 7
C     \ev
C     !USES:
      IMPLICIT NONE
C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "FFIELDS.h"
#include "ICEPLUME.h"

#ifdef ALLOW_PTRACERS
#include "PTRACERS_PARAMS.h"
#include "PTRACERS_START.h"
#include "PTRACERS_FIELDS.h"
#endif

#ifdef ALLOW_DIAGNOSTICS_VEC
#include "ICEPLUME_DV_FIELDS.h"
#endif

#ifdef ALLOW_SHELFICE
#include "SHELFICE.h"
#endif

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
      _RL myTime
      INTEGER myIter
      INTEGER myThid
      CHARACTER*(MAX_LEN_MBUF) msgBuf

#ifdef ALLOW_ICEPLUME

C     !LOCAL VARIABLES:
C     == Local variables ==
C     I,J,K,bi,bj  :: loop indices
C     msgBuf       :: Informational/error message buffer
C     sProf, tProf, uProf, vProf :: salt, pot. temperature and
C                            uVel and vVel profiles
C     ptrProf  :: ambient ptracer profile
C     ptrPlumeCum :: cumulative quantity of ptracer in plume
C     ptrPlume :: quantity of ptracer in plume outflow
C     eps5     :: for thermodynamics (see pkg icefront)
C     maxDepth :: vertical extent of domain (m)
C     plumeAreaInCell :: surface area of plume in contact with ice in that cell (m^2)
C     negSum, posSum :: sum of negative and positive contributions to the plume volume
C     posNegRatio    :: ratio of the above
C     wVelCell       :: vertical velocity component at cell centres
C     hVelCell_tangential :: horizontal velocity component at cell centres
C     meanVel :: ice tangental velocity
catnC rho_0 :: average density of seawater
catnC lhfusion :: latent heat of fusion for solid to liquid phase, Joule
catn              = ICEPLUMElatentHeat, so we use that.
C     rho_shelfice :: 917. kg/m3
C     secInday :: number of seconds in a day 
catn: introducing a new flag below to avoid confusion between where pos(neg) 
catn: mask is NorthSouth or EastWest oriented glacier front
C     GlacierFront_is_NorthSouth :: 1 for the glacier oriented NorthSouth
C                                :: 0 for EastWest, default is -9999

      INTEGER bi, bj
      INTEGER J, K, I
#ifdef ALLOW_DIAGNOSTICS_VEC
      INTEGER i_dv, j_dv, bi_dv, bj_dv
      INTEGER fld3D, loc_count, nFlds3D, numvalid, vec_id, vec_tag
#endif
      INTEGER GlacierFront_is_NorthSouth

      _RL eps5
      _RL plumeAreaInCell
      _RL negSum, posSum, posNegRatio
      _RL meltVol
      _RL posOut, saltCum, sOut, prFace
      INTEGER kp
      _RL wVelCell, hVelCell_tangential, meanVel
      _RL sw_temp, sw_ptmp
catn: Many loose rhos being defined.  Here we make rho_0 = rhoConst
catn: Second, make rho_shelfice a readable param to avoid hardcoding.
catn      _RL rho_0
catn      _RL rho_shelfice
      _RL secInday
      _RL tmpT, tmpS
      _RL tModMWplume
      external SW_TEMP
      external SW_PTMP

catn      PARAMETER(rho_0=1027.0D0)
      PARAMETER(secInday=86400.0D0)
catn      PARAMETER(rho_shelfice=917.0D0)
catn      PARAMETER(lhfusion=333.55E+3)
catn: be careful about overlap, i think these will all need them.
      _RL gT_iceplumeBG   (1-Olx:sNx+Olx,1-Oly:sNy+Oly,Nr)
      _RL gS_iceplumeBG   (1-Olx:sNx+Olx,1-Oly:sNy+Oly,Nr)
      _RL Fwflux3D        (1-Olx:sNx+Olx,1-Oly:sNy+Oly,Nr)
      _RL HeatFlux3D      (1-Olx:sNx+Olx,1-Oly:sNy+Oly,Nr)
      _RL rProfPlume3D    (1-Olx:sNx+Olx,1-Oly:sNy+Oly,Nr)
      _RL wProfPlume3D    (1-Olx:sNx+Olx,1-Oly:sNy+Oly,Nr)
      _RL tProfPlume3D    (1-Olx:sNx+Olx,1-Oly:sNy+Oly,Nr)
      _RL sProfPlume3D    (1-Olx:sNx+Olx,1-Oly:sNy+Oly,Nr)
      _RL mProfPlume3D    (1-Olx:sNx+Olx,1-Oly:sNy+Oly,Nr)
      _RL mProfAv3D       (1-Olx:sNx+Olx,1-Oly:sNy+Oly,Nr)
      INTEGER kSrf        (1-Olx:sNx+Olx,1-Oly:sNy+Oly)
catn no longer need runoff[Vel,Rad] if we read in Qsg, see below
c      _RL runoffVel       (1-Olx:sNx+Olx,1-Oly:sNy+Oly)
c      _RL runoffRad       (1-Olx:sNx+Olx,1-Oly:sNy+Oly)
#ifdef ICEPLUME_ALLOW_DETACHED_PLUME 
      _RL thetaProfPlume3D(1-Olx:sNx+Olx,1-Oly:sNy+Oly,Nr)
      _RL distanceProfPlume3D(1-Olx:sNx+Olx,1-Oly:sNy+Oly,Nr)
#endif /* ICEPLUME_ALLOW_DETACHED_PLUME */

#ifdef ALLOW_PTRACERS
      INTEGER iTracer
      _RL ptrPlume (PTRACERS_num)
      _RL ptrPlumeCum (PTRACERS_num)
      _RL ptrProf  (Nr,PTRACERS_num)
#endif

catn: brought out of do-loop, factor for heat tendency
catn: mass2rUnit=1/rhoConst (set in data), e.g., 1/1029. [m3/kg]
catn: HeatCapacity_Cp = 3994. [J/kg/degC]
cat:  so eps5 = 1/1029./3994. = 2.433193035422911e-07 [degC m3 / J]
      eps5 = mass2rUnit/HeatCapacity_Cp

#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_ENTER('ICEPLUME_CALC',myThid)
#endif
C -----------------------------------------
C Enter into loops
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)

catn: Initialize tendencies for plume
        DO k=1,Nr
         DO j=1-Oly,sNy+Oly
          DO i=1-Olx,sNx+Olx
           gT_iceplumeBG(i,j,k)=0. _d 0
           gS_iceplumeBG(i,j,k)=0. _d 0
           addMass3Dbg(i,j,k,bi,bj)=0. _d 0
           Fwflux3D(i,j,k)=0. _d 0
           Heatflux3D(i,j,k)=0. _d 0
           rProfPlume3D(i,j,k)=0. _d 0
           wProfPlume3D(i,j,k)=0. _d 0
           tProfPlume3D(i,j,k)=0. _d 0
           sProfPlume3D(i,j,k)=0. _d 0
           mProfPlume3D(i,j,k)=0. _d 0
           mProfAv3D(i,j,k)=0. _d 0
#ifdef ICEPLUME_ALLOW_DETACHED_PLUME 
           thetaProfPlume3D(i,j,k)=0. _d 0
           distanceProfPlume3D(i,j,k)=0. _d 0
#endif /* ICEPLUME_ALLOW_DETACHED_PLUME */
          ENDDO
         ENDDO
        ENDDO

C     Top wet level of each column: kSurfC is 1 in open water and the
C     level just below the ice base under pkg/shelfice. Plume quantities
C     are only computed (and only non-zero) for kSrf <= k.
        DO j=1-Oly,sNy+Oly
          DO i=1-Olx,sNx+Olx
           kSrf(i,j) = MAX( 1, MIN( kSurfC(i,j,bi,bj), Nr ) )
          ENDDO
        ENDDO


c when read in Qsg, initialize LOCAL runoff[Vel,Rad] array, no bibj
catn now removing runoff[Vel,Rad] local arrays if read Qsg
c        DO j=1-Oly,sNy+Oly
c         DO i=1-Olx,sNx+Olx
c          runoffVel(i,j)=0. _d 0
c          runoffRad(i,j)=0. _d 0
c         ENDDO
c        ENDDO

catn: When read in Qsg, add a loop outside to calculate the 2d fields
catn: runoff(Vel,Rad), which depends on geometry dLtangential
ckiki: use[sheet,cone,detatch,etc]plume not set here yet. 
ckiki: moved loop to later
c        DO j=1-OLy,sNy+OLy
c         DO i=1-OLx,sNx+OLx
c          IF ( plumeMask(I,J,bi,bj) .GT. 0 ) THEN
c            dLtangential = dyG(I,J,bi,bj)
c          ELSEIF ( plumeMask(I,J,bi,bj) .LT. 0 ) THEN
c            dLtangential = dxG(I,J,bi,bj)
c          ENDIF
c          IF(runoffQsg(I,J,bi,bj) .GT. 0. _d 0) THEN
c           runoffVel(i,j) = wVel_sg_0
c           IF ( useSheetPlume ) THEN
c            runoffRad(i,j)=runoffQsg(I,J,bi,bj)/
c     &                     (dLtangential*runoffVel(i,j))
c           ELSEIF (useConePlume .OR. useDetachPlume) THEN
c            runoffRad(i,j)=(twoRL*runoffQsg(I,J,bi,bj)/
c     &                     (runoffVel(i,j)*pi))**halfRL
catn to Kiki: please check this and put in correct physics
catn for truncated and buoy plume
c           ENDIF
c          ELSE
c           runoffVel(i,j) = 0. _d 0
c           runoffRad(i,j) = 0. _d 0
c          ENDIF
c         ENDDO
c        ENDDO

C Loop over interior points only: plumeMask/runoffQsg are read without
C an exchange (overlaps are 0), the background-melt velocities use i+1 /
C j+1, and addMass is only needed in the interior. Restricting the loop
C makes the result independent of the tile decomposition.
        DO J = 1,sNy
         DO I = 1,sNx
catn: check if we need overlap at all? if not we do without
catn        DO J = 1,sNy
catn         DO I = 1,sNx

#ifdef ALLOW_PTRACERS

C Clear local plume ptracer variables

          DO iTracer = 1,PTRACERS_num
           ptrPlume(iTracer)    = 0.D0
           ptrPlumeCum(iTracer) = 0.D0
          ENDDO

#endif /* ALLOW_PTRACERS */

catn initialize  
          GlacierFront_is_NorthSouth = -9999

C Check to see if there is ice in that cell. If not, skip to end.

          IF ( plumeMask(I,J,bi,bj) .NE. 0 ) THEN

C               PRINT *,'MWPRINT Found a plume or melt point at ',I,J

C Read from the plume mask which type of plume should be used in this cell.

C 1 = ice but no plume (melting only)
C 2 = sheet plume (Jenkins)
C 3 = half-conical plume (Morton/Slater)
C 4 = both sheet plume and half-conical plume (NOT YET IMPLEMENTED)
C 5 = detaching conical plume (Goldberg)
C 6 = truncated line plume (Jacksom et al., 2017, GRL) (Kiki)
C 7 = buoyancy driven sublayer, Wells and Worster 2008 (Kiki)

C POSITIVE values indicate ice front is orientated north-south
C NEGATIVE values indicate ice front is orientated east-west
catn: must use the sign as above to determine delta_[x,y]

C If there is subglacial discharge but no plume type defined, there will be no
C plume.

catn: This definition must go before the Q_sg, moved here from below.
catn: Also, we must avoid hardcoding delta_x and dLtangential to [x,y] of
catn: model grid because it depends on the geometry of the plume, which
catn: is including in the sign of the plumeMask.  So, since we are
catn: looping through I,J starting at the entry of this call, best to
catn: move defn of delta_[x,y] to above when we are checking plumeMask.

catn: In original code, it has been assummed that dyG is across the
catn: glacier and dxG is down the fjord toward open ocean. So
catn: we first define using that assumption:
C         Cell resolution
          IF ( plumeMask(I,J,bi,bj) .GT. 0 ) THEN
            GlacierFront_is_NorthSouth = 1
            dLnormal = dxG(I,J,bi,bj)
            dLtangential = dyG(I,J,bi,bj)
          ELSEIF ( plumeMask(I,J,bi,bj) .LT. 0 ) THEN
            GlacierFront_is_NorthSouth = 0
catn: but now we check, if the plume geometry is such that dxG is across
catn: the glacier and dyG along the fjord toward open ocean:
            dLnormal = dyG(I,J,bi,bj)
            dLtangential = dxG(I,J,bi,bj)
          ENDIF

          IF ( plumeMask(I,J,bi,bj) .EQ. -1
     &         .OR. plumeMask(I,J,bi,bj) .EQ. 1) THEN
            useSheetPlume = .FALSE.
            useConePlume = .FALSE.
            useDetachPlume = .FALSE.
            useTruncPlume = .FALSE.
            useBuoyPlume = .FALSE.
          ELSEIF ( plumeMask(I,J,bi,bj) .EQ. -2
     &         .OR. plumeMask(I,J,bi,bj) .EQ. 2) THEN
            useSheetPlume = .TRUE.
            useConePlume = .FALSE.
            useDetachPlume = .FALSE.
            useTruncPlume = .FALSE.
           useBuoyPlume = .FALSE.
          ELSEIF ( plumeMask(I,J,bi,bj) .EQ. -3
     &         .OR. plumeMask(I,J,bi,bj) .EQ. 3) THEN
            useSheetPlume = .FALSE.
            useConePlume = .TRUE.
            useDetachPlume = .FALSE.
            useTruncPlume = .FALSE.
            useBuoyPlume = .FALSE.
          ELSEIF ( plumeMask(I,J,bi,bj) .EQ. -4
     &         .OR. plumeMask(I,J,bi,bj) .EQ. 4) THEN
            useSheetPlume = .TRUE.
            useConePlume = .TRUE.
            useDetachPlume = .FALSE.
            useTruncPlume = .FALSE.
            useBuoyPlume = .FALSE.
          ELSEIF ( plumeMask(I,J,bi,bj) .EQ. -5
     &         .OR. plumeMask(I,J,bi,bj) .EQ. 5) THEN
#ifdef ICEPLUME_ALLOW_DETACHED_PLUME          
            useSheetPlume = .FALSE.
            useConePlume = .FALSE.
            useDetachPlume = .TRUE.
            useTruncPlume = .FALSE.
            useBuoyPlume = .FALSE.
#else /* ICEPLUME_ALLOW_DETACHED_PLUME */
            WRITE(msgBuf,'(2A)')
     &      'cannot use detaching plume without ',
     &      '#define ICEPLUME_ALLOW_DETACHED_PLUME'
            CALL PRINT_ERROR( msgBuf, myThid )
            CALL ALL_PROC_DIE( myThid )
            STOP 'ABNORMAL END: S/R ICEPLUME_CALC'
#endif
          ELSEIF ( plumemask(I,J,bi,bj) .EQ. -6
     &         .OR. plumemask(I,J,bi,bj) .EQ. 6) THEN
            useSheetPlume = .FALSE.
            useConePlume = .FALSE.
            useDetachPlume = .FALSE.
            useTruncPlume = .TRUE.
            useBuoyPlume = .FALSE.
          ELSEIF ( plumemask(I,J,bi,bj) .EQ. -7
     &         .OR. plumemask(I,J,bi,bj) .EQ. 7) THEN
            useSheetPlume = .FALSE.
            useConePlume = .FALSE.
            useDetachPlume = .FALSE.
            useTruncPlume = .FALSE.
            useBuoyPlume = .TRUE.
          ELSE

            WRITE(msgBuf,'(A,1PE12.4,A,2I6,A,2I4)')
     &      'ICEPLUME_CALC: plumeMask = ', plumeMask(I,J,bi,bj),
     &      ' at i,j =', I, J, ' bi,bj =', bi, bj
            CALL PRINT_ERROR( msgBuf, myThid )
            WRITE(msgBuf,'(A)')
     &      'ICEPLUME_CALC: plume mask value must be between -7 and 7'
            CALL PRINT_ERROR( msgBuf, myThid )
            CALL ALL_PROC_DIE( myThid )
            STOP 'ABNORMAL END: S/R ICEPLUME_CALC'

          ENDIF
C Reset all per-column plume state so that nothing (r_sg, L_sg, volFlux,
C volFluxDiff, plume profiles) can carry over from the previous column.
          Q_sg = 0. _d 0
          w_sg = 0. _d 0
          r_sg = 0. _d 0
          L_sg = 0. _d 0
          DO k=1,Nr+1
           volFlux(k) = 0. _d 0
           rProfPlume(k) = 0. _d 0
           wProfPlume(k) = 0. _d 0
           tProfPlume(k) = 0. _d 0
           sProfPlume(k) = 0. _d 0
           aProfPlume(k) = 0. _d 0
           mIntProfPlume(k) = 0. _d 0
#ifdef ICEPLUME_ALLOW_DETACHED_PLUME
           thetaProfPlume(k) = 0. _d 0
           distanceProfPlume(k) = 0. _d 0
#endif
          ENDDO
          DO k=1,Nr
           volFluxDiff(k) = 0. _d 0
          ENDDO

C Plume type +-1 is "melt only, no plume": discharge there has no plume
C model to carry it and would previously have used leftover r_sg from
C another column. Stop with a clear message instead.
          IF ( ( plumeMask(I,J,bi,bj) .EQ. -1
     &           .OR. plumeMask(I,J,bi,bj) .EQ. 1 )
     &         .AND. runoffQsg(I,J,bi,bj) .GT. 0. _d 0 ) THEN
            WRITE(msgBuf,'(A,1PE12.4,A,2I6,A,2I4)')
     &      'ICEPLUME_CALC: runoffQsg = ', runoffQsg(I,J,bi,bj),
     &      ' m3/s in a plumeMask=+-1 (melt-only) cell, i,j =', I, J,
     &      ' bi,bj =', bi, bj
            CALL PRINT_ERROR( msgBuf, myThid )
            WRITE(msgBuf,'(2A)') 'ICEPLUME_CALC: give the cell a ',
     &      'plume type (2,3,5,6,7) or move the discharge'
            CALL PRINT_ERROR( msgBuf, myThid )
            CALL ALL_PROC_DIE( myThid )
            STOP 'ABNORMAL END: S/R ICEPLUME_CALC'
          ENDIF

ckiki: inserted loop from above here
catn: In addition, since [w,r]_sg are assigned to runoff[Vel,Rad]
catn  immediately here, we remove these two local arrays altogether
          IF(runoffQsg(I,J,bi,bj) .GT. 0. _d 0) THEN
           Q_sg=runoffQsg(I,J,bi,bj)
           w_sg = wVel_sg_0
           IF ( useSheetPlume ) THEN
            r_sg=runoffQsg(I,J,bi,bj)/
     &                     (dLtangential*w_sg)
           ELSEIF ( useTruncPlume ) THEN
            r_sg=runoffQsg(I,J,bi,bj)/
     &                     (plumeLength(I,J,bi,bj)*w_sg)
            L_sg=plumeLength(i,j,bi,bj)   
           ELSEIF (useConePlume .OR. useDetachPlume) THEN
            r_sg=(twoRL*runoffQsg(I,J,bi,bj)/
     &                     (w_sg*pi))**halfRL
           ELSEIF ( useBuoyPlume ) THEN
C           sheet geometry (r_sg was never set for type 7, so the previous
C           column's value was used)
            r_sg=runoffQsg(I,J,bi,bj)/
     &                     (dLtangential*w_sg)
           ENDIF
          ELSE
           Q_sg = 0. _d 0
           w_sg = 0. _d 0
           r_sg = 0. _d 0
          ENDIF
C Read in the subglacial discharge for this place and time
C (only the lowermost cell in column is read)
catn: when read in Qsg, runoff[Vel,Rad] are local, no bibj
cccatn           w_sg = runoffVel(I,J,bi,bj)
cccatn           r_sg = runoffRad(I,J,bi,bj)
ccc           w_sg = runoffVel(I,J)
ccc           r_sg = runoffRad(I,J)
catn: changed 2 files into 1 since the only thing we need is Qsg

catn: There is a potential error for delta_y here because it is yet
catn initialized, and yet plays a role in Q_sg calc.  Then after that we
catn define dLtangential, and then compute (further down) volumeFlux based on
catn the redefined dLtangential.  atn has now moved defn of delta_[x,y] above.
catn Lastly, [x,y] has now been moved to normal/tangential to glacier front.
catn           IF ( useSheetPlume ) THEN
catn           Q_sg = w_sg*r_sg*dLtangential
catn           ELSEIF ( useConePlume ) THEN
catn           Q_sg = w_sg*((pi*r_sg**2.)/2.)
catn           ELSEIF ( useDetachPlume ) THEN
catn           Q_sg = w_sg*((pi*r_sg**2.)/2.)
catn           ELSEIF ( useTruncPlume ) THEN
catnckiki: needs to be implemented using Lp, for now only placeholder
catn           Q_sg = w_sg*r_sg
catn#ifdef ICEPLUME_ALLOW_SCHULZ22
catn           ELSEIF ( useBuoyPlume ) THEN
catn           Q_sg = 1.0 _d -10
catn#endif /* ICEPLUME_ALLOW_SCHULZ22 */
catn           ELSE
catn           Q_sg = 0. _d 0
catn           ENDIF
c          IF (runoffQsg(I,J,bi,bj) .GT. 0. _d 0 .AND. 
c     &               w_sg .GT. 0. _d 0) THEN
c            Q_sg=runoffQsg(I,J,bi,bj)
c            IF ( useSheetPlume ) THEN
c             r_sg=Q_sg/(dLtangential*w_sg)
c            ELSEIF ( useTruncPlume ) THEN
c             r_sg=Q_sg/(plumeLength(I,J,bi,bj)*w_sg)
c            L_sg=plumeLength(i,j,bi,bj)   
c            ELSEIF (useConePlume .OR. useDetachPlume) THEN
c             r_sg=(twoRL*Q_sg/(w_sg*pi))**halfRL
c            ELSEIF ( useBuoyPlume ) THEN
c             r_sg = Q_sg/(dLtangential*w_sg)
c            ENDIF
c          ELSE
c            Q_sg = 0. _d 0
c          ENDIF
c
C          PRINT *, 'MWPRINT At',I,J,' Q_sg = ',Q_sg


C Create variables with temperature, salinity
C and velocity profiles for that column

CMW           DO K = 1,Nr
            DO K=kSrf(I,J),Nr
C           Tracers
C           Pressure in dbar (1 dbar = 1e4 Pa), as SW_TEMP/SW_PTMP expect;
C           the previous factor 1e-6 gave MPa, i.e. 100x too small.
            prProf(k) = ABS(rC(k))*rhoConst*gravity*1.0 _d -4
            sProf(K)  = salt(I,J,K,bi,bj)         ! Salinity
            ptProf(K) = theta(I,J,K,bi,bj)        ! Potential Temperature
            tProf(k)  = 
     &      SW_TEMP(sProf(k),ptProf(k),prProf(k),0. _d 0) ! Temperature

C       IF ((kSrf(I,J).gt.1).and.
C     &        (plumeMask(I,J,bi,bj).eq.6)) THEN
C          WRITE(msgBuf,'(A,I5,A,F10.5)') "ptProf(k=",k,
C     &       " ) = ",ptProf(k)
C          CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
C     &                     SQUEEZE_RIGHT, myThid )
C          ENDIF

#ifdef ALLOW_PTRACERS
            DO iTracer = 1,PTRACERS_num
             ptrProf(k,iTracer) = pTracer(I,J,K,bi,bj,iTracer)
            ENDDO
#endif /* ALLOW_PTRACERS */

C           Velocities
            vProf(k) = ABS(vVel(I,J,K,bi,bj))          ! v velocity
            uProf(K) = ABS(uVel(I,J,K,bi,bj))          ! u Velocity

            delta_z(k) = drF(K)

           ENDDO

C     Above the top wet level (under an ice shelf) there is no water: fill
C     the ambient profiles with the top-wet-cell values so that nothing
C     left over from a previously solved column can be used (e.g. by the
C     ODE profile interpolation), and no plume flux exists there.
           DO K=1,kSrf(I,J)-1
            prProf(k) = prProf(kSrf(I,J))
            sProf(k)  = sProf(kSrf(I,J))
            ptProf(k) = ptProf(kSrf(I,J))
            tProf(k)  = tProf(kSrf(I,J))
            uProf(k)  = 0. _d 0
            vProf(k)  = 0. _d 0
            delta_z(k) = drF(K)
#ifdef ALLOW_PTRACERS
            DO iTracer = 1,PTRACERS_num
             ptrProf(k,iTracer) = ptrProf(kSrf(I,J),iTracer)
            ENDDO
#endif /* ALLOW_PTRACERS */
           ENDDO

C Vertical extent of domain
          maxDepth = rF(Nr+1)
C Depth of seabed at glacer front
          iceDepth = R_low(I,J,bi,bj)


          IF ( iceDepth .EQ. 0 ) THEN
            WRITE(msgBuf,'(A,I6,A,I6,A,2I4,A)')
     &      'ICEPLUME_CALC: plume specified in cell I = ', I,
     &      ', J = ', J, ' (bi,bj =', bi, bj,
     &      ') but depth of this cell = 0'
            CALL PRINT_ERROR( msgBuf, myThid )
            CALL ALL_PROC_DIE( myThid )
            STOP 'ABNORMAL END: S/R ICEPLUME_CALC III'
          ENDIF

C Find grid layer at depth of ice face (move this to initialisation)
          icedepthK = 0
          DO K=1,Nr+1
C          DO K=kSrf(I,J),Nr
           IF ( rF(K) .EQ. iceDepth ) iceDepthK = K
          ENDDO

C If a matching grid layer is not found, this may be because the bottom layer is a partial cell
C In this case, start in cell above partial cell
          IF ( iceDepthK .EQ. 0 ) THEN
           DO K=1,Nr
C           DO K=kSrf(I,J),Nr
            IF ( rF(K) .GT. iceDepth ) THEN
              IF ( rF(K+1) .LT. iceDepth ) THEN
               iceDepthK = K
              ENDIF
            ENDIF
           ENDDO
          ENDIF

C If we still cannot find the bottom cell

          IF ( iceDepthK .EQ. 0 ) THEN
             WRITE(msgBuf,'(A,I6,A,I6,A,2I4,A)')
     &       'ICEPLUME_CALC: unable to identify index of cell at '//
     &       'grounding line, I = ', I, ', J = ', J,
     &       ' (bi,bj =', bi, bj, ')'
             CALL PRINT_ERROR( msgBuf, myThid )
             CALL ALL_PROC_DIE( myThid )
             STOP 'ABNORMAL END: S/R ICEPLUME_CALC IV'
          ENDIF

C A column whose only wet cell is a partial bottom cell that is also the
C top wet cell (thin water column at an ice-shelf grounding line): the
C search above returns iceDepthK = kSrf and the plume loops would be
C empty, losing the discharge. Integrate across that single cell instead.
          IF ( iceDepthK .LE. kSrf(I,J) ) iceDepthK = kSrf(I,J) + 1

C --- If there is subglacial outflow in that column, then parameterise plume ---

C          IF ((kSrf(I,J).gt.1).and.
C     &        (plumeMask(I,J,bi,bj).eq.6)) THEN
C          WRITE(msgBuf,'(A,I5,A,I5)') "Before Qsg check: kSrf=",
C     &       kSrf(I,J)," iceDepthK=",iceDepthK
C         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
C     &                     SQUEEZE_RIGHT, myThid )
C         WRITE(msgBuf,'(A,F10.2,A,F10.2)') "Qsg=",
C     &       Q_sg,", runoffQsg(I,J,bi,bj)=",runoffQsg(I,J,bi,bj)
C          CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
C     &                     SQUEEZE_RIGHT, myThid )
C          WRITE(msgBuf,'(A,F10.2)') "maskC",maskC(i,j,1,bi,bj)
C          CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
C     &                     SQUEEZE_RIGHT, myThid )
C          ENDIF

           IF ( Q_sg.GT.0 ) THEN

C         IF ((kSrf(I,J).gt.1).and.
C     &        (plumeMask(I,J,bi,bj).eq.6)) THEN
C          WRITE(msgBuf,'(A,I5,A,I5)') "Before model: kSrf=",
C     &       kSrf(I,J)," iceDepthK=",iceDepthK
C          CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
C     &                     SQUEEZE_RIGHT, myThid )
C          ENDIF

C The plume can rise at most to the top of the water column: the
C ice-shelf base under pkg/shelfice, the surface otherwise.
            plumeTopK = kSrf(I,J)
            plumeTopZ = Ro_surf(I,J,bi,bj)
            plumeI  = I
            plumeJ  = J
            plumeBi = bi
            plumeBj = bj
C This routine calculates T, S and W and r profiles for plume
            CALL ICEPLUME_PLUME_MODEL (mythid)

C          IF ((kSrf(I,J).gt.1).and.
C     &        (plumeMask(I,J,bi,bj).eq.6)) THEN
C          WRITE(msgBuf,'(A,I5,A,I5)') "After model: kSrf=",
C     &       kSrf(I,J)," iceDepthK=",iceDepthK
C          CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
C     &                     SQUEEZE_RIGHT, myThid )
C          ENDIF

C Calculate vertical plume volume flux...

CMW            DO k=1,Nr
            DO k=kSrf(I,J),Nr
C ... after checking to see if we are above the base of the ice face...
             IF ( K .LT. iceDepthK ) THEN
C ... assuming specified plume horizontal extent (for sheet flow)...
              IF ( useSheetPlume ) THEN
               volFlux(k) = wProfPlume(k)*rProfPlume(k)*dLtangential
C ... or assuming half-conical form
              ELSEIF ( useConePlume ) THEN
               volFlux(k)=pi*(rProfPlume(k)**twoRL)*wProfPlume(k)/twoRL
C ... assuming detached-conical form
              ELSEIF ( useDetachPlume ) THEN
               volFlux(k)=rProfPlume(k)*wProfPlume(k)
C ... assuming truncated line plume 
              ELSEIF ( useTruncPlume ) THEN
               volFlux(k) =wProfPlume(k)*rProfPlume(k)
     &                         *plumeLength(I,J,bi,bj)
C ... assuming buoyancy driven sublayer plume (as sheet plume)                                
              ELSEIF ( useBuoyPlume ) THEN
               volFlux(k) = wProfPlume(k)*rProfPlume(k)*dLtangential
              ENDIF
             ELSE
              volFlux(k) = 0. _d 0
             ENDIF
            ENDDO

C A couple of corrections:
C - even if plume is still buoyant, it cannot flow through the fjord surface
C   or through the base of an ice shelf: no flux through the top face of
C   the top wet cell (or anywhere above it)
            DO k=1,kSrf(I,J)
             volFlux(k) = 0. _d 0
            ENDDO
C - the initial volume flux is equal to runoff [m3/s]
            volflux(iceDepthK) = Q_sg

C Calculate volume flux differential to give entrainment / extrainment, [m3/s]
C First clear volfluxdiff

CMW            DO K = 1,Nr
            DO K=kSrf(I,J),Nr
             volfluxdiff(K) = 0. _d 0
            ENDDO

CMW            DO k=1,iceDepthK-1
            DO K=kSrf(I,J),iceDepthK-1
             volFluxDiff(k) = volFlux(k+1) - volFlux(k)
            ENDDO

catn: The loop below indicates that when conserveMass is true, then we
catn: scale the volume flux total to reduce the effect of Q_sg into the
catn: domain.  An external test shows that before entering this loop, we
catn: have a non-zero sum(volFluxDiff), but after entering and getting
catn: scaled, sum(volFluxDiff) = 0 (within precision), and any posSum at
catn: any particular klev is scaled while all negSum were untouched.
            IF ( conserveMass ) THEN
C Scale output to compensate for entrainment lost in expanding of output layer
C i.e. so that there is no net flow over boundary

catn: Sign convention: if volFluxDiff(k) is less than 1, there is a net
catn: divergence == volume loss at the klev, and will contrib to negSum
             negSum = 0. _d 0
             posSum = 0. _d 0

CMW             DO K = 1,Nr
             DO K=kSrf(I,J),Nr
              IF ( volFluxDiff(K) .LT. 0 ) THEN
               negSum = volFluxDiff(K) + negSum
              ELSE
               posSum = volFluxDiff(K) + posSum
              ENDIF
             ENDDO

catn: The check of posSum not eq 0 implies there is at least ONE klev where
catn: there is a net CONVERGENCE into the cell in one or more klev.  So we
catn: take the ratio of -negSum/posSum to check for the total columm
catn: what the net convergence is.  If this ratio is one, negSum=posSum
catn: and there is no net water excess.  If posSum>negSum, there is
catn: excess water input that needs to be taken care of under the
catn: scenario where ocean volume cannot be changed.  So take the
catn: example negSum = 4 and posSum = 5, ratio is 4/5.  This ratio will
catn: now be mult across all volFluxDiff>0, such that after scaling the
catn: net column convergence is zero
             IF ( posSum .NE. 0 ) THEN
              posNegRatio = -negSum / posSum

CMW              DO K = 1,Nr
             DO K=kSrf(I,J),Nr
               IF ( volFluxDiff(K) .GT. 0 )
     &            volFluxDiff(K) = volFluxDiff(K) * posNegRatio
              ENDDO
             ENDIF

            ENDIF  !conserveMass

C Plume meltwater. volFlux grows with height by entrained ambient water
C AND by meltwater from the ice face (mIntProfPlume = height-integrated
C melt volume flux, m3/s). Previously all of a layer's volume increase
C was removed from the ambient cell, i.e. the meltwater volume was taken
C from the ocean at the ambient salinity and the plume melt removed salt
C (M*S_ambient). Book the meltwater as a freshwater source instead: only
C the ambient part is removed from the cell, so the plume adds Q_sg plus
C the meltwater volume. (Not done with conserveMass, whose purpose is to
C add no net volume.)
            IF ( .NOT.conserveMass ) THEN
             DO K=kSrf(I,J),iceDepthK-1
              IF ( volFluxDiff(K) .LT. 0. _d 0 ) THEN
               meltVol = mIntProfPlume(K) - mIntProfPlume(K+1)
               meltVol = MAX( 0. _d 0, MIN( meltVol, -volFluxDiff(K) ) )
               volFluxDiff(K) = volFluxDiff(K) + meltVol
              ENDIF
             ENDDO
            ENDIF

C Outflow properties consistent with the outflow volume. The outflow
C volume of layer k is volFlux(k+1) - volFlux(k) (the flux through its
C bottom face when the plume stops at its top face), but the code used the
C plume T/S at the top face k, i.e. after the plume had also entrained
C (and melted) inside layer k - water that is never removed from the
C ambient cell. That created salt/heat (e.g. 4.5e6 g/kg m3 per day for a
C 100 m3/s plume stopped by a 300 m ice shelf). With conserveMass=.FALSE.:
C  - salinity = (entrained salt + Q_sg*S_sg_0)/(outflow volume), exact
C    since meltwater adds no salt;
C  - temperature = plume temperature at the bottom face of the outflow
C    layer (the water actually crossing it).
            posOut  = 0. _d 0
            saltCum = Q_sg*S_sg_0
            DO K=kSrf(I,J),iceDepthK-1
             IF ( volFluxDiff(K) .GT. 0. _d 0 ) THEN
              posOut = posOut + volFluxDiff(K)
             ELSE
              saltCum = saltCum - volFluxDiff(K)*sProf(K)
             ENDIF
            ENDDO
            sOut = 0. _d 0
            IF ( posOut .GT. 0. _d 0 ) sOut = saltCum / posOut

#ifdef ALLOW_PTRACERS

C Add up total sum of each tracer in plume
CMW            DO K=1,iceDepthK-1
            DO K=kSrf(I,J),iceDepthK-1
             IF (volFLuxDiff(k) .LT. 0. ) THEN
              DO iTracer = 1,PTRACERS_num
              ptrPlumeCum(iTracer) 
     &          = ptrPlumeCum(iTracer)
     &             +(-volFluxDiff(k)*ptrProf(k,iTracer))
              ENDDO
             ENDIF
            ENDDO

C Add ptracers in runoff

            IF ( useInputPtracers ) THEN
             DO iTracer = 1,PTRACERS_num

              IF (ptracerMask(I,J,iTracer,bi,bj) .NE. 0 ) THEN

               ptrPlumeCum(iTracer) =
     &         ptrPlumeCum(iTracer) +
     &         ptracerMask(I,J,iTracer,bi,bj) * ! ptracerMask is now a nx by ny by n_ptracers matrix
     &         volFlux(iceDepthK)

              ENDIF
             ENDDO
            ENDIF

C Calculate concentration of tracer in outflow: total tracer carried by
C the plume divided by the total outflow volume (the previous code divided
C by the volFluxDiff of the last outflow layer only, which is wrong when
C the outflow spans more than one layer)
            posSum = 0. _d 0
            DO K=kSrf(I,J),iceDepthK-1
             IF ( volFluxDiff(k) .GT. 0. _d 0 )
     &          posSum = posSum + volFluxDiff(k)
            ENDDO
            IF ( posSum .GT. 0. _d 0 ) THEN
             DO iTracer = 1,PTRACERS_num
              ptrPlume(iTracer) = ptrPlumeCum(iTracer) / posSum
             ENDDO
            ENDIF


#endif /* ALLOW_PTRACERS */

           ELSE ! ( Q_sg .EQ. 0 )

C If no subglacial output, then there is no plume
            DO k = 1,Nr
             rProfPlume(K) = 0.D0
             wProfPlume(K) = 0.D0
             tProfPlume(K) = 0.D0
             sProfPlume(K) = 0.D0
             aProfPlume(K) = 0.D0
             mIntProfPlume(K) = 0.D0
#ifdef ICEPLUME_ALLOW_DETACHED_PLUME 
             thetaProfPlume(k) = 0.D0
             distanceProfPlume(k) = 0.D0
#endif
            ENDDO
           ENDIF ! ( Q_sg.NE.0 ) THEN

C Send outputs to 3D grid for diagnostics
           IF ( usePlumeDiagnostics ) THEN
            DO K = 1,Nr
             rProfPlume3D(I,J,K) = rProfPlume(k)
C             IF (rProfPlume(k) .gt. 0) THEN
C             print *,'MWPRINT Setting plume area at ',I,J,K,
C     &         ' to',rProfPlume(k)
C             ENDIF
             wProfPlume3D(I,J,K) = wProfPlume(k)
C             IF (wProfPlume(k) .gt. 0) THEN
C               print *,'MWPRINT Setting plume w at ',I,J,K,
C     &         ' to',wProfPlume(k)
C               ENDIF
             tProfPlume3D(I,J,K) = tProfPlume(k)
C             IF (tProfPlume(k) .gt. 0) THEN
C               print *,'MWPRINT Setting plume t at ',I,J,K,
C     &         ' to',tProfPlume(k)
C               ENDIF
             sProfPlume3D(I,J,K) = sProfPlume(k)
C             IF (sProfPlume(k) .gt. 0) THEN
C               print *,'MWPRINT Setting plume s at ',I,J,K,
C     &         ' to',sProfPlume(k)
C               ENDIF
#ifdef ICEPLUME_ALLOW_DETACHED_PLUME 
             thetaProfPlume3D(I,J,K) = thetaProfPlume(k)
             distanceProfPlume3D(I,J,K) = distanceProfPlume(k)
#endif
            ENDDO
           ENDIF

C-------- Calculate melt rates ----------------------------

CMW           DO K = 1,Nr
             DO K=kSrf(I,J),Nr

C Check to see if we are above the sea bed
            IF ( K .GE. iceDepthK ) THEN
C If not then there is no melting             
             mProfAv(k) = 0.D0
             mProfPlume(k) = 0.D0
             mProf(k)      = 0.D0
            ELSE !k le iceDepthK

C If there is a plume in that cell, then need to calculate plume melt rate [m/time]
C distinct to background melt rate. Plume melt rate is already encorporated in 
C the plrume model, and taken into account in the temperature and salinity of the
C plume outflow. It is useful though to have it available as a diagnostic.
             plumeAreaInCell = 0.0

catn: note that mProfPlume below is not impacted by scaling of volFluxDiff
             IF ( ( Q_sg .NE. 0 ) .AND. 
     &        (useConePlume .OR. useSheetPlume.OR.useDetachPlume
     &         .OR. useBuoyPlume
     &        .OR. useTruncPlume) ) THEN
catn: aProfPlume(k): integrated contact area, calc in iceplume_plume_model [m2]
              plumeAreaInCell = aProfPlume(k) - aProfPlume(k+1)
catn: if the area above is bigger than below, implying plume area expands
              IF (plumeAreaInCell .gt. 0.0) then
catn: mIntProfPlume: integrated melt rate [m/day]
               mProfPlume(k) =(mIntProfPlume(k)-mIntProfPlume(k+1))/ 
     &                           plumeAreaInCell
     &                           * secInday
              ELSE
               mProfPlume (k) = 0.0
              ENDIF

             ELSE  !Q_sg eq 0
C If there is no plume in that cell, set plume melt rate to zero
              mProfPlume(k) = 0.D0
             ENDIF

C Calculate the background melt rate (i.e. not generated by plumes). This will 
C then be used to update the temperature and salinity in the adjacent cells.
C Velocities are calculated at cell faces - find averages for cell centres.
C Does not include velocity perpendicular to ice - this differs depending on 
C orientation of ice front

catn: Compute only one horz tangential vel based on ice front orientation:
catn: If North-South (East-West), we use vvel (uvel)
             IF ( GlacierFront_is_NorthSouth .EQ. 1 ) THEN
              hVelCell_tangential = (ABS(vVel(I,J,K,bi,bj))
     &                 +ABS(vVel(I,J+1,K,bi,bj))) / twoRL
             ELSEIF ( GlacierFront_is_NorthSouth .EQ. 0 ) THEN
              hVelCell_tangential = (ABS(uVel(I,J,K,bi,bj))
     &                 +ABS(uVel(I+1,J,K,bi,bj))) / twoRL
             ELSE
catn: for adjoint
              hVelCell_tangential = 0. _d 0
             ENDIF

             IF ( K .LT. Nr ) THEN
              wVelCell = (ABS(wVel(I,J,K,bi,bj))
     &                 +ABS(wVel(I,J,K+1,bi,bj))) / twoRL
             ELSE
              wVelCell = ABS(wVel(I,J,K,bi,bj)) / twoRL
             ENDIF

catn: From above:
C POSITIVE values indicate ice front is orientated north-south
C NEGATIVE values indicate ice front is orientated east-west
catn: So there is a discrepancy between above and below.  The importance
catn: is to settle on one convention and remain consistent.  For now, we
catn: stay with convention established above, now improved with the new 
catn: flag GlacierFront_is_NorthSouth
catn            IF ( plumeMask(I,J,bi,bj) .LT. 0 ) THEN
catn-C Negative mask values indicate east-west ice front orientation
catn             IF ( GlacierFront_is_Northsouth .EQ. 0 ) THEN
catn              meanVel = ((wVelCell**twoRL)+(uVelCell**twoRL))**HalfRL
catn-            ELSEIF ( plumeMask(I,J,bi,bj) .GT. 0 ) THEN
catn-C Positive mask values indicate north-south ice front orientation
catn             ELSEIF ( GlacierFront_is_NorthSouth .EQ. 1 ) THEN
catn              meanVel = ((wVelCell**twoRL)+(vVelCell**twoRL))**HalfRL
catn             ENDIF
catn: From Tom: what we want here is meanVel tangential to the ice front.
              meanVel = ((wVelCell**twoRL)
     &                  +(hVelCell_tangential**twoRL))**HalfRL

catn 30.Nov.2022: to make this function more general, let us compute
catn meltrate output in m/s , and convert to m/day outside as needed
             tmpT = 0. _d 0
             tmpS = 0. _d 0
ckiki: assuming that submarine melting forms little meltwater plumes
ckiki: that we do not resolve for computational reasons
ckiki: we can mimic the water mass transformation within these plumes
ckiki: and not give the ambient fjord temperature to the ice front
ckiki: ie the 3 equation model, but a reduced temperature of 
ckiki: Tprime =0.9Tfjord -0.2, could so same for salinity but thts not
ckiki: sensitive
             tModMWplume = slopeTmod*tProf(k)+interceptTmod
             CALL ICEPLUME_MELTRATE(
     I            tModMWplume,sProf(k),meanVel,rC(k),
     O            mProf(k) , tmpT, tmpS)
             mProf(k)=mProf(k)*secInday

C Get average melt rate. This is useful for visualling melt patterns and 
C assessing overall melt rate of glacier. Unit: [m/time]

C the following should apply to both conical and sheet plume models
             IF ( ( Q_sg .NE. 0 ) ) THEN
              plumeAreaInCell = aProfPlume(k) - aProfPlume(k+1)
C              print *,'plumeAreaInCell',plumeAreaInCell,'k=',k
catn: This is problematic here that delta_y is entering as if
catn: this geometry is strictly applied for only ice front aligning
catn: along u-direction such that the cross-sectional dlength is delta_y
catn: Now we change delta_y to dLtangential
              IF ( plumeAreaInCell .LE. dLtangential*delta_z(k) ) THEN
               IF ( plumeAreaInCell .LE. 0 ) THEN
C If there is no plume in cell, then the melt rate is equal to the background melt rate.
                mprofAv(k) = mProf(K)
               ELSE
C If there is a plume in cell, calculate average melt rate
catn: Again the problematic assumption of geometry with delta_y below;
catn: We now change delta_y to dLtangential
                mProfAv(k) = (mProfPlume(k)*plumeAreaInCell
     &              +mProf(k)*(dLtangential*delta_z(k)-plumeAreaInCell))
     &              / (dLtangential * delta_z(k))
C          IF ((k.eq.50)) THEN
C          print *,'MWPRINT Plume melt should be ',
C     &         'set to',mProfAv(k),'at k=50'
C          ENDIF
C Scale down background melt rate to account for area occupied by plume
C (necessary so that tendency terms are not over estimated)
                mProf(k) = mProf(k)*(1-plumeAreaInCell/
     &                  (dLtangential*delta_z(k)))
               ENDIF
              ELSE
C if the plume contact area is larger than the cell area, we assume there is
C no background melting
               mProfAv(k) = mProfPlume(k)*plumeAreaInCell / 
     &                     (dLtangential*delta_z(k))
               mProf(k) = 0.
              ENDIF
             ELSE ! Q_sg eq 0, not coneplume or sheet plume

C If it is not a plume cell, then no plume melting.
              mProfPlume(k) = 0. _d 0
              mProfAv(k) = mProf(k)
C              IF ((k.eq.50)) THEN
C              PRINT *, 'Setting to 0 for some reason'
C              ENDIF
             ENDIF ! Q_sg
            ENDIF ! are we above sea bed?

C            IF ((k.eq.50).and.(i.eq.2).and.(j.eq.17)) THEN
C               print *,'MWPRINT Plume melt at k=50 is ',
C     &         ' to',mProfPlume(k)
C               ENDIF

C Send outputs to 3D grid for diagnostics
catn: Since we are saving these 3D diags, should know how to calc offline
catn: their contributions to addMass3D(plume,bg)
            IF ( usePlumeDiagnostics ) THEN
             mProfPlume3D(I,J,K) = mProfPlume(k)
             mProfAv3D(I,J,K) = mProfAv(k)
C             IF ((k.eq.50).and.(j.eq.17).and.(i.eq.2)) THEN
C               print *,'MWPRINT Setting plume melt at ',I,J,K,
C     &         ' to',mProfPlume(k),mProfPlume3D(I,J,K)
C               ENDIF
            ENDIF

C          IF ((k.eq.50).and.(i.eq.2).and.(j.eq.17)) THEN
C          print *,'MWPRINT Plume melt check 1 at ',2,17,50,
C     &         ' to',mProfPlume3D(2,17,50)
C          ENDIF

C ------------Tendencies------------------------------
C These are applied in cells where there is no plume. The idea is that in these areas there are likely to be 
C local subgrid convection cells. As such, it is most realistic to apply changes in T and S to ice edge cell. 
C Where there is a plume, products of melt are transported upwards so no local changes applied.
C The thermodynamics in this section are taken from pkg/icefront
C Tendencies are applied in S/R EXTERNAL_FORCING

C To convert from melt rate (m d^-1) to freshwater flux (kg m^-2 d^-1)
catn: m/day * kg/m3 = kg / day / m^2. Where does 94.22 come from? I also think the 
catn: quoted unit above is incorrect and fwflux is in kg/m2/s, same as in icefront.  
catn: rhoshelfice/secInday = 917./(24*3600) = 0.0106 = 1/94.22
catn: Sign convention from pkg/icefront: inward (negative, into the ocean) fwFlux 
catn: implies glacier melting due to outward (positive, leaving the ocean) heatflux
catn: The last step we need to take care of is that FwFlux(k) is actually computed 
catn: using rhoShelfIce, but outside of this addMass* is converted to Eta using 
catn: rUnit2mass = rhoConst.  How to reconcile this?  Do we mult FwFlux back by 
catn: rhoShelfIce to convert to m3/s, then rUnit2mass? If we look at how eccov4r5 
catn: pkg/shelfice is treating fwflux from background melt, they use rUnit2mass and 
catn: NOT rhoShelfIce.  In light of this, we REMOVE rhoShelfIce from the calculation 
catn: of FWflux and use rUnit2Mass=rhoConst.
catn            FwFlux(k) = -mProf(k)/94.22
catn            FwFlux(k) = -mProf(k)*rhoShelfIce/secInday
catn: this is using mProf, so it is background melt only
            FwFlux(k) = -mProf(k)*rUnit2Mass/secInday
catn  Now: m/day*day/sec*kg/m3=kg/m2/s
catn: Lastly, to convert this FwFlux (kg/m2/s) to addMass3Dbg [kg/s], mult by the vertical
catn: area (dLtangential*dz).  Note we have already taken care at this point to ensure 
catn: dLtangential is along the glacier front.  In addition, we have already taken care to
catn: reduce this addmass3dbg to zero if the whole cell is plume.
catn      addMass3Dbg(I,j,k,bi,bj)=Fwflux3D(I,J,k,bi,bj)*dLtangential*delta_z(k)*rUnit2mass 
catn: pay attention to sign! addMass is minus FwFlux!! (following pkg/shelfice of eccov4r5)
catn: kg/m2/s*m2=kg/s
            addMass3Dbg(I,J,k,bi,bj)=-FwFlux(k)*dLtangential*delta_z(k)

C Heat required to melt that much ice (W m^-2)
catn: 333.55 x 10^3 is the Latent heat of fusion: When melting 1kg of ice, 333.55 kJ of 
catn: energy is absorbed with no temperature change.  While melting, the heat energy 
catn: needed to change a substance from solid to liquid at atm pressure is the latent 
catn: heat of fusion. The liquid phase has a higher internal energy than the solid phase.
catn: From statement above, in an ocean-ice system, when there is an dheff ice melt, implying 
catn: -dheff*rhoi/rhosw=fwflux (negative) to the ocean, the ocean first needs to put energy 
catn: into raising the ice temperature to freezing, then put more energy into solid->liquid 
catn: phase change. Thus, from the ocean perspective, it loses dheat to raise ice temp plus 
catn: dheat_Latent convert solid->liquid.  So a +dheff gives -fwflx (neg, inward, into ocean),
catn: the ice system gains energy and the ocean loses energy and cools down 
catn: (outward, postive=leaving ocean). dheat_lost_by_ocean is opposite sign of fwflux. 
catn: kg/s/m2*J/kg=J/s/m2=W/m2
catn            heatflux(k) = -FwFlux(k)*333.55E+3
catn This is based on mprof, so it is background fluxes
            heatflux(k) = -FwFlux(k)*ICEPLUMElatentHeat

C Create local (no overlap) arrays of heat and freshwater flux from background melting
catn FwFlux kg/s/m2 , HeatFlux: W/m2
catn: More problematic assumption about I vs J and geometry of glacier front
catn  Now addressed with GlacierFront_is_NorthSouth flag
            IF (GlacierFront_is_NorthSouth .eq. 1) THEN
              IF ( ( J .GT. 0 ) .AND. ( J .LT. sNy+1 ) ) THEN
               Fwflux3D(I,J,k) = FwFlux(k)
               HeatFlux3D(I,J,k) = HeatFlux(k)
              ENDIF
            ELSEIF (GlacierFront_is_NorthSouth .eq. 0) THEN
              IF ( ( I .GT. 0 ) .AND. ( I .LT. sNx+1 ) ) THEN
               Fwflux3D(I,J,k) = FwFlux(k)
               HeatFlux3D(I,J,k) = HeatFlux(k)
              ENDIF
            ENDIF

C     Compute tendencies (as for pkg/icefront)
catn: In ini_parms.F: mass2rUnit = recip_rhoConst, e.g., 1/999.8=1.0002e-3, HeatCapacity_Cp=3994
catn: Sign convention: using pkg/icefront as a guide: ocean will take the negative of the 
catn: heatflux and positive of fwflux above, i.e., opposite sign of what was computed above.  
catn: J/s/m2 * degC.m3/J = degC/s.m ; kg/s/m2 m3/kg g/kg = g/kg/s.m
            iceplumeBG_TendT(I,J,K,bi,bj) = - HeatFlux3D(I,J,K)*eps5
            iceplumeBG_TendS(I,J,K,bi,bj) = FWFlux3D(I,J,K)*
     &                mass2rUnit * sProf(k)
C     Scale by icefrontlength, which is the ratio of the horizontal length
C     of the ice front in each model grid cell divided by the grid cell area.
C     (icefrontlength = dy / dxdy = 1 / dx)
catn: More problematic assumption of geometry of delta_x and delta_y
catn: Solved now by rename to dLtangential and dLnormal
catn: Scale factor: dLtangential/(dLnormal*dLtantengial) = 1/dLnormal
catn degC/s.m/m = degC/s; g/kg/s.m/m=g/kg/s
            iceplumeBG_TendT(I,j,K,bi,bj) = 
     &            iceplumeBG_TendT(I,j,K,bi,bj)
     &            * 1./dLnormal
            iceplumeBG_TendS(I,j,K,bi,bj) = 
     &            iceplumeBG_TendS(I,j,K,bi,bj)
     &            * 1./dLnormal
           ENDDO !the big k-loop

catn: Note: currently we are not taking care of the mass added from FwFlux (background melt).
catn: To fully take into account all fw added to the ocean, we will need to either add it
catn: to addMass array, or create yet another 3D field specific to iceplume where we store
catn: it there.  The latter, though not desirable as we are carrying more 3d fields around,
catn: is very good for diagnostic purpose, as we can for sure compute the ratio of 
catn: background melt vs plume-induced entrainment vs Qsg

C The plume transport is undertaken using the addMass terms.
C addMass terms for volume (kg/s), pot. temperature, salinity and ptracers are input 
C into the correct locations in 3D arrays
catn: Need to first compute addMass specific to plume and store it independently IFF we
catn: want to use separate processes to treat Tendencies at different places in the code
C No plume source/sink above the top wet level (inside the ice shelf)
           DO K=1,kSrf(I,J)-1
            temp_AddMass3Dplume(I,J,K,bi,bj) = 0. _d 0
            salt_AddMass3Dplume(I,J,K,bi,bj) = 0. _d 0
            addMass3Dplume(I,J,K,bi,bj) = 0. _d 0
#ifdef ALLOW_PTRACERS
            DO iTracer = 1,PTRACERS_num
             ptr_AddMass3D(I,J,K,bi,bj,iTracer) = 0. _d 0
            ENDDO
#endif /* ALLOW_PTRACERS */
           ENDDO
           IF ( Q_sg.NE.0. _d 0 ) THEN
C Find temperature and salinity of plume outflow
CMW             DO k = 1,Nr        !second kloop
             DO K=kSrf(I,J),Nr
catn: volFluxDiff entering here: when there is net convergence at this k
               IF ( volFluxDiff(k) .GT. 0. _d 0 ) THEN
                IF ( conserveMass ) THEN
                 temp_AddMass3Dplume(I,J,K,bi,bj) = ! convert to pot. temp
     &            SW_PTMP(sProfPlume(k),tProfPlume(k),prProf(k),0. _d 0)
                 salt_AddMass3Dplume(I,J,K,bi,bj) = sProfPlume(k)
                ELSE
                 kp = MIN( k+1, iceDepthK )
                 prFace = ABS(rF(kp))*rhoConst*gravity*1.0 _d -4
                 temp_AddMass3Dplume(I,J,K,bi,bj) =
     &            SW_PTMP(sOut,tProfPlume(kp),prFace,0. _d 0)
                 salt_AddMass3Dplume(I,J,K,bi,bj) = sOut
                ENDIF
#ifdef ALLOW_PTRACERS
                DO iTracer = 1,PTRACERS_num
                 ptr_AddMass3D(I,J,K,bi,bj,iTracer)
     &                             = ptrPlume(iTracer)
                ENDDO
#endif /* ALLOW_PTRACERS */
               ELSE 
                temp_AddMass3Dplume(I,J,K,bi,bj) = ptProf(k)
                salt_AddMass3Dplume(I,J,K,bi,bj) = sProf(k)
#ifdef ALLOW_PTRACERS
                DO iTracer = 1,PTRACERS_num
                 ptr_AddMass3D(I,J,K,bi,bj,iTracer) 
     &                             = ptrProf(k,iTracer)
                ENDDO
#endif /* ALLOW_PTRACERS */
               ENDIF

catn: Now also filling out array for background melt, although not needed for 
catn: tendency because we already take care of it in here with tend[T,S]

C Convert m3/s into kg/s
catn: volFluxDiff is now converted to addMass, so already taken care of massflux into ocean at depth
catn: Need to avoid hardcoded 1000 (missing decimal), will use rhoConstFresh
catn: Note: in integr_continuity, we convert addMass to volume using rhoConst=1029, so we will need to 
catn: use that here to avoid a difference in "mass" effect between what is computed in here and outside
catn             addMass(I,j,k,bi,bj) = volFLuxDiff(k)*1000 
catn              addMass(I,j,k,bi,bj) = volFLuxDiff(k)*rhoConstFresh 
catn: rUnit2mass = rhoConst, this is how mitgcm converts between volume and mass
              addMass3Dplume(I,j,k,bi,bj) = volFLuxDiff(k)*rUnit2mass
             ENDDO      !kloop

           ELSE !Q_sg eq 0

C             DO K = 1,Nr
             DO K=kSrf(I,J),Nr
              temp_AddMass3Dplume(I,J,K,bi,bj) = ptProf(k)
              salt_AddMass3Dplume(I,J,K,bi,bj) = sProf(k)
#ifdef ALLOW_PTRACERS
              DO iTracer = 1,PTRACERS_num
               ptr_AddMass3D(I,J,K,bi,bj,iTracer) = ptrProf(k,iTracer)
              ENDDO
#endif /* ALLOW_PTRACERS */
catn: Cannot reset a global field in here, now change to process-specfic field
catn              addMass(I,j,k,bi,bj) = 0.D0
              addMass3Dplume(I,j,k,bi,bj) = 0. _d 0
             ENDDO

           ENDIF

         ENDIF ! plumeMask .NE. 0

C DO J loop
         ENDDO

C DO I loop
         ENDDO

catn: here we calc fields that live outside this routine
#ifdef ALLOW_ADDFLUID
C Background (non-plume) melt: its freshwater mass belongs with the
C freshwater (salinity) tendency. If applyIcePlumeBGTendS = .FALSE. the
C meltwater is not added at all (previously the mass was still added,
C bringing in water at the local salinity, i.e. creating salt).
         IF ( .NOT.applyIcePlumeBGTendS ) THEN
          DO K=1,Nr
           DO J = 1-Oly,sNy+Oly
            DO I = 1-Olx,sNx+Olx
             addMass3Dbg(I,J,K,bi,bj) = 0. _d 0
            ENDDO
           ENDDO
          ENDDO
         ENDIF
         DO K=1,Nr
          DO J = 1,sNy
           DO I = 1,sNx
             addMass(I,J,K,bi,bj) = addMass(I,J,K,bi,bj)
     &         + addMass3Dplume(I,J,K,bi,bj) + addMass3Dbg(I,J,K,bi,bj)
           ENDDO
          ENDDO
         ENDDO
#endif
C-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

C       print *,'MWPRINT Plume melt check 2 at ',2,17,50,
C    &         ' to',mProfPlume3D(2,17,50,bi,bj)

#ifdef ALLOW_DIAGNOSTICS_VEC

      DO k=1,Nr
       DO j=1-Oly,sNy+Oly
        DO i=1-Olx,sNx+Olx
          icefrntA_dv(i,j,k,bi,bj)= mProfAv3D(i,j,k)
          icefrntM_dv(i,j,k,bi,bj)= mProfPlume3D(i,j,k)
          icefrntR_dv(i,j,k,bi,bj)= rProfPlume3D(i,j,k)
          icefrntS_dv(i,j,k,bi,bj)= sProfPlume3D(i,j,k)
          icefrntT_dv(i,j,k,bi,bj)= tProfPlume3D(i,j,k)
          icefrntW_dv(i,j,k,bi,bj)= wProfPlume3D(i,j,k)
        ENDDO
       ENDDO
      ENDDO

#endif /* DIAGNOSTICS_VEC */

C Save plume values
#ifdef ALLOW_DIAGNOSTICS

         IF ( useDiagnostics .AND. usePlumeDiagnostics ) THEN

C Transfer to local (no bi,bj indices) and interior only arrays
C (can not seem to get the diagnostics to output properly with the
C other options!)
        DO K = 1,Nr
         DO J = 1-Oly,sNy+Oly
          DO I = 1-Olx,sNx+Olx
            gT_iceplumeBG(i,j,k)=iceplumeBG_tendT(i,j,k,bi,bj)
            gS_iceplumeBG(i,j,k)=iceplumeBG_tendS(i,j,k,bi,bj)
           ENDDO !I
          ENDDO  !J
         ENDDO   !K

C Output diagnostics
          DO K=1,Nr
C Here, the first k is the layer in the output field in which to save the data
C and the second k is layer from which this data is taken in the original field
c flag 2 for with overlap, 3 for without overlap; 1--3: bi,bj done outside
          CALL DIAGNOSTICS_FILL_RS(wProfPlume3D,'icefrntW',
     &         k,k,2,bi,bj,myThid)
          CALL DIAGNOSTICS_FILL_RS(tProfPlume3D,'icefrntT',
     &         k,k,2,bi,bj,myThid)
          CALL DIAGNOSTICS_FILL_RS(sProfPlume3D,'icefrntS',
     &         k,k,2,bi,bj,myThid)
          CALL DIAGNOSTICS_FILL_RS(rProfPlume3D,'icefrntR',
     &         k,k,2,bi,bj,myThid)
C          print *,'MWPRINT Outputting mProfPlume3D at ',2,17,k,
C     &         ' to',mProfPlume3D(2,17,k)
          CALL DIAGNOSTICS_FILL_RS(mProfPlume3D,'icefrntM',
     &         k,k,2,bi,bj,myThid)
          CALL DIAGNOSTICS_FILL_RS(mProfAv3D,'icefrntA',
     &         k,k,2,bi,bj,myThid)
#ifdef ICEPLUME_ALLOW_DETACHED_PLUME
          CALL DIAGNOSTICS_FILL_RS(thetaProfPlume3D,'PlumAngl',
     &         k,k,2,bi,bj,myThid)
          CALL DIAGNOSTICS_FILL_RS(distanceProfPlume3D,'PlumDist',
     &         k,k,2,bi,bj,myThid)
#endif
          CALL DIAGNOSTICS_FILL_RS(gT_iceplumeBG,'IP_gTbg ',
     &         k,k,2,bi,bj,myThid)
          CALL DIAGNOSTICS_FILL_RS(gS_iceplumeBG,'IP_gSbg ',
     &         k,k,2,bi,bj,myThid)
          CALL DIAGNOSTICS_FILL(addMass3Dplume,'IPmasspl',
     &         k,k,1,bi,bj,myThid)
          CALL DIAGNOSTICS_FILL(addMass3Dbg,'IPmassbg',
     &         k,k,1,bi,bj,myThid)
          ENDDO

catn since none of below needs bi,bj, they need reset here.
C Clear local arrays otherwise results replicate on other tiles
        DO K = 1,Nr
         DO J = 1-Oly,sNy+Oly
          DO I = 1-Olx,sNx+Olx
C         DO J = 1,sNy
C          DO I = 1,sNx
           wProfPlume3d(I,J,K) = 0.d0
           tProfPlume3d(I,J,K) = 0.d0
           sProfPlume3d(I,J,K) = 0.d0
           rProfPlume3d(I,J,K) = 0.d0
           mProfPlume3d(I,J,K) = 0.d0
           mProfAv3d(I,J,K) = 0.d0
#ifdef ICEPLUME_ALLOW_DETACHED_PLUME
           thetaProfPlume3d(I,J,K) = 0.d0
           distanceProfPlume3d(I,J,K) = 0.d0
#endif
           gT_iceplumeBG(i,j,k)=0. _d 0
           gS_iceplumeBG(i,j,k)=0. _d 0
           rProfPlume3D(i,j,k)=0. _d 0
           wProfPlume3D(i,j,k)=0. _d 0
           tProfPlume3D(i,j,k)=0. _d 0
           sProfPlume3D(i,j,k)=0. _d 0
           mProfPlume3D(i,j,k)=0. _d 0
           mProfAv3D(i,j,k)=0. _d 0
#ifdef ICEPLUME_ALLOW_DETACHED_PLUME 
           thetaProfPlume3D(i,j,k)=0. _d 0
           distanceProfPlume3D(i,j,k)=0. _d 0
#endif /* ICEPLUME_ALLOW_DETACHED_PLUME */

          ENDDO
         ENDDO
        ENDDO

         ENDIF

#endif /* ALLOW_DIAGNOSTICS */ 

C     end bi/bj-loops
       ENDDO
      ENDDO

C Plume columns are only solved in the tile interior, so addMass was left
C without the plume sources in the overlaps; the core uses addMass in the
C overlaps (e.g. mom_fluxform.F, addMass(i-1)/addMass(j-1)), so a plume in
C the first/last row or column of a tile made the velocities on the tile
C edge inconsistent and the MPI run lost volume (-3.2 m3/s of a 61 m3/s
C source in a 4-rank test). Exchange addMass after adding the sources.
C (The plume T/S/ptracer source fields are only used in the interior;
C exchanging them changed nothing, bit for bit, so they are not.)
#ifdef ALLOW_ADDFLUID
      CALL EXCH_XYZ_RL( addMass, myThid )
#endif

#endif /* ALLOW_ICEPLUME */

      RETURN
      END



