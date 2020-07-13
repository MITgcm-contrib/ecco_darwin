#include "CPP_OPTIONS.h"
#include "PTRACERS_OPTIONS.h"
#include "DARWIN_OPTIONS.h"

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_DARWIN

#ifdef ALLOW_CARBON

CBOP
C !ROUTINE: DIC_BUDGETSALT

C !INTERFACE: ==========================================================
      SUBROUTINE DIC_BUDGETSALT( PTR_DIC , PTR_ALK, PTR_PO4, PTR_SIL,
     O           deltaSalt,
     I           bi,bj,imin,imax,jmin,jmax,
     I           myIter,myTime,myThid)

C !DESCRIPTION:
C  Calculate the carbon air-sea flux terms              
C  following external_forcing_dic.F (OCMIP run) from Mick            

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "DYNVARS.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "FFIELDS.h"
#include "DARWIN_SIZE.h"
#include "DARWIN_IO.h"
#include "DARWIN_FLUX.h"
#ifdef USE_EXFWIND
#include "EXF_FIELDS.h"
#endif

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
C  myIter               :: current timestep
C  myTime               :: current time
c  PTR_DIC              :: DIC tracer field
      INTEGER myIter, myThid
      _RL myTime
      _RL  PTR_DIC(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL  PTR_ALK(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL  PTR_PO4(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL  PTR_SIL(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER iMin,iMax,jMin,jMax, bi, bj

C !LOCAL VARIABLES: ====================================================
       INTEGER I,J, kLev, it
C Number of iterations for pCO2 solvers...
C Solubility relation coefficients
      _RL SchmidtNoDIC(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL pCO2sat(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL Kwexch(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL pisvel(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
C local variables for carbon chem
      _RL surfdic(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL surfalk(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL surfphos(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL surfsi(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL surfsalt(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL surftemp(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL budgetSaltPert(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
#ifdef ALLOW_OLD_VIRTUALFLUX
      _RL VirtualFlux(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
#endif
C local variables for CO2_FLUX_BUDGET
      _RL deltaSalt(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
CEOP

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      kLev=1

cc if coupled to atmsopheric model, use the
cc Co2 value passed from the coupler
c#ifndef USE_ATMOSCO2
cC PRE-INDUSTRIAL STEADY STATE pCO2 = 278.0 ppmv
c       DO j=1-OLy,sNy+OLy
c        DO i=1-OLx,sNx+OLx
c           AtmospCO2(i,j,bi,bj)=278.0 _d -6
c        ENDDO
c       ENDDO
c#endif
C =================================================================
C determine inorganic carbon chem coefficients
        DO j=jmin,jmax
         DO i=imin,imax
c put bounds on tracers so pH solver doesn't blow up
             surfdic(i,j) = 
     &           max(100. _d 0 , min(4000. _d 0, PTR_DIC(i,j)))*1e-3
     &                          * maskC(i,j,kLev,bi,bj)
             surfalk(i,j) = 
     &           max(100. _d 0 , min(4000. _d 0, PTR_ALK(i,j)))*1e-3
     &                          * maskC(i,j,kLev,bi,bj)
             surfphos(i,j)  = 
     &           max(1. _d -10, min(10. _d 0, PTR_PO4(i,j)))*1e-3
     &                          * maskC(i,j,kLev,bi,bj)
             surfsi(i,j)   = 
     &           max(1. _d -8, min(500. _d 0, PTR_SIL(i,j)))*1e-3
     &                          * maskC(i,j,kLev,bi,bj)
             surfsalt(i,j) = 
     &           max(4. _d 0, min(50. _d 0, salt(i,j,kLev,bi,bj)))
             surftemp(i,j) = 
     &           max(-4. _d 0, min(39. _d 0, theta(i,j,kLev,bi,bj)))
             budgetSaltPert(i,j) = surfsalt(i,j) + budgetPert
         ENDDO
        ENDDO

         CALL CARBON_COEFFS(
     I                       surftemp,budgetSaltPert,
     I                       bi,bj,iMin,iMax,jMin,jMax,myThid)
C====================================================================

       DO j=jmin,jmax
        DO i=imin,imax
C Compute AtmosP and Kwexch_Pre which are re-used for flux of O2

#ifdef USE_PLOAD
C Convert anomalous pressure pLoad (in Pa) from atmospheric model
C to total pressure (in Atm)
C Note: it is assumed the reference atmospheric pressure is 1Atm=1013mb
C       rather than the actual ref. pressure from Atm. model so that on
C       average AtmosP is about 1 Atm.
                AtmosP(i,j,bi,bj)= 1. _d 0 + pLoad(i,j,bi,bj)/Pa2Atm
#endif

C Pre-compute part of exchange coefficient: pisvel*(1-fice)
C Schmidt number is accounted for later
#ifdef USE_EXFWIND
                pisvel(i,j)=0.337 _d 0 *wspeed(i,j,bi,bj)**2/3.6 _d 5
cBX linear piston velocity after Krakauer et al. (2006), Eq. 3
cBX  using <k> = 20, n=0.5, and <u^n> = 2.6747 (as determined from 2010
cBX  EXFwspee field from cube92 run)
cDc              pisvel(i,j)=20 _d 0 *(wspeed(i,j,bi,bj)**0.5
cDc     &              /2.6747 _d 0) /3.6 _d 5
#else
              pisvel(i,j)=0.337 _d 0 *wind(i,j,bi,bj)**2/3.6 _d 5
#endif
              Kwexch_Pre(i,j,bi,bj) = pisvel(i,j)
     &                              * (1. _d 0 - FIce(i,j,bi,bj))

        ENDDO
       ENDDO

c pCO2 solver...
C$TAF LOOP = parallel
       DO j=jmin,jmax
C$TAF LOOP = parallel
        DO i=imin,imax

          IF ( maskC(i,j,kLev,bi,bj).NE.0. _d 0 ) THEN
            CALL CALC_PCO2_APPROX(
     I        surftemp(i,j),budgetSaltPert(i,j),
     I        surfdic(i,j), surfphos(i,j),
     I        surfsi(i,j),surfalk(i,j),
     I        ak1(i,j,bi,bj),ak2(i,j,bi,bj),
     I        ak1p(i,j,bi,bj),ak2p(i,j,bi,bj),ak3p(i,j,bi,bj),
     I        aks(i,j,bi,bj),akb(i,j,bi,bj),akw(i,j,bi,bj),
     I        aksi(i,j,bi,bj),akf(i,j,bi,bj),
     I        ak0(i,j,bi,bj), fugf(i,j,bi,bj),
     I        ff(i,j,bi,bj),
     I        bt(i,j,bi,bj),st(i,j,bi,bj),ft(i,j,bi,bj),
     U        pH(i,j,bi,bj),pCO2(i,j,bi,bj),CO3(i,j,bi,bj),
     I        myThid )
          ELSE
              pH(i,j,bi,bj) = 0. _d 0
 	      pCO2(i,j,bi,bj) = 0. _d 0
	      CO3(i,j,bi,bj) = 0. _d 0 
          ENDIF
        ENDDO
       ENDDO

C update tendency      
         DO j=jmin,jmax
          DO i=imin,imax
           if(budgetTStep1.EQ.0) then
C if first timestep
C this is problematic at restart; clean-up later
            deltaSalt(i,j) = 0. _d 0
            previousSalt(i,j,bi,bj) = surfsalt(i,j)
           else
C current value - value from last timestep
            deltaSalt(i,j) = surfsalt(i,j) - previousSalt(i,j,bi,bj)
            previousSalt(i,j,bi,bj) = surfsalt(i,j)
           endif
          ENDDO
         ENDDO

        RETURN
        END
#endif  /*ALLOW_CARBON*/

#endif  /*DARWIN*/
#endif  /*ALLOW_PTRACERS*/
c ==================================================================
