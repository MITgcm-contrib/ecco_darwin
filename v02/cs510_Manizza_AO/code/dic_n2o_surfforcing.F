#include "CPP_OPTIONS.h"
#include "PTRACERS_OPTIONS.h"
#include "DARWIN_OPTIONS.h"

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_DARWIN

#ifdef ALLOW_CARBON

CBOP
C !ROUTINE: O2_SURFFORCING

C !INTERFACE: ==========================================================
      SUBROUTINE DIC_N2O_SURFFORCING( PTR_N2O, SGN2O,
     I           bi,bj,iMin,iMax,jMin,jMax,
     I           myIter, myTime, myThid )

C !DESCRIPTION:
C Calculate the oxygen air-sea flux terms               

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "DYNVARS.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "FFIELDS.h"
#include "DARWIN_FLUX.h"

c  !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
C  myIter               :: current timestep
C  myTime               :: current time
C  PTR_O2               :: oxygen tracer field
      _RL myTime
      _RL  PTR_N2O(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      INTEGER iMin,iMax,jMin,jMax, bi, bj
      INTEGER myIter, myThid

c  !OUTPUT PARAMETERS: ===================================================
C  SGO2                  :: air-sea exchange of oxygen
      _RL  SGN2O(1-OLx:sNx+OLx,1-OLy:sNy+OLy)

C !LOCAL VARIABLES: ===================================================
C I, J, K - Loop counters
      INTEGER I,J,K
C Solubility relation coefficients
      _RL SchmidtNoN2O(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL N2Osat(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL Kwexch(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL FlxN2O(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL AtmosN2O(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL  aTT 
      _RL  aTK
      _RL  aTS
      _RL  aTS2
      _RL  aTS3
      _RL  aTS4
      _RL  aTS5
      _RL  o2s
      _RL  ttemp
      _RL  stemp
      _RL  oCnew
CEOP

C
C NEW COEFFICIENTS for N2O calculation


      _RL a1,a2,a3,b1,b2,b3
      _RL part1, part2, part3,sal,tk
      _RL soln2o,pN2Oatm

         a1 = -62.7062
         a2 =  97.3066
         a3 =  24.1406
         b1 = -0.05842
         b2 =  0.033193
         b3 = -0.0051313

C Setting atmospheric pressure of N2O
C units are in ppb or nanoatm
C example 320 ppb = 320E-9

         pN2Oatm = 320E-9


      K=1

C calculate SCHMIDT NO. for O2
        DO j=jmin,jmax
          DO i=imin,imax
            IF (maskC(i,j,k,bi,bj).NE.0.) THEN
              ttemp = theta(i,j,k,bi,bj)
              stemp = salt(i,j,k,bi,bj)

              SchmidtNoN2O(i,j) = 
     &            sn2o1 
     &          + sn2o2 * ttemp
     &          + sn2o3 * ttemp*ttemp  
     &          + sn2o4 * ttemp*ttemp*ttemp

C Determine surface flux of O2 
C exchange coeff accounting for ice cover and Schmidt no.
C Kwexch_Pre= pisvel*(1-fice): previously computed in dic_surfforcing.F

              Kwexch(i,j) = Kwexch_Pre(i,j,bi,bj)
     &                    / sqrt(SchmidtNoN2O(i,j)/660.0 _d 0)


C
C Convert celsius into kelvin
C
             tk = ttemp+273.15


C     MM Coefficients of solubilty factors of N2O

            soln2o=(exp(a1+a2*(100/tk)+a3*log(tk/100)
     &     +stemp*(b1+b2*(tk/100)+b3*((tk/100)**2))))


C convert into mol/m3/atm
             soln2o = soln2o*1E3


c Determine flux, inc. correction for local atmos surface pressure
cQQ PTR_O2?
C The flux is calculated in mmol N2O

              FlxN2O(i,j) = Kwexch(i,j)*
     &                     ((pN2Oatm*soln2o*1E3)
     &                      - PTR_N2O(i,j,1))


C              write(*,*),'===================================='
C              write(*,*),'FluxN2O GAS EXCH. = ', FluxN2O(i,j)
C              write(*,*),'pN2Oatm_N2O = ',pN2Oatm
C              write(*,*),'PTR_N2O = ',PTR_N2O(i,j,1)
C              write(*,*),'pN2Oatm*soln2o*1E3 = ',pN2Oatm*soln2o*1E3
C              write(*,*),'====================================' 




            ELSE
              FlxN2O(i,j) = 0. _d 0
            ENDIF


          END DO
        END DO

C update surface tendencies      
        DO j=jmin,jmax
          DO i=imin,imax
           SGN2O(i,j)= FlxN2O(i,j)
     &         *recip_drF(K) * recip_hFacC(i,j,K,bi,bj)
           
C               write(*,*),'SGNO2 GAS EXCH. = ', SGN2O(i,j)

 
          ENDDO
         ENDDO


        RETURN
        END
#endif  /*ALLOW_CARBON*/

#endif  /*DARWIN*/
#endif  /*ALLOW_PTRACERS*/
c ==================================================================


