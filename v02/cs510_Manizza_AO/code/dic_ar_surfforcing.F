#include "CPP_OPTIONS.h"
#include "PTRACERS_OPTIONS.h"
#include "DARWIN_OPTIONS.h"

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_DARWIN

#ifdef ALLOW_CARBON

CBOP
C !ROUTINE: AR_SURFFORCING

C !INTERFACE: ==========================================================
      SUBROUTINE DIC_AR_SURFFORCING( PTR_AR, SGAR,
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
      _RL  PTR_AR(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      INTEGER iMin,iMax,jMin,jMax, bi, bj
      INTEGER myIter, myThid

c  !OUTPUT PARAMETERS: ===================================================
C  SGAR                  :: air-sea exchange ofargon 
      _RL  SGAR(1-OLx:sNx+OLx,1-OLy:sNy+OLy)

C !LOCAL VARIABLES: ===================================================
C I, J, K - Loop counters
      INTEGER I,J,K
C Solubility relation coefficients
      _RL SchmidtNoAr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL ARsat(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL Kwexch(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL FlxAR(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
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
      _RL tempS

C NEW VARIABLES
C FOR ARGON 

      _RL rhosw
      _RL sumAr
      _RL Arsol

C  Density of seawater (kg/m3)

      rhosw = 1024.5


CEOP


      K=1

C calculate SCHMIDT NO. for Argon 
        DO j=jmin,jmax
          DO i=imin,imax
            IF (maskC(i,j,k,bi,bj).NE.0.) THEN

             ttemp = theta(i,j,k,bi,bj)
             stemp = salt(i,j,k,bi,bj)





C              write(*,*),'TEMP ARGON',  ttemp
C              write(*,*),'SALT ARGON', stemp
c              write(*,*),'sar1 = ',sar1
C              write(*,*),'sar2 = ',sar2
C              write(*,*),'sar3 = ',sar3
C              write(*,*),'sar4 = ',sar4
 
              

              SchmidtNoAr(i,j) = 
     &            sar1 
     &          + sar2 * ttemp
     &          + sar3 * ttemp*ttemp  
     &          + sar4 * ttemp*ttemp*ttemp


C              write(*,*),' SchmidtNoAr = ',  SchmidtNoAr(i,j)

C Determine surface flux of Argon 
C exchange coeff accounting for ice cover and Schmidt no.
C Kwexch_Pre= pisvel*(1-fice): previously computed in dic_surfforcing.F

              Kwexch(i,j) = Kwexch_Pre(i,j,bi,bj)
     &                    / sqrt(SchmidtNoAr(i,j)/660.0 _d 0)


C convert T to scaled temperature

            tempS = log((298.15 - ttemp)/(273.15 + ttemp));


C determine saturation Ar 
C using Hamme and Emerson


          sumAr=  A0ar+(A1ar*tempS)+(A2ar*tempS**2)+
     &   (A3ar*(tempS**3))+stemp*(B0ar+(B1ar*tempS)+
     &   (B2ar*(tempS**2)))

C          write(*,*), 'sumAr = ', sumAr


C  Solubilty factor computed in micromol/Kg
C  and conversion into mmol/m3

           Arsol=EXP(sumAr)*(rhosw/1e3)
C conversion into mol/m3
C           Arsol=Arsol/1e3

      Arsat(i,j)=Arsol



C Determine flux, inc. correction for local atmos surface pressure
              FlxAR(i,j) = Kwexch(i,j)*
     &                     (AtmosP(i,j,bi,bj)*Arsat(i,j) 
     &                      - PTR_AR(i,j,K)) 


C             write(*,*),'*****************************'
c             write(*,*),'PTR_AR = ',PTR_AR(i,j,K)
c             write(*,*),'Arsat = ',Arsat(i,j)
c             write(*,*),'AtmosP  = ', AtmosP(i,j,bi,bj)
c             write(*,*),'FluxAR GAS EXCH = ', FluxAR(i,j)
c             write(*,*),'*****************************'


            ELSE
              FlxAR(i,j) = 0. _d 0
            ENDIF


          END DO
        END DO

C update surface tendencies      
        DO j=jmin,jmax
          DO i=imin,imax
           SGAR(i,j)= FlxAR(i,j)
     &         *recip_drF(K) * recip_hFacC(i,j,K,bi,bj)

C              write(*,*),'SGAR = ', SGAR(i,j)
          ENDDO
         ENDDO


        RETURN
        END
#endif  /*ALLOW_CARBON*/

#endif  /*DARWIN*/
#endif  /*ALLOW_PTRACERS*/
c ==================================================================


