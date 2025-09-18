C TODO
C - have to call fechem after model?

#include "DARWIN_OPTIONS.h"
#ifdef ALLOW_SEAICE
#include "SEAICE_OPTIONS.h"
#endif

CBOP
C !ROUTINE: DARWIN_FORCING
C !INTERFACE: ==========================================================
      SUBROUTINE DARWIN_FORCING( Ptrdummy,
     I                        bi, bj, iMin, iMax, jMin, jMax,
     I                        myIter, myTime, myThid )

C !DESCRIPTION:

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "GRID.h"
#include "EEPARAMS.h"
#ifdef ALLOW_EXCH2
#include "W2_EXCH2_SIZE.h"
#include "W2_EXCH2_TOPOLOGY.h"
#include "W2_EXCH2_PARAMS.h"
#endif
#include "PARAMS.h"
#include "DYNVARS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "PTRACERS_FIELDS.h"
#ifdef ALLOW_SEAICE
#include "SEAICE_SIZE.h"
#include "SEAICE.h"
#endif
#ifdef ALLOW_RADTRANS
#include "RADTRANS_SIZE.h"
#endif
#ifdef ALLOW_DARWIN
#include "GCHEM.h"
#include "DARWIN_SIZE.h"
#include "DARWIN_INDICES.h"
#include "DARWIN_DIAGS.h"
#include "DARWIN_PARAMS.h"
#include "DARWIN_TRAITS.h"
#include "DARWIN_FIELDS.h"
#include "DARWIN_EXF_FIELDS.h"
#endif
C#ifdef ALLOW_WAVE_FORCING 
#include "EXF_OPTIONS.h"
#include "EXF_FIELDS.h"
C#endif
#ifdef DARWIN_NUTRIENT_RUNOFF
#include "EXF_INTERP_SIZE.h"
#include "DARWIN_EXF_PARAMS.h"
#endif

C !INPUT PARAMETERS: ===================================================
C  Ptrdummy             :: dummy for darwin2 compatibility
C  myThid               :: thread number
      _RL Ptrdummy(*)
      _RL myTime
      INTEGER iMin, iMax, jMin, jMax, bi, bj, myIter, myThid
CEOP

#ifdef ALLOW_DARWIN

C!LOCAL VARIABLES: ====================================================
C  i,j                  :: loop indices
C  k                    :: vertical level
      LOGICAL  DIAGNOSTICS_IS_ON
      EXTERNAL DIAGNOSTICS_IS_ON
      INTEGER i,j,k,kdn,iTr,l,isub
      CHARACTER*8 diagname
      _RL dTsub(Nr)
      _RL midTime, subTime
      _RL gPtr(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, Nr, nDarwin)
      _RL gPtrTemp(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, Nr, nDarwin)
      _RL PAR(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, Nr, nlam)
      _RL diags(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr, darwin_nDiag)
      _RL gDIC(1-OLx:sNx+OLx, 1-OLy:sNy+OLy)
      _RL gALK(1-OLx:sNx+OLx, 1-OLy:sNy+OLy)
      _RL gDICsurfForc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL gALKsurfForc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL gNO3surfForc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL gNO2surfForc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL gNH4surfForc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)  
      _RL gPO4surfForc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL gFeTsurfForc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL gSiO2surfForc(1-OLx:sNx+OLx, 1-OLy:sNy+OLy)
      _RL gO2surfForc(1-OLx:sNx+OLx, 1-OLy:sNy+OLy)
      _RL gO2(1-OLx:sNx+OLx, 1-OLy:sNy+OLy)
      _RL surfFe(1-OLx:sNx+OLx, 1-OLy:sNy+OLy)
      _RL freeFe(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, Nr)
      _RL diffFe1(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, Nr)
      _RL diffFe2(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, Nr)
      _RL diffFe(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, Nr)
      _RL scv,scav_poc
      _RL flx, POCl
      _RL ptr(nDarwin), gtr(nDarwin), PARl(nlam)
      _RL chlout(nPhoto)
      _RL diagsl(darwin_nDiag)
      _RL photoTempFunc(nplank)
      _RL hetTempFunc(nplank)
      _RL grazTempFunc(nplank)
      _RL reminTempFunc
      _RL mortTempFunc
      _RL mort2TempFunc
      _RL uptakeTempFunc
      _RL tmp
      _RL sedFe(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL cNO3(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL cNO2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL cNH4(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)      
      _RL cPO4(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL cSi(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL cFe(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL sNO3(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL sNO2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL sNH4(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL sPO4(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL sSi(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL sFe(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL consDIC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL consDIC_PIC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL respirDIC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL reminDIC_DOC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL reminDIC_POC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL disscDIC_PIC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL cALK(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL sALK(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL cO2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL sO2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)

C macroalgae model variables
      _RL hor_vel
      _RL day_yr  
      _RL grate_madum,mort_madum
      _RL gTdum, gHdum, gQdum, gEdum
      _RL UptakeNout, dEout, mortout, growthout
      _RL excess_growth
      _RL excess_mort
      _RL prev_qmag
      _RL prev_length
      _RL length_below
      _RL biomass_up
      _RL biomass_down
      INTEGER k_bottom      
      _RL grate_ma(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)  ! kelp growth rate
      _RL mortrate_ma(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)  ! mortality rate      
C terms controlling growth of Kelp      
      _RL gT_ma(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)   ! gT, temperature
      _RL gQ_ma(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)   ! gQ, nitrogen quota
      _RL gH_ma(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)   ! gH, wave motion 
      _RL gE_ma(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)   ! gE, irradiance
      _RL UptakeN(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL exudation_mag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL mortality_mag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL growth_mag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL maglength(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)

      INTEGER curr_date(4), first_day(4), numdays(4)
	  
      DO k=1,Nr
       dTsub(k) = PTRACERS_dTLev(k)/nsubtime
      ENDDO

C     time at middle of sub-timestep
      midTime = myTime - deltaTclock + .5*deltaTclock/nsubtime
C     time at end of sub-timestep
      subTime = myTime - deltaTclock + deltaTclock/nsubtime
      DO isub=1,nsubtime

C === reset tendencies =================================================
      DO itr=1,nDarwin
      DO k=1,Nr
      DO j=jMin,jMax
      DO i=iMin,iMax
        gPtr(i,j,k,iTr) = 0.0 _d 0
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO k=1,Nr
      DO j=jMin,jMax
      DO i=iMin,iMax
        diffFe1(i,j,k) = 0.0 _d 0
        diffFe2(i,j,k) = 0.0 _d 0
        diffFe(i,j,k) = 0.0 _d 0
      ENDDO
      ENDDO
      ENDDO

#ifdef DARWIN_ALLOW_MACROALGAE   
C initialize total macroalgae length  
C mag_tot_length(i,j,bi,bj) is iterative through k so do not include 
C total length until the end of k loop  
      DO j=jMin,jMax
      DO i=iMin,iMax
      mag_tot_length(i,j,bi,bj) = mag_cum_length(i,j,bi,bj)       
      mag_cum_length(i,j,bi,bj) = 0.0 _d 0 
      ENDDO
      ENDDO
#endif
      
C === light ============================================================
C     Initialize Chl from balanced-growth Chl:C if requested
C     and check Chl:C bounds.
C     Note: myIter has already been incremented
c      IF (myIter-1 .EQ. darwin_chlIter0) THEN
c        CALL DARWIN_INIT_CHL(bi, bj, subTime, myIter, myThid)
c      ENDIF
c
#ifdef ALLOW_SEAICE
      IF (DARWIN_useSEAICE) THEN
       DO j=jMin,jMax
        DO i=iMin,iMax
         iceFrac(i,j,bi,bj) = AREA(i,j,bi,bj)
        ENDDO
       ENDDO
      ENDIF
#endif

      CALL TIMER_START('DARWIN_LIGHT [DARWIN_FORCING]',myThid)
#ifdef ALLOW_RADTRANS
      CALL DARWIN_LIGHT_RADTRANS(PAR,subTime,bi,bj,iMin,iMax,jMin,jMax,
     &                        subTime,myIter,myThid)
#else
      CALL DARWIN_LIGHT(PAR, midTime, bi, bj, iMin, iMax, jMin, jMax, 
     & subTime, myIter, myThid)
#endif
      CALL TIMER_STOP ('DARWIN_LIGHT [DARWIN_FORCING]',myThid)

C === dic ==============================================================
#ifdef DARWIN_ALLOW_CARBON
C carbon air-sea interaction
      CALL TIMER_START('DARWIN_SURFFORCING [DARWIN_FORCING]',myThid)
      CALL DARWIN_SURFFORCING(
     O                    gDIC, gALK, gO2,
     I                    bi,bj,iMin,iMax,jMin,jMax,
     I                    subTime,myIter,myThid)
      CALL TIMER_STOP ('DARWIN_SURFFORCING [DARWIN_FORCING]',myThid)
      DO j=jMin,jMax
      DO i=iMin,iMax
        gPtr(i,j,1,iDIC) = gPtr(i,j,1,iDIC) + gDIC(i,j)
        gPtr(i,j,1,iALK) = gPtr(i,j,1,iALK) + gALK(i,j)
        gPtr(i,j,1,iO2)  = gPtr(i,j,1,iO2)  + gO2(i,j)
C        gDICsurfForc(i,j) = surfaceForcingPTr(i,j,bi,bj,iDIC) 
C     &   * recip_drF(1)*recip_hFacC(i,j,1,bi,bj)	
      ENDDO
      ENDDO
#endif

C surface fluxes from E/P/runoff
      DO j=jMin,jMax
      DO i=iMin,iMax
        gDICsurfForc(i,j) = surfaceForcingPTr(i,j,bi,bj,iDIC)
     &   * recip_drF(1)*recip_hFacC(i,j,1,bi,bj)
        gALKsurfForc(i,j) = surfaceForcingPTr(i,j,bi,bj,iALK)
     &   * recip_drF(1)*recip_hFacC(i,j,1,bi,bj)
        gNO3surfForc(i,j) = surfaceForcingPTr(i,j,bi,bj,iNO3)
     &   * recip_drF(1)*recip_hFacC(i,j,1,bi,bj)
        gNO2surfForc(i,j) = surfaceForcingPTr(i,j,bi,bj,iNO2)
     &   * recip_drF(1)*recip_hFacC(i,j,1,bi,bj)
        gNH4surfForc(i,j) = surfaceForcingPTr(i,j,bi,bj,iNH4)
     &   * recip_drF(1)*recip_hFacC(i,j,1,bi,bj)
        gPO4surfForc(i,j) = surfaceForcingPTr(i,j,bi,bj,iPO4)
     &   * recip_drF(1)*recip_hFacC(i,j,1,bi,bj)
        gFeTsurfForc(i,j) = surfaceForcingPTr(i,j,bi,bj,iFeT)
     &   * recip_drF(1)*recip_hFacC(i,j,1,bi,bj)
        gSiO2surfForc(i,j) = surfaceForcingPTr(i,j,bi,bj,iSiO2)
     &   * recip_drF(1)*recip_hFacC(i,j,1,bi,bj)
        gO2surfForc(i,j) = surfaceForcingPTr(i,j,bi,bj,iO2)
     &   * recip_drF(1)*recip_hFacC(i,j,1,bi,bj)
      ENDDO
      ENDDO
      
C === nutrient runoff ==================================================
#ifdef DARWIN_NUTRIENT_RUNOFF
      CALL DARWIN_ADD_SURFFORC(
     U                          gPtr(1-OLx,1-OLy,1,iDOC),
     I                          DOCrunoff(1-OLx,1-OLy,bi,bj), 1 _d 0,
     I                          bi, bj, iMin, iMax, jMin, jMax,
     I                          myIter, myTime, myThid )

      CALL DARWIN_ADD_SURFFORC(
     U                          gPtr(1-OLx,1-OLy,1,iDON),
     I                          DONrunoff(1-OLx,1-OLy,bi,bj), 1 _d 0,
     I                          bi, bj, iMin, iMax, jMin, jMax,
     I                          myIter, myTime, myThid )

      CALL DARWIN_ADD_SURFFORC(
     U                          gPtr(1-OLx,1-OLy,1,iDOP),
     I                          DOPrunoff(1-OLx,1-OLy,bi,bj), 1 _d 0,
     I                          bi, bj, iMin, iMax, jMin, jMax,
     I                          myIter, myTime, myThid )

      CALL DARWIN_ADD_SURFFORC(
     U                          gPtr(1-OLx,1-OLy,1,iNO3),
     I                          DINrunoff(1-OLx,1-OLy,bi,bj),
     I                          R_NO3_DIN_runoff,
     I                          bi, bj, iMin, iMax, jMin, jMax,
     I                          myIter, myTime, myThid )

      CALL DARWIN_ADD_SURFFORC(
     U                          gPtr(1-OLx,1-OLy,1,iNO2),
     I                          DINrunoff(1-OLx,1-OLy,bi,bj),
     I                          R_NO2_DIN_runoff,
     I                          bi, bj, iMin, iMax, jMin, jMax,
     I                          myIter, myTime, myThid )

      CALL DARWIN_ADD_SURFFORC(
     U                          gPtr(1-OLx,1-OLy,1,iNH4),
     I                          DINrunoff(1-OLx,1-OLy,bi,bj),
     I                          R_NH4_DIN_runoff,
     I                          bi, bj, iMin, iMax, jMin, jMax,
     I                          myIter, myTime, myThid )

      CALL DARWIN_ADD_SURFFORC(
     U                          gPtr(1-OLx,1-OLy,1,iNO3),
     I                          NO3runoff(1-OLx,1-OLy,bi,bj), 1 _d 0,
     I                          bi, bj, iMin, iMax, jMin, jMax,
     I                          myIter, myTime, myThid )

      CALL DARWIN_ADD_SURFFORC(
     U                          gPtr(1-OLx,1-OLy,1,iNO2),
     I                          NO2runoff(1-OLx,1-OLy,bi,bj), 1 _d 0,
     I                          bi, bj, iMin, iMax, jMin, jMax,
     I                          myIter, myTime, myThid )

      CALL DARWIN_ADD_SURFFORC(
     U                          gPtr(1-OLx,1-OLy,1,iNH4),
     I                          NH4runoff(1-OLx,1-OLy,bi,bj), 1 _d 0,
     I                          bi, bj, iMin, iMax, jMin, jMax,
     I                          myIter, myTime, myThid )

      CALL DARWIN_ADD_SURFFORC(
     U                          gPtr(1-OLx,1-OLy,1,iPO4),
     I                          IPrunoff(1-OLx,1-OLy,bi,bj),
     I                          R_DIP_IP_runoff,
     I                          bi, bj, iMin, iMax, jMin, jMax,
     I                          myIter, myTime, myThid )

      CALL DARWIN_ADD_SURFFORC(
     U                          gPtr(1-OLx,1-OLy,1,iSiO2),
     I                          DSirunoff(1-OLx,1-OLy,bi,bj), 1 _d 0,
     I                          bi, bj, iMin, iMax, jMin, jMax,
     I                          myIter, myTime, myThid )

      CALL DARWIN_ADD_SURFFORC(
     U                          gPtr(1-OLx,1-OLy,1,iPOC),
     I                          POCrunoff(1-OLx,1-OLy,bi,bj), 1 _d 0,
     I                          bi, bj, iMin, iMax, jMin, jMax,
     I                          myIter, myTime, myThid )

      CALL DARWIN_ADD_SURFFORC(
     U                          gPtr(1-OLx,1-OLy,1,iPON),
     I                          PONrunoff(1-OLx,1-OLy,bi,bj), 1 _d 0,
     I                          bi, bj, iMin, iMax, jMin, jMax,
     I                          myIter, myTime, myThid )

      CALL DARWIN_ADD_SURFFORC(
     U                          gPtr(1-OLx,1-OLy,1,iPOP),
     I                          POPrunoff(1-OLx,1-OLy,bi,bj), 1 _d 0,
     I                          bi, bj, iMin, iMax, jMin, jMax,
     I                          myIter, myTime, myThid )

      CALL DARWIN_ADD_SURFFORC(
     U                          gPtr(1-OLx,1-OLy,1,iDIC),
     I                          DICrunoff(1-OLx,1-OLy,bi,bj), 1 _d 0,
     I                          bi, bj, iMin, iMax, jMin, jMax,
     I                          myIter, myTime, myThid )

#ifdef DARWIN_ALLOW_CARBON
      CALL DARWIN_ADD_SURFFORC(
     U                          gPtr(1-OLx,1-OLy,1,iALK),
     I                          DICrunoff(1-OLx,1-OLy,bi,bj),
     I                          R_ALK_DIC_runoff,
     I                          bi, bj, iMin, iMax, jMin, jMax,
     I                          myIter, myTime, myThid )
#endif

      CALL DARWIN_ADD_SURFFORC(
     U                          gPtr(1-OLx,1-OLy,1,iFeT),
     I                          IPrunoff(1-OLx,1-OLy,bi,bj),
     I                          R_DIP_IP_runoff*R_DFe_DIP_runoff,
     I                          bi, bj, iMin, iMax, jMin, jMax,
     I                          myIter, myTime, myThid )

      CALL DARWIN_ADD_SURFFORC(
     U                          gPtr(1-OLx,1-OLy,1,iDOFe),
     I                          DOPrunoff(1-OLx,1-OLy,bi,bj),
     I                          R_DOFe_DOP_runoff,
     I                          bi, bj, iMin, iMax, jMin, jMax,
     I                          myIter, myTime, myThid )

      CALL DARWIN_ADD_SURFFORC(
     U                          gPtr(1-OLx,1-OLy,1,iPOFe),
     I                          POPrunoff(1-OLx,1-OLy,bi,bj),
     I                          R_POFe_POP_runoff,
     I                          bi, bj, iMin, iMax, jMin, jMax,
     I                          myIter, myTime, myThid )
     
#ifdef DARWIN_ALLOW_TSS
      CALL DARWIN_ADD_SURFFORC(
     U                          gPtr(1-OLx,1-OLy,1,iTSS),
     I                          TSSrunoff(1-OLx,1-OLy,bi,bj),
     I                          bi, bj, iMin, iMax, jMin, jMax,
     I                          myIter, myTime, myThid )
#endif

#endif

C === iron =============================================================
      CALL TIMER_START('DARWIN_FE_CHEM [DARWIN_FORCING]',myThid)
      CALL DARWIN_FE_CHEM(
     U                 Ptracer(1-OLx,1-OLy,1,bi,bj,iFeT),
     O                 diffFe1(1-OLx,1-OLy,1),
     O                 freeFe(1-OLx, 1-OLy, 1),
     I                 bi, bj, iMin, iMax, jMin, jMax, myThid)
      CALL TIMER_STOP ('DARWIN_FE_CHEM [DARWIN_FORCING]',myThid)

C     iron dust input
      CALL TIMER_START('DARWIN_FE_DUST [DARWIN_FORCING]',myThid)
      DO j=jMin,jMax
      DO i=iMin,iMax
      IF (hFacC(i,j,1,bi,bj) .NE. 0.) THEN
        gPtr(i,j,1,iFeT) =
     &  gPtr(i,j,1,iFeT) + alpfe * inputFe(i,j,bi,bj) /
     &                           drF(1) / hFacC(i,j,1,bi,bj)
      ENDIF
      ENDDO
      ENDDO
      CALL TIMER_STOP ('DARWIN_FE_DUST [DARWIN_FORCING]',myThid)

C     scavenging
      CALL TIMER_START('DARWIN_FE_SCAV [DARWIN_FORCING]',myThid)
      DO k=1,Nr
      DO j=jMin,jMax
      DO i=iMin,iMax
#ifdef DARWIN_PART_SCAV_POP
        scav_poc = MAX(0 _d 0, Ptracer(i,j,k,bi,bj,iPOP))/scav_R_POPPOC
        scv = scav_rat*scav_inter*(scav_poc**scav_exp)
#elif defined(DARWIN_PART_SCAV)
        scav_poc = MAX(0 _d 0, Ptracer(i,j,k,bi,bj,iPOC))
        scv = scav_rat*scav_inter*(scav_poc**scav_exp)
#else
        scv = scav
#endif
        gPtr(i,j,k,iFeT) = gPtr(i,j,k,iFeT) -
     &                     scv*MAX(0 _d 0, freefe(i,j,k))
      ENDDO
      ENDDO
      ENDDO
      CALL TIMER_STOP ('DARWIN_FE_SCAV [DARWIN_FORCING]',myThid)

C     iron sediment source (in bottom grid cell above kMaxFeSed)
      CALL TIMER_START('DARWIN_FE_SED [DARWIN_FORCING]',myThid)
      DO k = kMinFeSed, kMaxFeSed
       kdn = MIN(Nr, k+1)
       DO j=jMin,jMax
        DO i=iMin,iMax
         IF (hFacC(i,j,k,bi,bj) .GT. 0. .AND.
     &       (hFacC(i,j,kdn,bi,bj) .EQ. 0. .OR. k .EQ. Nr)) THEN
#ifdef DARWIN_IRON_SED_SOURCE_VARIABLE
# ifdef DARWIN_IRON_SED_SOURCE_POP
          flx = fesedflux_pcm*wp_sink*R_CP_fesed*
     &            MAX(0 _d 0, Ptracer(i,j,k-1,bi,bj,iPOP))
#else
          POCl = MAX(0 _d 0, Ptracer(i,j,k,bi,bj,iPOC))
          flx = fesedflux_pcm*wc_sink*POCl - fesedflux_min
          flx = MAX(0. _d 0, flx)
# endif
#else
          flx = fesedflux
#endif
	  sedFe(i,j,k) = flx/(drF(k)*hFacC(i,j,k,bi,bj))
          gPtr(i,j,k,iFeT) = gPtr(i,j,k,iFeT) +
     &                             flx/(drF(k)*hFacC(i,j,k,bi,bj))
         ENDIF
        ENDDO
       ENDDO
      ENDDO
      CALL TIMER_STOP ('DARWIN_FE_SED [DARWIN_FORCING]',myThid)

C === plankton =========================================================
      CALL TIMER_START('DARWIN_PLANKTON [DARWIN_FORCING]',myThid)
      do k=1,Nr
      do j=jMin,jMax
      do i=iMin,iMax
      if (hFacC(i,j,k,bi,bj) .GT. 0.) then
        CALL DARWIN_TEMPFUNC(Theta(i,j,k,bi,bj),
     &         photoTempFunc, hetTempFunc, grazTempFunc,
     &         reminTempFunc, mortTempFunc, mort2TempFunc,
     &         uptakeTempFunc, myThid)
        DO iTr=1,nDarwin
         ptr(iTr) = Ptracer(i, j, k, bi, bj, iTr)
         gtr(iTr) = gPtr(i, j, k, iTr)
        ENDDO
        DO l=1,nlam
         PARl(l) = MAX(0 _d 0, PAR(i, j, k, l))
        ENDDO
        CALL DARWIN_PLANKTON(
     I                 ptr,
     U                 gtr,
     O                 chlout, diagsl,
     I                 PARl, photoTempFunc,
     I                 hetTempFunc,
     I                 grazTempFunc,
     I                 reminTempFunc,
     I                 mortTempFunc, mort2TempFunc,
     I                 uptakeTempFunc,
#ifdef DARWIN_DEBUG
     I                 myxgloballo+(bi-1)*sNx+i-1,
     I                 myygloballo+(bj-1)*sNy+j-1,
     I                 k, dTsub(k),
#endif
     I                 subTime, myIter, myThid )


#ifdef DARWIN_ALLOW_MACROALGAE
C === macrolagae =========================================================
      CALL TIMER_START('DARWIN_MACROALGAE [DARWIN_FORCING]',myThid)
C macroalgae code                     
        k_bottom = kLowC(i,j,bi,bj)
C       #################Remove after tests ################################## 
        k_bottom = 3 
C        ######################################################################
C compute horizontal velocity  
        hor_vel=sqrt(
     &  (0.5*(uvel(i,j,k,bi,bj)+uvel(i+1,j,k,bi,bj)))**2 +
     &  (0.5*(vvel(i,j,k,bi,bj)+vvel(i,j+1,k,bi,bj)))**2)  
C Vertical transfer of growth and mortality from previous timestep     
        excess_growth = growth_transfer(i,j,k,bi,bj)
        excess_mort = mort_transfer(i,j,k,bi,bj)
C Reset transfer vector        
        growth_transfer(i,j,k,bi,bj) = 0.0 _d 0
        mort_transfer(i,j,k,bi,bj) = 0.0 _d 0
C Set quota to value of quota from the bellow grid cell if transfer
        IF(Qmag_transfer(i,j,k,bi,bj).GT.0.0)THEN
         Ptracer(i, j, k, bi, bj, iMagQ) = Qmag_transfer(i,j,k,bi,bj)       
         ptr(iMagQ) = Qmag_transfer(i,j,k,bi,bj)
         Qmag_transfer(i,j,k,bi,bj) = 0.0 _d 0
        ENDIF
C Check if there is biomass in the grid cell above        
        biomass_up = 0.0
        IF(k.GT.1)THEN
          biomass_up = Ptracer(i, j, k-1, bi, bj, iMagB)
        ENDIF
C Save biomass from the lowest grid cell        
        biomass_down = Ptracer(i, j, k_bottom, bi, bj, iMagB)
C Compute frond length (if length above max grid cell thickness then keep it at grid cell thickness)
C Save length from previous time step        
        prev_length = mag_length(i,j,k,bi,bj)
        IF(k.EQ.k_bottom)THEN  
C ----- Subsurface ------
          IF(mag_length(i,j,k-1,bi,bj).GT.0.0)THEN
            mag_length(i,j,k,bi,bj) = drF(k)*hFacC(i,j,k,bi,bj)
          ELSE
            mag_length(i,j,k,bi,bj) = min(
     &       (ptr(iMagB)*
     &       ptr(iMagQ))*
     &       (1.0/mp_spp_length_subsurface),
     &       drF(k)*hFacC(i,j,k,bi,bj))
          ENDIF       
        ELSEIF(k.GT.1.AND.k.LT.k_bottom)THEN
C ----- Water column ------
          IF(mag_length(i,j,k-1,bi,bj).GT.0.0)THEN
            mag_length(i,j,k,bi,bj) = drF(k)*hFacC(i,j,k,bi,bj)
          ELSE
            mag_length(i,j,k,bi,bj) = min(
     &       (ptr(iMagB)*
     &       ptr(iMagQ))*
     &       (1.0/mp_spp_length_watercolumn),
     &       drF(k)*hFacC(i,j,k,bi,bj))
          ENDIF           
        ELSEIF(k.EQ.1)THEN
C ----- Canopy ------
C Allows for the frond to grow longer than the thickness of the grid cell 
C for mimicking horizontal growth of canopy at the surface
C but cannot exceed max length   
           length_below = mag_tot_length(i,j,bi,bj) - 
     &       mag_length(i,j,k,bi,bj)                 
           mag_length(i,j,k,bi,bj) = min(
     &       (ptr(iMagB)*
     &       ptr(iMagQ))*
     &       (1.0/mp_spp_length_canopy),
     &       max(0.0,(mp_spp_maxlength - length_below)))            
        ELSE
         mag_length(i,j,k,bi,bj) = 0.0 _d 0        
        ENDIF
C Compute total length (here is cumulative so it reaches total frond length only at kmax/kbottom)           
        mag_cum_length(i,j,bi,bj) = mag_cum_length(i,j,bi,bj) +
     &   mag_length(i,j,k,bi,bj)  
C Reset quota to zero if biomass is exhausted 
        IF(prev_length.GT.0.0.AND.mag_length(i,j,k,bi,bj).EQ.0.0)THEN
          Ptracer(i, j, k, bi, bj, iMagQ) = 0.0       
          ptr(iMagQ) = 0.0
        ENDIF
                                      
        CALL DARWIN_MACROALGAE(
     I   k,k_bottom,deltaTclock,
     I   drF(k),hFacC(i,j,k,bi,bj),  
     I   theta(i,j,k,bi,bj),hor_vel,parl,
     I   WvHeight(i,j,bi,bj),WvPeriod(i,j,bi,bj),
     I   ptr,
     I   mag_cum_length(i,j,bi,bj),
     I   mag_length(i,j,k,bi,bj),
     I   mag_tot_length(i,j,bi,bj),
     I   biomass_up,biomass_down,
     U   gtr,
     U   excess_growth,excess_mort,
     O   prev_qmag,
     O   gTdum, gQdum, gEdum, grate_madum,
     O   mort_madum,UptakeNout,
     O   mortout, dEout, growthout)
                  
C Transfer vertical and horizontal growth of macroalgae if carrying capacity reached [gDW/m3/s]       
        IF (mp_spp_vertical) THEN                                    
C Vertical transfer for next timestep
          IF (k.GT.1) THEN   
            growth_transfer(i,j,k-1,bi,bj) =  
     &       growth_transfer(i,j,k-1,bi,bj) + excess_growth
            mort_transfer(i,j,k-1,bi,bj) =  
     &       mort_transfer(i,j,k-1,bi,bj) + excess_mort           
           IF(prev_qmag.GT.0.0)THEN
             Qmag_transfer(i,j,k-1,bi,bj) = prev_qmag
           ELSE
             Qmag_transfer(i,j,k-1,bi,bj) = 0.0
           ENDIF                 
          ENDIF          
        ENDIF   
          
C       growth and mortality rates [s-1] and limitation terms [0-1]
        grate_ma(i,j,k) = grate_madum
        gT_ma(i,j,k) = gTdum
        gQ_ma(i,j,k) = gQdum 
        gE_ma(i,j,k) = gEdum
        mortrate_ma(i,j,k) = mort_madum
C       loss/gain terms from macroalgae
C       uptakeN and exudation in [mmolN/m3/s]
C       mortality in [mmolC/m3/s]
        UptakeN(i,j,k) = UptakeNout
        exudation_mag(i,j,k) = dEout 
        mortality_mag(i,j,k) = mortout
        growth_mag(i,j,k) = growthout
        
C end macroalgae code
      CALL TIMER_STOP('DARWIN_MACROALGAE [DARWIN_FORCING]',myThid)
#endif

        DO l=1,darwin_nDiag
         diags(i, j, k, l) = diagsl(l)
        ENDDO
#ifdef DARWIN_ALLOW_CONS
        DARWIN_Nfix(i,j,k,bi,bj) = diagsl(iNfix)
        DARWIN_Ndenit(i,j,k,bi,bj) = diagsl(iDenitN)
#endif
#ifndef DARWIN_ALLOW_CHLQUOTA
#ifdef ALLOW_RADTRANS
        DO l=1,nPhoto
         chlPrev(i, j, k, bi, bj, l) = chlout(l)
        ENDDO
#else
        tmp = 0 _d 0
        DO l=1,nPhoto
         tmp = tmp + chlout(l)
        ENDDO
        chlPrev(i, j, k, bi, bj) = tmp
#endif
#endif
        DO iTr=1,nDarwin
         gPtr(i, j, k, iTr) = gtr(iTr)
        ENDDO
        
      endif
      enddo
      enddo
      enddo
      CALL TIMER_STOP('DARWIN_PLANKTON [DARWIN_FORCING]',myThid)

C === sinking ==========================================================
      CALL TIMER_START('DARWIN_SINKING [DARWIN_FORCING]',myThid)
      CALL DARWIN_SINKING( Ptracer,gPtr,bi,bj,iMin,iMax,jMin,jMax,
     &                  subTime,myIter,myThid )
      CALL TIMER_STOP ('DARWIN_SINKING [DARWIN_FORCING]',myThid)

C === apply tendencies to tracers ======================================
      CALL TIMER_START('DARWIN_STEP [DARWIN_FORCING]',myThid)
      DO iTr=1,nDarwin
      DO k=1,Nr
      DO j=jMin,jMax
      DO i=iMin,iMax
        pTracer(i,j,k,bi,bj,iTr)=pTracer(i,j,k,bi,bj,iTr)
     &                      +dTsub(k)*gPtr(i,j,k,iTr)*maskInC(i,j,bi,bj)  
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      CALL TIMER_STOP ('DARWIN_STEP [DARWIN_FORCING]',myThid)

C === iron =============================================================
C     re-apply free iron limit to FeT
      CALL TIMER_START('DARWIN_FE_CHEM [DARWIN_FORCING]',myThid)
      CALL DARWIN_FE_CHEM(
     U                 pTracer(1-OLx,1-OLy,1,bi,bj,iFeT),
     O                 diffFe2(1-OLx,1-OLy,1),
     O                 freeFe(1-OLx,1-OLy,1),
     I                 bi, bj, iMin, iMax, jMin, jMax, myThid)
      CALL TIMER_STOP ('DARWIN_FE_CHEM [DARWIN_FORCING]',myThid)

      DO k=1,Nr
      DO j=jMin,jMax
      DO i=iMin,iMax
        diffFe(i,j,k) = hFacC(i,j,k,bi,bj) 
     &   * ((diffFe1(i,j,k) + diffFe2(i,j,k)) / deltaTClock)
      ENDDO
      ENDDO
      ENDDO
C === diagnostics ======================================================
#ifdef ALLOW_DIAGNOSTICS
      IF (useDIAGNOSTICS) THEN
       CALL TIMER_START('DIAGS_FILL [DARWIN_FORCING]',myThid)
       DO l = 1, nlam
        WRITE(diagname, '(A,I3.3)') 'PAR', l
        CALL DIAGNOSTICS_FILL(PAR(1-OLx,1-OLy,1,l),diagname,0,Nr,2,
     &        bi,bj,myThid)
       ENDDO
       IF (DIAGNOSTICS_IS_ON('PAR     ', myThid)) THEN
        DO l=2,nlam
         DO k=1,Nr
          DO j=1,sNy
           DO i=1,sNx
            PAR(i,j,k,1) = PAR(i,j,k,1) + PAR(i,j,k,l)
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        WRITE(diagname, '(A)') 'PAR'
        CALL DIAGNOSTICS_FILL(PAR,diagname,0,Nr,2,bi,bj,myThid)
       ENDIF
       
       CALL DIAGNOSTICS_FILL(diags(1-OLx,1-OLy,1,iPP),   'PP      ',
     &       0,Nr,2,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(diags(1-OLx,1-OLy,1,iNfix), 'Nfix    ',
     &       0,Nr,2,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(diags(1-OLx,1-OLy,1,iDenit),'Denit   ',
     &       0,Nr,2,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(diags(1-OLx,1-OLy,1,iDenitN),'DenitN  ',
     &       0,Nr,2,bi,bj,myThid)

       DO k=1,Nr
        DO j=1,sNy
         DO i=1,sNx
          cNO3(i,j,k) = hFacC(i,j,k,bi,bj)
     &     * diags(i,j,k,iConsNO3)
          cNO2(i,j,k) = hFacC(i,j,k,bi,bj)
     &     * diags(i,j,k,iConsNO2)   
          cNH4(i,j,k) = hFacC(i,j,k,bi,bj)
     &     * diags(i,j,k,iConsNH4)
          cPO4(i,j,k) = hFacC(i,j,k,bi,bj)
     &     * diags(i,j,k,iConsPO4)
          cSi(i,j,k) = hFacC(i,j,k,bi,bj)
     &     * diags(i,j,k,iConsSi)
          cFe(i,j,k) = hFacC(i,j,k,bi,bj)
     &     * diags(i,j,k,iConsFe)
          consDIC(i,j,k) = hFacC(i,j,k,bi,bj)
     &     * diags(i,j,k,iConsDIC)
          consDIC_PIC(i,j,k) = hFacC(i,j,k,bi,bj)
     &     * diags(i,j,k,iConsDIC_PIC)
          respirDIC(i,j,k) = hFacC(i,j,k,bi,bj)
     &     * diags(i,j,k,iRespirDIC)
          reminDIC_DOC(i,j,k) = hFacC(i,j,k,bi,bj)
     &     * diags(i,j,k,iReminDIC_DOC)
          reminDIC_POC(i,j,k) = hFacC(i,j,k,bi,bj)
     &     * diags(i,j,k,iReminDIC_POC)
          disscDIC_PIC(i,j,k) = hFacC(i,j,k,bi,bj)
     &     * diags(i,j,k,iDisscDIC_PIC)    
          cALK(i,j,k) = hFacC(i,j,k,bi,bj)
     &     * diags(i,j,k,iConsALK)
          cO2(i,j,k) = hFacC(i,j,k,bi,bj)
     &     * diags(i,j,k,iConsO2)     

         ENDDO
        ENDDO
       ENDDO

C      nutrient consumption diagnostics
       CALL DIAGNOSTICS_FILL(cNO3(1-OLx,1-OLy,1),
     &  'C_NO3   ',0,Nr,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(cNO2(1-OLx,1-OLy,1),
     &  'C_NO2   ',0,Nr,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(cNH4(1-OLx,1-OLy,1),
     &  'C_NH4   ',0,Nr,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(cPO4(1-OLx,1-OLy,1),
     &  'C_PO4   ',0,Nr,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(cSi(1-OLx,1-OLy,1),
     &  'C_Si    ',0,Nr,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(cFe(1-OLx,1-OLy,1),
     & 'C_Fe    ',0,Nr,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(consDIC(1-OLx,1-OLy,1),
     & 'cDIC    ',0,Nr,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(consDIC_PIC(1-OLx,1-OLy,1),
     & 'cDIC_PIC',0,Nr,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(respirDIC(1-OLx,1-OLy,1),
     & 'respDIC ',0,Nr,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(reminDIC_DOC(1-OLx,1-OLy,1),
     & 'rDIC_DOC',0,Nr,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(reminDIC_POC(1-OLx,1-OLy,1),
     & 'rDIC_POC',0,Nr,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(disscDIC_PIC(1-OLx,1-OLy,1),
     & 'dDIC_PIC',0,Nr,1,bi,bj,myThid)     
       CALL DIAGNOSTICS_FILL(cALK(1-OLx,1-OLy,1),
     & 'C_ALK   ',0,Nr,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(cO2(1-OLx,1-OLy,1),
     & 'C_O2    ',0,Nr,1,bi,bj,myThid)

C      compute 'source' diagnostics: tendencies without consumption;
C      add full tendencies to consumption and store result back into
C      diags(iCons*) for convenience
       DO k=1,Nr
        DO j=1,sNy
         DO i=1,sNx
          sALK(i,j,k) = hFacC(i,j,k,bi,bj) * (diags(i,j,k,iConsALK)
     &     + gPtr(i,j,k,iALK))
          sNO3(i,j,k) = hFacC(i,j,k,bi,bj) * (diags(i,j,k,iConsNO3)
     &     + gPtr(i,j,k,iNO3))
          sNO2(i,j,k) = hFacC(i,j,k,bi,bj) * (diags(i,j,k,iConsNO2)
     &     + gPtr(i,j,k,iNO2))
          sNH4(i,j,k) = hFacC(i,j,k,bi,bj) * (diags(i,j,k,iConsNH4)
     &     + gPtr(i,j,k,iNH4))
          sPO4(i,j,k) = hFacC(i,j,k,bi,bj) * ((diags(i,j,k,iConsPO4))
     &     + gPtr(i,j,k,iPO4))
          sSi(i,j,k) = hFacC(i,j,k,bi,bj) * (diags(i,j,k,iConsSi)
     &     + gPtr(i,j,k,iSiO2))
          sFe(i,j,k) = hFacC(i,j,k,bi,bj) * (diags(i,j,k,iConsFe)
     &     + gPtr(i,j,k,iFeT))
          sO2(i,j,k) = hFacC(i,j,k,bi,bj) * (diags(i,j,k,iConsO2)
     &     + gPtr(i,j,k,iO2))
 
         ENDDO
        ENDDO
       ENDDO

       CALL DIAGNOSTICS_FILL(sALK(1-OLx,1-OLy,1),
     &  'S_ALK   ',0,Nr,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(sNO3(1-OLx,1-OLy,1),
     &  'S_NO3   ',0,Nr,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(sNO2(1-OLx,1-OLy,1),
     &  'S_NO2   ',0,Nr,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(sNH4(1-OLx,1-OLy,1),
     &  'S_NH4   ',0,Nr,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(sPO4(1-OLx,1-OLy,1),
     &  'S_PO4   ',0,Nr,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(sSi(1-OLx,1-OLy,1),
     &  'S_Si    ',0,Nr,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(sFe(1-OLx,1-OLy,1),
     &  'S_Fe    ',0,Nr,1,bi,bj,myThid)
      CALL DIAGNOSTICS_FILL(sO2(1-OLx,1-OLy,1),
     &  'S_O2    ',0,Nr,1,bi,bj,myThid)

#ifdef DARWIN_ALLOW_MACROALGAE     
C macroalgae diagnostics
      CALL DIAGNOSTICS_FILL(gT_ma(1-OLx,1-OLy,1),
     &  'gT_mag  ',0,Nr,1,bi,bj,myThid)

      CALL DIAGNOSTICS_FILL(gQ_ma(1-OLx,1-OLy,1),
     &  'gQ_mag  ',0,Nr,1,bi,bj,myThid)

      CALL DIAGNOSTICS_FILL(gE_ma(1-OLx,1-OLy,1),
     &  'gE_mag  ',0,Nr,1,bi,bj,myThid)

      CALL DIAGNOSTICS_FILL(gH_ma(1-OLx,1-OLy,1),
     &  'gH_mag  ',0,Nr,1,bi,bj,myThid)

      CALL DIAGNOSTICS_FILL(grate_ma(1-OLx,1-OLy,1),
     &  'mu_mag  ',0,Nr,1,bi,bj,myThid)

      CALL DIAGNOSTICS_FILL(mortrate_ma(1-OLx,1-OLy,1),
     &  'mrtrte_m',0,Nr,1,bi,bj,myThid)
     
      CALL DIAGNOSTICS_FILL(UptakeN(1-OLx,1-OLy,1),
     &  'UptkN   ',0,Nr,1,bi,bj,myThid)   
     
      CALL DIAGNOSTICS_FILL(exudation_mag(1-OLx,1-OLy,1),
     &  'exud_mag',0,Nr,1,bi,bj,myThid)  
     
      CALL DIAGNOSTICS_FILL(mortality_mag(1-OLx,1-OLy,1),
     &  'mort_mag',0,Nr,1,bi,bj,myThid)  
     
      CALL DIAGNOSTICS_FILL(growth_mag(1-OLx,1-OLy,1),
     &  'grth_mag',0,Nr,1,bi,bj,myThid)  
     
      CALL DIAGNOSTICS_FILL(mag_length,
     &  'lgth_mag',0,Nr,1,bi,bj,myThid) 
     
      CALL DIAGNOSTICS_FILL(mag_tot_length,
     &  'tot_lgth',0,1,2,bi,bj,myThid)
     
#endif

      DO iTr=1,nDarwin
       DO k=1,Nr
        DO j=1,sNy
         DO i=1,sNx
          gPtrTemp(i,j,k,iTr) = hFacC(i,j,k,bi,bj)*
     &     gPtr(i,j,k,iTr)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

       DO iTr=1,nDarwin
        diagname = 'gDAR'//PTRACERS_ioLabel(iTr)
        CALL DIAGNOSTICS_FILL(gPtrTemp(1-OLx,1-OLy,1,iTr), diagname,
     &   0,Nr,2,bi,bj,myThid)
     &  
       ENDDO

#ifdef DARWIN_DIAG_PERTYPE
       DO iTr=1,nplank
        WRITE(diagname, '(A,I4.4)') 'PP', iTr
        CALL DIAGNOSTICS_FILL(diags(1-OLx,1-OLy,1,iPPplank+iTr-1),
     &        diagname,0,Nr,2,bi,bj,myThid)
        WRITE(diagname, '(A,I4.4)') 'GR', iTr
        CALL DIAGNOSTICS_FILL(diags(1-OLx,1-OLy,1,iGRplank+iTr-1),
     &        diagname,0,Nr,2,bi,bj,myThid)
        WRITE(diagname, '(A,I4.4)') 'GrGn', iTr
        CALL DIAGNOSTICS_FILL(diags(1-OLx,1-OLy,1,iGrGn+iTr-1),
     &        diagname,0,Nr,2,bi,bj,myThid)
       ENDDO
#endif
       IF (DIAGNOSTICS_IS_ON('sfcSolFe', myThid)) THEN
        DO j=jMin,jMax
         DO i=iMin,iMax
          IF (hFacC(i,j,1,bi,bj) .EQ. 0.) THEN
           surfFe(i,j) = 0 _d 0
          ELSE
           surfFe(i,j) = alpfe * inputFe(i,j,bi,bj)
          ENDIF
         ENDDO
        ENDDO
        CALL DIAGNOSTICS_FILL(surfFe,'sfcSolFe',0,1,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(sedFe,'sedFe   ',0,Nr,1,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(diffFe,'freeFe  ',0,Nr,1,bi,bj,myThid)
       ENDIF

#ifdef DARWIN_ALLOW_CARBON
       CALL DIAGNOSTICS_FILL(pH,     'pH      ',0,Nr,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(pCO2,   'pCO2    ',0,Nr,1,bi,bj,myThid)

       CALL DIAGNOSTICS_FILL(FluxCO2,'fluxCO2 ',0,1,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(gDIC,   'gDICsurf',0,1,2,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(gDICsurfForc,'gDICEpr ',0,1,2,
     &  bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(gALKsurfForc,'gALKEpr ',0,1,2,
     &  bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(gO2surfForc, 'gO2Epr  ',0,1,2,
     &  bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(gNO3surfForc,'gNO3Epr ',0,1,2,
     &  bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(gNO2surfForc,'gNO2Epr ',0,1,2,
     &  bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(gNH4surfForc,'gNH4Epr ',0,1,2,
     &  bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(gPO4surfForc,'gPO4Epr ',0,1,2,
     &  bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(gFeTsurfForc,'gFeTEpr ',0,1,2,
     &  bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(gSiO2surfForc,'gSiO2Epr',0,1,2,
     &  bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(gALK,   'gALKsurf',0,1,2,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(gO2,    'gO2surf ',0,1,2,bi,bj,myThid)
#endif
       CALL TIMER_STOP ('DIAGS_FILL [DARWIN_FORCING]',myThid)
      ENDIF
#endif

      midTime = midTime + deltaTclock/nsubtime
      subTime = subTime + deltaTclock/nsubtime
C     isub
      ENDDO

#endif /* ALLOW_DARWIN */

      RETURN
      END
