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
#ifdef DARWIN_ALLOW_CARBON
       _RL  t
       _RL  s
       _RL  Tlocal
       _RL  Slocal
       _RL  ta
       _RL  pt
       _RL  sit
       _RL  tk
       _RL  tk100
       _RL  tk1002
       _RL  dlogtk
       _RL  sqrtis
       _RL  sqrts
       _RL  s15
       _RL  scl
       _RL  x1
       _RL  x2
       _RL  s2
       _RL  xacc
       _RL  invtk
       _RL  is
       _RL  is2
       _RL  bdepth
       _RL  cdepth
       _RL  pressc
       _RL  Ksp_T_Calc
       _RL  xvalue
       _RL  zdum
       _RL  tmpa1
       _RL  tmpa2
       _RL  tmpa3
       _RL  logKspc
       _RL  dv
       _RL  dk
       _RL  pfactor
       _RL  bigR
       
       _RL dicl
       _RL docl
       _RL picl
       _RL alkl
       _RL o2l
       _RL po4l
       _RL sil    

       _RL pCO2SolverTemp
       _RL pCO2SolverSal
       _RL pCO2SolverDic
       _RL pCO2SolverPo4
       _RL pCO2SolverSi
       _RL pCO2SolverAlk
#endif

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
     &               subTime, myIter, myThid)
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

#ifdef DARWIN_ALLOW_CARBON
c compute Ksp as a function of temperature and pressure

             Tlocal = theta(i,j,k,bi,bj)
             Slocal = salt(i,j,k,bi,bj)

            bdepth = 0.0d0
            cdepth = 0.0d0
            pressc = 1.0d0

            do l = 1,k
             cdepth = bdepth + 0.5d0*drF(l)
             bdepth = bdepth + drF(l)
             pressc = 1.0d0 + 0.1d0*cdepth
            end do

            if (maskC(i,j,k,bi,bj).NE.0. _d 0) then 
             t = Tlocal
             s = max(4. _d 0, Slocal)

             tk = 273.15 + t
             tk100 = tk/100.0
             tk1002=tk100*tk100
             invtk=1.0/tk
             dlogtk=log(tk)
             is=19.924*s/(1000.-1.005*s)
             is2=is*is
             sqrtis=sqrt(is)
             s2=s*s
             sqrts=sqrt(s)
             s15=s**1.5
             scl=s/1.80655

c f = k0(1-pH2O)*correction term for non-ideality
c Weiss & Price (1980, Mar. Chem., 8, 347-359; Eq 13 with table 6 values)
             ff(i,j,bi,bj) = exp(-162.8301 + 218.2968/tk100 +
     &        90.9241*log(tk100)-1.47696*tk1002 +
     &        s*(0.025695-0.025225*tk100 +
     &        0.0049867*tk1002))

c K0 from Weiss 1974
             ak0(i,j,bi,bj) = exp(93.4517/tk100-60.2409 + 
     &        23.3585*log(tk100) +
     &        s*(0.023517-0.023656*tk100 + 
     &        0.0047036*tk1002))

c k1 = [H][HCO3]/[H2CO3]
c k2 = [H][CO3]/[HCO3]
c Millero p.664 (1995) using Mehrbach et al. data on seawater scale 
             ak1(i,j,bi,bj)=10**(-1*(3670.7*invtk - 
     &        62.008+9.7944*dlogtk -
     &        0.0118*s+0.000116*s2))
           
             ak2(i,j,bi,bj)=10**(-1*(1394.7*invtk+4.777 - 
     &        0.0184*s+0.000118*s2))

c NOW PRESSURE DEPENDENCE:
c Following Takahashi (1981) GEOSECS report - quoting Culberson and 
c Pytkowicz (1968)
c pressc = pressure in bars
             ak1(i,j,bi,bj) = ak1(i,j,bi,bj) *
     &        exp((24.2-0.085*t)*(pressc-1.0)/(83.143*tk))

c FIRST GO FOR K2: According to GEOSECS (1982) report
c            ak2(i,j,bi,bj) = ak2(i,j,bi,bj) *
c    &        exp((26.4-0.040*t)*(pressc-1.0)/(83.143*tk))

c SECOND GO FOR K2: corrected coeff according to CO2sys documentation
c E. Lewis and D. Wallace (1998) ORNL/CDIAC-105
             ak2(i,j,bi,bj) = ak2(i,j,bi,bj) *
     &        exp((16.4-0.040*t)*(pressc-1.0)/(83.143*tk))

c kb = [H][BO2]/[HBO2]
c Millero p.669 (1995) using data from dickson (1990)
             akb(i,j,bi,bj)=exp((-8966.90-2890.53*sqrts-77.942*s +
     &        1.728*s15-0.0996*s2)*invtk +
     &        (148.0248+137.1942*sqrts+1.62142*s) +
     &        (-24.4344-25.085*sqrts-0.2474*s) *
     &        dlogtk+0.053105*sqrts*tk)

c Mick and Karsten - Dec 04
c ADDING pressure dependence based on Millero (1995), p675
c with additional info from CO2sys documentation (E. Lewis and 
c D. Wallace, 1998 - see endnotes for commentary on Millero, 95)
             bigR = 83.145 
             dv = -29.48+0.1622*t+2.608d-3*t*t
             dk = -2.84d-3
             pfactor = - (dv/(bigR*tk))*pressc 
     &        + (0.5*dk/(bigR*tk))*pressc*pressc

             akb(i,j,bi,bj) = akb(i,j,bi,bj)*exp(pfactor)

c k1p = [H][H2PO4]/[H3PO4]
c DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
             ak1p(i,j,bi,bj) = exp(-4576.752*invtk+115.525 - 
     &        18.453*dlogtk +
     &        (-106.736*invtk+0.69171)*sqrts +
     &        (-0.65643*invtk-0.01844)*s)

c k2p = [H][HPO4]/[H2PO4]
c DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))
             ak2p(i,j,bi,bj) = exp(-8814.715*invtk+172.0883 - 
     &        27.927*dlogtk +
     &        (-160.340*invtk+1.3566)*sqrts +
     &        (0.37335*invtk-0.05778)*s)

c k3p = [H][PO4]/[HPO4]
c DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
             ak3p(i,j,bi,bj) = exp(-3070.75*invtk-18.141 + 
     &        (17.27039*invtk+2.81197) *
     &        sqrts+(-44.99486*invtk-0.09984)*s)

c ksi = [H][SiO(OH)3]/[Si(OH)4]
c Millero p.671 (1995) using data from Yao and Millero (1995)
             aksi(i,j,bi,bj) = exp(-8904.2*invtk+117.385 - 
     &        19.334*dlogtk +
     &        (-458.79*invtk+3.5913)*sqrtis +
     &        (188.74*invtk-1.5998)*is +
     &        (-12.1652*invtk+0.07871)*is2 +
     &        log(1.0-0.001005*s))

c kw = [H][OH]
c Millero p.670 (1995) using composite data
             akw(i,j,bi,bj) = exp(-13847.26*invtk+148.9652 - 
     &        23.6521*dlogtk + 
     &        (118.67*invtk-5.977+1.0495*dlogtk) *
     &        sqrts-0.01615*s)

c ks = [H][SO4]/[HSO4]
c dickson (1990, J. chem. Thermodynamics 22, 113)
             aks(i,j,bi,bj)=exp(-4276.1*invtk+141.328 - 
     &        23.093*dlogtk +
     &        (-13856*invtk+324.57-47.986*dlogtk)*sqrtis +
     &        (35474*invtk-771.54+114.723*dlogtk)*is -
     &        2698*invtk*is**1.5+1776*invtk*is2 +
     &        log(1.0-0.001005*s))

c kf = [H][F]/[HF]
c dickson and Riley (1979) -- change pH scale to total
             akf(i,j,bi,bj)=exp(1590.2*invtk-12.641+1.525*sqrtis +
     &        log(1.0-0.001005*s) +
     &        log(1.0+(0.1400/96.062)*(scl)/aks(i,j,bi,bj)))

c Calculate concentrations for borate, sulfate, and fluoride
c Uppstrom (1974)
             bt(i,j,bi,bj) = 0.000232*scl/10.811

c Morris & Riley (1966)
             st(i,j,bi,bj) = 0.14*scl/96.062

c Riley (1965)
             ft(i,j,bi,bj) = 0.000067*scl/18.9984

c solubility product for calcite
C Following Mucci (1983) - from Zeebe/Wolf-Gladrow equic.m
             tmpa1 = -171.9065-(0.077993*tk)+(2839.319/tk) +
     &        (71.595*log10(tk))
         
             tmpa2 = +(-0.77712+(0.0028426*tk)+(178.34/tk))*sqrts
             tmpa3 = -(0.07711*s)+(0.0041249*s15)
             logKspc = tmpa1+tmpa2+tmpa3
             Ksp_T_Calc = 10.0**logKspc

c alternative pressure dependence from Ingle (1975)
             zdum = (pressc*10.0d0-10.0d0)/10.0d0
             xvalue = ((48.8d0-0.53d0*t)*zdum +
     &        (-0.00588d0+0.0001845d0*t)*zdum*zdum) /
     &        (188.93d0*(t+273.15d0))

             Ksp_TP(i,j,bi,bj) = Ksp_T_Calc*10**(xvalue)

            else
             ff(i,j,bi,bj)=0.d0
             ak0(i,j,bi,bj)= 0.d0
             ak1(i,j,bi,bj)= 0.d0
             ak2(i,j,bi,bj)= 0.d0
             akb(i,j,bi,bj)= 0.d0
             ak1p(i,j,bi,bj) = 0.d0
             ak2p(i,j,bi,bj) = 0.d0
             ak3p(i,j,bi,bj) = 0.d0
             aksi(i,j,bi,bj) = 0.d0
             akw(i,j,bi,bj) = 0.d0
             aks(i,j,bi,bj)= 0.d0
             akf(i,j,bi,bj)= 0.d0
             bt(i,j,bi,bj) = 0.d0
             st(i,j,bi,bj) = 0.d0
             ft(i,j,bi,bj) = 0.d0
             Ksp_TP(i,j,bi,bj) = 0.d0
            endif

c compute CO3 for DARWIN_PLANKTON and sediment fluxes
            if (maskC(i,j,k,bi,bj).NE.0. _d 0) then

             pCO2SolverTemp = max(-4. _d 0, min(39. _d 0, Tlocal))
             pCO2SolverSal = max(4. _d 0, min(50. _d 0, Slocal))
              
c check bounds for PCO2 solver
             dicl  = max(Ptracer(i,j,k,bi,bj,iDIC  ),0. _d 0)
             docl  = max(Ptracer(i,j,k,bi,bj,iDOC  ),0. _d 0)
             pocl  = max(Ptracer(i,j,k,bi,bj,iPOC  ),0. _d 0)
             picl  = max(Ptracer(i,j,k,bi,bj,iPIC  ),0. _d 0)
             alkl  = max(Ptracer(i,j,k,bi,bj,iALK  ),0. _d 0)
             o2l   = max(Ptracer(i,j,k,bi,bj,iO2   ),0. _d 0)
             po4l  = max(Ptracer(i,j,k,bi,bj,iPO4  ),0. _d 0)
             sil   = max(Ptracer(i,j,k,bi,bj,iSiO2 ),0. _d 0)

 	         pCO2SolverDic = max(400. _d 0, min(4000. _d 0,dicl))
             pCO2SolverAlk = max(400. _d 0, min(4000. _d 0,alkl))
             pCO2SolverPo4 = max(1. _d -10, min(10. _d 0,po4l))
             pCO2SolverSi = max(1. _d -8, min(500. _d 0,sil))
             
             print *,'pCO2SolverDic =',pCO2SolverDic
             print *,'pCO2SolverAlk =',pCO2SolverAlk
             print *,'pCO2SolverPo4 =',pCO2SolverPo4
             print *,'pCO2SolverSi =',pCO2SolverSi
c convert to mol m^-3
             pCO2SolverDic = pCO2SolverDic*1.0 _d -3
             pCO2SolverAlk = pCO2SolverAlk*1.0 _d -3
             pCO2SolverPo4 = pCO2SolverPo4*1.0 _d -3
             pCO2SolverSi = pCO2SolverSi*1.0 _d -3 
            
 
             CALL DARWIN_CALC_PCO2_APPROX(
     I        pCO2SolverTemp,pCO2SolverSal,
     I        pCO2SolverDic,pCO2SolverPo4,
     I        pCO2SolverSi,pCO2SolverAlk,
     I        ak1(i,j,bi,bj),ak2(i,j,bi,bj),
     I        ak1p(i,j,bi,bj),ak2p(i,j,bi,bj),ak3p(i,j,bi,bj),
     I        aks(i,j,bi,bj),akb(i,j,bi,bj),akw(i,j,bi,bj),
     I        aksi(i,j,bi,bj),akf(i,j,bi,bj),
     I        ak0(i,j,bi,bj),fugf(i,j,bi,bj),
     I        ff(i,j,bi,bj),
     I        bt(i,j,bi,bj),st(i,j,bi,bj),ft(i,j,bi,bj),
     U        pH(i,j,k,bi,bj),pCO2(i,j,k,bi,bj),CO3(i,j,k,bi,bj),
     I        myThid)
	    
             else

              pH(i,j,k,bi,bj) = 0. _d 0
              pCO2(i,j,k,bi,bj) = 0. _d 0
              CO3(i,j,k,bi,bj) = 0. _d 0
        
	     endif

      print *,'darwin_forcing Ksp_TP(i,j,bi,bj) = ',Ksp_TP(i,j,bi,bj)
      
      calcium(i,j,k,bi,bj) = 1.028 _d -2*Slocal/35. _d 0
     &        * 1. _d 3
      print *,'darwin_forcing picl  = ',picl
      print *,'darwin_forcing calcium(i,j,k,bi,bj) =',
     &   calcium(i,j,k,bi,bj) 
      print *,'darwin_forcing CO3(i,j,bi,bj)=',CO3(i,j,k,bi,bj)
      
#endif
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
#ifdef DARWIN_ALLOW_CARBON
     I                 picl,
     I                 Ksp_TP(i,j,bi,bj), 
     I                 CO3(i,j,k,bi,bj),
     I                 calcium(i,j,k,bi,bj),
     O                 disscPIC(i, j, k, bi, bj),
     O                 omegaC(i, j, k, bi, bj),
#endif   
#ifdef DARWIN_DEBUG
     I                 myxgloballo+(bi-1)*sNx+i-1,
     I                 myygloballo+(bj-1)*sNy+j-1,
     I                 k, dTsub(k),
#endif
     I                 subTime, myIter, myThid )
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

#ifdef DARWIN_ALLOW_CARBON
        print *,'darwin_forcing disscPIC(i, j, k, bi, bj) = ',
     &    disscPIC(i, j, k, bi, bj)
        print *,'darwin_forcing omegaC(i, j, k, bi, bj) = ', 
     &    omegaC(i, j, k, bi, bj)
#endif
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
       CALL DIAGNOSTICS_FILL(disscPIC(1-OLx,1-OLy,1,bi,bj),'disscPIC',
     &  0,Nr,2,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(omegaC(1-OLx,1-OLy,1,bi,bj),'omegaC  ',
     &  0,Nr,2,bi,bj,myThid)
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
