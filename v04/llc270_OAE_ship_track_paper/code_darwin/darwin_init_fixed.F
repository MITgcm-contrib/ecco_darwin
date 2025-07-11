C $Header: /home/ubuntu/mnt/e9_copy/MITgcm_contrib/ecco_darwin/v4_llc270/code_darwin/darwin_init_fixed.F,v 1.7 2020/01/08 16:05:31 dcarroll Exp $
C $Name:  $

#include "CPP_OPTIONS.h"
#include "PTRACERS_OPTIONS.h"
#include "DARWIN_OPTIONS.h"

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_DARWIN

c===============================================================================
C===============================================================================
CStartofinterface
      SUBROUTINE DARWIN_INIT_FIXED(myThid)                        
C     =============== Global data ==========================================
C     === Global variables ===
         implicit none
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "GCHEM.h"
#include "DARWIN_SIZE.h"
#include "DARWIN.h"
#include "DARWIN_FLUX.h"

        INTEGER myThid
C============== Local variables ============================================
C     msgBuf    - Informational/error meesage buffer
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      _RL pday
      INTEGER i,j,k,bi,bj,nz
      INTEGER tmp
      INTEGER prec
      CHARACTER*(MAX_LEN_MBUF) fn
C     /--------------------------------------------------------------\
C     |  initialise common block biochemical parameters               |
C     \--------------------------------------------------------------/

      WRITE(msgBuf,'(A)')
     &'// ======================================================='
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') '// Darwin loading parameters'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)')
     &'// ======================================================='
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )


c ANNA set fixed params for WAVEBANDS
#ifdef WAVEBANDS
        call wavebands_init_fixed(myThid)
#endif 

c c define 1 day in seconds
         pday = 86400.0 _d 0


c attenuation coefficients
c ANNA only if not wavebands
#ifndef WAVEBANDS
      k0=        4. _d -2                   !atten coefficient water(m^-1)
#ifdef GEIDER
      kc=        4. _d -2                   !atten coefficient chl ((mmol chl/m3)-1)
#else
      kc=        4. _d -2*16. _d 0*1. _d 0        !atten coefficient phy((uM m)-1)
#endif
#endif
c ANNA endif


c par parameters
      parfrac=  0.4 _d 0                    !fraction Qsw that is PAR
      parconv= 1. _d 0/0.2174 _d 0             !conversion from W/m2 to uEin/m2/s
c for chl
      chlpmax=40. _d 0                      ! mg Chl/mmolP
      chlpmin=16. _d 0                      ! mg Chl/mmolP
      istar=90. _d 0                        ! w/m2
c iron related
C ECCO-Darwin JAMES paper value
      alpfe=     0.9273 _d 0                  !solubility of Fe dust
      scav=      0.40 _d 0/(360. _d 0*86400. _d 0) !iron chem scavenging rate (s-1)
      ligand_tot=1. _d -3                   !total ligand (uM)
      ligand_stab=2. _d 5                   !ligand stability rate ratio
cswd  freefemax = .4 _d -3                  ! max free fe
      freefemax = .1 _d -3                  ! max free fe

#ifdef IRON_SED_SOURCE
c iron sediment source
cswd  depthfesed=1000.0 _d 0       !depth above which to add sediment source
      depthfesed=7000.0 _d 0       !depth above which to add sediment source
      fesedflux =1.0 _d 0 * 1.0 _d -3 / (86400.0 _d 0) !iron flux (mmol/m2/s)
      fesedflux_pcm =0.68 _d 0 * 1.0 _d -3 !iron flux (mmol/m3/s) per
c                                            mmol POC/m3/s
#endif
#ifdef PART_SCAV
cswd  scav_rat=0.005 _d 0 /(86400.0 _d 0)
C ECCO-Darwin JAMES paper value
      scav_rat= 9.3208 _d 0 * 0.005 _d 0 /(86400.0 _d 0)
      scav_inter=0.079 _d 0
      scav_exp=0.58 _d 0
#endif

c depth for denitrification to start
      depthdenit=185.0 _d 0

c critical oxygen for O2/NO3 remineralization
      O2crit = 6.0 _d 0    !(Lipschultz et al 1990, DSR 37, 1513-1541)
c ratio of no3 to p in denitrification process
      denit_np = 120.0 _d 0

c
c oxidation rates for ammonium and nitrite
c i.e. Knita ...    NH4 -> NO2
c i.e. Knitb ...    NO2 -> NO3
       Knita = 1.0 _d 0/(.50 _d 0*pday)
       Knitb = 1.0 _d 0/(10.0 _d 0*pday)
c critical light level (muEin/m2/s) after which oxidation starts
       PAR0 = 10. _d 0
c
#ifndef GEIDER
c set growth days ...small or big organism?
C ECCO-Darwin JAMES paper value
       Smallgrow = 0.6918 _d 0
       Biggrow   = 0.4226 _d 0
       Smallgrowrange = 0. _d 0
       Biggrowrange = 0. _d 0
       diaz_growfac = 2. _d 0
#endif
c set mort days ...small or big organism?
       Smallmort = 10. _d 0
       Bigmort   = 10. _d 0
       Smallmortrange = 0. _d 0
       Bigmortrange = 0. _d 0
c set export fraction ...small or big organism?
cswd   Smallexport = 0.2 _d 0
       Smallexport = 0.05 _d 0
cswd   Bigexport   = 0.5 _d 0
       Bigexport   = 0.3 _d 0
c set sinking rates (m/s)... small or big organism?
       SmallSink = 0.0 _d 0/pday
       BigSink =  0.5 _d 0/pday   !0.5 _d 0/pday
c set parameters for light function for phyto growth
#ifndef GEIDER
       smallksatpar =    0.12 _d 0 ! 0.8 _d 0
       smallksatparstd = 0.20 _d 0  ! 0.3 _d 0
cswd   smallkinhib =     6.0 _d 0   ! 2.0 _d 0
       smallkinhib =     3.0 _d 0
       smallkinhibstd =  0.10 _d 0  ! 0.5 _d 0
       Bigksatpar =      0.12 _d 0  ! 0.35 _d 0
       Bigksatparstd =   0.06 _d 0  ! 0.1 _d 0
       Bigkinhib =       1.0 _d 0   ! 0.5 _d 0
       Bigkinhibstd =    0.05 _d 0   ! 0.1 _d 0
#endif
#ifdef  GEIDER
c for Pcm -- should be growth rates, but using old variables
c note these are in terms of days - converted to 1/s later
C ECCO-Darwin JAMES paper value
       Smallgrow = 0.6918 _d 0
       Biggrow   = 0.4226 _d 0
       Smallgrowrange = 0. _d 0
       Biggrowrange = 0. _d 0
       diaz_growfac = 2. _d 0
c
       smallchl2cmax = 0.2 _d 0      !mg Chl (mmol C)
       smallchl2cmaxrange = 0.3 _d 0  !mg Chl (mmol C)
       Bigchl2cmax = 0.5 _d 0         !mg Chl (mmol C)
       Bigchl2cmaxrange = 0.3 _d 0    !mg Chl (mmol C)

c ANNA_Q units for alpha are same as expected: mmol C (mg chla)-1 (uEin)-1 (m)2
c       smallalphachl =  1. _d -6     !mmol C (uEin/m-2)-1 (mg Chl)-1 
c       smallalphachlrange = 1. _d -6 !mmol C (uEin/m-2)-1 (mg Chl)-1
c       Bigalphachl = 6. _d -7        !mmol C (uEin/m-2)-1 (mg Chl)-1
c       Bigalphachlrange = 4. _d -7   !mmol C (uEin/m-2)-1 (mg Chl)-1
c ANNA mQyield vals are from alphachl / aphy_chl which for now is 0.02
c ANNA ranges for mQyield are same as alphachl but reduced by factor 100
       smallmQyield = 5. _d -5        !mmol C (uEin)-1
       smallmQyieldrange = 1. _d -4   !mmol C (uEin)-1
       BigmQyield = 3. _d -5          !mmol C (uEin)-1
       BigmQyieldrange = 4. _d -5     !mmol C (uEin)-1

c ANNA value of aphy_chl_ave = 0.02 - its the mean of all spectras used as input data
       aphy_chl_ave = 0.02 _d 0      !m2 (mg chla)-1 (ie. x chla gives absorption m-1)

c inhib for Prochl?
C       inhibcoef_geid_val = 1.2 _d 0 !DUMMY VAL
       inhibcoef_geid_val = 0 _d 0 !DUMMY VAL
#ifdef DYNAMIC_CHL
       acclimtimescl = 1./(60. _d 0 *60. _d 0 *24. _d 0 * 20. _d 0)
#endif
#endif
c

c set temperature function
       tempcoeff1 = 1. _d 0/3. _d 0
       tempcoeff2_small = 0.001 _d 0
       tempcoeff2_big   = 0.0003 _d 0
       tempcoeff3 = 1.04 _d 0
       tempmax = 30. _d 0      ! 32. _d 0
       temprange = 32. _d 0    ! 30. _d 0
       tempnorm = 0.3 _d 0  ! 1. _d 0
       tempdecay = 4. _d 0
c set phosphate half stauration constants .. small or big organism
       SmallPsat=0.015 _d 0
       BigPsat=0.035 _d 0
       ProcPsat=0.01 _d 0
       UniDzPsat=0.012 _d 0
       SmallPsatrange=0.02 _d 0
       BigPsatrange=0.02 _d 0
       ProcPsatrange=0.005 _d 0
       UniDzPsatrange=0.02 _d 0
c set NH4/NO2 frac, so that NH4/NO2 can be preferred nitrogen source
       ksatNH4fac=.50 _d 0
       ksatNO2fac=1.0 _d 0
c set prochl lower half-sat (used only for mutants)
       prochlPsat=.85 _d 0
c ammonia and nitrite inhibition
       sig1 = 4.6 _d 0
       sig2 = 4.6 _d 0
       sig3 = 4.6 _d 0
       ngrowfac = 1. _d 0
       ilight = 2. _d 0
c set si half sat
       val_ksatsi=1. _d 0
c set nutrient ratios for phyto
       val_R_SiP_diatom=16.0 _d 0  ! 32 for Fanny's runs
       val_R_NP=16.0 _d 0
       val_RFeP=1.0 _d -3
       val_R_NP_diaz=40.0 _d  0
       val_RFeP_diaz=30.0 _d 0 * val_RFeP
       val_R_PC=120.0 _d 0
       val_R_PICPOC=0.8 _d 0
#ifdef OLD_GRAZE
c grazing hlaf saturation
       kgrazesat = 0.1 _d 0
c set grazing rates .. small or big organism?
       GrazeFast = 1.0 _d 0/(5.0 _d 0*pday)
       GrazeSlow = 1.0 _d 0/(30.0 _d 0*pday)
c set grazing effeciency
       GrazeEffsmall=0.6 _d 0
       GrazeEffbig  =0.2 _d 0
c set grazing of diatom factor
       diatomgraz = 0.9 _d 0
       coccograz = 0.7 _d 0
       olargegraz = 0.7 _d 0
#else
c grazing hlaf saturation
c      kgrazesat = 0.1 _d 0
cswd   kgrazesat = 0.1 _d 0
       kgrazesat = 0.085 _d 0
c      phygrazmin = 1 _d -5
       phygrazmin = 1 _d -10
c set grazing rates .. small or big organism?
c      GrazeFast = 1.0 _d 0/(5.0 _d 0*pday)
cswd   GrazeFast = 1.0 _d 0/(2.0 _d 0*pday)
       GrazeFast = 1.25 _d 0 * 1.0 _d 0/(2.0 _d 0*pday)
c      GrazeSlow = 1.0 _d 0/(30.0 _d 0*pday)
cswd   GrazeSlow = 1.0 _d 0/(7.0 _d 0*pday)
       GrazeSlow = 1.25 _d 0 * 1.0 _d 0/(7.0 _d 0*pday)
c set grazing effeciency
       GrazeEfflow= 0.2 _d 0
       GrazeEffmod= 0.5 _d 0
       GrazeEffhi = 0.7 _d 0
c set palatibility
       palathi = 1.0 _d 0
       palatlo = 0.2 _d 0
c set palatibilty diatom factor
cswd   diatomgraz = 0.7 _d 0
C ECCO-Darwin JAMES paper value
       diatomgraz = 0.8456 _d 0
cswd   coccograz = 0.6 _d 0
       coccograz = 0.85 _d 0
cswd   olargegraz = 1.0 _d 0
       olargegraz = 0.90 _d 0
c set faction graz to POM
cswd   ExGrazfracbig = 0.8 _d 0
       ExGrazfracbig = 0.6 _d 0
cswd   ExGrazfracsmall = 0.8 _d 0
       ExGrazfracsmall = 0.3 _d 0
#endif
c set zoo mortality
       ZoomortSmall = 1.0 _d 0/(30.0 _d 0*pday)
       ZoomortBig = 1.0 _d 0/(30.0 _d 0*pday)
c set zoo exportfrac
cswd   ZooexfacSmall = 0.2 _d 0
       ZooexfacSmall = 0.1 _d 0
cswd   ZooexfacBig = 0.7 _d 0
       ZooexfacBig = 0.5 _d 0
c minimum phyto (below which grazing and mortality doesn't happen)
c      phymin = 1 _d -10
c      phymin = 1 _d -50
       phymin = 1 _d -20
c DOM remin rates
cswd     Kdop = 1.0 _d 0/(100.0 _d 0*pday)
         Kdop = 2.0 _d 0 * 1.0 _d 0/(100.0 _d 0*pday)
cswd     Kdon = 1.0 _d 0/(100.0 _d 0*pday)
         Kdon = 2.0 _d 0 * 1.0 _d 0/(100.0 _d 0*pday)
cswd     KdoFe = 1.0 _d 0/(100.0 _d 0*pday)
         KdoFe = 2.0 _d 0 * 1.0 _d 0/(100.0 _d 0*pday)
c Particulate detritus remin rates
c z* = wx_sink/Kremin_X
c for e-folding length scale, z* = 300 m
c choose Kremin_X = 1/30 day-1, wx_sink = 10 m day-1
cswd     Kpremin_P = 1.0 _d 0/(50.0 _d 0*pday)
         Kpremin_P = 2.0 _d 0/(50.0 _d 0*pday)
         Kpremin_N = Kpremin_P
         Kpremin_Fe = Kpremin_P
cswd     Kpremin_Si = 1.0 _d 0/(300.0 _d 0*pday)
         Kpremin_Si = 2.0 _d 0/(300.0 _d 0*pday)
c sinking rate for particulate matter (m/s)
         wp_sink = 10.0 _d 0/pday
         wn_sink = wp_sink
         wfe_sink = wp_sink
         wsi_sink = wp_sink

#ifdef ALLOW_CARBON
         R_OP = 170 _d 0
         Kdoc = 1.0 _d 0/(100.0 _d 0*pday)
         Kpremin_C = 1.0 _d 0/(50.0 _d 0*pday)
         Kdissc =  1.0 _d 0/(300.0 _d 0*pday)
         wc_sink = wp_sink
         wpic_sink = 15.0 _d 0/pday
         permil      = 1. _d 0 / 1024.5 _d 0
         Pa2Atm      = 1.01325 _d 5
#endif

C make sure we have reserved enough space in Ptracers
      IF ( nCompZooMax .LT. 4 ) THEN
        WRITE(msgBuf,'(A,A,I3)')
     &    'DARWIN_INIT_FIXED: ERROR: 4 zooplankton components, but ',
     &    'nCompZooMax = ', nCompZooMax
        CALL PRINT_ERROR( msgBuf , 1)
        STOP 'ABNORMAL END: S/R DARWIN_INIT_FIXED'
      ENDIF
      DO nz = 1,nzmax
        iZooP (nz)  = iZoo +                   (nz-1)*strideTypeZoo
        iZooN (nz)  = iZoo + 1*strideCompZoo + (nz-1)*strideTypeZoo
        iZooFe(nz)  = iZoo + 2*strideCompZoo + (nz-1)*strideTypeZoo
        iZooSi(nz)  = iZoo + 3*strideCompZoo + (nz-1)*strideTypeZoo
      ENDDO
#ifdef ALLOW_CARBON
       DO nz = 1,nzmax
        iZooC (nz)  = iZoC +                   (nz-1)
      ENDDO
#endif

#ifdef CO2_FLUX_BUDGET
       budgetTStep1 = 0. _d 0
       budgetPert = 1.0 _d -2
#endif

#ifdef DAR_DIAG_DIVER
c only look at grid point with a minimum biomass
        diver_thresh0=1 _d -12
c diver1 - if any type greater than
        diver_thresh1=1 _d -8
c diver2 - if more than this proportion of total biomass
        diver_thresh2=1 _d -3
c diver3 - fraction of biomass to count
        diver_thresh3=.999 _d 0
c diver4 - fraction of maximum species
        diver_thresh4=1 _d -5
#endif

c set up diagnostics
#ifdef ALLOW_MNC
      IF ( useMNC ) THEN
        CALL DARWIN_MNC_INIT( myThid )
#ifdef ALLOW_CARBON
        CALL DIC_MNC_INIT( myThid )
#endif
      ENDIF
#endif /* ALLOW_MNC */

COJ set up diagnostics
#ifdef ALLOW_DIAGNOSTICS
      IF ( useDIAGNOSTICS ) THEN
        CALL DARWIN_DIAGNOSTICS_INIT( myThid )
#ifdef ALLOW_CARBON
        CALL DIC_DIAGNOSTICS_INIT( myThid )
#endif
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */
COJ



        RETURN
        END
C============================================================================
#endif
#endif
