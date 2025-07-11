C $Header: /home/ubuntu/mnt/e9_copy/MITgcm_contrib/ecco_darwin/v4_llc270/code_darwin/DARWIN_IO.h,v 1.3 2019/08/26 05:33:09 dcarroll Exp $
C $Name:  $

#include "DARWIN_OPTIONS.h"

c DARWIN_IO.h
c --------------------------------------------------------------------
c  Description:  input and diagnostic output for DARWIN pkg
c          Stephanie Dutkiewicz: Spring 2006
c NOTE: additional diagnostics still needed
c --------------------------------------------------------------------      

C--   COMMON /DARWIN_FILENAMES/
C  darwin_iceFile       :: file name of seaice fraction
C  darwin_ironFile      :: file name of aeolian iron flux
C  darwin_PARFile       :: file name of Photosynthetically Active Radiation at surface
C  darwin_nutWVelFile   :: file name of nutrient wvel
C  darwin_waterabsorbFile    :: file name of water absorption/scattering spectra
C  darwin_phytoabsorbFile    :: file name of phyto absorption/scattering spectra
C  darwin_particleabsorbFile :: file name of particle absorption/scattering spectra
C  darwin_surfacespecFile    :: file name of spectrum of incident light
C  darwin_acdomFile          :: file name of CDOM absorption spectrum
C  darwin_forcingPeriod :: perioDARWIN forcing parameter specific for DARWIN (seconds)
C  darwin_forcingCycle  :: perioDARWIN forcing parameter specific for DARWIN (seconds)
C  darwin_seed          :: seed for the random number generator
C  darwin_dicSurfFluxFile    :: file name of dic surface flux 
C  darwin_alkSurfFluxFile    :: file name of alk surface flux
C  darwin_caSurfFluxFile    :: file name of calcium surface flux
C  darwin_BBLFile            :: file name of BBL file
 
      COMMON /DARWIN_FILENAMES/
     &        darwin_iceFile,
     &        darwin_ironFile,
     &        darwin_PARFile,
     &        darwin_nutWVelFile,
     &        darwin_PO4_relaxFile, darwin_NO3_relaxFile,
     &        darwin_FeT_relaxFile, darwin_Si_relaxFile,
     &        darwin_relaxscale,
     &        darwin_PO4_fluxFile, darwin_NO3_FluxFile,
     &        darwin_FeT_fluxFile, darwin_Si_fluxFile,
     &        darwin_oasim_edFile, darwin_oasim_esFile,
     &        darwin_waterabsorbFile, darwin_phytoabsorbFile,
     &        darwin_surfacespecFile, darwin_acdomFile,
     &        darwin_particleabsorbFile,
     &        darwin_forcingPeriod, darwin_forcingCycle,
     &        darwin_seed,
     &        darwin_dicSurfFluxFile,
     &        darwin_alkSurfFluxFile,
     &        darwin_caSurfFluxFile,
     &        darwin_BBLFile

      CHARACTER*(MAX_LEN_FNAM) darwin_iceFile
      CHARACTER*(MAX_LEN_FNAM) darwin_ironFile
      CHARACTER*(MAX_LEN_FNAM) darwin_PARFile
      CHARACTER*(MAX_LEN_FNAM) darwin_NutWVelFile
      CHARACTER*(MAX_LEN_FNAM) darwin_PO4_relaxFile
      CHARACTER*(MAX_LEN_FNAM) darwin_NO3_relaxFile
      CHARACTER*(MAX_LEN_FNAM) darwin_FeT_relaxFile
      CHARACTER*(MAX_LEN_FNAM) darwin_Si_relaxFile
      CHARACTER*(MAX_LEN_FNAM) darwin_PO4_fluxFile
      CHARACTER*(MAX_LEN_FNAM) darwin_NO3_fluxFile
      CHARACTER*(MAX_LEN_FNAM) darwin_FeT_fluxFile
      CHARACTER*(MAX_LEN_FNAM) darwin_Si_fluxFile
      CHARACTER*(MAX_LEN_FNAM) darwin_oasim_edFile
      CHARACTER*(MAX_LEN_FNAM) darwin_oasim_esFile
      CHARACTER*(MAX_LEN_FNAM) darwin_waterabsorbFile
      CHARACTER*(MAX_LEN_FNAM) darwin_phytoabsorbFile
      CHARACTER*(MAX_LEN_FNAM) darwin_particleabsorbFile
      CHARACTER*(MAX_LEN_FNAM) darwin_surfacespecFile
      CHARACTER*(MAX_LEN_FNAM) darwin_acdomFile
      CHARACTER*(MAX_LEN_FNAM) darwin_dicSurfFluxFile
      CHARACTER*(MAX_LEN_FNAM) darwin_alkSurfFluxFile
      CHARACTER*(MAX_LEN_FNAM) darwin_caSurfFluxFile
      CHARACTER*(MAX_LEN_FNAM) darwin_BBLFile

      _RL darwin_relaxscale
      _RL     darwin_forcingPeriod
      _RL     darwin_forcingCycle
      INTEGER darwin_seed

c INPUT:
c    fice           - ice fraction
c    inputFe        - aeolian input of iron
c    sur_par        - surface PAR
      COMMON /DARWIN_INPUT/
     &       fice, inputFe, sur_par
#ifdef NUT_SUPPLY
     &      ,nut_wvel
#endif
#ifdef RELAX_NUTS
     &      ,po4_obs, no3_obs, fet_obs, si_obs
#endif
#ifdef FLUX_NUTS
     &      ,po4_flx, no3_flx, fet_flx, si_flx
#endif
#ifdef ADKINS_SURF_FLUX
     &      ,dicSurf_flx, alkSurf_flx, caSurf_flx
#endif
c
       _RL   fice(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx, nSy)
       _RL   inputFe(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx, nSy)
       _RL   sur_par(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx, nSy)
#ifdef NUT_SUPPLY
       _RL   nut_wvel(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx, nSy)
#endif
#ifdef RELAX_NUTS
       _RL   po4_obs(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx, nSy)
       _RL   no3_obs(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx, nSy)
       _RL   fet_obs(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx, nSy)
       _RL   si_obs(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx, nSy)
#endif
#ifdef FLUX_NUTS
       _RL   po4_flx(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx, nSy)
       _RL   no3_flx(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx, nSy)
       _RL   fet_flx(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx, nSy)
       _RL   si_flx(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx, nSy)
#endif
#ifdef ADKINS_SURF_FLUX  
       _RL   dicSurf_flx(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
       _RL   alkSurf_flx(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
       _RL   caSurf_flx(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
#endif
#ifdef ALLOW_SED_DISS_FLUX
       _RL   BBLThickness(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
#endif
c
c OUPUT DIAGNOSTICS
c    PPave        - average primary production
c    Nfixave      - average N fixation
c    Zoograzave   - average zooplankton consumption
c    Parave       - average PAR
c    Chlave       - average chlorophyll
c    npzd_timeave  - time for averaging

      COMMON /DARWIN_OUPUT/
     &      PPave, Nfixave, Zoograzave,
     &      PARave, Chlave, Denitave, 
c ANNA_TAVE
#ifdef WAVES_DIAG_PCHL
     &      Pchlave,
#endif
#ifdef DAR_DIAG_ACDOM
     &      aCDOMave,
#endif
c ANNA end TAVE
#ifdef DAR_DIAG_RSTAR
     &      Rstarave, RNstarave,
#endif
#ifdef DAR_DIAG_DIVER
     &      Diver1ave, Diver2ave, Diver3ave, Diver4ave,
#endif
#ifdef DAR_DIAG_GROW
     &      Growave, Growsqave,
#endif
#ifdef ALLOW_DIAZ
#ifdef DAR_DIAG_NFIXP
     &      NfixPave,
#endif
#endif
     &      dar_timeave
c  
       _RL  PPave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx, nSy)
       _RL  Nfixave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx, nSy)
       _RL  Zoograzave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx, nSy)
       _RL  PARave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx, nSy)
       _RL  Chlave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx, nSy)
       _RL  Denitave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx, nSy)
c ANNA_TAVE
#ifdef WAVES_DIAG_PCHL
       _RL  Pchlave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx, nSy,npmax)
#endif
#ifdef DAR_DIAG_ACDOM
       _RL  aCDOMave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx,nSy)
#endif
c ANNA end TAVE
#ifdef DAR_DIAG_RSTAR
       _RL  Rstarave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx, nSy,npmax)
       _RL  RNstarave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx, nSy,npmax)
#endif
#ifdef DAR_DIAG_DIVER
       _RL  Diver1ave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx, nSy)
       _RL  Diver2ave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx, nSy)
       _RL  Diver3ave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx, nSy)
       _RL  Diver4ave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx, nSy)
#endif
#ifdef DAR_DIAG_GROW
       _RL Growave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx, nSy,npmax)
       _RL Growsqave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx, nSy,npmax)
#endif
#ifdef ALLOW_DIAZ
#ifdef DAR_DIAG_NFIXP
       _RL  NfixPave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx, nSy,npmax)
#endif
#endif
       _RL  DAR_timeave(nSx,nSy)

#ifdef CHECK_CONS
      COMMON /DARWIN_CHECK_CONS/
     &      DAR_cons_unit1, DAR_cons_unit2,
     &      DAR_cons_unit3, DAR_cons_unit4,
     &      DAR_cons_unit5, DAR_cons_unit6,
     &      DAR_cons_unit7
       INTEGER DAR_cons_unit1
       INTEGER DAR_cons_unit2
       INTEGER DAR_cons_unit3
       INTEGER DAR_cons_unit4
       INTEGER DAR_cons_unit5
       INTEGER DAR_cons_unit6
       INTEGER DAR_cons_unit7
#endif

#ifdef ALLOW_PAR_DAY
C     PARday    :: array for accumulating/storing daily-averaged PAR
      COMMON /DARWIN_PAR_DAY/ PARday,
     &      darwin_PARavPeriod, darwin_PARnav
      _RL PARday(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx,nSy,2)
      _RL darwin_PARavPeriod
      INTEGER darwin_PARnav
#endif

#ifdef ALLOW_CARBON
C--   COMMON /DIC_FILENAMES/
C  DIC_windFile    :: file name of wind speeds
C  DIC_atmospFile  :: file name of atmospheric pressure
C  dic_pCO2          :: Atmospheric pCO2 to be rad in data.dic
C  dic_int*          :: place holder to read in a integer number, set at run time

      COMMON /DIC_FILENAMES/
     &        DIC_windFile, DIC_atmospFile, 
     &        dic_pCO2, dic_int1, dic_int2, dic_int3, dic_int4
      CHARACTER*(MAX_LEN_FNAM) DIC_windFile
      CHARACTER*(MAX_LEN_FNAM) DIC_atmospFile
      _RL dic_pCO2
      INTEGER dic_int1
      INTEGER dic_int2
      INTEGER dic_int3
      INTEGER dic_int4
#endif
