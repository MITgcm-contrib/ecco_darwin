#include "DARWIN_OPTIONS.h"
#include "EXF_OPTIONS.h"

CBOP
C !ROUTINE: DARWIN_EXF_READPARMS

C !INTERFACE: ==========================================================
      SUBROUTINE DARWIN_EXF_READPARMS( iUnit, oUnit, myThid )

C !DESCRIPTION:
C     Initialize DARWIN parameters, read in data.darwin

C !USES: ===============================================================
      IMPLICIT NONE
#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#ifdef ALLOW_EXCH2
#include "W2_EXCH2_SIZE.h"
#include "W2_EXCH2_TOPOLOGY.h"
#endif
#ifdef ALLOW_GCHEM
#include "GCHEM.h"
#endif
#ifdef ALLOW_EXF
#include "EXF_PARAM.h"
#include "EXF_INTERP_SIZE.h"
#include "EXF_INTERP_PARAM.h"
#endif
#ifdef ALLOW_RADTRANS
#include "RADTRANS_SIZE.h"
#include "RADTRANS_PARAMS.h"
#endif
#include "DARWIN_SIZE.h"
#include "DARWIN_PARAMS.h"
#include "DARWIN_EXF_PARAMS.h"
#include "DARWIN_EXF_FIELDS.h"

C !INPUT PARAMETERS: ===================================================
C  iUnit   :: unit number for reading
C  oUnit   :: unit number for writing
C  myThid  :: thread number
      INTEGER iUnit, oUnit, myThid

C !OUTPUT PARAMETERS: ==================================================
C  none
CEOP

#ifdef ALLOW_DARWIN

C !LOCAL VARIABLES: ====================================================
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      LOGICAL darwin_haveIce
      LOGICAL darwin_haveWind
      LOGICAL darwin_havepCO2
      LOGICAL darwin_haveDOCrunoff
      LOGICAL darwin_havetDOCrunoff
      LOGICAL darwin_haveDONrunoff
      LOGICAL darwin_haveDOPrunoff
      LOGICAL darwin_haveDSirunoff
      LOGICAL darwin_haveDICrunoff
      LOGICAL darwin_haveALKrunoff
      INTEGER j, errCount

      NAMELIST/darwin_forcing_params/
     &    darwin_useEXFwind,
     &    darwin_useQsw,
     &    darwin_useSEAICE,
C
     &    PARfile,
     &    PARperiod,
     &    PARRepCycle,
     &    PARStartTime,
     &    PARstartdate1,
     &    PARstartdate2,
     &    PARconst,
     &    PAR_exfremo_intercept,
     &    PAR_exfremo_slope,
     &    PARmask,
     &    darwin_inscal_PAR,
C
     &    ironfile,
     &    ironperiod,
     &    ironRepCycle,
     &    ironStartTime,
     &    ironstartdate1,
     &    ironstartdate2,
     &    ironconst,
     &    iron_exfremo_intercept,
     &    iron_exfremo_slope,
     &    ironmask,
     &    darwin_inscal_iron,
C
     &    icefile,
     &    iceperiod,
     &    iceRepCycle,
     &    iceStartTime,
     &    icestartdate1,
     &    icestartdate2,
     &    iceconst,
     &    ice_exfremo_intercept,
     &    ice_exfremo_slope,
     &    icemask,
     &    darwin_inscal_ice,
C
     &    windfile,
     &    windperiod,
     &    windRepCycle,
     &    windStartTime,
     &    windstartdate1,
     &    windstartdate2,
     &    windconst,
     &    wind_exfremo_intercept,
     &    wind_exfremo_slope,
     &    windmask,
     &    darwin_inscal_wind,
C
     &    pCO2file,
     &    pCO2period,
     &    pCO2RepCycle,
     &    pCO2StartTime,
     &    pCO2startdate1,
     &    pCO2startdate2,
     &    pCO2const,
     &    pCO2_exfremo_intercept,
     &    pCO2_exfremo_slope,
     &    pCO2mask,
     &    darwin_inscal_pCO2,
C
     &    DOCrunofffile,
     &    DOCrunoffperiod,
     &    DOCrunoffRepCycle,
     &    DOCrunoffStartTime,
     &    DOCrunoffstartdate1,
     &    DOCrunoffstartdate2,
     &    DOCrunoffconst,
     &    DOCrunoff_exfremo_intercept,
     &    DOCrunoff_exfremo_slope,
     &    DOCrunoffmask,
     &    darwin_inscal_DOCrunoff,
C
     &    tDOCrunofffile,
     &    tDOCrunoffperiod,
     &    tDOCrunoffRepCycle,
     &    tDOCrunoffStartTime,
     &    tDOCrunoffstartdate1,
     &    tDOCrunoffstartdate2,
     &    tDOCrunoffconst,
     &    tDOCrunoff_exfremo_intercept,
     &    tDOCrunoff_exfremo_slope,
     &    tDOCrunoffmask,
     &    darwin_inscal_tDOCrunoff,
C
     &    DONrunofffile,
     &    DONrunoffperiod,
     &    DONrunoffRepCycle,
     &    DONrunoffStartTime,
     &    DONrunoffstartdate1,
     &    DONrunoffstartdate2,
     &    DONrunoffconst,
     &    DONrunoff_exfremo_intercept,
     &    DONrunoff_exfremo_slope,
     &    DONrunoffmask,
     &    darwin_inscal_DONrunoff,
C
     &    DOPrunofffile,
     &    DOPrunoffperiod,
     &    DOPrunoffRepCycle,
     &    DOPrunoffStartTime,
     &    DOPrunoffstartdate1,
     &    DOPrunoffstartdate2,
     &    DOPrunoffconst,
     &    DOPrunoff_exfremo_intercept,
     &    DOPrunoff_exfremo_slope,
     &    DOPrunoffmask,
     &    darwin_inscal_DOPrunoff,
C
     &    DSirunofffile,
     &    DSirunoffperiod,
     &    DSirunoffRepCycle,
     &    DSirunoffStartTime,
     &    DSirunoffstartdate1,
     &    DSirunoffstartdate2,
     &    DSirunoffconst,
     &    DSirunoff_exfremo_intercept,
     &    DSirunoff_exfremo_slope,
     &    DSirunoffmask,
     &    darwin_inscal_DSirunoff,
C
     &    DICrunofffile,
     &    DICrunoffperiod,
     &    DICrunoffRepCycle,
     &    DICrunoffStartTime,
     &    DICrunoffstartdate1,
     &    DICrunoffstartdate2,
     &    DICrunoffconst,
     &    DICrunoff_exfremo_intercept,
     &    DICrunoff_exfremo_slope,
     &    DICrunoffmask,
     &    darwin_inscal_DICrunoff,
C
     &    ALKrunofffile,
     &    ALKrunoffperiod,
     &    ALKrunoffRepCycle,
     &    ALKrunoffStartTime,
     &    ALKrunoffstartdate1,
     &    ALKrunoffstartdate2,
     &    ALKrunoffconst,
     &    ALKrunoff_exfremo_intercept,
     &    ALKrunoff_exfremo_slope,
     &    ALKrunoffmask,
     &    darwin_inscal_ALKrunoff

#ifdef USE_EXF_INTERPOLATION
      NAMELIST/darwin_interp_params/
C
     &    PAR_lon0,
     &    PAR_lat0,
     &    PAR_nlon,
     &    PAR_nlat,
     &    PAR_lon_inc,
     &    PAR_interpMethod,
     &    PAR_lat_inc,
C
     &    iron_lon0,
     &    iron_lat0,
     &    iron_nlon,
     &    iron_nlat,
     &    iron_lon_inc,
     &    iron_interpMethod,
     &    iron_lat_inc,
C
     &    ice_lon0,
     &    ice_lat0,
     &    ice_nlon,
     &    ice_nlat,
     &    ice_lon_inc,
     &    ice_interpMethod,
     &    ice_lat_inc,
C
     &    wind_lon0,
     &    wind_lat0,
     &    wind_nlon,
     &    wind_nlat,
     &    wind_lon_inc,
     &    wind_interpMethod,
     &    wind_lat_inc,
C
     &    pCO2_lon0,
     &    pCO2_lat0,
     &    pCO2_nlon,
     &    pCO2_nlat,
     &    pCO2_lon_inc,
     &    pCO2_interpMethod,
     &    pCO2_lat_inc,
C
     &    DOCrunoff_lon0,
     &    DOCrunoff_lat0,
     &    DOCrunoff_nlon,
     &    DOCrunoff_nlat,
     &    DOCrunoff_lon_inc,
     &    DOCrunoff_interpMethod,
     &    DOCrunoff_lat_inc,
C
     &    tDOCrunoff_lon0,
     &    tDOCrunoff_lat0,
     &    tDOCrunoff_nlon,
     &    tDOCrunoff_nlat,
     &    tDOCrunoff_lon_inc,
     &    tDOCrunoff_interpMethod,
     &    tDOCrunoff_lat_inc,
C
     &    DONrunoff_lon0,
     &    DONrunoff_lat0,
     &    DONrunoff_nlon,
     &    DONrunoff_nlat,
     &    DONrunoff_lon_inc,
     &    DONrunoff_interpMethod,
     &    DONrunoff_lat_inc,
C
     &    DOPrunoff_lon0,
     &    DOPrunoff_lat0,
     &    DOPrunoff_nlon,
     &    DOPrunoff_nlat,
     &    DOPrunoff_lon_inc,
     &    DOPrunoff_interpMethod,
     &    DOPrunoff_lat_inc,
C
     &    DSirunoff_lon0,
     &    DSirunoff_lat0,
     &    DSirunoff_nlon,
     &    DSirunoff_nlat,
     &    DSirunoff_lon_inc,
     &    DSirunoff_interpMethod,
     &    DSirunoff_lat_inc,
C
     &    DICrunoff_lon0,
     &    DICrunoff_lat0,
     &    DICrunoff_nlon,
     &    DICrunoff_nlat,
     &    DICrunoff_lon_inc,
     &    DICrunoff_interpMethod,
     &    DICrunoff_lat_inc,
C
     &    ALKrunoff_lon0,
     &    ALKrunoff_lat0,
     &    ALKrunoff_nlon,
     &    ALKrunoff_nlat,
     &    ALKrunoff_lon_inc,
     &    ALKrunoff_interpMethod,
     &    ALKrunoff_lat_inc
#endif

C Open and read the data.darwin file

      darwin_useEXFwind = .FALSE.
      darwin_useQsw = .FALSE.

      PARfile = ' '
      PARperiod             = 0.0 _d 0
      PARRepCycle           = repeatPeriod
      PARStartTime          = UNSET_RL
      PARstartdate1         = 0
      PARstartdate2         = 0
      PARconst              = UNSET_RL
      PAR_exfremo_intercept = 0.0 _d 0
      PAR_exfremo_slope     = 0.0 _d 0
      PARmask = 'c'

      darwin_inscal_PAR    =  1. _d 0

      ironfile = ' '
      ironperiod             = 0.0 _d 0
      ironRepCycle           = repeatPeriod
      ironStartTime          = UNSET_RL
      ironstartdate1         = 0
      ironstartdate2         = 0
      ironconst              = 0.0 _d 0
      iron_exfremo_intercept = 0.0 _d 0
      iron_exfremo_slope     = 0.0 _d 0
      ironmask = 'c'

      darwin_inscal_iron    =  1. _d 0

      icefile = ' '
      iceperiod             = 0.0 _d 0
      iceRepCycle           = repeatPeriod
      iceStartTime          = UNSET_RL
      icestartdate1         = 0
      icestartdate2         = 0
      iceconst              = UNSET_RL
      ice_exfremo_intercept = 0.0 _d 0
      ice_exfremo_slope     = 0.0 _d 0
      icemask = 'c'

      darwin_inscal_ice    =  1. _d 0

      windfile = ' '
      windperiod             = 0.0 _d 0
      windRepCycle           = repeatPeriod
      windStartTime          = UNSET_RL
      windstartdate1         = 0
      windstartdate2         = 0
      windconst              = UNSET_RL
      wind_exfremo_intercept = 0.0 _d 0
      wind_exfremo_slope     = 0.0 _d 0
      windmask = 'c'

      darwin_inscal_wind    =  1. _d 0

      pCO2file = ' '
      pCO2period             = 0.0 _d 0
      pCO2RepCycle           = repeatPeriod
      pCO2StartTime          = UNSET_RL
      pCO2startdate1         = 0
      pCO2startdate2         = 0
      pCO2const              = UNSET_RL
      pCO2_exfremo_intercept = 0.0 _d 0
      pCO2_exfremo_slope     = 0.0 _d 0
      pCO2mask = 'c'

      darwin_inscal_pCO2    =  1. _d 0

      DOCrunofffile = ' '
      DOCrunoffperiod             = 0.0 _d 0
      DOCrunoffRepCycle           = repeatPeriod
      DOCrunoffStartTime          = UNSET_RL
      DOCrunoffstartdate1         = 0
      DOCrunoffstartdate2         = 0
      DOCrunoffconst              = UNSET_RL
      DOCrunoff_exfremo_intercept = 0.0 _d 0
      DOCrunoff_exfremo_slope     = 0.0 _d 0
      DOCrunoffmask = 'c'

      darwin_inscal_DOCrunoff    =  1. _d 0

      tDOCrunofffile = ' '
      tDOCrunoffperiod             = 0.0 _d 0
      tDOCrunoffRepCycle           = repeatPeriod
      tDOCrunoffStartTime          = UNSET_RL
      tDOCrunoffstartdate1         = 0
      tDOCrunoffstartdate2         = 0
      tDOCrunoffconst              = UNSET_RL
      tDOCrunoff_exfremo_intercept = 0.0 _d 0
      tDOCrunoff_exfremo_slope     = 0.0 _d 0
      tDOCrunoffmask = 'c'

      darwin_inscal_tDOCrunoff    =  1. _d 0

      DONrunofffile = ' '
      DONrunoffperiod             = 0.0 _d 0
      DONrunoffRepCycle           = repeatPeriod
      DONrunoffStartTime          = UNSET_RL
      DONrunoffstartdate1         = 0
      DONrunoffstartdate2         = 0
      DONrunoffconst              = UNSET_RL
      DONrunoff_exfremo_intercept = 0.0 _d 0
      DONrunoff_exfremo_slope     = 0.0 _d 0
      DONrunoffmask = 'c'

      darwin_inscal_DONrunoff    =  1. _d 0

      DOPrunofffile = ' '
      DOPrunoffperiod             = 0.0 _d 0
      DOPrunoffRepCycle           = repeatPeriod
      DOPrunoffStartTime          = UNSET_RL
      DOPrunoffstartdate1         = 0
      DOPrunoffstartdate2         = 0
      DOPrunoffconst              = UNSET_RL
      DOPrunoff_exfremo_intercept = 0.0 _d 0
      DOPrunoff_exfremo_slope     = 0.0 _d 0
      DOPrunoffmask = 'c'

      darwin_inscal_DOPrunoff    =  1. _d 0

      DSirunofffile = ' '
      DSirunoffperiod             = 0.0 _d 0
      DSirunoffRepCycle           = repeatPeriod
      DSirunoffStartTime          = UNSET_RL
      DSirunoffstartdate1         = 0
      DSirunoffstartdate2         = 0
      DSirunoffconst              = UNSET_RL
      DSirunoff_exfremo_intercept = 0.0 _d 0
      DSirunoff_exfremo_slope     = 0.0 _d 0
      DSirunoffmask = 'c'

      darwin_inscal_DSirunoff    =  1. _d 0

      DICrunofffile = ' '
      DICrunoffperiod             = 0.0 _d 0
      DICrunoffRepCycle           = repeatPeriod
      DICrunoffStartTime          = UNSET_RL
      DICrunoffstartdate1         = 0
      DICrunoffstartdate2         = 0
      DICrunoffconst              = UNSET_RL
      DICrunoff_exfremo_intercept = 0.0 _d 0
      DICrunoff_exfremo_slope     = 0.0 _d 0
      DICrunoffmask = 'c'

      darwin_inscal_DICrunoff    =  1. _d 0

      ALKrunofffile = ' '
      ALKrunoffperiod             = 0.0 _d 0
      ALKrunoffRepCycle           = repeatPeriod
      ALKrunoffStartTime          = UNSET_RL
      ALKrunoffstartdate1         = 0
      ALKrunoffstartdate2         = 0
      ALKrunoffconst              = UNSET_RL
      ALKrunoff_exfremo_intercept = 0.0 _d 0
      ALKrunoff_exfremo_slope     = 0.0 _d 0
      ALKrunoffmask = 'c'

      darwin_inscal_ALKrunoff    =  1. _d 0

#ifdef USE_EXF_INTERPOLATION
      PAR_lon0 = inp_lon0
      PAR_lat0 = inp_lat0
      PAR_nlon = inp_gNx
      PAR_nlat = inp_gNy
      PAR_lon_inc = inp_dLon
      PAR_interpMethod  = 1

      iron_lon0 = inp_lon0
      iron_lat0 = inp_lat0
      iron_nlon = inp_gNx
      iron_nlat = inp_gNy
      iron_lon_inc = inp_dLon
      iron_interpMethod  = 1

      ice_lon0 = inp_lon0
      ice_lat0 = inp_lat0
      ice_nlon = inp_gNx
      ice_nlat = inp_gNy
      ice_lon_inc = inp_dLon
      ice_interpMethod  = 1

      wind_lon0 = inp_lon0
      wind_lat0 = inp_lat0
      wind_nlon = inp_gNx
      wind_nlat = inp_gNy
      wind_lon_inc = inp_dLon
      wind_interpMethod  = 1

      pCO2_lon0 = inp_lon0
      pCO2_lat0 = inp_lat0
      pCO2_nlon = inp_gNx
      pCO2_nlat = inp_gNy
      pCO2_lon_inc = inp_dLon
      pCO2_interpMethod  = 1

      DOCrunoff_lon0 = inp_lon0
      DOCrunoff_lat0 = inp_lat0
      DOCrunoff_nlon = inp_gNx
      DOCrunoff_nlat = inp_gNy
      DOCrunoff_lon_inc = inp_dLon
      DOCrunoff_interpMethod  = 1

      tDOCrunoff_lon0 = inp_lon0
      tDOCrunoff_lat0 = inp_lat0
      tDOCrunoff_nlon = inp_gNx
      tDOCrunoff_nlat = inp_gNy
      tDOCrunoff_lon_inc = inp_dLon
      tDOCrunoff_interpMethod  = 1

      DONrunoff_lon0 = inp_lon0
      DONrunoff_lat0 = inp_lat0
      DONrunoff_nlon = inp_gNx
      DONrunoff_nlat = inp_gNy
      DONrunoff_lon_inc = inp_dLon
      DONrunoff_interpMethod  = 1

      DOPrunoff_lon0 = inp_lon0
      DOPrunoff_lat0 = inp_lat0
      DOPrunoff_nlon = inp_gNx
      DOPrunoff_nlat = inp_gNy
      DOPrunoff_lon_inc = inp_dLon
      DOPrunoff_interpMethod  = 1

      DSirunoff_lon0 = inp_lon0
      DSirunoff_lat0 = inp_lat0
      DSirunoff_nlon = inp_gNx
      DSirunoff_nlat = inp_gNy
      DSirunoff_lon_inc = inp_dLon
      DSirunoff_interpMethod  = 1

      DICrunoff_lon0 = inp_lon0
      DICrunoff_lat0 = inp_lat0
      DICrunoff_nlon = inp_gNx
      DICrunoff_nlat = inp_gNy
      DICrunoff_lon_inc = inp_dLon
      DICrunoff_interpMethod  = 1

      ALKrunoff_lon0 = inp_lon0
      ALKrunoff_lat0 = inp_lat0
      ALKrunoff_nlon = inp_gNx
      ALKrunoff_nlat = inp_gNy
      ALKrunoff_lon_inc = inp_dLon
      ALKrunoff_interpMethod  = 1

      DO j=1,MAX_LAT_INC
        PAR_lat_inc(j) = inp_dLat(j)
        iron_lat_inc(j) = inp_dLat(j)
        ice_lat_inc(j) = inp_dLat(j)
        wind_lat_inc(j) = inp_dLat(j)
        pCO2_lat_inc(j) = inp_dLat(j)
        DOCrunoff_lat_inc(j) = inp_dLat(j)
        tDOCrunoff_lat_inc(j) = inp_dLat(j)
        DONrunoff_lat_inc(j) = inp_dLat(j)
        DOPrunoff_lat_inc(j) = inp_dLat(j)
        DSirunoff_lat_inc(j) = inp_dLat(j)
        DICrunoff_lat_inc(j) = inp_dLat(j)
        ALKrunoff_lat_inc(j) = inp_dLat(j)
      ENDDO
#endif /* USE_EXF_INTERPOLATION */

C ======================================================================

      PARconst = UNSET_RL
      iceconst = UNSET_RL
      windconst = UNSET_RL

      READ(UNIT=iUnit, NML=darwin_forcing_params)
#ifdef USE_EXF_INTERPOLATION
      READ(UNIT=iUnit, NML=darwin_interp_params)
#endif

C ======================================================================
C set alternate forcing sources

      DARWIN_haveSurfPAR = PARFile.NE.' ' .OR. PARconst.NE.UNSET_RL
      IF (PARconst .EQ. UNSET_RL) THEN
        PARconst = 0 _d 0
      ENDIF

      DARWIN_haveIce = iceFile.NE.' ' .OR. iceconst.NE.UNSET_RL
      IF (iceconst.EQ.UNSET_RL) THEN
        iceconst = 0 _d 0
      ENDIF

      DARWIN_haveWind = windFile.NE.' ' .OR. windconst.NE.UNSET_RL
      IF (windconst.EQ.UNSET_RL) THEN
        windconst = 5 _d 0
      ENDIF

      DARWIN_havepCO2 = pCO2File.NE.' ' .OR. pCO2const.NE.UNSET_RL
      IF (pCO2const.EQ.UNSET_RL) THEN
        pCO2const = 278 _d -6
      ENDIF

      DARWIN_haveDOCrunoff = DOCrunofffile.NE.' ' .OR.
     & DOCrunoffconst.NE.UNSET_RL
      IF (DOCrunoffconst.EQ.UNSET_RL) THEN
        DOCrunoffconst = 0 _d 0
      ENDIF

      DARWIN_havetDOCrunoff = tDOCrunofffile.NE.' ' .OR.
     & tDOCrunoffconst.NE.UNSET_RL
      IF (tDOCrunoffconst.EQ.UNSET_RL) THEN
        tDOCrunoffconst = 0 _d 0
      ENDIF

      DARWIN_haveDONrunoff = DONrunofffile.NE.' ' .OR.
     & DONrunoffconst.NE.UNSET_RL
      IF (DONrunoffconst.EQ.UNSET_RL) THEN
        DONrunoffconst = 0 _d 0
      ENDIF

      DARWIN_haveDOPrunoff = DOPrunofffile.NE.' ' .OR.
     & DOPrunoffconst.NE.UNSET_RL
      IF (DOPrunoffconst.EQ.UNSET_RL) THEN
        DOPrunoffconst = 0 _d 0
      ENDIF

      DARWIN_haveDSirunoff = DSirunofffile.NE.' ' .OR.
     & DSirunoffconst.NE.UNSET_RL
      IF (DSirunoffconst.EQ.UNSET_RL) THEN
        DSirunoffconst = 0 _d 0
      ENDIF

      DARWIN_haveDICrunoff = DICrunofffile.NE.' ' .OR.
     & DICrunoffconst.NE.UNSET_RL
      IF (DICrunoffconst.EQ.UNSET_RL) THEN
        DICrunoffconst = 0 _d 0
      ENDIF

      DARWIN_haveALKrunoff = ALKrunofffile.NE.' ' .OR.
     & ALKrunoffconst.NE.UNSET_RL
      IF (ALKrunoffconst.EQ.UNSET_RL) THEN
        ALKrunoffconst = 0 _d 0
      ENDIF

      IF (DARWIN_haveSurfPAR .AND. DARWIN_useQsw) THEN
        WRITE(msgBuf,'(2A)') 'DARWIN_EXF_READPARMS: ',
     &    'cannot set both PARFile/PARconst and darwin_useQsw'
       CALL PRINT_ERROR( msgBuf, myThid )
       STOP'ABNORMAL END: DARWIN_EXF_READPARMS: inconsistent light forc'
      ENDIF

      IF (DARWIN_haveIce .AND. DARWIN_useSEAICE) THEN
        WRITE(msgBuf,'(2A)') 'DARWIN_EXF_READPARMS: ',
     &    'cannot set both iceFile/iceconst and darwin_useSEAICE'
       CALL PRINT_ERROR( msgBuf, myThid )
       STOP 'ABNORMAL END: DARWIN_EXF_READPARMS: inconsistent ice forc'
      ENDIF

#ifdef ALLOW_SEAICE
      IF (useSEAICE .AND. .NOT. DARWIN_useSEAICE) THEN
        WRITE(msgBuf,'(2A)') '** WARNING ** DARWIN_EXF_READPARMS: ',
     &    'iceFile or iceconst set, not using ice area from seaice pkg'
        CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
      ENDIF
#endif

      IF (DARWIN_haveWind .AND. DARWIN_useEXFwind) THEN
        WRITE(msgBuf,'(2A)') 'DARWIN_EXF_READPARMS: ',
     &    'cannot set both windFile/windconst and darwin_useEXFwind'
       CALL PRINT_ERROR( msgBuf, myThid )
       STOP 'ABNORMAL END: DARWIN_EXF_READPARMS: inconsistent wind forc'
      ENDIF

#ifdef DARWIN_ALLOW_CARBON
      IF (.NOT.DARWIN_haveWind .AND. .NOT.DARWIN_useEXFwind) THEN
        WRITE(msgBuf,'(2A,F7.3,A)')
     &    '** WARNING ** DARWIN_EXF_READPARMS: ',
     &    'default constant wind speed ',windconst,' m/s is used'
        CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
      ENDIF

      IF (.NOT.DARWIN_havepCO2) THEN
        WRITE(msgBuf,'(2A,E12.6,A)')
     &    '** WARNING ** DARWIN_EXF_READPARMS: ',
     &    'default pCO2 of ',pCO2const,' atm is used'
        CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
      ENDIF
#endif

#ifdef ALLOW_RADTRANS
C     darwin iceFile is used for carbon air-sea exchanges while
C     RT_iceFile is used for light.
C     If one is set but not the other we stop.
C     If they are both set but different, we warn (see further below).
      IF (DARWIN_haveIce .NEQV. RT_haveIce) THEN
       WRITE(msgBuf,'(A,A)') 'DARWIN_EXF_READPARMS: ',
     &  'One of darwin and radtrans does not have ice forcing set.'
       CALL PRINT_ERROR( msgBuf, myThid )
       WRITE(msgBuf,'(A,A)') 'DARWIN_EXF_READPARMS: ',
     &  'If you really want this, set iceconst or RT_iceconst to zero.'
       CALL PRINT_ERROR( msgBuf, myThid )
       STOP 'ABNORMAL END: DARWIN_EXF_READPARMS: ice inconsistency'
      ENDIF
#endif

#ifdef ALLOW_RADTRANS
C     warn if darwin and radtrans ice forcing are different
      IF (iceFile.NE.RT_iceFile .OR. iceconst.NE.RT_iceconst) THEN
       WRITE(msgBuf,'(A,A)') '** WARNING ** DARWIN_EXF_READPARMS: ',
     &  'iceFile or iceconst different from RT_iceFile or RT_iceconst.'
       CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
       WRITE(msgBuf,'(A,A)') '** WARNING ** DARWIN_EXF_READPARMS: ',
     &  'darwin settings will be used for carbon surface forcing and'
       CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
       WRITE(msgBuf,'(A,A)') '** WARNING ** DARWIN_EXF_READPARMS: ',
     &  'radtrans settings for light.'
       CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      ENDIF
#endif

      errCount = 0
      IF ( useExfYearlyFields ) THEN
       IF ( PARRepCycle.NE.0. ) THEN
        WRITE(msgBuf,'(2A)') 'DARWIN_EXF_READPARMS: The use of ',
     &       'useExfYearlyFields AND PARRepCycle is not implemented'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
       ENDIF
       IF ( ironRepCycle.NE.0. ) THEN
        WRITE(msgBuf,'(2A)') 'DARWIN_EXF_READPARMS: The use of ',
     &       'useExfYearlyFields AND ironRepCycle is not implemented'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
       ENDIF
       IF ( iceRepCycle.NE.0. ) THEN
        WRITE(msgBuf,'(2A)') 'DARWIN_EXF_READPARMS: The use of ',
     &       'useExfYearlyFields AND iceRepCycle is not implemented'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
       ENDIF
       IF ( windRepCycle.NE.0. ) THEN
        WRITE(msgBuf,'(2A)') 'DARWIN_EXF_READPARMS: The use of ',
     &       'useExfYearlyFields AND windRepCycle is not implemented'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
       ENDIF
       IF ( pCO2RepCycle.NE.0. ) THEN
        WRITE(msgBuf,'(2A)') 'DARWIN_EXF_READPARMS: The use of ',
     &       'useExfYearlyFields AND pCO2RepCycle is not implemented'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
       ENDIF
       IF ( DOCrunoffRepCycle.NE.0. ) THEN
        WRITE(msgBuf,'(2A)') 'DARWIN_EXF_READPARMS: The use of ',
     &   'useExfYearlyFields AND',
     &   'DOCrunoffRepCycle is not implemented'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF
       IF ( tDOCrunoffRepCycle.NE.0. ) THEN
        WRITE(msgBuf,'(2A)') 'DARWIN_EXF_READPARMS: The use of ',
     &   'useExfYearlyFields AND',
     &   'tDOCrunoffRepCycle is not implemented'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
       ENDIF
       IF ( DONrunoffRepCycle.NE.0. ) THEN
        WRITE(msgBuf,'(2A)') 'DARWIN_EXF_READPARMS: The use of ',
     &   'useExfYearlyFields AND',
     &   'DONrunoffRepCycle is not implemented'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
       ENDIF
       IF ( DOPrunoffRepCycle.NE.0. ) THEN
        WRITE(msgBuf,'(2A)') 'DARWIN_EXF_READPARMS: The use of ',
     &   'useExfYearlyFields AND',
     &   'DOPrunoffRepCycle is not implemented'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
       ENDIF
       IF ( DSirunoffRepCycle.NE.0. ) THEN
        WRITE(msgBuf,'(2A)') 'DARWIN_EXF_READPARMS: The use of ',
     &   'useExfYearlyFields AND',
     &   'DSirunoffRepCycle is not implemented'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
       ENDIF
       IF ( DICrunoffRepCycle.NE.0. ) THEN
        WRITE(msgBuf,'(2A)') 'DARWIN_EXF_READPARMS: The use of ',
     &   'useExfYearlyFields AND',
     &   'DICrunoffRepCycle is not implemented'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
       ENDIF
       IF ( ALKrunoffRepCycle.NE.0. ) THEN
        WRITE(msgBuf,'(2A)') 'DARWIN_EXF_READPARMS: The use of ',
     &   'useExfYearlyFields AND',
     &   'ALKrunoffRepCycle is not implemented'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
       ENDIF
      ENDIF
      IF ( errCount.GE.1 ) THEN
        WRITE(msgBuf,'(A,I3,A)')
     &     'DARWIN_EXF_READPARMS: detected', errCount,' fatal error(s)'
        CALL PRINT_ERROR( msgBuf, myThid )
        CALL ALL_PROC_DIE( 0 )
        STOP 'ABNORMAL END: S/R DARWIN_EXF_READPARMS'
      ENDIF

      IF (oUnit .GE. 0) THEN
       WRITE(UNIT=oUnit, NML=DARWIN_FORCING_PARAMS)
#ifdef USE_EXF_INTERPOLATION
       WRITE(UNIT=oUnit, NML=DARWIN_INTERP_PARAMS)
#endif
      ENDIF

#endif /* ALLOW_DARWIN */

      RETURN
      END
