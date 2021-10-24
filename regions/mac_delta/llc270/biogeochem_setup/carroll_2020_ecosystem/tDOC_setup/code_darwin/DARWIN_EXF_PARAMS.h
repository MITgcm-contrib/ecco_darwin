#ifdef ALLOW_DARWIN

CBOP
C     !ROUTINE: DARWIN_EXF_PARAMS.h
C     !INTERFACE:
C #include DARWIN_EXF_PARAMS.h

C     !DESCRIPTION:
C Contains parameters for reading forcing for darwin package through exf
C
C Requires: EXF_OPTIONS.h
C Requires: EXF_INTERP_SIZE.h

      COMMON/darwin_forcing_exf_params_l/
     &    darwin_loadFieldsEarly
      LOGICAL darwin_loadFieldsEarly

C PAR forcing parameters for exf

      _RL PARStartTime

      COMMON/darwin_forcing_PAR_c/
     &    PARmask
      COMMON/darwin_forcing_PAR_i/
     &    PARstartdate1, PARstartdate2
      COMMON/darwin_forcing_PAR_r/
     &    PARStartTime,
     &    PARperiod, PARRepCycle, PARconst,
     &    PAR_exfremo_intercept, PAR_exfremo_slope,
     &    darwin_inscal_PAR
      CHARACTER*1 PARmask
      INTEGER PARstartdate1
      INTEGER PARstartdate2
      _RL PARperiod
      _RL PARRepCycle
      _RL PARconst
      _RL PAR_exfremo_intercept
      _RL PAR_exfremo_slope
      _RL darwin_inscal_PAR

#ifdef USE_EXF_INTERPOLATION
      COMMON/darwin_interp_PAR_i/
     &    PAR_nlon, PAR_nlat, PAR_interpMethod
      COMMON/darwin_interp_PAR_r/
     &    PAR_lon0, PAR_lat0, PAR_lon_inc, PAR_lat_inc
      INTEGER PAR_interpMethod, PAR_nlon, PAR_nlat
      _RL  PAR_lon0
      _RL  PAR_lat0
      _RL  PAR_lon_inc
      _RL  PAR_lat_inc(MAX_LAT_INC)
#endif

C iron forcing parameters for exf

      _RL ironStartTime

      COMMON/darwin_forcing_iron_c/
     &    ironmask
      COMMON/darwin_forcing_iron_i/
     &    ironstartdate1, ironstartdate2
      COMMON/darwin_forcing_iron_r/
     &    ironStartTime,
     &    ironperiod, ironRepCycle, ironconst,
     &    iron_exfremo_intercept, iron_exfremo_slope,
     &    darwin_inscal_iron
      CHARACTER*1 ironmask
      INTEGER ironstartdate1
      INTEGER ironstartdate2
      _RL ironperiod
      _RL ironRepCycle
      _RL ironconst
      _RL iron_exfremo_intercept
      _RL iron_exfremo_slope
      _RL darwin_inscal_iron

#ifdef USE_EXF_INTERPOLATION
      COMMON/darwin_interp_iron_i/
     &    iron_nlon, iron_nlat, iron_interpMethod
      COMMON/darwin_interp_iron_r/
     &    iron_lon0, iron_lat0, iron_lon_inc, iron_lat_inc
      INTEGER iron_interpMethod, iron_nlon, iron_nlat
      _RL  iron_lon0
      _RL  iron_lat0
      _RL  iron_lon_inc
      _RL  iron_lat_inc(MAX_LAT_INC)
#endif

C ice forcing parameters for exf

      _RL iceStartTime

      COMMON/darwin_forcing_ice_c/
     &    icemask
      COMMON/darwin_forcing_ice_i/
     &    icestartdate1, icestartdate2
      COMMON/darwin_forcing_ice_r/
     &    iceStartTime,
     &    iceperiod, iceRepCycle, iceconst,
     &    ice_exfremo_intercept, ice_exfremo_slope,
     &    darwin_inscal_ice
      CHARACTER*1 icemask
      INTEGER icestartdate1
      INTEGER icestartdate2
      _RL iceperiod
      _RL iceRepCycle
      _RL iceconst
      _RL ice_exfremo_intercept
      _RL ice_exfremo_slope
      _RL darwin_inscal_ice

#ifdef USE_EXF_INTERPOLATION
      COMMON/darwin_interp_ice_i/
     &    ice_nlon, ice_nlat, ice_interpMethod
      COMMON/darwin_interp_ice_r/
     &    ice_lon0, ice_lat0, ice_lon_inc, ice_lat_inc
      INTEGER ice_interpMethod, ice_nlon, ice_nlat
      _RL  ice_lon0
      _RL  ice_lat0
      _RL  ice_lon_inc
      _RL  ice_lat_inc(MAX_LAT_INC)
#endif

C wind forcing parameters for exf

      _RL windStartTime

      COMMON/darwin_forcing_wind_c/
     &    windmask
      COMMON/darwin_forcing_wind_i/
     &    windstartdate1, windstartdate2
      COMMON/darwin_forcing_wind_r/
     &    windStartTime,
     &    windperiod, windRepCycle, windconst,
     &    wind_exfremo_intercept, wind_exfremo_slope,
     &    darwin_inscal_wind
      CHARACTER*1 windmask
      INTEGER windstartdate1
      INTEGER windstartdate2
      _RL windperiod
      _RL windRepCycle
      _RL windconst
      _RL wind_exfremo_intercept
      _RL wind_exfremo_slope
      _RL darwin_inscal_wind

#ifdef USE_EXF_INTERPOLATION
      COMMON/darwin_interp_wind_i/
     &    wind_nlon, wind_nlat, wind_interpMethod
      COMMON/darwin_interp_wind_r/
     &    wind_lon0, wind_lat0, wind_lon_inc, wind_lat_inc
      INTEGER wind_interpMethod, wind_nlon, wind_nlat
      _RL  wind_lon0
      _RL  wind_lat0
      _RL  wind_lon_inc
      _RL  wind_lat_inc(MAX_LAT_INC)
#endif

C pCO2 forcing parameters for exf

      _RL pCO2StartTime

      COMMON/darwin_forcing_pCO2_c/
     &    pCO2mask
      COMMON/darwin_forcing_pCO2_i/
     &    pCO2startdate1, pCO2startdate2
      COMMON/darwin_forcing_pCO2_r/
     &    pCO2StartTime,
     &    pCO2period, pCO2RepCycle, pCO2const,
     &    pCO2_exfremo_intercept, pCO2_exfremo_slope,
     &    darwin_inscal_pCO2
      CHARACTER*1 pCO2mask
      INTEGER pCO2startdate1
      INTEGER pCO2startdate2
      _RL pCO2period
      _RL pCO2RepCycle
      _RL pCO2const
      _RL pCO2_exfremo_intercept
      _RL pCO2_exfremo_slope
      _RL darwin_inscal_pCO2

#ifdef USE_EXF_INTERPOLATION
      COMMON/darwin_interp_pCO2_i/
     &    pCO2_nlon, pCO2_nlat, pCO2_interpMethod
      COMMON/darwin_interp_pCO2_r/
     &    pCO2_lon0, pCO2_lat0, pCO2_lon_inc, pCO2_lat_inc
      INTEGER pCO2_interpMethod, pCO2_nlon, pCO2_nlat
      _RL  pCO2_lon0
      _RL  pCO2_lat0
      _RL  pCO2_lon_inc
      _RL  pCO2_lat_inc(MAX_LAT_INC)
#endif

C DOC forcing parameters for exf

      _RL DOCrunoffStartTime

      COMMON/darwin_forcing_DOCrunoff_c/
     &    DOCrunoffmask
      COMMON/darwin_forcing_DOCrunoff_i/
     &    DOCrunoffstartdate1, DOCrunoffstartdate2
      COMMON/darwin_forcing_DOCrunoff_r/
     &    DOCrunoffStartTime,
     &    DOCrunoffperiod, DOCrunoffRepCycle, DOCrunoffconst,
     &    DOCrunoff_exfremo_intercept, DOCrunoff_exfremo_slope,
     &    darwin_inscal_DOCrunoff
      CHARACTER*1 DOCrunoffmask
      INTEGER DOCrunoffstartdate1
      INTEGER DOCrunoffstartdate2
      _RL DOCrunoffperiod
      _RL DOCrunoffRepCycle
      _RL DOCrunoffconst
      _RL DOCrunoff_exfremo_intercept
      _RL DOCrunoff_exfremo_slope
      _RL darwin_inscal_DOCrunoff

#ifdef USE_EXF_INTERPOLATION
      COMMON/darwin_interp_DOCrunoff_i/
     &    DOCrunoff_nlon, DOCrunoff_nlat, DOCrunoff_interpMethod
      COMMON/darwin_interp_DOCrunoff_r/
     &    DOCrunoff_lon0, DOCrunoff_lat0, DOCrunoff_lon_inc,
     %    DOCrunoff_lat_inc
      INTEGER DOCrunoff_interpMethod, DOCrunoff_nlon, DOCrunoff_nlat
      _RL  DOCrunoff_lon0
      _RL  DOCrunoff_lat0
      _RL  DOCrunoff_lon_inc
      _RL  DOCrunoff_lat_inc(MAX_LAT_INC)
#endif

C tDOC forcing parameters for exf

      _RL tDOCrunoffStartTime

      COMMON/darwin_forcing_tDOCrunoff_c/
     &    tDOCrunoffmask
      COMMON/darwin_forcing_tDOCrunoff_i/
     &    tDOCrunoffstartdate1, tDOCrunoffstartdate2
      COMMON/darwin_forcing_tDOCrunoff_r/
     &    tDOCrunoffStartTime,
     &    tDOCrunoffperiod, tDOCrunoffRepCycle, tDOCrunoffconst,
     &    tDOCrunoff_exfremo_intercept, tDOCrunoff_exfremo_slope,
     &    darwin_inscal_tDOCrunoff
      CHARACTER*1 tDOCrunoffmask
      INTEGER tDOCrunoffstartdate1
      INTEGER tDOCrunoffstartdate2
      _RL tDOCrunoffperiod
      _RL tDOCrunoffRepCycle
      _RL tDOCrunoffconst
      _RL tDOCrunoff_exfremo_intercept
      _RL tDOCrunoff_exfremo_slope
      _RL darwin_inscal_tDOCrunoff

#ifdef USE_EXF_INTERPOLATION
      COMMON/darwin_interp_tDOCrunoff_i/
     &    tDOCrunoff_nlon, tDOCrunoff_nlat, tDOCrunoff_interpMethod
      COMMON/darwin_interp_tDOCrunoff_r/
     &    tDOCrunoff_lon0, tDOCrunoff_lat0, tDOCrunoff_lon_inc,
     %    tDOCrunoff_lat_inc
      INTEGER tDOCrunoff_interpMethod, tDOCrunoff_nlon, tDOCrunoff_nlat
      _RL  tDOCrunoff_lon0
      _RL  tDOCrunoff_lat0
      _RL  tDOCrunoff_lon_inc
      _RL  tDOCrunoff_lat_inc(MAX_LAT_INC)
#endif

C DON forcing parameters for exf

      _RL DONrunoffStartTime

      COMMON/darwin_forcing_DONrunoff_c/
     &    DONrunoffmask
      COMMON/darwin_forcing_DONrunoff_i/
     &    DONrunoffstartdate1, DONrunoffstartdate2
      COMMON/darwin_forcing_DONrunoff_r/
     &    DONrunoffStartTime,
     &    DONrunoffperiod, DONrunoffRepCycle, DONrunoffconst,
     &    DONrunoff_exfremo_intercept, DONrunoff_exfremo_slope,
     &    darwin_inscal_DONrunoff
      CHARACTER*1 DONrunoffmask
      INTEGER DONrunoffstartdate1
      INTEGER DONrunoffstartdate2
      _RL DONrunoffperiod
      _RL DONrunoffRepCycle
      _RL DONrunoffconst
      _RL DONrunoff_exfremo_intercept
      _RL DONrunoff_exfremo_slope
      _RL darwin_inscal_DONrunoff

#ifdef USE_EXF_INTERPOLATION
      COMMON/darwin_interp_DONrunoff_i/
     &    DONrunoff_nlon, DONrunoff_nlat, DONrunoff_interpMethod
      COMMON/darwin_interp_DONrunoff_r/
     &    DONrunoff_lon0, DONrunoff_lat0, DONrunoff_lon_inc,
     &    DONrunoff_lat_inc
      INTEGER DONrunoff_interpMethod, DONrunoff_nlon, DONrunoff_nlat
      _RL  DONrunoff_lon0
      _RL  DONrunoff_lat0
      _RL  DONrunoff_lon_inc
      _RL  DONrunoff_lat_inc(MAX_LAT_INC)
#endif

C DOP forcing parameters for exf

      _RL DOPrunoffStartTime

      COMMON/darwin_forcing_DOPrunoff_c/
     &    DOPrunoffmask
      COMMON/darwin_forcing_DOPrunoff_i/
     &    DOPrunoffstartdate1, DOPrunoffstartdate2
      COMMON/darwin_forcing_DOPrunoff_r/
     &    DOPrunoffStartTime,
     &    DOPrunoffperiod, DOPrunoffRepCycle, DOPrunoffconst,
     &    DOPrunoff_exfremo_intercept, DOPrunoff_exfremo_slope,
     &    darwin_inscal_DOPrunoff
      CHARACTER*1 DOPrunoffmask
      INTEGER DOPrunoffstartdate1
      INTEGER DOPrunoffstartdate2
      _RL DOPrunoffperiod
      _RL DOPrunoffRepCycle
      _RL DOPrunoffconst
      _RL DOPrunoff_exfremo_intercept
      _RL DOPrunoff_exfremo_slope
      _RL darwin_inscal_DOPrunoff

#ifdef USE_EXF_INTERPOLATION
      COMMON/darwin_interp_DOPrunoff_i/
     &    DOPrunoff_nlon, DOPrunoff_nlat, DOPrunoff_interpMethod
      COMMON/darwin_interp_DOPrunoff_r/
     &    DOPrunoff_lon0, DOPrunoff_lat0, DOPrunoff_lon_inc,
     &    DOPrunoff_lat_inc
      INTEGER DOPrunoff_interpMethod, DOPrunoff_nlon, DOPrunoff_nlat
      _RL  DOPrunoff_lon0
      _RL  DOPrunoff_lat0
      _RL  DOPrunoff_lon_inc
      _RL  DOPrunoff_lat_inc(MAX_LAT_INC)
#endif

C DSi forcing parameters for exf

      _RL DSirunoffStartTime

      COMMON/darwin_forcing_DSirunoff_c/
     &    DSirunoffmask
      COMMON/darwin_forcing_DSirunoff_i/
     &    DSirunoffstartdate1, DSirunoffstartdate2
      COMMON/darwin_forcing_DSirunoff_r/
     &    DSirunoffStartTime,
     &    DSirunoffperiod, DSirunoffRepCycle, DSirunoffconst,
     &    DSirunoff_exfremo_intercept, DSirunoff_exfremo_slope,
     &    darwin_inscal_DSirunoff
      CHARACTER*1 DSirunoffmask
      INTEGER DSirunoffstartdate1
      INTEGER DSirunoffstartdate2
      _RL DSirunoffperiod
      _RL DSirunoffRepCycle
      _RL DSirunoffconst
      _RL DSirunoff_exfremo_intercept
      _RL DSirunoff_exfremo_slope
      _RL darwin_inscal_DSirunoff

#ifdef USE_EXF_INTERPOLATION
      COMMON/darwin_interp_DSirunoff_i/
     &    DSirunoff_nlon, DSirunoff_nlat, DSirunoff_interpMethod
      COMMON/darwin_interp_DSirunoff_r/
     &    DSirunoff_lon0, DSirunoff_lat0, DSirunoff_lon_inc,
     &    DSirunoff_lat_inc
      INTEGER DSirunoff_interpMethod, DSirunoff_nlon, DSirunoff_nlat
      _RL  DSirunoff_lon0
      _RL  DSirunoff_lat0
      _RL  DSirunoff_lon_inc
      _RL  DSirunoff_lat_inc(MAX_LAT_INC)
#endif

C DIC forcing parameters for exf

      _RL DICrunoffStartTime

      COMMON/darwin_forcing_DICrunoff_c/
     &    DICrunoffmask
      COMMON/darwin_forcing_DICrunoff_i/
     &    DICrunoffstartdate1, DICrunoffstartdate2
      COMMON/darwin_forcing_DICrunoff_r/
     &    DICrunoffStartTime,
     &    DICrunoffperiod, DICrunoffRepCycle, DICrunoffconst,
     &    DICrunoff_exfremo_intercept, DICrunoff_exfremo_slope,
     &    darwin_inscal_DICrunoff
      CHARACTER*1 DICrunoffmask
      INTEGER DICrunoffstartdate1
      INTEGER DICrunoffstartdate2
      _RL DICrunoffperiod
      _RL DICrunoffRepCycle
      _RL DICrunoffconst
      _RL DICrunoff_exfremo_intercept
      _RL DICrunoff_exfremo_slope
      _RL darwin_inscal_DICrunoff

#ifdef USE_EXF_INTERPOLATION
      COMMON/darwin_interp_DICrunoff_i/
     &    DICrunoff_nlon, DICrunoff_nlat, DICrunoff_interpMethod
      COMMON/darwin_interp_DICrunoff_r/
     &    DICrunoff_lon0, DICrunoff_lat0, DICrunoff_lon_inc,
     &    DICrunoff_lat_inc
      INTEGER DICrunoff_interpMethod, DICrunoff_nlon, DICrunoff_nlat
      _RL  DICrunoff_lon0
      _RL  DICrunoff_lat0
      _RL  DICrunoff_lon_inc
      _RL  DICrunoff_lat_inc(MAX_LAT_INC)
#endif

C ALK forcing parameters for exf

      _RL ALKrunoffStartTime

      COMMON/darwin_forcing_ALKrunoff_c/
     &    ALKrunoffmask
      COMMON/darwin_forcing_ALKrunoff_i/
     &    ALKrunoffstartdate1, ALKrunoffstartdate2
      COMMON/darwin_forcing_ALKrunoff_r/
     &    ALKrunoffStartTime,
     &    ALKrunoffperiod, ALKrunoffRepCycle, ALKrunoffconst,
     &    ALKrunoff_exfremo_intercept, ALKrunoff_exfremo_slope,
     &    darwin_inscal_ALKrunoff
      CHARACTER*1 ALKrunoffmask
      INTEGER ALKrunoffstartdate1
      INTEGER ALKrunoffstartdate2
      _RL ALKrunoffperiod
      _RL ALKrunoffRepCycle
      _RL ALKrunoffconst
      _RL ALKrunoff_exfremo_intercept
      _RL ALKrunoff_exfremo_slope
      _RL darwin_inscal_ALKrunoff

#ifdef USE_EXF_INTERPOLATION
      COMMON/darwin_interp_ALKrunoff_i/
     &    ALKrunoff_nlon, ALKrunoff_nlat, ALKrunoff_interpMethod
      COMMON/darwin_interp_ALKrunoff_r/
     &    ALKrunoff_lon0, ALKrunoff_lat0, ALKrunoff_lon_inc,
     &    ALKrunoff_lat_inc
      INTEGER ALKrunoff_interpMethod, ALKrunoff_nlon, ALKrunoff_nlat
      _RL  ALKrunoff_lon0
      _RL  ALKrunoff_lat0
      _RL  ALKrunoff_lon_inc
      _RL  ALKrunoff_lat_inc(MAX_LAT_INC)
#endif
CEOP

#endif /* ALLOW_DARWIN */
