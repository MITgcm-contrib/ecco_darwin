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
     &    PAR_lon0, PAR_lat0, PAR_lon_inc,
     &    PAR_lat_inc
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
     &    iron_lon0, iron_lat0, iron_lon_inc,
     &    iron_lat_inc
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
     &    ice_lon0, ice_lat0, ice_lon_inc,
     &    ice_lat_inc
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
     &    wind_lon0, wind_lat0, wind_lon_inc,
     &    wind_lat_inc
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
     &    pCO2_lon0, pCO2_lat0, pCO2_lon_inc,
     &    pCO2_lat_inc
      INTEGER pCO2_interpMethod, pCO2_nlon, pCO2_nlat
      _RL  pCO2_lon0
      _RL  pCO2_lat0
      _RL  pCO2_lon_inc
      _RL  pCO2_lat_inc(MAX_LAT_INC)
#endif

C ventHe3 forcing parameters for exf

      _RL ventHe3StartTime

      COMMON/darwin_forcing_ventHe3_c/
     &    ventHe3mask
      COMMON/darwin_forcing_ventHe3_i/
     &    ventHe3startdate1, ventHe3startdate2
      COMMON/darwin_forcing_ventHe3_r/
     &    ventHe3StartTime,
     &    ventHe3period, ventHe3RepCycle, ventHe3const,
     &    ventHe3_exfremo_intercept, ventHe3_exfremo_slope,
     &    darwin_inscal_ventHe3
      CHARACTER*1 ventHe3mask
      INTEGER ventHe3startdate1
      INTEGER ventHe3startdate2
      _RL ventHe3period
      _RL ventHe3RepCycle
      _RL ventHe3const
      _RL ventHe3_exfremo_intercept
      _RL ventHe3_exfremo_slope
      _RL darwin_inscal_ventHe3

#ifdef USE_EXF_INTERPOLATION
      COMMON/darwin_interp_ventHe3_i/
     &    ventHe3_nlon, ventHe3_nlat, ventHe3_interpMethod
      COMMON/darwin_interp_ventHe3_r/
     &    ventHe3_lon0, ventHe3_lat0, ventHe3_lon_inc,
     &    ventHe3_lat_inc
      INTEGER ventHe3_interpMethod, ventHe3_nlon, ventHe3_nlat
      _RL  ventHe3_lon0
      _RL  ventHe3_lat0
      _RL  ventHe3_lon_inc
      _RL  ventHe3_lat_inc(MAX_LAT_INC)
#endif

C DOCrunoff forcing parameters for exf

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
     &    DOCrunoff_lat_inc
      INTEGER DOCrunoff_interpMethod, DOCrunoff_nlon, DOCrunoff_nlat
      _RL  DOCrunoff_lon0
      _RL  DOCrunoff_lat0
      _RL  DOCrunoff_lon_inc
      _RL  DOCrunoff_lat_inc(MAX_LAT_INC)
#endif

C rDOC forcing parameters for exf

      _RL rDOCrunoffStartTime

      COMMON/darwin_forcing_rDOCrunoff_c/
     &    rDOCrunoffmask
      COMMON/darwin_forcing_rDOCrunoff_i/
     &    rDOCrunoffstartdate1, rDOCrunoffstartdate2
      COMMON/darwin_forcing_rDOCrunoff_r/
     &    rDOCrunoffStartTime,
     &    rDOCrunoffperiod, rDOCrunoffRepCycle, rDOCrunoffconst,
     &    rDOCrunoff_exfremo_intercept, rDOCrunoff_exfremo_slope,
     &    darwin_inscal_rDOCrunoff
      CHARACTER*1 rDOCrunoffmask
      INTEGER rDOCrunoffstartdate1
      INTEGER rDOCrunoffstartdate2
      _RL rDOCrunoffperiod
      _RL rDOCrunoffRepCycle
      _RL rDOCrunoffconst
      _RL rDOCrunoff_exfremo_intercept
      _RL rDOCrunoff_exfremo_slope
      _RL darwin_inscal_rDOCrunoff

#ifdef USE_EXF_INTERPOLATION
      COMMON/darwin_interp_rDOCrunoff_i/
     &    rDOCrunoff_nlon, rDOCrunoff_nlat, rDOCrunoff_interpMethod
      COMMON/darwin_interp_rDOCrunoff_r/
     &    rDOCrunoff_lon0, rDOCrunoff_lat0, rDOCrunoff_lon_inc,
     %    rDOCrunoff_lat_inc
      INTEGER rDOCrunoff_interpMethod, rDOCrunoff_nlon, rDOCrunoff_nlat
      _RL  rDOCrunoff_lon0
      _RL  rDOCrunoff_lat0
      _RL  rDOCrunoff_lon_inc
      _RL  rDOCrunoff_lat_inc(MAX_LAT_INC)
#endif

C CDOM forcing parameters for exf

      _RL CDOMrunoffStartTime

      COMMON/darwin_forcing_CDOMrunoff_c/
     &    CDOMrunoffmask
      COMMON/darwin_forcing_CDOMrunoff_i/
     &    CDOMrunoffstartdate1, CDOMrunoffstartdate2
      COMMON/darwin_forcing_CDOMrunoff_r/
     &    CDOMrunoffStartTime,
     &    CDOMrunoffperiod, CDOMrunoffRepCycle, CDOMrunoffconst,
     &    CDOMrunoff_exfremo_intercept, CDOMrunoff_exfremo_slope,
     &    darwin_inscal_CDOMrunoff
      CHARACTER*1 CDOMrunoffmask
      INTEGER CDOMrunoffstartdate1
      INTEGER CDOMrunoffstartdate2
      _RL CDOMrunoffperiod
      _RL CDOMrunoffRepCycle
      _RL CDOMrunoffconst
      _RL CDOMrunoff_exfremo_intercept
      _RL CDOMrunoff_exfremo_slope
      _RL darwin_inscal_CDOMrunoff

#ifdef USE_EXF_INTERPOLATION
      COMMON/darwin_interp_CDOMrunoff_i/
     &    CDOMrunoff_nlon, CDOMrunoff_nlat, CDOMrunoff_interpMethod
      COMMON/darwin_interp_CDOMrunoff_r/
     &    CDOMrunoff_lon0, CDOMrunoff_lat0, CDOMrunoff_lon_inc,
     %    CDOMrunoff_lat_inc
      INTEGER CDOMrunoff_interpMethod, CDOMrunoff_nlon, CDOMrunoff_nlat
      _RL  CDOMrunoff_lon0
      _RL  CDOMrunoff_lat0
      _RL  CDOMrunoff_lon_inc
      _RL  CDOMrunoff_lat_inc(MAX_LAT_INC)
#endif

C DONrunoff forcing parameters for exf

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

C DOPrunoff forcing parameters for exf

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

C NO3runoff forcing parameters for exf

      _RL NO3runoffStartTime

      COMMON/darwin_forcing_NO3runoff_c/
     &    NO3runoffmask
      COMMON/darwin_forcing_NO3runoff_i/
     &    NO3runoffstartdate1, NO3runoffstartdate2
      COMMON/darwin_forcing_NO3runoff_r/
     &    NO3runoffStartTime,
     &    NO3runoffperiod, NO3runoffRepCycle, NO3runoffconst,
     &    NO3runoff_exfremo_intercept, NO3runoff_exfremo_slope,
     &    darwin_inscal_NO3runoff
      CHARACTER*1 NO3runoffmask
      INTEGER NO3runoffstartdate1
      INTEGER NO3runoffstartdate2
      _RL NO3runoffperiod
      _RL NO3runoffRepCycle
      _RL NO3runoffconst
      _RL NO3runoff_exfremo_intercept
      _RL NO3runoff_exfremo_slope
      _RL darwin_inscal_NO3runoff

#ifdef USE_EXF_INTERPOLATION
      COMMON/darwin_interp_NO3runoff_i/
     &    NO3runoff_nlon, NO3runoff_nlat, NO3runoff_interpMethod
      COMMON/darwin_interp_NO3runoff_r/
     &    NO3runoff_lon0, NO3runoff_lat0, NO3runoff_lon_inc,
     &    NO3runoff_lat_inc
      INTEGER NO3runoff_interpMethod, NO3runoff_nlon, NO3runoff_nlat
      _RL  NO3runoff_lon0
      _RL  NO3runoff_lat0
      _RL  NO3runoff_lon_inc
      _RL  NO3runoff_lat_inc(MAX_LAT_INC)
#endif

C NO2runoff forcing parameters for exf

      _RL NO2runoffStartTime

      COMMON/darwin_forcing_NO2runoff_c/
     &    NO2runoffmask
      COMMON/darwin_forcing_NO2runoff_i/
     &    NO2runoffstartdate1, NO2runoffstartdate2
      COMMON/darwin_forcing_NO2runoff_r/
     &    NO2runoffStartTime,
     &    NO2runoffperiod, NO2runoffRepCycle, NO2runoffconst,
     &    NO2runoff_exfremo_intercept, NO2runoff_exfremo_slope,
     &    darwin_inscal_NO2runoff
      CHARACTER*1 NO2runoffmask
      INTEGER NO2runoffstartdate1
      INTEGER NO2runoffstartdate2
      _RL NO2runoffperiod
      _RL NO2runoffRepCycle
      _RL NO2runoffconst
      _RL NO2runoff_exfremo_intercept
      _RL NO2runoff_exfremo_slope
      _RL darwin_inscal_NO2runoff

#ifdef USE_EXF_INTERPOLATION
      COMMON/darwin_interp_NO2runoff_i/
     &    NO2runoff_nlon, NO2runoff_nlat, NO2runoff_interpMethod
      COMMON/darwin_interp_NO2runoff_r/
     &    NO2runoff_lon0, NO2runoff_lat0, NO2runoff_lon_inc,
     &    NO2runoff_lat_inc
      INTEGER NO2runoff_interpMethod, NO2runoff_nlon, NO2runoff_nlat
      _RL  NO2runoff_lon0
      _RL  NO2runoff_lat0
      _RL  NO2runoff_lon_inc
      _RL  NO2runoff_lat_inc(MAX_LAT_INC)
#endif

C NH4runoff forcing parameters for exf

      _RL NH4runoffStartTime

      COMMON/darwin_forcing_NH4runoff_c/
     &    NH4runoffmask
      COMMON/darwin_forcing_NH4runoff_i/
     &    NH4runoffstartdate1, NH4runoffstartdate2
      COMMON/darwin_forcing_NH4runoff_r/
     &    NH4runoffStartTime,
     &    NH4runoffperiod, NH4runoffRepCycle, NH4runoffconst,
     &    NH4runoff_exfremo_intercept, NH4runoff_exfremo_slope,
     &    darwin_inscal_NH4runoff
      CHARACTER*1 NH4runoffmask
      INTEGER NH4runoffstartdate1
      INTEGER NH4runoffstartdate2
      _RL NH4runoffperiod
      _RL NH4runoffRepCycle
      _RL NH4runoffconst
      _RL NH4runoff_exfremo_intercept
      _RL NH4runoff_exfremo_slope
      _RL darwin_inscal_NH4runoff

#ifdef USE_EXF_INTERPOLATION
      COMMON/darwin_interp_NH4runoff_i/
     &    NH4runoff_nlon, NH4runoff_nlat, NH4runoff_interpMethod
      COMMON/darwin_interp_NH4runoff_r/
     &    NH4runoff_lon0, NH4runoff_lat0, NH4runoff_lon_inc,
     &    NH4runoff_lat_inc
      INTEGER NH4runoff_interpMethod, NH4runoff_nlon, NH4runoff_nlat
      _RL  NH4runoff_lon0
      _RL  NH4runoff_lat0
      _RL  NH4runoff_lon_inc
      _RL  NH4runoff_lat_inc(MAX_LAT_INC)
#endif

C PO4runoff forcing parameters for exf

      _RL PO4runoffStartTime

      COMMON/darwin_forcing_PO4runoff_c/
     &    PO4runoffmask
      COMMON/darwin_forcing_PO4runoff_i/
     &    PO4runoffstartdate1, PO4runoffstartdate2
      COMMON/darwin_forcing_PO4runoff_r/
     &    PO4runoffStartTime,
     &    PO4runoffperiod, PO4runoffRepCycle, PO4runoffconst,
     &    PO4runoff_exfremo_intercept, PO4runoff_exfremo_slope,
     &    darwin_inscal_PO4runoff
      CHARACTER*1 PO4runoffmask
      INTEGER PO4runoffstartdate1
      INTEGER PO4runoffstartdate2
      _RL PO4runoffperiod
      _RL PO4runoffRepCycle
      _RL PO4runoffconst
      _RL PO4runoff_exfremo_intercept
      _RL PO4runoff_exfremo_slope
      _RL darwin_inscal_PO4runoff

#ifdef USE_EXF_INTERPOLATION
      COMMON/darwin_interp_PO4runoff_i/
     &    PO4runoff_nlon, PO4runoff_nlat, PO4runoff_interpMethod
      COMMON/darwin_interp_PO4runoff_r/
     &    PO4runoff_lon0, PO4runoff_lat0, PO4runoff_lon_inc,
     &    PO4runoff_lat_inc
      INTEGER PO4runoff_interpMethod, PO4runoff_nlon, PO4runoff_nlat
      _RL  PO4runoff_lon0
      _RL  PO4runoff_lat0
      _RL  PO4runoff_lon_inc
      _RL  PO4runoff_lat_inc(MAX_LAT_INC)
#endif

C DSirunoff forcing parameters for exf

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

C POCrunoff forcing parameters for exf

      _RL POCrunoffStartTime

      COMMON/darwin_forcing_POCrunoff_c/
     &    POCrunoffmask
      COMMON/darwin_forcing_POCrunoff_i/
     &    POCrunoffstartdate1, POCrunoffstartdate2
      COMMON/darwin_forcing_POCrunoff_r/
     &    POCrunoffStartTime,
     &    POCrunoffperiod, POCrunoffRepCycle, POCrunoffconst,
     &    POCrunoff_exfremo_intercept, POCrunoff_exfremo_slope,
     &    darwin_inscal_POCrunoff
      CHARACTER*1 POCrunoffmask
      INTEGER POCrunoffstartdate1
      INTEGER POCrunoffstartdate2
      _RL POCrunoffperiod
      _RL POCrunoffRepCycle
      _RL POCrunoffconst
      _RL POCrunoff_exfremo_intercept
      _RL POCrunoff_exfremo_slope
      _RL darwin_inscal_POCrunoff

#ifdef USE_EXF_INTERPOLATION
      COMMON/darwin_interp_POCrunoff_i/
     &    POCrunoff_nlon, POCrunoff_nlat, POCrunoff_interpMethod
      COMMON/darwin_interp_POCrunoff_r/
     &    POCrunoff_lon0, POCrunoff_lat0, POCrunoff_lon_inc,
     &    POCrunoff_lat_inc
      INTEGER POCrunoff_interpMethod, POCrunoff_nlon, POCrunoff_nlat
      _RL  POCrunoff_lon0
      _RL  POCrunoff_lat0
      _RL  POCrunoff_lon_inc
      _RL  POCrunoff_lat_inc(MAX_LAT_INC)
#endif

C POPrunoff forcing parameters for exf

      _RL POPrunoffStartTime

      COMMON/darwin_forcing_POPrunoff_c/
     &    POPrunoffmask
      COMMON/darwin_forcing_POPrunoff_i/
     &    POPrunoffstartdate1, POPrunoffstartdate2
      COMMON/darwin_forcing_POPrunoff_r/
     &    POPrunoffStartTime,
     &    POPrunoffperiod, POPrunoffRepCycle, POPrunoffconst,
     &    POPrunoff_exfremo_intercept, POPrunoff_exfremo_slope,
     &    darwin_inscal_POPrunoff
      CHARACTER*1 POPrunoffmask
      INTEGER POPrunoffstartdate1
      INTEGER POPrunoffstartdate2
      _RL POPrunoffperiod
      _RL POPrunoffRepCycle
      _RL POPrunoffconst
      _RL POPrunoff_exfremo_intercept
      _RL POPrunoff_exfremo_slope
      _RL darwin_inscal_POPrunoff

#ifdef USE_EXF_INTERPOLATION
      COMMON/darwin_interp_POPrunoff_i/
     &    POPrunoff_nlon, POPrunoff_nlat, POPrunoff_interpMethod
      COMMON/darwin_interp_POPrunoff_r/
     &    POPrunoff_lon0, POPrunoff_lat0, POPrunoff_lon_inc,
     &    POPrunoff_lat_inc
      INTEGER POPrunoff_interpMethod, POPrunoff_nlon, POPrunoff_nlat
      _RL  POPrunoff_lon0
      _RL  POPrunoff_lat0
      _RL  POPrunoff_lon_inc
      _RL  POPrunoff_lat_inc(MAX_LAT_INC)
#endif

C PONrunoff forcing parameters for exf

      _RL PONrunoffStartTime

      COMMON/darwin_forcing_PONrunoff_c/
     &    PONrunoffmask
      COMMON/darwin_forcing_PONrunoff_i/
     &    PONrunoffstartdate1, PONrunoffstartdate2
      COMMON/darwin_forcing_PONrunoff_r/
     &    PONrunoffStartTime,
     &    PONrunoffperiod, PONrunoffRepCycle, PONrunoffconst,
     &    PONrunoff_exfremo_intercept, PONrunoff_exfremo_slope,
     &    darwin_inscal_PONrunoff
      CHARACTER*1 PONrunoffmask
      INTEGER PONrunoffstartdate1
      INTEGER PONrunoffstartdate2
      _RL PONrunoffperiod
      _RL PONrunoffRepCycle
      _RL PONrunoffconst
      _RL PONrunoff_exfremo_intercept
      _RL PONrunoff_exfremo_slope
      _RL darwin_inscal_PONrunoff

#ifdef USE_EXF_INTERPOLATION
      COMMON/darwin_interp_PONrunoff_i/
     &    PONrunoff_nlon, PONrunoff_nlat, PONrunoff_interpMethod
      COMMON/darwin_interp_PONrunoff_r/
     &    PONrunoff_lon0, PONrunoff_lat0, PONrunoff_lon_inc,
     &    PONrunoff_lat_inc
      INTEGER PONrunoff_interpMethod, PONrunoff_nlon, PONrunoff_nlat
      _RL  PONrunoff_lon0
      _RL  PONrunoff_lat0
      _RL  PONrunoff_lon_inc
      _RL  PONrunoff_lat_inc(MAX_LAT_INC)
#endif

C DICrunoff forcing parameters for exf

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

C ALKrunoff forcing parameters for exf

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

C--   COMMON/DARWIN_RUNOFF_PARAMS_r/ Parameters for runoff nutrient input
C     R_ALK_DIC_runoff  :: [eq/mol C]     ALK:DIC ratio in runoff
C     R_NO3_DIN_runoff  :: [mol N/mol N]  NO3 fraction of DIN in runoff
C     R_NO2_DIN_runoff  :: [mol N/mol N]  NO2 fraction of DIN in runoff
C     R_NH4_DIN_runoff  :: [mol N/mol N]  NH4 fraction of DIN in runoff
C     R_DIP_IP_runoff   :: [mol P/mol P]  ratio of dissolved to total inorganic P in runoff
C     R_DFe_DIP_runoff  :: [mol Fe/mol P] Fe:P ratio for inorganic P in runoff
C     R_DOFe_DOP_runoff :: [mol Fe/mol P] Fe:P ratio for dissolved organic matter in runoff
C     R_POFe_POP_runoff :: [mol Fe/mol P] Fe:P ratio for particulate organic matter in runoff
      COMMON/DARWIN_RUNOFF_PARAMS_R/
     &    R_ALK_DIC_runoff,
     &    R_NO3_DIN_runoff, R_NO2_DIN_runoff, R_NH4_DIN_runoff,
     &    R_DIP_IP_runoff, R_DFe_DIP_runoff,
     &    R_DOFe_DOP_runoff,
     &    R_POFe_POP_runoff
      _RL R_ALK_DIC_runoff
      _RL R_NO3_DIN_runoff
      _RL R_NO2_DIN_runoff
      _RL R_NH4_DIN_runoff
      _RL R_DIP_IP_runoff
      _RL R_DFe_DIP_runoff
      _RL R_DOFe_DOP_runoff
      _RL R_POFe_POP_runoff
CEOP

#endif /* ALLOW_DARWIN */
