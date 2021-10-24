#include "DARWIN_OPTIONS.h"
#include "EXF_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: DARWIN_EXF_INIT_VARIA

C     !INTERFACE:
      SUBROUTINE DARWIN_EXF_INIT_VARIA( myThid )

C     !DESCRIPTION:
C     Initialize DARWIN constants

C     !USES:
      IMPLICIT NONE
#include "EEPARAMS.h"
#include "SIZE.h"
#include "GRID.h"
#include "PARAMS.h"
#ifdef ALLOW_EXF
#include "EXF_INTERP_SIZE.h"
#endif
#include "DARWIN_SIZE.h"
#include "DARWIN_EXF_PARAMS.h"
#include "DARWIN_EXF_FIELDS.h"

C     !INPUT PARAMETERS:
      INTEGER myThid
CEOP

#ifdef ALLOW_DARWIN
#ifdef ALLOW_EXF

#ifdef ALLOW_RADTRANS
      INTEGER l
#endif

C     these load constant forcing fields

      call exf_init_fld (
     &     'PAR', PARfile, PARmask,
     &     PARperiod, darwin_inscal_PAR, PARconst,
     &     surfPAR, PAR0, PAR1,
#ifdef USE_EXF_INTERPOLATION
     &     PAR_lon0, PAR_lon_inc,
     &     PAR_lat0, PAR_lat_inc,
     &     PAR_nlon, PAR_nlat, xC, yC, PAR_interpMethod,
#endif
     &     mythid )

      call exf_init_fld (
     &     'iron', ironfile, ironmask,
     &     ironperiod, darwin_inscal_iron, ironconst,
     &     inputFe, iron0, iron1,
#ifdef USE_EXF_INTERPOLATION
     &     iron_lon0, iron_lon_inc,
     &     iron_lat0, iron_lat_inc,
     &     iron_nlon, iron_nlat, xC, yC, iron_interpMethod,
#endif
     &     mythid )

      call exf_init_fld (
     &     'ice', icefile, icemask,
     &     iceperiod, darwin_inscal_ice, iceconst,
     &     iceFrac, ice0, ice1,
#ifdef USE_EXF_INTERPOLATION
     &     ice_lon0, ice_lon_inc,
     &     ice_lat0, ice_lat_inc,
     &     ice_nlon, ice_nlat, xC, yC, ice_interpMethod,
#endif
     &     mythid )

      call exf_init_fld (
     &     'wind', windfile, windmask,
     &     windperiod, darwin_inscal_wind, windconst,
     &     windSpeed, wind0, wind1,
#ifdef USE_EXF_INTERPOLATION
     &     wind_lon0, wind_lon_inc,
     &     wind_lat0, wind_lat_inc,
     &     wind_nlon, wind_nlat, xC, yC, wind_interpMethod,
#endif
     &     mythid )

      call exf_init_fld (
     &     'pCO2', pCO2file, pCO2mask,
     &     pCO2period, darwin_inscal_pCO2, pCO2const,
     &     atmospCO2, pCO20, pCO21,
#ifdef USE_EXF_INTERPOLATION
     &     pCO2_lon0, pCO2_lon_inc,
     &     pCO2_lat0, pCO2_lat_inc,
     &     pCO2_nlon, pCO2_nlat, xC, yC, pCO2_interpMethod,
#endif
     &     mythid )

      call exf_init_fld (
     &     'DOCrunoff', DOCrunofffile, DOCrunoffmask,
     &     DOCrunoffperiod, darwin_inscal_DOCrunoff,
     &     DOCrunoffconst,
     &     DOCrunoff, DOCrunoff0, DOCrunoff1,
#ifdef USE_EXF_INTERPOLATION
     &     DOCrunoff_lon0, DOCrunoff_lon_inc,
     &     DOCrunoff_lat0, DOCrunoff_lat_inc,
     &     DOCrunoff_nlon, DOCrunoff_nlat, xC, yC,
     &     DOCrunoff_interpMethod,
#endif
     &     mythid )

      call exf_init_fld (
     &     'tDOCrunoff', tDOCrunofffile, tDOCrunoffmask,
     &     tDOCrunoffperiod, darwin_inscal_tDOCrunoff,
     &     tDOCrunoffconst,
     &     tDOCrunoff, tDOCrunoff0, tDOCrunoff1,
#ifdef USE_EXF_INTERPOLATION
     &     tDOCrunoff_lon0, tDOCrunoff_lon_inc,
     &     tDOCrunoff_lat0, tDOCrunoff_lat_inc,
     &     tDOCrunoff_nlon, tDOCrunoff_nlat, xC, yC,
     &     tDOCrunoff_interpMethod,
#endif
     &     mythid )

      call exf_init_fld (
     &     'DONrunoff', DONrunofffile, DONrunoffmask,
     &     DONrunoffperiod, darwin_inscal_DONrunoff,
     &     DONrunoffconst,
     &     DONrunoff, DONrunoff0, DONrunoff1,
#ifdef USE_EXF_INTERPOLATION
     &     DONrunoff_lon0, DONrunoff_lon_inc,
     &     DONrunoff_lat0, DONrunoff_lat_inc,
     &     DONrunoff_nlon, DONrunoff_nlat, xC, yC,
     &     DONrunoff_interpMethod,
#endif
     &     mythid )

      call exf_init_fld (
     &     'DOPrunoff', DOPrunofffile, DOPrunoffmask,
     &     DOPrunoffperiod, darwin_inscal_DOPrunoff,
     &     DOPrunoffconst,
     &     DOPrunoff, DOPrunoff0, DOPrunoff1,
#ifdef USE_EXF_INTERPOLATION
     &     DOPrunoff_lon0, DOPrunoff_lon_inc,
     &     DOPrunoff_lat0, DOPrunoff_lat_inc,
     &     DOPrunoff_nlon, DOPrunoff_nlat, xC, yC,
     &     DOPrunoff_interpMethod,
#endif
     &     mythid )

      call exf_init_fld (
     &     'DSirunoff', DSirunofffile, DSirunoffmask,
     &     DSirunoffperiod, darwin_inscal_DSirunoff,
     &     DSirunoffconst,
     &     DSirunoff, DSirunoff0, DSirunoff1,
#ifdef USE_EXF_INTERPOLATION
     &     DSirunoff_lon0, DSirunoff_lon_inc,
     &     DSirunoff_lat0, DSirunoff_lat_inc,
     &     DSirunoff_nlon, DSirunoff_nlat, xC, yC,
     &     DSirunoff_interpMethod,
#endif
     &     mythid )

      call exf_init_fld (
     &     'DICrunoff', DICrunofffile, DICrunoffmask,
     &     DICrunoffperiod, darwin_inscal_DICrunoff,
     &     DICrunoffconst,
     &     DICrunoff, DICrunoff0, DICrunoff1,
#ifdef USE_EXF_INTERPOLATION
     &     DICrunoff_lon0, DICrunoff_lon_inc,
     &     DICrunoff_lat0, DICrunoff_lat_inc,
     &     DICrunoff_nlon, DICrunoff_nlat, xC, yC,
     &     DICrunoff_interpMethod,
#endif
     &     mythid )

      call exf_init_fld (
     &     'ALKrunoff', ALKrunofffile, ALKrunoffmask,
     &     ALKrunoffperiod, darwin_inscal_ALKrunoff,
     &     ALKrunoffconst,
     &     ALKrunoff, ALKrunoff0, ALKrunoff1,
#ifdef USE_EXF_INTERPOLATION
     &     ALKrunoff_lon0, ALKrunoff_lon_inc,
     &     ALKrunoff_lat0, ALKrunoff_lat_inc,
     &     ALKrunoff_nlon, ALKrunoff_nlat, xC, yC,
     &     ALKrunoff_interpMethod,
#endif
     &     mythid )

      IF (darwin_loadFieldsEarly) THEN
C       these load time-varying forcing fields
C       (usually done in darwin_fields_load, but may be needed for chl)
        call exf_set_fld(
     &     'PAR', PARfile, PARmask,
     &     PARStartTime, PARperiod, PARRepCycle,
     &     darwin_inscal_PAR,
     &     PAR_exfremo_intercept, PAR_exfremo_slope,
     &     surfPAR, PAR0, PAR1,
#ifdef USE_EXF_INTERPOLATION
     &     PAR_lon0, PAR_lon_inc,
     &     PAR_lat0, PAR_lat_inc,
     &     PAR_nlon, PAR_nlat, xC, yC, PAR_interpMethod,
#endif
     &     starttime, nIter0, mythid )

        call exf_set_fld(
     &     'iron', ironfile, ironmask,
     &     ironStartTime, ironperiod, ironRepCycle,
     &     darwin_inscal_iron,
     &     iron_exfremo_intercept, iron_exfremo_slope,
     &     inputFe, iron0, iron1,
#ifdef USE_EXF_INTERPOLATION
     &     iron_lon0, iron_lon_inc,
     &     iron_lat0, iron_lat_inc,
     &     iron_nlon, iron_nlat, xC, yC, iron_interpMethod,
#endif
     &     starttime, nIter0, mythid )

        call exf_set_fld(
     &     'ice', icefile, icemask,
     &     iceStartTime, iceperiod, iceRepCycle,
     &     darwin_inscal_ice,
     &     ice_exfremo_intercept, ice_exfremo_slope,
     &     iceFrac, ice0, ice1,
#ifdef USE_EXF_INTERPOLATION
     &     ice_lon0, ice_lon_inc,
     &     ice_lat0, ice_lat_inc,
     &     ice_nlon, ice_nlat, xC, yC, ice_interpMethod,
#endif
     &     starttime, nIter0, mythid )

        call exf_set_fld(
     &     'wind', windfile, windmask,
     &     windStartTime, windperiod, windRepCycle,
     &     darwin_inscal_wind,
     &     wind_exfremo_intercept, wind_exfremo_slope,
     &     windSpeed, wind0, wind1,
#ifdef USE_EXF_INTERPOLATION
     &     wind_lon0, wind_lon_inc,
     &     wind_lat0, wind_lat_inc,
     &     wind_nlon, wind_nlat, xC, yC, wind_interpMethod,
#endif
     &     starttime, nIter0, mythid )

        call exf_set_fld(
     &     'pCO2', pCO2file, pCO2mask,
     &     pCO2StartTime, pCO2period, pCO2RepCycle,
     &     darwin_inscal_pCO2,
     &     pCO2_exfremo_intercept, pCO2_exfremo_slope,
     &     atmospCO2, pCO20, pCO21,
#ifdef USE_EXF_INTERPOLATION
     &     pCO2_lon0, pCO2_lon_inc,
     &     pCO2_lat0, pCO2_lat_inc,
     &     pCO2_nlon, pCO2_nlat, xC, yC, pCO2_interpMethod,
#endif
     &     starttime, nIter0, mythid )

        call exf_set_fld(
     &     'DOCrunoff', DOCrunofffile, DOCrunoffmask,
     &     DOCrunoffStartTime, DOCrunoffperiod, DOCrunoffRepCycle,
     &     darwin_inscal_DOCrunoff,
     &     DOCrunoff_exfremo_intercept, DOCrunoff_exfremo_slope,
     &     DOCrunoff, DOCrunoff0, DOCrunoff1,
#ifdef USE_EXF_INTERPOLATION
     &     DOCrunoff_lon0, DOCrunoff_lon_inc,
     &     DOCrunoff_lat0, DOCrunoff_lat_inc,
     &     DOCrunoff_nlon, DOCrunoff_nlat, xC, yC,
     &     DOCrunoff_interpMethod,
#endif
     &     starttime, nIter0, mythid )

        call exf_set_fld(
     &     'tDOCrunoff', tDOCrunofffile, tDOCrunoffmask,
     &     tDOCrunoffStartTime, tDOCrunoffperiod, tDOCrunoffRepCycle,
     &     darwin_inscal_tDOCrunoff,
     &     tDOCrunoff_exfremo_intercept, tDOCrunoff_exfremo_slope,
     &     tDOCrunoff, tDOCrunoff0, tDOCrunoff1,
#ifdef USE_EXF_INTERPOLATION
     &     tDOCrunoff_lon0, tDOCrunoff_lon_inc,
     &     tDOCrunoff_lat0, tDOCrunoff_lat_inc,
     &     tDOCrunoff_nlon, tDOCrunoff_nlat, xC, yC,
     &     tDOCrunoff_interpMethod,
#endif
     &     starttime, nIter0, mythid )

        call exf_set_fld(
     &     'DONrunoff', DONrunofffile, DONrunoffmask,
     &     DONrunoffStartTime, DONrunoffperiod, DONrunoffRepCycle,
     &     darwin_inscal_DONrunoff,
     &     DONrunoff_exfremo_intercept, DONrunoff_exfremo_slope,
     &     DONrunoff, DONrunoff0, DONrunoff1,
#ifdef USE_EXF_INTERPOLATION
     &     DONrunoff_lon0, DONrunoff_lon_inc,
     &     DONrunoff_lat0, DONrunoff_lat_inc,
     &     DONrunoff_nlon, DONrunoff_nlat, xC, yC,
     &     DONrunoff_interpMethod,
#endif
     &     starttime, nIter0, mythid )

        call exf_set_fld(
     &     'DOPrunoff', DOPrunofffile, DOPrunoffmask,
     &     DOPrunoffStartTime, DOPrunoffperiod, DOPrunoffRepCycle,
     &     darwin_inscal_DOPrunoff,
     &     DOPrunoff_exfremo_intercept, DOPrunoff_exfremo_slope,
     &     DOPrunoff, DOPrunoff0, DOPrunoff1,
#ifdef USE_EXF_INTERPOLATION
     &     DOPrunoff_lon0, DOPrunoff_lon_inc,
     &     DOPrunoff_lat0, DOPrunoff_lat_inc,
     &     DOPrunoff_nlon, DOPrunoff_nlat, xC, yC,
     &     DOPrunoff_interpMethod,
#endif
     &     starttime, nIter0, mythid )

        call exf_set_fld(
     &     'DSirunoff', DSirunofffile, DSirunoffmask,
     &     DSirunoffStartTime, DSirunoffperiod, DSirunoffRepCycle,
     &     darwin_inscal_DSirunoff,
     &     DSirunoff_exfremo_intercept, DSirunoff_exfremo_slope,
     &     DSirunoff, DSirunoff0, DSirunoff1,
#ifdef USE_EXF_INTERPOLATION
     &     DSirunoff_lon0, DSirunoff_lon_inc,
     &     DSirunoff_lat0, DSirunoff_lat_inc,
     &     DSirunoff_nlon, DSirunoff_nlat, xC, yC,
     &     DSirunoff_interpMethod,
#endif
     &     starttime, nIter0, mythid )

        call exf_set_fld(
     &     'DICrunoff', DICrunofffile, DICrunoffmask,
     &     DICrunoffStartTime, DICrunoffperiod, DICrunoffRepCycle,
     &     darwin_inscal_DICrunoff,
     &     DICrunoff_exfremo_intercept, DICrunoff_exfremo_slope,
     &     DICrunoff, DICrunoff0, DICrunoff1,
#ifdef USE_EXF_INTERPOLATION
     &     DICrunoff_lon0, DICrunoff_lon_inc,
     &     DICrunoff_lat0, DICrunoff_lat_inc,
     &     DICrunoff_nlon, DICrunoff_nlat, xC, yC,
     &     DICrunoff_interpMethod,
#endif
     &     starttime, nIter0, mythid )

        call exf_set_fld(
     &     'ALKrunoff', ALKrunofffile, ALKrunoffmask,
     &     ALKrunoffStartTime, ALKrunoffperiod, ALKrunoffRepCycle,
     &     darwin_inscal_ALKrunoff,
     &     ALKrunoff_exfremo_intercept, ALKrunoff_exfremo_slope,
     &     ALKrunoff, ALKrunoff0, ALKrunoff1,
#ifdef USE_EXF_INTERPOLATION
     &     ALKrunoff_lon0, ALKrunoff_lon_inc,
     &     ALKrunoff_lat0, ALKrunoff_lat_inc,
     &     ALKrunoff_nlon, ALKrunoff_nlat, xC, yC,
     &     ALKrunoff_interpMethod,
#endif
     &     starttime, nIter0, mythid )

        _EXCH_XY_RL( surfPAR, myThid )

        _EXCH_XY_RL( inputFe, myThid )

        _EXCH_XY_RL( iceFrac, myThid )

        _EXCH_XY_RL( windSpeed, myThid )

        _EXCH_XY_RL( atmospCO2, myThid )

        _EXCH_XY_RL( DOCrunoff, myThid )

        _EXCH_XY_RL( tDOCrunoff, myThid )

        _EXCH_XY_RL( DONrunoff, myThid )

        _EXCH_XY_RL( DOPrunoff, myThid )

        _EXCH_XY_RL( DSirunoff, myThid )

        _EXCH_XY_RL( DICrunoff, myThid )

        _EXCH_XY_RL( ALKrunoff, myThid )

      ENDIF

#endif /* ALLOW_EXF */
#endif /* ALLOW_DARWIN */

      RETURN
      END
