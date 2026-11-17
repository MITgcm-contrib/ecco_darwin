#include "RADTRANS_OPTIONS.h"
#include "EXF_OPTIONS.h"

CBOP
C !ROUTINE: RADTRANS_READPARMS

C !INTERFACE: ==========================================================
      SUBROUTINE RADTRANS_READPARMS( myThid )

C !DESCRIPTION:
C     Initialize radtrans parameters, read in data.radtrans

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
#include "RADTRANS_SIZE.h"
#include "RADTRANS_PARAMS.h"
#include "RADTRANS_EXF_PARAMS.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      INTEGER myThid

C !OUTPUT PARAMETERS: ==================================================
C  none
CEOP

#ifdef ALLOW_RADTRANS

C !LOCAL VARIABLES: ====================================================
C  iUnit   :: unit number for I/O
C  msgBuf  :: message buffer
      integer iUnit, oUnit, l, errCount, m
      _RL gmW(7), gmK(7), wgt, wsum
      CHARACTER*(MAX_LEN_MBUF) msgBuf

      NAMELIST/RADTRANS_FORCING_PARAMS/
C
     &    RT_Edfile,
     &    RT_Esfile,
     &    RT_E_mask,
     &    RT_E_period,
     &    RT_E_RepCycle,
     &    RT_E_startTime,
     &    RT_E_startdate1,
     &    RT_E_startdate2,
     &    RT_Ed_const,
     &    RT_Ed_exfremo_intercept,
     &    RT_Ed_exfremo_slope,
     &    RT_inscal_Ed,
     &    RT_Es_const,
     &    RT_Es_exfremo_intercept,
     &    RT_Es_exfremo_slope,
     &    RT_inscal_Es,
C
     &    RT_icefile,
     &    RT_iceperiod,
     &    RT_iceRepCycle,
     &    RT_iceStartTime,
     &    RT_icestartdate1,
     &    RT_icestartdate2,
     &    RT_iceconst,
     &    RT_ice_exfremo_intercept,
     &    RT_ice_exfremo_slope,
     &    RT_icemask,
     &    RT_inscal_ice

#ifdef USE_EXF_INTERPOLATION
      NAMELIST/RADTRANS_INTERP_PARAMS/
C
     &    RT_E_lon0,
     &    RT_E_lat0,
     &    RT_E_nlon,
     &    RT_E_nlat,
     &    RT_E_lon_inc,
     &    RT_E_interpMethod,
     &    RT_E_lat_inc,
C
     &    RT_ice_lon0,
     &    RT_ice_lat0,
     &    RT_ice_nlon,
     &    RT_ice_nlat,
     &    RT_ice_lon_inc,
     &    RT_ice_interpMethod,
     &    RT_ice_lat_inc
#endif

      NAMELIST/RADTRANS_PARAMS/
     &    RT_refract_water,
     &    RT_rmud_max,
     &    RT_wbEdges,
     &    RT_wbRefWLs,
     &    RT_kmax,
     &    RT_useOASIMrmud,
     &    RT_useMeanCosSolz,
     &    RT_useNoonSolz,
     &    RT_sfcIrrThresh,
     &    RT_oasimWgt,
     &    RT_iceTransmit, RT_iceNcat, RT_iceMultDimSnow,
     &    RT_iceI0, RT_iceKappa, RT_snowKappaDry, RT_snowKappaWet,
     &    RT_snowI0Wet, RT_snowSSL, RT_snowMin,
     &    RT_pondFrac, RT_pondAlb, RT_pondI0, RT_iceSpectral

      NAMELIST/RADTRANS_DEPENDENT/
     &    RT_wbWidths, RT_kapShape

C     Grenfell & Maykut (1977) spectral extinction shape (see below)
      DATA gmW / 400. _d 0, 500. _d 0, 550. _d 0, 600. _d 0,
     &           650. _d 0, 700. _d 0, 750. _d 0 /
      DATA gmK / 1.1 _d 0, 1.1 _d 0, 1.2 _d 0, 1.5 _d 0,
     &           2.0 _d 0, 3.3 _d 0, 5.5 _d 0 /

      IF ( .NOT.useRADTRANS ) THEN
C-    pkg RADTRANS is not used
        _BEGIN_MASTER(myThid)
C-    Track pkg activation status:
C     print a (weak) warning if data.radtrans is found
         CALL PACKAGES_UNUSED_MSG( 'useRADTRANS', ' ', ' ' )
        _END_MASTER(myThid)
        RETURN
      ENDIF

C Open and read the data.radtrans file

      _BEGIN_MASTER(myThid)

      DO l = 1, nlam
       RT_Edfile(l) = ' '
       RT_Esfile(l) = ' '
       RT_inscal_Ed(l)            = 1.0 _d 0
       RT_inscal_Es(l)            = 1.0 _d 0
       RT_Ed_const(l)             = 0.0 _d 0
       RT_Es_const(l)             = 0.0 _d 0
       RT_Ed_exfremo_intercept(l) = 0.0 _d 0
       RT_Es_exfremo_intercept(l) = 0.0 _d 0
       RT_Ed_exfremo_slope(l)     = 0.0 _d 0
       RT_Es_exfremo_slope(l)     = 0.0 _d 0
       RT_oasimWgt(l)             = 1.0 _d 0
      ENDDO
      RT_E_period     = 0.0 _d 0
      RT_E_RepCycle   = repeatPeriod
      RT_E_StartTime  = UNSET_RL
      RT_E_startdate1 = 0
      RT_E_startdate2 = 0
      RT_E_mask = 'c'


      RT_icefile = ' '
      RT_iceperiod             = 0.0 _d 0
      RT_iceRepCycle           = repeatPeriod
      RT_iceStartTime          = UNSET_RL
      RT_icestartdate1         = 0
      RT_icestartdate2         = 0
      RT_iceconst              = UNSET_RL
      RT_ice_exfremo_intercept = 0.0 _d 0
      RT_ice_exfremo_slope     = 0.0 _d 0
      RT_icemask = 'c'

      RT_inscal_ice    =  1.0 _d 0

#ifdef USE_EXF_INTERPOLATION
      RT_E_lon0 = inp_lon0
      RT_E_lat0 = inp_lat0
      RT_E_nlon = inp_gNx
      RT_E_nlat = inp_gNy
      RT_E_lon_inc = inp_dLon
      RT_E_interpMethod  = 1

      RT_ice_lon0 = inp_lon0
      RT_ice_lat0 = inp_lat0
      RT_ice_nlon = inp_gNx
      RT_ice_nlat = inp_gNy
      RT_ice_lon_inc = inp_dLon
      RT_ice_interpMethod  = 1

      DO l=1,MAX_LAT_INC
        RT_E_lat_inc(l) = inp_dLat(l)
        RT_ice_lat_inc(l) = inp_dLat(l)
      ENDDO
#endif /* USE_EXF_INTERPOLATION */

C     these are taken from Ackleson, et al. 1994 (JGR)
      RT_rd = 1.5 _d 0
      RT_ru = 3.0 _d 0
      RT_rmus = 1.0/0.83 _d 0
      RT_rmuu = 1.0/0.4 _d 0

      RT_refract_water = 1.341 _d 0        ! refractive index of seawater
      RT_rmud_max = 1.5 _d 0
      DO l = 1,nlam
        RT_wbEdges(l) = -1 _d 0
        RT_wbRefWLs(l) = -1 _d 0
      ENDDO
      RT_kmax = Nr
      RT_useOASIMrmud = useOASIM
      RT_useMeanCosSolz = .FALSE.
      RT_useNoonSolz = .FALSE.
      RT_sfcIrrThresh = 1 _d -4
C     under-ice light (defaults = recommended values; README.md)
      RT_iceTransmit    = .FALSE.
C     1 = mean thickness, consistent with SEAICE_multDim=1 thermodynamics
C     and within x2 of point benchmarks; 7 (Castellani et al. 2017) is an
C     upper-bound sensitivity: the flat [0,2h] pdf puts 1/7 of the ice at
C     ~h/7 with ~hs/7 snow and raises snow-covered transmittance ~4x
      RT_iceNcat        = 1
      RT_iceMultDimSnow = .TRUE.
C     Briegleb et al. (2004) CCSM3 visible i0
      RT_iceI0          = 0.70 _d 0
C     Lebrun et al. (2023) tuned PAR extinction; Light et al. (2008) 0.6-0.8
      RT_iceKappa       = 1.0 _d 0
C     Lebrun et al. (2023) snow: dry i0=1, kappa 7; wet 3-cm SSL i0=0.45, kappa 5
      RT_snowKappaDry   = 7.0 _d 0
      RT_snowKappaWet   = 5.0 _d 0
      RT_snowI0Wet      = 0.45 _d 0
      RT_snowSSL        = 0.03 _d 0
      RT_snowMin        = 0.02 _d 0
C     ponds (Lebrun et al. 2023 ISP: albedo 0.27, i0 0.75); fraction from
C     Nicolaus et al. (2012), FYI 0.42 / MYI 0.23 seasonal maxima
      RT_pondFrac       = 0.2 _d 0
      RT_pondAlb        = 0.27 _d 0
      RT_pondI0         = 0.75 _d 0
      RT_iceSpectral    = .TRUE.

      WRITE(msgBuf,'(A)') ' RADTRANS_READPARMS: opening data.radtrans'
      CALL PRINT_MESSAGE(msgBuf, standardMessageUnit,
     &                   SQUEEZE_RIGHT , 1)

      CALL OPEN_COPY_DATA_FILE(
     I                   'data.radtrans', 'RADTRANS_READPARMS',
     O                   iUnit,
     I                   myThid )

      READ(UNIT=iUnit, NML=radtrans_forcing_params)
#ifdef USE_EXF_INTERPOLATION
      READ(UNIT=iUnit, NML=radtrans_interp_params)
#endif
      READ(UNIT=iUnit, NML=radtrans_params)
#ifdef SINGLE_DISK_IO
      CLOSE(iUnit)
#else
      CLOSE(iUnit,STATUS='DELETE')
#endif /* SINGLE_DISK_IO */

C ======================================================================

      RT_haveIce = RT_iceFile.NE.' ' .OR. RT_iceconst.NE.UNSET_RL

#ifdef ALLOW_SEAICE
      RT_useSEAICE = .not. RT_haveIce
      IF (RT_useSEAICE) THEN
        WRITE(msgBuf,'(2A)') '** WARNING ** RT_READPARMS: ',
     &  'RT_iceFile and RT_iceconst unset, using area from seaice pkg'
        CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
      ELSE
        WRITE(msgBuf,'(2A)') '** WARNING ** RT_READPARMS: ',
     &  'RT_iceFile or RT_iceconst set, not using area from seaice pkg'
        CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
      ENDIF
#endif

      IF (RT_iceconst.EQ.UNSET_RL) THEN
        RT_iceconst = 0 _d 0
      ENDIF

C ======================================================================
C band widths used to convert total-waveband irradiance to per nm
      RT_wbTotalWidth = 0.0
      DO l=1,nlam
        RT_wbWidths(l) = RT_wbEdges(l+1) -
     &                  RT_wbEdges(l)
        RT_wbTotalWidth = RT_wbTotalWidth + RT_wbWidths(l)
      ENDDO
      IF (RT_wbTotalWidth.LE.0) then
        WRITE(msgBuf,'(2A)') 'RADTRANS_READPARMS: ',
     &  'please provide wavebabnds in RT_wbEdges.'
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R RADTRANS_READPARMS'
      ENDIF
      DO l = 1,nlam
        IF (RT_wbRefWLs(l) .LT. 0.0) THEN
          RT_wbRefWLs(l) = 0.5*(RT_wbEdges(l) +
     &                          RT_wbEdges(l+1))
        ENDIF
      ENDDO

C spectral shape of ice/snow extinction: Grenfell & Maykut (1977, Fig. 4)
C white interior ice, read from the figure (600/700 nm interpolated);
C interpolated to RT_wbRefWLs and normalised to a band-width-weighted
C mean of 1 over the radtrans wavebands, so RT_iceKappa etc stay PAR means
      wsum = 0.0
      wgt  = 0.0
      DO l = 1,nlam
       RT_kapShape(l) = 1.0 _d 0
       IF (RT_iceSpectral) THEN
        IF (RT_wbRefWLs(l).LE.gmW(1)) THEN
         RT_kapShape(l) = gmK(1)
        ELSEIF (RT_wbRefWLs(l).GE.gmW(7)) THEN
         RT_kapShape(l) = gmK(7)
        ELSE
         DO m = 1,6
          IF (RT_wbRefWLs(l).GE.gmW(m) .AND.
     &        RT_wbRefWLs(l).LE.gmW(m+1)) THEN
           RT_kapShape(l) = gmK(m) + (gmK(m+1)-gmK(m))
     &        *(RT_wbRefWLs(l)-gmW(m))/(gmW(m+1)-gmW(m))
          ENDIF
         ENDDO
        ENDIF
       ENDIF
       wsum = wsum + RT_wbWidths(l)*RT_kapShape(l)
       wgt  = wgt  + RT_wbWidths(l)
      ENDDO
      DO l = 1,nlam
       RT_kapShape(l) = RT_kapShape(l)*wgt/wsum
      ENDDO
      IF (RT_iceNcat.LT.1) THEN
        WRITE(msgBuf,'(2A)') 'RADTRANS_READPARMS: ',
     &  'RT_iceNcat must be >= 1'
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R RADTRANS_READPARMS'
      ENDIF

C ----------------------------------------------------------------------
      IF ( myProcId.EQ.0 .AND. myThid.EQ.1 ) THEN
        CALL MDSFINDUNIT( oUnit, mythid )
        OPEN(oUnit,file='radtrans_params.txt',status='unknown')
        WRITE(UNIT=oUnit, NML=RADTRANS_FORCING_PARAMS)
#ifdef USE_EXF_INTERPOLATION
        WRITE(UNIT=oUnit, NML=RADTRANS_INTERP_PARAMS)
#endif
        WRITE(UNIT=oUnit, NML=RADTRANS_PARAMS)
        WRITE(UNIT=oUnit, NML=RADTRANS_DEPENDENT)
        CLOSE(oUnit)
      ENDIF

      WRITE(msgBuf,'(A)') ' ==================================='
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )

      _END_MASTER(myThid)

C Everyone else must wait for the parameters to be loaded
      _BARRIER

#endif /* ALLOW_RADTRANS */

      RETURN
      END

