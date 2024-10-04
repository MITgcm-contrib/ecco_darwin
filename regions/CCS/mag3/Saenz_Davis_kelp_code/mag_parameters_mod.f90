
! ======================================================================
! Global Paramaters, Types, & Variables
! ======================================================================
module mag_parameters_mod

  use mag_kinds_mod, only : r8,i4

  implicit none
  public
  save

  integer, parameter, public ::                 &
       breakage_Duarte_Ferreira = 0             , &   ! fancy breakage
       breakage_Rodrigues       = 1             , &   ! death rate scales with swh
       nx             = 2160                    , &   ! grid is fixed
       ny             = 4320                    , &    ! grid is fixed
       chunk_size     = 40                           !

  real(kind=r8), parameter, public ::        &
      c0     =    0.0_r8                   , &
      c1     =    1.0_r8                   , &
      c2     =    2.0_r8                   , &
      c3     =    3.0_r8                   , &
      c4     =    4.0_r8                   , &
      c10    =   10.0_r8                   , &
      c1000  = 1000.0_r8                   , &
      p001   =    0.001_r8                 , &
      p5     =    0.5_r8                   , &
      pi     =    3.14159265358979323846_r8


  !---------------------------------------------------------------------
  !  Unit Conversion
  !---------------------------------------------------------------------

  real(kind=r8), parameter, public :: &
	  mw_n      = 14.00672_r8,  & ! molecular weight N
      sphr      = 3600.0_r8,    & ! number of seconds in an hour
      spd       = 86400.0_r8,   & ! number of seconds in a day
      dpy       = 365.0_r8,     & ! number of days in a year
      spy       = dpy*spd,      & ! number of seconds in a year
      hrps      = c1 / sphr,    & ! number of hours in a second
      dps       = c1 / spd,     & ! number of days in a second
      ypd       = c1 / dpy,     & ! number of years in a day
      yps       = c1 / spy,     & ! number of years in a second
      cmperm    = 100.0_r8,     & ! cm per meter
      mpercm    = .01_r8          ! meters per cm

  !---------------------------------------------------------------------
  !  Physical Constants
  !---------------------------------------------------------------------

  real(kind=r8), parameter, public :: &
      vonkar    =   0.4_r8,            & ! von Karman constant
      T0_Kelvin = 273.15_r8,           & ! freezing T of fresh water (K)
      K_Boltz   =   8.617330350e-5_r8, & ! Boltzmann constant (eV/K)
      rho_sw    =   1.026_r8,          & ! density of salt water (g/cm^3)
      epsC      =   1.0e-8_r8,         & ! small C concentration (mmol C/m^3)
      epsTinv   =   3.17e-8_r8,        & ! small inverse time scale (1/year) (1/sec)
      molw_Fe   =  55.845_r8,          & ! molecular weight of iron (gFe / mol Fe)
      R13C_std  =   1.0_r8,            & ! actual 13C/12C PDB standard ratio (Craig, 1957) = 1123.72e-5_r8
      R14C_std =    1.0_r8               ! actual 14C/12C NOSAMS standard ratio = 11.76e-13_r8

  !real, parameter :: d_a(7) = (/ -1.36471e-1, 4.68181e-2, 8.07004e-1, -7.45353e-3, -2.94418e-3, 3.43570e-5, 3.48658e-5/)

  !---------------------------------------------------------------------
  !  Parameters Enumeration
  !---------------------------------------------------------------------

  integer(kind=i4), parameter, public :: &
      npar = 37          , &
	  mp_spp_Vmax = 1   , &
	  mp_spp_Ks_NO3 = 2   , &
	  mp_spp_kcap = 3   , &
	  mp_spp_Gmax_cap = 4   , &
	  mp_spp_PARs = 5   , &
	  mp_spp_PARc = 6   , &
	  mp_spp_Q0 = 7   , &
	  mp_spp_Qmin = 8   , &
	  mp_spp_Qmax = 9   , &
	  mp_spp_BtoSA = 10   , &
	  mp_spp_line_sep = 11   , &
	  mp_spp_kcap_rate = 12   , &
	  mp_spp_Topt1 = 13   , &
	  mp_spp_K1 = 14   , &
	  mp_spp_Topt2 = 15   , &
	  mp_spp_K2 = 16   , &
	  mp_spp_CD = 17   , &
	  mp_spp_dry_sa = 18   , &
	  mp_spp_dry_wet = 19   , &
	  mp_spp_E = 20   , &
	  mp_spp_seed = 21   , &
	  mp_spp_death = 22   , &
	  mp_harvest_type = 23   , &
      mp_harvest_schedule = 24, &
      mp_harvest_avg_period = 25, &
   	  mp_harvest_kg = 26     , &
	  mp_harvest_f = 27     , &
	  mp_harvest_nmax = 28     , &
	  mp_breakage_type = 29     , &
	  mp_dte = 30            , & ! matlab datenum od time step
	  mp_dt_mag = 31            , &! time step length (days)
      mp_N_flux_limit = 32      , &
      mp_farm_depth = 33        , &
      mp_wave_mort_factor = 34      , &
      mp_N_uptake_type = 35      , &
      mp_growth_lim_type = 36      , &
      mp_Q_lim_type = 37
  CONTAINS

!	real(kind=d8) PURE FUNCTION sind(deg)
!		sind = sin(deg*3.141592653589793_d8/180._d8)

!	real(kind=d8) PURE FUNCTION asind(n)
!		asind = asin(n)/3.141592653589793_d8*180._d8

!	real(kind=d8) PURE FUNCTION cosd(deg)
!		cosd = cos(deg*3.141592653589793_d8/180._d8)

!	real(kind=d8) PURE FUNCTION acosd(n)
!		acosd = acos(n)/3.141592653589793_d8*180._d8

    real(kind=r8) PURE FUNCTION daylength_r8(lat,lon,alt,dte)
!    real(kind=r8) FUNCTION daylength(lat,lon,alt,dte)

      use mag_kinds_mod, only : r8
      implicit none

      real(kind=r8), INTENT(IN) :: lat,lon,alt,dte
      real(kind=r8) :: dte_2000_days, n2000, Js, M, C, lambda, delta, h, omega ! Jt

      ! main function that computes daylength and noon time
      ! https://en.wikipedia.org/wiki/Sunrise_equation

      ! number of days since Jan 1st, 2000 12:00 UT
      !dte_2000_days = 730490.0_r8
      !n2000 = dte - dte_2000_days + 68.184_r8/86400._r8
      n2000 = dte - 0.5 + 68.184_r8/86400._r8

      ! mean solar moon
      Js = n2000 - lon/360._r8

      ! solar mean anomaly
      M = mod(357.5291_r8 + 0.98560028_r8*Js,360._r8)

      ! center
      C = 1.9148_r8*sind(M) + 0.0200_r8*sind(2.0_r8*M) + 0.0003_r8*sind(3.0_r8*M)

      ! ecliptic longitude
      lambda = mod(M + C + 180._r8 + 102.9372_r8,360._r8)

      ! solar transit -- don't need for only daylength
      !Jt = 2451545.5_r8 + Js + 0.0053_r8*sind(M) - 0.0069_r8*sind(2.0_r8*lambda)

      ! Sun declination
      delta = asind(sind(lambda)*sind(23.44_r8))

      ! hour angle (day expressed in geometric degrees)
      h = (sind(-0.83_r8 - 2.076_r8*sqrt(alt)/60._r8) - sind(lat)*sind(delta))/(cosd(lat)*cosd(delta))


      ! to avoid meaningless complex angles: forces omega to 0 or 12h
      if (h < -1) then
        omega = 180._r8
      elseif (h > 1) then
        omega = 0._r8
      else
        omega = acosd(h)
      endif

      !print *,dte_2000_days,n2000,Js,M,C,lambda,delta,h,omega,lat,lon,dte

      daylength_r8 = omega/180._r8

    end FUNCTION daylength_r8

    real(kind=d8) PURE FUNCTION daylength(lat,lon,alt,dte)
!    real(kind=r8) FUNCTION daylength(lat,lon,alt,dte)

      use mag_kinds_mod, only : r8,d8
      implicit none

      real(kind=r8), INTENT(IN) :: lat,lon,alt,dte
      real(kind=d8) :: dte_2000_days, n2000, Js, M, C, lambda, delta, h, omega ! Jt

      ! main function that computes daylength and noon time
      ! https://en.wikipedia.org/wiki/Sunrise_equation

      ! number of days since Jan 1st, 2000 12:00 UT
      !dte_2000_days = 730490.0_d8
      !n2000 = dte - dte_2000_days + 68.184_d8/86400._d8
      n2000 = dte - 0.5 + 68.184_d8/86400._d8

      ! mean solar moon
      Js = n2000 - lon/360._d8

      ! solar mean anomaly
      M = mod(357.5291_d8 + 0.98560028_d8*Js,360._d8)

      ! center
      C = 1.9148_d8*sind(M) + 0.0200_d8*sind(2.0_d8*M) + 0.0003_d8*sind(3.0_d8*M)

      ! ecliptic longitude
      lambda = mod(M + C + 180._d8 + 102.9372_d8,360._d8)

      ! solar transit -- don't need for only daylength
      !Jt = 2451545.5_d8 + Js + 0.0053_d8*sind(M) - 0.0069_d8*sind(2.0_d8*lambda)

      ! Sun declination
      delta = asind(sind(lambda)*sind(23.44_d8))

      ! hour angle (day expressed in geometric degrees)
      h = (sind(-0.83_d8 - 2.076_d8*sqrt(alt)/60._d8) - sind(lat)*sind(delta))/(cosd(lat)*cosd(delta))


      ! to avoid meaningless complex angles: forces omega to 0 or 12h
      if (h < -1) then
        omega = 180._d8
      elseif (h > 1) then
        omega = 0._d8
      else
        omega = acosd(h)
      endif

      !print *,dte_2000_days,n2000,Js,M,C,lambda,delta,h,omega,lat,lon,dte

      daylength = omega/180._d8

    end FUNCTION daylength


    real(kind=r8) PURE FUNCTION lambda_NO3(magu,Tw,CD,VmaxNO3,KsNO3,NO3)

      use mag_kinds_mod, only : r8, i4
      implicit none

      real(kind=r8), INTENT(IN) :: magu,Tw,CD,VmaxNO3,KsNO3,NO3

      integer(kind=i4), parameter :: n_length = 25
      real(kind=r8), parameter :: pi = 3.14159265358979323846_r8
      real(kind=r8), parameter :: visc = 1.e-6_r8 * 86400._r8
      real(kind=r8), parameter :: Dm = (18_r8*3.65e-11_r8 + 9.72e-10_r8) * 86400._r8
      real(kind=r8), dimension(n_length) :: vval  ! two v's b/c val seems to be a keyword in fortran

      integer(kind=i4) :: n
      real(kind=r8) :: DBL,oscillatory,flow,beta, magu_m, Tw_s,NO3_u

      ! unit conversions, minima
      magu_m = max(1.0_r8,magu*86400.0_r8) ! m/day
      Tw_s = max(0.01_r8,Tw)/86400.0_r8  ! whaa...?
      NO3_u = NO3*1000.0_r8 ! converting from uM to umol/m3


      DBL = 10._r8 * (visc / (sqrt(CD) * abs(magu_m)))

      ! 1. Oscillatory Flow

      do n = 1,n_length
        vval(n) = (1-exp((-Dm * n**2 * pi**2 *Tw_s)/(2._r8*DBL**2)))/(n**2 * pi**2)
      enddo
      oscillatory = ((4._r8*DBL)/Tw_s) * sum(vval)

      ! 2. Uni-directional Flow
      flow = Dm / DBL

      beta = flow + oscillatory

      lambda_NO3 = c1 + (VmaxNO3 / (beta*KsNO3)) - (NO3_u/KsNO3);

    end FUNCTION lambda_NO3


    real(kind=r8) PURE FUNCTION temp_lim(sst,Topt1,K1,Topt2,K2)

      use mag_kinds_mod, only : r8
      implicit none

      real(kind=r8), INTENT(IN) :: sst,Topt1,K1,Topt2,K2

      if (sst >= Topt1) then
        if (sst <= Topt2) then
          temp_lim = 1._r8
        else
          temp_lim = exp(-K2*(sst-Topt2)**2)
        endif
      else
        temp_lim = exp(-K1*(sst-Topt1)**2)
      endif
      temp_lim = max(0._r8,temp_lim)
      temp_lim = min(1._r8,temp_lim)

    end FUNCTION temp_lim


end module mag_parameters_mod

