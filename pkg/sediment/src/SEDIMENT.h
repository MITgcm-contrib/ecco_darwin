CBOP
C !ROUTINE: SEDIMENT.h
C !INTERFACE:
C #include "SEDIMENT.h"

C !DESCRIPTION:
C *==================================================================*
C | Header for pkg "sediment": suspended sediment classes carried as
C | passive tracers (pkg/ptracers, concentration in kg/m^3), with
C | settling, deposition to and erosion from a bed store (kg/m^2),
C | and a sediment load on the added mass (addMass, e.g. a river).
C | Optional ROMS-CSTMS-like bed (sedNlay > 0): bed layers with an
C | active layer, a bulk critical stress that consolidates / swells
C | toward an equilibrium profile, and mixed-bed erosion (Sherwood et
C | al. 2018, GMD 11; Moriarty et al. 2017, Biogeosciences 14), and
C | an optional wave contribution to the bed stress (files, e.g. from
C | WAVEWATCH III, or a parametric fetch- and depth-limited wind sea).
C | The bed elevation is not changed (no morphodynamic feedback).
C *==================================================================*
CEOP

C     sedMaxCl :: maximum number of sediment classes
      INTEGER sedMaxCl
      PARAMETER ( sedMaxCl = 4 )
C     sedMaxLay :: maximum number of bed layers
      INTEGER sedMaxLay
      PARAMETER ( sedMaxLay = 10 )

C--   Parameters (namelist SEDIMENT_PARM01 in data.sediment)
C     sedNcl      :: number of sediment classes in use (<= sedMaxCl)
C     sedPtrIdx   :: pTracer index holding each class (kg/m^3)
C     sedWs       :: settling velocity [m/s] (> 0 downward)
C     sedTauCE    :: critical bed shear stress for erosion [Pa]
C     sedEroM     :: erosion rate constant [kg/m^2/s]:
C                    E = sedEroM*(tau_b/sedTauCE - 1) when tau_b > sedTauCE
C     sedTauCD    :: critical stress for deposition [Pa] (Krone):
C                    deposited fraction of the settling flux out of the
C                    bottom layer is 1 - tau_b/sedTauCD ; <= 0: always 1
C     sedRiverC   :: concentration of the added fluid (addMass) [kg/m^3]
C     sedBulkRho  :: dry bulk density of the deposit [kg/m^3] (only for
C                    the equivalent-thickness diagnostic and monitor)
C     sedBedFile  :: initial bed mass file [kg/m^2] (' ': clean bed)
C     sedRiverCFile :: optional 2-D map of the concentration of the added
C                    fluid [kg/m^3] (e.g. one value per river at its
C                    source cell); replaces sedRiverC where given
C     sedDryDepth :: no erosion from columns shallower than this [m];
C                    deposition ramps from 1 at 2*sedDryDepth to 0 here
C     sedDryDeposit :: .TRUE.: deposit from drying columns (films) too
C     sedRefDepth :: the bed stress uses the mean velocity of the lowest
C                    sedRefDepth [m] of the column (not just the bottom
C                    cell, which can be a thin partial cell damped by
C                    drag); <= 0: bottom wet level only
C     sedMonFreq  :: frequency of the %SEDIMENT_MON monitor [s]
C     sedConsTime :: fluff layer: deposits enter an easily eroded fluff
C                    layer that consolidates into the bed with this
C                    e-folding time [s]; <= 0 (default): no fluff layer
C     sedTauCEf   :: critical stress for erosion of the fluff [Pa]
C     sedEroMf    :: erosion rate constant of the fluff [kg/m^2/s];
C                    the consolidated bed only erodes once the fluff of
C                    the column is gone
      INTEGER sedNcl
      INTEGER sedPtrIdx(sedMaxCl)
      _RL sedWs(sedMaxCl), sedTauCE(sedMaxCl)
      _RL sedEroM(sedMaxCl), sedTauCD(sedMaxCl)
      _RL sedRiverC(sedMaxCl), sedBulkRho(sedMaxCl)
      _RL sedDryDepth, sedRefDepth, sedMonFreq, sedConsTime
      _RL sedTauCEf(sedMaxCl), sedEroMf(sedMaxCl)
      LOGICAL sedDryDeposit
      COMMON /SEDIMENT_PARAMS_L/ sedDryDeposit
      CHARACTER*(MAX_LEN_FNAM) sedBedFile(sedMaxCl)
      COMMON /SEDIMENT_PARAMS_I/ sedNcl, sedPtrIdx
      COMMON /SEDIMENT_PARAMS_R/
     &     sedWs, sedTauCE, sedEroM, sedTauCD, sedRiverC, sedBulkRho,
     &     sedDryDepth, sedRefDepth, sedMonFreq, sedConsTime,
     &     sedTauCEf, sedEroMf
      CHARACTER*(MAX_LEN_FNAM) sedRiverCFile(sedMaxCl)
      COMMON /SEDIMENT_PARAMS_C/ sedBedFile, sedRiverCFile
C     sedRivC     :: concentration of the added fluid in each cell
C                    (sedRiverC, or the sedRiverCFile map) [kg/m^3]
      _RL sedRivC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy,sedMaxCl)
      COMMON /SEDIMENT_RIVC/ sedRivC

C--   Layered bed (namelist SEDIMENT_PARM02, optional)
C     sedNlay     :: number of bed layers; 0 (default): one bed store per
C                    class (with the optional fluff layer above). Layer 1
C                    is the active layer, the only one that erodes
C     sedCohesive :: .TRUE. for cohesive (mud) classes
C     sedPorosity :: bed porosity phi: E_i = sedEroM*(1-phi)*f_i
C                    *(tau_b/tau_ce,i - 1), f_i the mass fraction of
C                    class i in the active layer (Moriarty et al. Eq. 2)
C     sedFcLo, sedFcHi :: cohesive mass fraction of the active layer at
C                    which the bed starts / ends turning cohesive: the
C                    cohesive-behaviour factor Pc goes linearly 0 -> 1,
C                    tau_ce,i = max( Pc*tau_cb + (1-Pc)*sedTauCE(i),
C                    sedTauCE(i) ) (Sherwood et al. Eq. 6); sedTauCE is
C                    then the particle critical stress of each class
C     sedTcbNew   :: bulk critical stress of fresh deposit [Pa]
C     sedTcbOffset, sedTcbSlope :: equilibrium bulk critical stress
C                    tau_cb,eq = exp( (ln z_m - offset)/slope ) [Pa] at
C                    mass depth z_m [kg/m^2] (Sherwood et al. Eq. 4),
C                    limited to sedTcbMin .. sedTcbMax
C     sedTcons, sedTswell :: time scales [s] on which tau_cb relaxes
C                    toward tau_cb,eq from below (consolidation) and
C                    from above (swelling)
C     sedZaCoef, sedZaMin :: active layer thickness
C                    za = max( sedZaMin, sedZaCoef*(tau_b - tau_ce) )
C                    [m, m/Pa] (Harris and Wiberg 1997)
C     sedLayMass  :: mass per unit area [kg/m^2] above which material
C                    moved below the active layer starts a new layer
      INTEGER sedNlay
      LOGICAL sedCohesive(sedMaxCl)
      _RL sedPorosity, sedFcLo, sedFcHi, sedTcbNew
      _RL sedTcbOffset, sedTcbSlope, sedTcbMin, sedTcbMax
      _RL sedTcons, sedTswell, sedZaCoef, sedZaMin, sedLayMass
      COMMON /SEDIMENT_BED_I/ sedNlay
      COMMON /SEDIMENT_BED_L/ sedCohesive
      COMMON /SEDIMENT_BED_R/
     &     sedPorosity, sedFcLo, sedFcHi, sedTcbNew,
     &     sedTcbOffset, sedTcbSlope, sedTcbMin, sedTcbMax,
     &     sedTcons, sedTswell, sedZaCoef, sedZaMin, sedLayMass

C--   Waves (namelist SEDIMENT_PARM03, optional)
C     sedWaveMode :: 0 (default): no waves; 1: bottom orbital velocity
C                    [m/s], period [s] and (optional) direction [deg,
C                    toward, counter-clockwise from +x] from files, e.g.
C                    WAVEWATCH III UBR and a mean period; 2: parametric
C                    wind sea (Young and Verhagen 1996) from sedWindSpeed,
C                    sedWindDir and sedFetch, limited by the local depth;
C                    3: significant height [m], period [s] and direction
C                    from files (e.g. WAVEWATCH III HTSGW, PERPW, DIRPW
C                    interpolated to the model grid), the orbital velocity
C                    then from the local depth as in mode 2
C     sedWaveUbFile, sedWavePerFile, sedWaveDirFile :: files of mode 1,
C     sedWaveHsFile :: significant wave height file of mode 3 (with
C                    sedWavePerFile and sedWaveDirFile),
C                    sedWaveNrec records every sedWaveDt [s], cyclic,
C                    linear in time between records
C     sedWindSpeed :: 10 m wind speed [m/s] (mode 2)
C     sedWindDir   :: direction the wind blows toward [deg ccw from +x]
C     sedFetch     :: fetch [m] (mode 2)
C     sedWaveGamma :: depth-limited breaking, Hs <= sedWaveGamma*depth
C     sedZ0        :: bed roughness length for the wave friction factor
C                    f_w = 1.39*(A/sedZ0)**-0.52 [m] (Soulsby 1997)
C     sedWaveIce   :: .TRUE.: orbital velocity times (1 - sea-ice area)
C     Wave and current stresses combine as tau_m = tau_c*(1 + 1.2*
C     (tau_w/(tau_c+tau_w))**3.2), tau_b = |tau_m + tau_w| (vector sum at
C     the wave-current angle) (Soulsby 1995); tau_b drives erosion,
C     deposition (sedTauCD) and the active layer
      INTEGER sedWaveMode, sedWaveNrec
      _RL sedWaveDt, sedWindSpeed, sedWindDir, sedFetch
      _RL sedWaveGamma, sedZ0
      LOGICAL sedWaveIce
      CHARACTER*(MAX_LEN_FNAM) sedWaveUbFile, sedWavePerFile
      CHARACTER*(MAX_LEN_FNAM) sedWaveDirFile, sedWaveHsFile
      COMMON /SEDIMENT_WAVE_I/ sedWaveMode, sedWaveNrec
      COMMON /SEDIMENT_WAVE_R/ sedWaveDt, sedWindSpeed, sedWindDir,
     &     sedFetch, sedWaveGamma, sedZ0
      COMMON /SEDIMENT_WAVE_L/ sedWaveIce
      COMMON /SEDIMENT_WAVE_C/ sedWaveUbFile, sedWavePerFile,
     &     sedWaveDirFile, sedWaveHsFile

C--   Fields
C     sedBed   :: bed sediment mass per unit area [kg/m^2] (total,
C                 fluff included)
C     sedFluff :: part of sedBed in the fluff layer [kg/m^2]
C     sedTauB  :: bed shear stress of the latest step [Pa] (maximum
C                 wave + current stress with waves)
C     sedLay   :: mass of each class in each bed layer [kg/m^2]
C                 (sedNlay > 0; sedBed is their sum)
C     sedTcb   :: bulk critical stress of each bed layer [Pa]
C     sedZa    :: active layer thickness of the latest step [m]
C     sedTauW  :: wave bed stress [Pa]
C     sedUbw, sedPerW, sedDirW, sedHsW :: bottom orbital velocity [m/s],
C                 wave period [s], direction [rad, toward, ccw from +x]
C                 and significant height [m] (mode 2; 0 with files)
C     sedWrec  :: file records (1-based) held in sedWub0/1 etc. (modes 1,
C                 3; in mode 3 sedWub0/1 hold the significant height)
C     sedCumRiv, sedCumDep, sedCumEro :: cumulative global river input,
C                 deposition and erosion of each class [kg]
C     sedErrMax:: largest |change of suspended + bed mass - river input|
C                 made by SEDIMENT_STEP in one step, all classes [kg]
      _RL sedBed (1-OLx:sNx+OLx,1-OLy:sNy+OLy,sedMaxCl,nSx,nSy)
      _RL sedFluff(1-OLx:sNx+OLx,1-OLy:sNy+OLy,sedMaxCl,nSx,nSy)
      _RL sedTauB(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL sedCumRiv(sedMaxCl), sedCumDep(sedMaxCl)
      _RL sedCumEro(sedMaxCl)
      _RL sedErrMax
      COMMON /SEDIMENT_FIELDS_RL/
     &     sedBed, sedFluff, sedTauB, sedCumRiv, sedCumDep, sedCumEro,
     &     sedErrMax
      _RL sedLay (1-OLx:sNx+OLx,1-OLy:sNy+OLy,sedMaxLay,sedMaxCl,
     &            nSx,nSy)
      _RL sedTcb (1-OLx:sNx+OLx,1-OLy:sNy+OLy,sedMaxLay,nSx,nSy)
      _RL sedZa  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /SEDIMENT_BED_FIELDS/ sedLay, sedTcb, sedZa
      _RL sedTauW(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL sedUbw (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL sedPerW(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL sedDirW(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL sedHsW (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL sedWub0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL sedWub1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL sedWpe0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL sedWpe1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL sedWdi0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL sedWdi1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /SEDIMENT_WAVE_FIELDS/ sedTauW, sedUbw, sedPerW, sedDirW,
     &     sedHsW, sedWub0, sedWub1, sedWpe0, sedWpe1, sedWdi0, sedWdi1
      INTEGER sedWrec
      COMMON /SEDIMENT_WAVE_REC/ sedWrec
