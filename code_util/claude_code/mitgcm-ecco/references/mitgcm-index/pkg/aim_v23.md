# pkg/aim_v23

AIM intermediate-complexity atmospheric physics (SPEEDY v23: convection, clouds, radiation, surface fluxes) for atmosphere set-ups.

**pkg_depend:** +atm_common  (`+` requires, `-` excludes)
**runtime switch:** `useAIM_V23`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.aimphys`
**manual:** `doc/examples/examples.rst`, `doc/outp_pkgs/outp_pkgs.rst`

## Namelist parameters
### AIM_PARAMS
- `aim_useFMsurfBC` — select surface B.C. from Franco Molteni  _[ifdef ALLOW_AIM]_
- `aim_useMMsurfFc` — select Monthly Mean surface forcing (e.g., NCEP)  _[ifdef ALLOW_AIM]_
- `aim_surfForc_TimePeriod` — Length of forcing time period (e.g. 1 month)  _[ifdef ALLOW_AIM]_
- `aim_surfForc_NppCycle` — Number of time period per Cycle (e.g. 12)  _[ifdef ALLOW_AIM]_
- `aim_surfForc_TransRatio` — transition ratio from one month to the next  _[ifdef ALLOW_AIM]_
- `aim_surfPotTemp` — surf.Temp input file is in Pot.Temp (aim_useMMsurfFc)  _[ifdef ALLOW_AIM]_
- `aim_energPrecip` — account for energy of precipitation (snow & rain temp)  _[ifdef ALLOW_AIM]_
- `aim_splitSIOsFx` — compute separately Sea-Ice & Ocean surf. Flux (also land SW & LW) ; default=F as in original version  _[ifdef ALLOW_AIM]_
- `aim_MMsufx` — sufix for all Monthly Mean surface forcing files  _[ifdef ALLOW_AIM]_
- `aim_MMsufxLength` — Length of sufix (Monthly Mean surf. forcing files)  _[ifdef ALLOW_AIM]_
- `aim_LandFile` — file name for Land fraction  _[ifdef ALLOW_AIM]_
- `aim_albFile` — file name for Albedo input file   (F.M. surfBC)  _[ifdef ALLOW_AIM]_
- `aim_vegFile` — file name for vegetation fraction (F.M. surfBC)  _[ifdef ALLOW_AIM]_
- `aim_sstFile` — file name for  Sea.Surf.Temp      (F.M. surfBC)  _[ifdef ALLOW_AIM]_
- `aim_lstFile` — file name for Land.Surf.Temp      (F.M. surfBC)  _[ifdef ALLOW_AIM]_
- `aim_oiceFile` — file name for Sea Ice fraction    (F.M. surfBC)  _[ifdef ALLOW_AIM]_
- `aim_snowFile` — file name for Snow depth          (F.M. surfBC)  _[ifdef ALLOW_AIM]_
- `aim_swcFile` — file name for Soil Water content  (F.M. surfBC)  _[ifdef ALLOW_AIM]_
- `aim_qfxFile` — file name for ocean q-flux  _[ifdef ALLOW_AIM]_
- `aim_dragStrato` — stratospheric-drag damping time scale (s)  _[ifdef ALLOW_AIM]_
- `aim_selectOceAlbedo` — select ocean albedo scheme:  =0: constant (default)  _[ifdef ALLOW_AIM]_
- `aim_select_pCO2` — select AIM CO2 formulation:  _[ifdef ALLOW_AIM]_
- `aim_abs_pCO2` — pCO2 dependence coeff. of CO2 band LW absortion  _[ifdef ALLOW_AIM]_
- `aim_fixed_pCO2`  _[ifdef ALLOW_AIM]_
- `atmpCO2init`  _[ifdef ALLOW_AIM]_
- `aim_clrSkyDiag` — compute clear-sky radiation for diagnostics  _[ifdef ALLOW_AIM]_
- `aim_taveFreq`  _[ifdef ALLOW_AIM]_
- `aim_diagFreq` — Frequency^-1 for diagnostic output (s)  _[ifdef ALLOW_AIM]_
- `aim_tendFreq` — Frequency^-1 for tendencies output (s)  _[ifdef ALLOW_AIM]_
- `aim_timeave_mnc`  _[ifdef ALLOW_AIM]_
- `aim_snapshot_mnc`  _[ifdef ALLOW_AIM]_
- `aim_pickup_write_mnc`  _[ifdef ALLOW_AIM]_
- `aim_pickup_read_mnc`  _[ifdef ALLOW_AIM]_
### AIM_PAR_FOR
- `SOLC`  _[ifdef ALLOW_AIM]_
- `ALBSEA`  _[ifdef ALLOW_AIM]_
- `ALBICE`  _[ifdef ALLOW_AIM]_
- `ALBSN`  _[ifdef ALLOW_AIM]_
- `SDALB`  _[ifdef ALLOW_AIM]_
- `SWCAP`  _[ifdef ALLOW_AIM]_
- `SWWIL`  _[ifdef ALLOW_AIM]_
- `hSnowWetness` — snow depth (m) corresponding to maximum wetness  _[ifdef ALLOW_AIM]_
- `OBLIQ` — Obliquity (in degree) used with ALLOW_INSOLATION  _[ifdef ALLOW_AIM]_
### AIM_PAR_SFL
- `FWIND0`  _[ifdef ALLOW_AIM]_
- `FTEMP0`  _[ifdef ALLOW_AIM]_
- `FHUM0`  _[ifdef ALLOW_AIM]_
- `CDL`  _[ifdef ALLOW_AIM]_
- `CDS`  _[ifdef ALLOW_AIM]_
- `CHL`  _[ifdef ALLOW_AIM]_
- `CHS`  _[ifdef ALLOW_AIM]_
- `VGUST`  _[ifdef ALLOW_AIM]_
- `CTDAY`  _[ifdef ALLOW_AIM]_
- `DTHETA`  _[ifdef ALLOW_AIM]_
- `dTstab`  _[ifdef ALLOW_AIM]_
- `FSTAB`  _[ifdef ALLOW_AIM]_
- `HDRAG`  _[ifdef ALLOW_AIM]_
- `FHDRAG`  _[ifdef ALLOW_AIM]_
### AIM_PAR_CNV
- `PSMIN`  _[ifdef ALLOW_AIM]_
- `TRCNV`  _[ifdef ALLOW_AIM]_
- `QBL`  _[ifdef ALLOW_AIM]_
- `RHBL`  _[ifdef ALLOW_AIM]_
- `RHIL`  _[ifdef ALLOW_AIM]_
- `ENTMAX`  _[ifdef ALLOW_AIM]_
- `SMF`  _[ifdef ALLOW_AIM]_
### AIM_PAR_LSC
- `TRLSC`  _[ifdef ALLOW_AIM]_
- `RHLSC`  _[ifdef ALLOW_AIM]_
- `DRHLSC`  _[ifdef ALLOW_AIM]_
- `QSMAX`  _[ifdef ALLOW_AIM]_
### AIM_PAR_RAD
- `RHCL1`  _[ifdef ALLOW_AIM]_
- `RHCL2`  _[ifdef ALLOW_AIM]_
- `QACL1`  _[ifdef ALLOW_AIM]_
- `QACL2`  _[ifdef ALLOW_AIM]_
- `ALBCL`  _[ifdef ALLOW_AIM]_
- `EPSSW`  _[ifdef ALLOW_AIM]_
- `EPSLW`  _[ifdef ALLOW_AIM]_
- `EMISFC`  _[ifdef ALLOW_AIM]_
- `ABSDRY`  _[ifdef ALLOW_AIM]_
- `ABSAER`  _[ifdef ALLOW_AIM]_
- `ABSWV1`  _[ifdef ALLOW_AIM]_
- `ABSWV2`  _[ifdef ALLOW_AIM]_
- `ABSCL1`  _[ifdef ALLOW_AIM]_
- `ABSCL2`  _[ifdef ALLOW_AIM]_
- `ABLWIN`  _[ifdef ALLOW_AIM]_
- `ABLCO2`  _[ifdef ALLOW_AIM]_
- `ABLWV1`  _[ifdef ALLOW_AIM]_
- `ABLWV2`  _[ifdef ALLOW_AIM]_
- `ABLCL1`  _[ifdef ALLOW_AIM]_
- `ABLCL2`  _[ifdef ALLOW_AIM]_
### AIM_PAR_VDI
- `TRSHC`  _[ifdef ALLOW_AIM]_
- `TRVDI`  _[ifdef ALLOW_AIM]_
- `TRVDS`  _[ifdef ALLOW_AIM]_
- `RHGRAD`  _[ifdef ALLOW_AIM]_
- `SEGRAD`  _[ifdef ALLOW_AIM]_

## CPP options (defaults as shipped)
- `ALLOW_DEW_ON_LAND` (undef, AIM_OPTIONS.h) — allow dew to form on land (=negative evaporation)
- `ALLOW_INSOLATION` (undef, AIM_OPTIONS.h) — calculate top-atmosphere insolation using orbital parameters (obliquity, eccentricity ...) provided as run-time params
- `ALLOW_CLOUD_3D` (undef, AIM_OPTIONS.h) — allow 3D cloud fraction for computation of radiation
- `ALLOW_AIM_CO2` (undef, AIM_OPTIONS.h) — allow CO2 concentration
- `ALLOW_CLR_SKY_DIAG` (define, AIM_OPTIONS.h) — allow Clear-Sky diagnostic:
- `_KD2KA` (define, AIM_OPTIONS.h) — Macro mapping dynamics vertical indexing (KD) to AIM vertical indexing (KA). ( dynamics puts K=1 at bottom of atmos., AIM puts K=1 at top of atmos. )

## Headers
- `AIM2DYN.h` — AIM output fields in dynamics conforming arrays
- `AIM_CO2.h` — AIM CO2 fields.
- `AIM_FFIELDS.h` — AIM (surface) forcing fields.
- `AIM_GRID.h` — --   COMMON /AIM_GRID_R/: AIM surface and grid-related arrays WVSurf  : weights for vertical interpolation down to the surface (replace WVI(NLEV,2) in
- `AIM_OPTIONS.h` — CPP options file for AIM package
- `AIM_PARAMS.h` — Header file for AIM package parameters e.g.: output/input file & parameters; forcing & interface parameters;
- `AIM_SIZE.h` — MITgcm declaration of grid size. Latitudinal extent is one less than MITgcm ( i.e. NY-1) because MITgcm has dummy layer of land at northern most edge.
- `com_cnvcon.h` — --   COMMON /CNVCON/: Convection constants (initial. in INPHYS) PSMIN  = minimum (norm.) sfc. pressure for the occurrence of convection TRCNV  = time 
- `com_forcing.h` — -- Note: Variables which do not need to stay in common block (local var) are declare locally in each S/R that use them (commented with "cL"); Some var
- `com_forcon.h` — --   COMMON /FORCON/: Constants for forcing fields (initial. in INPHYS) SOLC   = Solar constant (area averaged) in W/m^2 ALBSEA = Albedo over sea ALBI
- `com_lsccon.h` — --   COMMON /LSCCON/: Constants for large-scale condendation (initial. in INPHYS) TRLSC  = Relaxation time (in hours) for specific humidity RHLSC  = M
- `com_physcon.h` — --   COMMON /PHYCON/: Physical constants (initial. in INPHYS) P0    = reference pressure                 [Pa=N/m2] GG    = gravity accel.             
- `com_physvar.h` — -- Note: Variables which do not need to stay in common block (local var) are declare locally in each S/R that use them (commented with "cL"); Some var
- `com_radcon.h` — --   COMMON /RADCON/: Radiation constants (initial. in INPHYS) RHCL1  = relative hum. corresponding to cloud cover = 0 RHCL2  = relative hum. correspo
- `com_radvar.h` — 2nd part of original file "com_radcon.h": contains temp. variables used within radiation scheme and passed as arguments to SOL_OZ, RADSW & RADLW (orig
- `com_sflcon.h` — --   COMMON /SFLCON/: Constants for surface fluxes (initial. in INPHYS) FWIND0 = ratio of near-sfc wind to lowest-level wind FTEMP0 = weight for near-
- `com_vdicon.h` — --   COMMON /VDICON/: Constants for vertical diffusion and shallow convection (initial. in INPHYS) TRSHC  = relaxation time (in hours) for shallow con
- `phy_const.h` — --   Constants for physical parametrization routines:

## Routines (42)
`aim_aim2dyn.F`, `aim_aim2land.F`, `aim_aim2sioce.F`, `aim_diagnostics.F`, `aim_diagnostics_init.F`, `aim_do_co2.F`, `aim_do_physics.F`, `aim_dyn2aim.F`, `aim_fields_load.F`, `aim_initialise.F`, `aim_land2aim.F`, `aim_land_impl.F`, `aim_mnc_init.F`, `aim_readparms.F`, `aim_sice2aim.F`, `aim_sice_impl.F`, `aim_surf_bc.F`, `aim_tendency_apply.F`, `aim_write_local.F`, `aim_write_phys.F`, `aim_write_tave.F`, `phy_convmf.F`, `phy_driver.F`, `phy_inphys.F`, `phy_lscond.F`, `phy_radiat.F`, `phy_shtorh.F`, `phy_snow_precip.F`, `phy_suflux_land.F`, `phy_suflux_ocean.F`, `phy_suflux_post.F`, `phy_suflux_prep.F`, `phy_suflux_sice.F`, `phy_vdifsc.F`

## Called from outside the package
- `AIM_TENDENCY_APPLY_S` ← `model/src/apply_forcing.F:854`
- `AIM_TENDENCY_APPLY_T` ← `model/src/apply_forcing.F:484`
- `AIM_TENDENCY_APPLY_U` ← `model/src/apply_forcing.F:106`
- `AIM_TENDENCY_APPLY_V` ← `model/src/apply_forcing.F:296`
- `AIM_DO_PHYSICS` ← `model/src/do_atmospheric_phys.F:146`
- `AIM_TENDENCY_APPLY_S` ← `model/src/external_forcing.F:684`
- `AIM_TENDENCY_APPLY_T` ← `model/src/external_forcing.F:362`
- `AIM_TENDENCY_APPLY_U` ← `model/src/external_forcing.F:69`
- `AIM_TENDENCY_APPLY_V` ← `model/src/external_forcing.F:209`
- `AIM_FIELDS_LOAD` ← `model/src/load_fields_driver.F:256`
- `AIM_INITIALISE` ← `model/src/packages_init_fixed.F:561`

## Verification experiments compiling it (4)
`aim.5l_Equatorial_Channel` `aim.5l_LatLon` `aim.5l_cs` `cpl_aim+ocn`
