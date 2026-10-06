# pkg/dic

Simple dissolved inorganic carbon biogeochemistry (DIC, ALK, PO4, DOP, O2, Fe) with carbonate chemistry and air-sea CO2 flux (OCMIP-style).

**pkg_depend:** +gchem  (`+` requires, `-` excludes)
**runtime switch:** `useDIC`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.dic`
**manual:** `doc/examples/global_oce_biogeo/global_oce_biogeo.rst`
**adjoint support files:** dic_ad_check_lev1_dir.h, dic_ad_check_lev2_dir.h, dic_ad_check_lev3_dir.h, dic_ad_check_lev4_dir.h, dic_ad_diff.list

## Namelist parameters
### ABIOTIC_PARMS
- `permil`
- `Pa2Atm`
- `selectBTconst` — estimates borate concentration from salinity:
- `selectFTconst` — estimates fluoride concentration from salinity:
- `selectHFconst` — sets the first dissociation constant for hydrogen fluoride:
- `selectK1K2const` — sets the 1rst & 2nd dissociation constants of carbonic acid:
- `selectPHsolver` — sets the pH solver to use:
- `useCalciteSaturation` — Dissolve calcium carbonate only below saturation horizon (needs DIC_CALCITE_SAT to be defined)
- `calcOmegaCalciteFreq` — Frequency at which 3D calcite saturation state, omegaC, is updated (s).
- `nIterCO3` — Number of iterations of the Follows 3D pH solver to calculate deep carbonate ion concenetration (no effect when using the Munhoven/SolveSapHe solvers).
- `selectCalciteDissolution` — flag to control calcite dissolution rate method: =0 : Constant dissolution rate; =1 : Follows (default) ; =2 : Keir (1980) Geochem. Cosmochem. Acta. ; =3 : Naviaux et al. 2019, Marine Chemistry
- `WsinkPIC` — sinking speed (m/s) of particulate inorganic carbon for calculation of calcite dissolution through the watercolumn
- `selectCalciteBottomRemin` — to either remineralize in bottom or top layer if flux reaches bottom layer; =0 : bottom, =1 : top
- `calciteDissolRate` — Rate constant (%) for calcite dissolution from Keir (1980) Geochem. Cosmochem. Acta.
- `calciteDissolExp` — Rate exponent for calcite dissolution from Keir (1980) Geochem. Cosmochem. Acta.
- `zca` — Scale depth for CaCO3 remineralization [m]
### BIOTIC_PARMS
- `DOPfraction` — fraction of new production going to DOP  _[ifdef DIC_BIOTIC]_
- `KDOPRemin` — DOP remineralization rate [1/s]  _[ifdef DIC_BIOTIC]_
- `KRemin` — remineralization power law coeffient  _[ifdef DIC_BIOTIC]_
- `zcrit` — Minimum Depth (m) over which biological activity is computed  _[ifdef DIC_BIOTIC]_
- `O2crit` — critical oxygen level [mol/m3]  _[ifdef DIC_BIOTIC]_
- `R_OP` — stochiometric ratios of nutrients  _[ifdef DIC_BIOTIC]_
- `R_CP` — stochiometric ratios of nutrients  _[ifdef DIC_BIOTIC]_
- `R_NP` — stochiometric ratios of nutrients (assumption of stoichometry of plankton and particulate  and dissolved organic matter)  _[ifdef DIC_BIOTIC]_
- `R_FeP` — stochiometric ratios of nutrients (assumption of stoichometry of plankton and particulate  and dissolved organic matter)  _[ifdef DIC_BIOTIC]_
- `parfrac` — fraction of Qsw that is PAR  _[ifdef DIC_BIOTIC]_
- `k0` — light attentuation coefficient of water [1/m]  _[ifdef DIC_BIOTIC]_
- `lit0` — half saturation constants for phosphate [mol P/m3], iron [mol Fe/m3] and light [W/m2]  _[ifdef DIC_BIOTIC]_
- `KPO4` — half saturation constants for phosphate [mol P/m3], iron [mol Fe/m3] and light [W/m2]  _[ifdef DIC_BIOTIC]_
- `KFE` — half saturation constants for phosphate [mol P/m3], iron [mol Fe/m3] and light [W/m2]  _[ifdef DIC_BIOTIC]_
- `kchl` — light attentuation coefficient of chlorophyll [m2/mg]  _[ifdef DIC_BIOTIC]_
- `alpfe` — solubility of aeolian fe [fraction]  _[ifdef DIC_BIOTIC]_
- `fesedflux_pcm` — ratio of sediment iron to sinking organic matter  _[ifdef DIC_BIOTIC]_
- `FeIntSec` — Sediment Fe flux, intersect value in: Fe_flux = fesedflux_pcm*pflux + FeIntSec  _[ifdef DIC_BIOTIC]_
- `freefemax` — max soluble free iron [mol/m3]  _[ifdef DIC_BIOTIC]_
- `KScav` — iron scavenging rate [1/s]  _[ifdef DIC_BIOTIC]_
- `ligand_stab` — ligand-free iron stability constant [m3/mol]  _[ifdef DIC_BIOTIC]_
- `ligand_tot` — uniform, invariant total free ligand conc [mol/m3]  _[ifdef DIC_BIOTIC]_
- `alphaUniform` — read in alphaUniform to fill in 2d array alpha  _[ifdef DIC_BIOTIC]_
- `rainRatioUniform` — read in rainRatioUniform to fill in 2d array rain_ratio  _[ifdef DIC_BIOTIC]_
### DIC_FORCING
- `DIC_windFile` — file name of wind speeds
- `DIC_atmospFile` — file name of atmospheric pressure
- `DIC_silicaFile` — file name of surface silica
- `DIC_deepSilicaFile` — file name of 3D silica fields
- `DIC_iceFile` — file name of seaice fraction
- `DIC_parFile` — file name of photosynthetically available radiation (PAR)
- `DIC_chlaFile` — file name of chlorophyll climatology
- `DIC_ironFile` — file name of aeolian iron flux
- `DIC_forcingPeriod` — periodic forcing parameter specific for dic (seconds)
- `DIC_forcingCycle` — periodic forcing parameter specific for dic (seconds)
- `dic_int1`
- `dic_int2` — number pCO2 entries to read from file
- `dic_int3` — start timestep
- `dic_int4` — timestep between file entries
- `dic_pCO2` — atmospheric pCO2 to be read from data.dic

## CPP options (defaults as shipped)
- `CARBONCHEM_SOLVESAPHE` (undef, DIC_OPTIONS.h) — Compile Munhoven (2013) "Solvesaphe" package for pH/pCO2 can still select Follows et al (2006) solver in data.dic, but will use solvesaphe dissociation coefficient options.
- `CARBONCHEM_TOTALPHSCALE` (undef, DIC_OPTIONS.h) — consistent with other coefficients (currently on the seawater scale). NOTE: Has NO effect when CARBONCHEM_SOLVESAPHE is defined (different coeffs are used).
- `DIC_BIOTIC` (define, DIC_OPTIONS.h) — BIOTIC OPTIONS
- `ALLOW_O2` (define, DIC_OPTIONS.h)
- `ALLOW_FE` (undef, DIC_OPTIONS.h)
- `READ_PAR` (undef, DIC_OPTIONS.h)
- `MINFE` (undef, DIC_OPTIONS.h)
- `DIC_NO_NEG` (undef, DIC_OPTIONS.h)
- `DIC_BOUNDS` (undef, DIC_OPTIONS.h)
- `USE_QSW` (undef, DIC_OPTIONS.h) — these all need to be defined for coupling to atmospheric model:
- `USE_QSW_UNDERICE` (undef, DIC_OPTIONS.h)
- `USE_PLOAD` (undef, DIC_OPTIONS.h)
- `ALLOW_OLD_VIRTUALFLUX` (undef, DIC_OPTIONS.h) — use surface salinity forcing (scaled by mean surf value) for DIC & ALK forcing
- `WATERVAP_BUG` (undef, DIC_OPTIONS.h) — put back bugs related to Water-Vapour in carbonate chemistry & air-sea fluxes
- `DIC_CALCITE_SAT` (undef, DIC_OPTIONS.h) — dissolution only below saturation horizon following method by Karsten Friis
- `LIGHT_CHL` (undef, DIC_OPTIONS.h) — Include self-shading effect by phytoplankton
- `SEDFE` (undef, DIC_OPTIONS.h) — Include iron sediment source using DOP flux
- `DIC_AD_SAFE` (undef, DIC_OPTIONS.h) — For Adjoint built

## Headers
- `DIC_ATMOS.h` — 
- `DIC_COST.h` — control variables QQ      _RL alpha(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy) QQ      _RL rain_ratio(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy) QQ      _RL KScav
- `DIC_CTRL.h` — DIC_XX.h Control of Biological Carbon Variables
- `DIC_LOAD.h` — --   COMMON /DIC_LOAD/ DIC_ldRec     :: time-record currently loaded (in temp arrays *[1]) chlinput      :: chlorophyll climatology input field [mg/m3
- `DIC_OPTIONS.h` — BOP
- `DIC_VARS.h` — Abiotic Carbon Variables
- `dic_ad_check_lev1_dir.h` — common CARBON_NEEDS ADJ STORE pH                 = comlev1, key = ikey_dynamics ADJ STORE pCO2               = comlev1, key = ikey_dynamics ADJ STORE 
- `dic_ad_check_lev2_dir.h` — common CARBON_NEEDS ADJ STORE pH                 = tapelev2, key = ilev_2 CADJ STORE pCO2               = tapelev2, key = ilev_2 CADJ STORE fIce      
- `dic_ad_check_lev3_dir.h` — common CARBON_NEEDS ADJ STORE pH                 = tapelev3, key = ilev_3 CADJ STORE pCO2               = tapelev3, key = ilev_3 CADJ STORE fIce      
- `dic_ad_check_lev4_dir.h` — common CARBON_NEEDS ADJ STORE pH                 = tapelev4, key = ilev_4 CADJ STORE pCO2               = tapelev4, key = ilev_4 CADJ STORE fIce      

## Routines (44)
`alk_surfforcing.F`, `bio_export.F`, `calcite_saturation.F`, `car_flux.F`, `car_flux_omega_top.F`, `carbon_chem.F`, `dic_atmos.F`, `dic_biotic_diags.F`, `dic_biotic_forcing.F`, `dic_biotic_init.F`, `dic_cost.F`, `dic_diagnostics_init.F`, `dic_fields_load.F`, `dic_fields_update.F`, `dic_ini_atmos.F`, `dic_ini_forcing.F`, `dic_init_fixed.F`, `dic_init_varia.F`, `dic_mnc_init.F`, `dic_read_co2_pickup.F`, `dic_read_pickup.F`, `dic_readparms.F`, `dic_set_control.F`, `dic_solvesaphe.F`, `dic_store_fluxco2.F`, `dic_surfforcing.F`, `dic_surfforcing_init.F`, `dic_tr_register.F`, `dic_write_pickup.F`, `fe_chem.F`, `insol.F`, `o2_surfforcing.F`, `phos_flux.F`

## Called from outside the package
- `CALC_PCO2_APPROX` ← `pkg/bling/bling_airseaflux.F:211`
- `CALC_PCO2_SOLVESAPHE` ← `pkg/bling/bling_airseaflux.F:197`
- `CARBON_COEFFS` ← `pkg/bling/bling_airseaflux.F:159`
- `DIC_COEFFS_SURF` ← `pkg/bling/bling_airseaflux.F:152`
- `AHINI_FOR_AT` ← `pkg/bling/bling_carbonate_init.F:239`
- `CALC_PCO2_APPROX` ← `pkg/bling/bling_carbonate_init.F:266`
- `CALC_PCO2_SOLVESAPHE` ← `pkg/bling/bling_carbonate_init.F:251`
- `CARBON_COEFFS_PRESSURE_DEP` ← `pkg/bling/bling_carbonate_init.F:218`
- `DIC_COEFFS_DEEP` ← `pkg/bling/bling_carbonate_init.F:210`
- `DIC_COEFFS_SURF` ← `pkg/bling/bling_carbonate_init.F:205`
- `CALC_PCO2_APPROX` ← `pkg/bling/bling_carbonate_sys.F:212`
- `CALC_PCO2_SOLVESAPHE` ← `pkg/bling/bling_carbonate_sys.F:189`
- `CARBON_COEFFS_PRESSURE_DEP` ← `pkg/bling/bling_carbonate_sys.F:134`
- `DIC_COEFFS_DEEP` ← `pkg/bling/bling_carbonate_sys.F:121`
- `DIC_COEFFS_SURF` ← `pkg/bling/bling_carbonate_sys.F:116`
- `ANW_INFSUP` ← `pkg/bling/bling_solvesaphe.F:1541,1800`
- `EQUATION_AT` ← `pkg/bling/bling_solvesaphe.F:1400,1450,1598,1851`
- `SOLVE_AT_FAST` ← `pkg/bling/bling_solvesaphe.F:322`
- `SOLVE_AT_GENERAL` ← `pkg/bling/bling_solvesaphe.F:304`
- `SOLVE_AT_GENERAL_SEC` ← `pkg/bling/bling_solvesaphe.F:313`
- `DIC_FIELDS_LOAD` ← `pkg/gchem/gchem_fields_load.F:35`
- `DIC_ATMOS` ← `pkg/gchem/gchem_forcing_sep.F:290,304`
- `DIC_BIOTIC_FORCING` ← `pkg/gchem/gchem_forcing_sep.F:150,160,168`
- `DIC_COST` ← `pkg/gchem/gchem_forcing_sep.F:309`
- `DIC_STORE_FLUXCO2` ← `pkg/gchem/gchem_forcing_sep.F:306`
- `DIC_INIT_FIXED` ← `pkg/gchem/gchem_init_fixed.F:44`
- `DIC_INIT_VARIA` ← `pkg/gchem/gchem_init_vari.F:69`
- `DIC_SURFFORCING_INIT` ← `pkg/gchem/gchem_init_vari.F:97`
- `DIC_READPARMS` ← `pkg/gchem/gchem_readparms.F:158`
- `DIC_TR_REGISTER` ← `pkg/gchem/gchem_tr_register.F:53`
- `DIC_WRITE_PICKUP` ← `pkg/gchem/gchem_write_pickup.F:43,60`

## Verification experiments compiling it (3)
`so_box_biogeo` `tutorial_dic_adjoffline` `tutorial_global_oce_biogeo`
