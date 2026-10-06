# pkg/bling

BLING biogeochemistry (Biogeochemistry with Light, Iron, Nutrients and Gases; Galbraith et al.) via gchem/ptracers.

**pkg_depend:** +gchem -dic  (`+` requires, `-` excludes)
**runtime switch:** `useBLING`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.bling`
**manual:** `doc/examples/examples.rst`
**adjoint support files:** bling_ad_check_lev1_dir.h, bling_ad_check_lev2_dir.h, bling_ad_check_lev3_dir.h, bling_ad_check_lev4_dir.h, bling_ad_diff.list

## Namelist parameters
### ABIOTIC_PARMS
- `permil`
- `Pa2Atm`
- `epsln`
- `selectBTconst` — estimates borate concentration from salinity:  _[ifdef CARBONCHEM_SOLVESAPHE]_
- `selectFTconst` — estimates fluoride concentration from salinity:  _[ifdef CARBONCHEM_SOLVESAPHE]_
- `selectHFconst` — sets the first dissociation constant for hydrogen fluoride:  _[ifdef CARBONCHEM_SOLVESAPHE]_
- `selectK1K2const` — sets the 1rst & 2nd dissociation constants of carbonic acid:  _[ifdef CARBONCHEM_SOLVESAPHE]_
- `selectPHsolver` — sets the pH solver to use:  _[ifdef CARBONCHEM_SOLVESAPHE]_
### BIOTIC_PARMS
- `CtoN`
- `CtoP`
- `NtoP`
- `HtoC`
- `NO3toN`
- `O2toN`
- `O2toP`
- `CatoN`
- `CatoP`
- `masstoN`
- `Pc_0_diaz`  _[ifndef USE_BLING_V1]_
- `alpha_photo`  _[ifndef USE_BLING_V1]_
- `gamma_DON`  _[ifndef USE_BLING_V1]_
- `k_Fe_diaz`  _[ifndef USE_BLING_V1]_
- `k_NO3`  _[ifndef USE_BLING_V1]_
- `k_PtoN`  _[ifndef USE_BLING_V1]_
- `k_FetoN`  _[ifndef USE_BLING_V1]_
- `k_NO3_sm`  _[ifndef USE_BLING_V1]_
- `k_NO3_lg`  _[ifndef USE_BLING_V1]_
- `k_PO4_sm`  _[ifndef USE_BLING_V1]_
- `k_PO4_lg`  _[ifndef USE_BLING_V1]_
- `k_Fe_sm`  _[ifndef USE_BLING_V1]_
- `k_Fe_lg`  _[ifndef USE_BLING_V1]_
- `PtoN_min`  _[ifndef USE_BLING_V1]_
- `PtoN_max`  _[ifndef USE_BLING_V1]_
- `FetoN_min`  _[ifndef USE_BLING_V1]_
- `FetoN_max`  _[ifndef USE_BLING_V1]_
- `kappa_eppley_diaz`  _[ifndef USE_BLING_V1]_
- `phi_dvm`  _[ifndef USE_BLING_V1]_
- `sigma_dvm`  _[ifndef USE_BLING_V1]_
- `k_Si`  _[ifndef USE_BLING_V1 & ifdef USE_SIBLING]_
- `gamma_Si_0`  _[ifndef USE_BLING_V1 & ifdef USE_SIBLING]_
- `kappa_remin_Si`  _[ifndef USE_BLING_V1 & ifdef USE_SIBLING]_
- `wsink_Si`  _[ifndef USE_BLING_V1 & ifdef USE_SIBLING]_
- `SitoN_uptake_min`  _[ifndef USE_BLING_V1 & ifdef USE_SIBLING]_
- `SitoN_uptake_max`  _[ifndef USE_BLING_V1 & ifdef USE_SIBLING]_
- `SitoN_uptake_scale`  _[ifndef USE_BLING_V1 & ifdef USE_SIBLING]_
- `SitoN_uptake_exp`  _[ifndef USE_BLING_V1 & ifdef USE_SIBLING]_
- `q_SitoN_diss`  _[ifndef USE_BLING_V1 & ifdef USE_SIBLING]_
- `alpha_max`  _[NOT(ifndef USE_BLING_V1)]_
- `alpha_min`  _[NOT(ifndef USE_BLING_V1)]_
- `gamma_biomass`  _[NOT(ifndef USE_BLING_V1)]_
- `k_FetoP`  _[NOT(ifndef USE_BLING_V1)]_
- `FetoP_max`  _[NOT(ifndef USE_BLING_V1)]_
- `Fe_lim_min`  _[NOT(ifndef USE_BLING_V1)]_
- `pivotal`
- `Pc_0`
- `lambda_0`
- `resp_frac`
- `chl_min`
- `theta_Fe_max_hi`
- `theta_Fe_max_lo`
- `gamma_irr_mem`
- `gamma_DOP`
- `gamma_POM`
- `k_Fe`
- `k_O2`
- `k_PO4`
- `kFe_eq_lig_max`
- `kFe_eq_lig_min`
- `kFe_eq_lig_Femin`
- `kFe_eq_lig_irr`
- `kFe_org`
- `kFe_inorg`
- `FetoC_sed`
- `remin_min`
- `oxic_min`
- `ligand`
- `kappa_eppley`
- `kappa_remin`
- `ca_remin_depth`
- `phi_DOM`
- `phi_sm`
- `phi_lg`
- `wsink0z`
- `wsink0`
- `wsinkacc`
- `parfrac`
- `alpfe`
- `k0`
- `MLmix_max`
- `chlsat_locTimWindow`
- `river_conc_po4`
- `river_dom_to_nut`
### BLING_FORCING
- `bling_windFile` — file name of wind speeds
- `bling_atmospFile` — file name of atmospheric pressure
- `bling_iceFile` — file name of sea ice fraction
- `bling_ironFile` — file name of aeolian iron flux
- `bling_silicaFile` — file name of surface silica
- `bling_psmFile` — file name of init small phyto biomass
- `bling_plgFile` — file name of init lg phyto biomass
- `bling_PdiazFile` — file name of init diaz biomass
- `bling_forcingPeriod` — period of forcing for biogeochemistry (seconds)
- `bling_forcingCycle` — periodic forcing parameter for biogeochemistry
- `bling_Pc_2dFile`
- `bling_Pc_2d_diazFile`
- `bling_k0_2dFile` — File containing a 2D spatial field of light attenuation coefficient (k0_2d, in m^-1). This coefficient regulates underwater light availability in the BLING model. If not specified, a constant value k0 (default= 0.04 m^-1) is applied for entire domain.
- `bling_alpha_photo2dFile`
- `bling_phi_DOM2dFile`
- `bling_k_Fe2dFile`
- `bling_k_Fe_diaz2dFile`
- `bling_gamma_POM2dFile`
- `bling_wsink0_2dFile`
- `bling_phi_sm2dFile`
- `bling_phi_lg2dFile`
- `bling_pCO2` — Atmospheric pCO2 to be read in data.bling
- `river_conc_po4`
- `river_dom_to_nut`
- `apco2file`  _[ifdef ALLOW_EXF]_
- `apco2startdate1`  _[ifdef ALLOW_EXF]_
- `apco2startdate2`  _[ifdef ALLOW_EXF]_
- `apco2RepCycle`  _[ifdef ALLOW_EXF]_
- `apco2period`  _[ifdef ALLOW_EXF]_
- `apco2StartTime`  _[ifdef ALLOW_EXF]_
- `exf_inscal_apco2`  _[ifdef ALLOW_EXF]_
- `exf_outscal_apco2`  _[ifdef ALLOW_EXF]_
- `apco2const`  _[ifdef ALLOW_EXF]_
- `apco2_exfremo_intercept`  _[ifdef ALLOW_EXF]_
- `apco2_exfremo_slope`  _[ifdef ALLOW_EXF]_
- `apco2_lon0`  _[ifdef ALLOW_EXF & ifdef USE_EXF_INTERPOLATION]_
- `apco2_lon_inc`  _[ifdef ALLOW_EXF & ifdef USE_EXF_INTERPOLATION]_
- `apco2_lat0`  _[ifdef ALLOW_EXF & ifdef USE_EXF_INTERPOLATION]_
- `apco2_lat_inc`  _[ifdef ALLOW_EXF & ifdef USE_EXF_INTERPOLATION]_
- `apco2_nlon`  _[ifdef ALLOW_EXF & ifdef USE_EXF_INTERPOLATION]_
- `apco2_nlat`  _[ifdef ALLOW_EXF & ifdef USE_EXF_INTERPOLATION]_
- `apco2_interpMethod`  _[ifdef ALLOW_EXF & ifdef USE_EXF_INTERPOLATION]_

## CPP options (defaults as shipped)
- `USE_BLING_V1` (undef, BLING_OPTIONS.h) — of BLING with 8 tracers and 3 phyto classes. For the original 6-tracer model of Galbraith et al (2010), define USE_BLING_V1 - but note the different order of tracers in data.ptracers
- `USE_SIBLING` (undef, BLING_OPTIONS.h) — Options for BLING+Nitrogen code: SiBLING: add a 9th tracer for silica
- `USE_BLING_DVM` (undef, BLING_OPTIONS.h) — apply remineralization from diel vertical migration
- `ADVECT_PHYTO` (undef, BLING_OPTIONS.h) — active tracer for total phytoplankton biomass
- `BLING_NO_NEG` (define, BLING_OPTIONS.h) — Prevents negative values in nutrient fields
- `MIN_NUT_LIM` (define, BLING_OPTIONS.h) — Use Liebig function instead of geometric mean of the nutrient limitations to calculate maximum phyto growth rate
- `SIZE_DEP_LIM` (undef, BLING_OPTIONS.h) — Allow different phytoplankton groups to have different growth rates and nutrient/light limitations. Parameters implemented have yet to be tuned
- `ML_MEAN_LIGHT` (undef, BLING_OPTIONS.h) — Assume that phytoplankton in the mixed layer experience the average light over the mixed layer (as in original BLING model)
- `ML_MEAN_PHYTO` (define, BLING_OPTIONS.h) — Assume that phytoplankton are homogenized in the mixed layer
- `BLING_USE_THRESHOLD_MLD` (undef, BLING_OPTIONS.h) — Calculate MLD using a threshold criterion. If undefined, MLD is calculated using the second derivative of rho(z)
- `USE_QSW` (undef, BLING_OPTIONS.h) — Determine PAR from shortwave radiation Qsw; otherwise determined from date and latitude (Do not define if not using pkg/exf)
- `PHYTO_SELF_SHADING` (undef, BLING_OPTIONS.h) — Light absorption scheme from Manizza et al. (2005), with self shading from phytoplankton
- `BLING_ADJOINT_SAFE` (define, BLING_OPTIONS.h) — Simplify some parts of the code that are problematic when using the adjoint
- `USE_BLING_DVM` (undef, BLING_OPTIONS.h) — For adjoint safe, do not call bling_dvm
- `CARBONCHEM_SOLVESAPHE` (undef, BLING_OPTIONS.h) — Compile "Solvesaphe" package (Munhoven 2013) for pH/pCO2 can still select Follows et al (2006) solver in data.bling, but will use solvesaphe dissociation coefficient options
- `CARBONCHEM_TOTALPHSCALE` (undef, BLING_OPTIONS.h) — consistent with other coefficients (currently on the seawater scale). NOTE: Has NO effect when CARBONCHEM_SOLVESAPHE is defined (different coeffs are used).
- `NEW_FRAC_EXP` (define, BLING_OPTIONS.h) — When calculating the fraction of sinking organic matter, use model biomass diagnostics.

## Headers
- `BLING_LOAD.h` — --   COMMON /BLING_LOAD/ BLING_ldRec     :: time-record currently loaded (in temp arrays *[1])
- `BLING_OPTIONS.h` — Package-specific Options & Macros go here
- `BLING_VARS.h` — Carbon chemistry variables
- `bling_ad_check_lev1_dir.h` — ADJ STORE pH                = comlev1, key = ikey_dynamics, kind=isbyte ADJ STORE fice              = comlev1, key = ikey_dynamics, kind=isbyte ADJ ST
- `bling_ad_check_lev2_dir.h` — ADJ STORE pH                = tapelev2, key = ilev_2 CADJ STORE fice              = tapelev2, key = ilev_2 CADJ STORE atmosP            = tapelev2, ke
- `bling_ad_check_lev3_dir.h` — ADJ STORE pH                = tapelev3, key = ilev_3 CADJ STORE fice              = tapelev3, key = ilev_3 CADJ STORE atmosP            = tapelev3, ke
- `bling_ad_check_lev4_dir.h` — ADJ STORE pH                = tapelev4, key = ilev_4 CADJ STORE fice              = tapelev4, key = ilev_4 CADJ STORE atmosP            = tapelev4, ke

## Routines (19)
`bling_airseaflux.F`, `bling_bio.F`, `bling_bio_nitrogen.F`, `bling_carbonate_init.F`, `bling_carbonate_sys.F`, `bling_diagnostics_init.F`, `bling_fields_load.F`, `bling_ini_forcing.F`, `bling_init_fixed.F`, `bling_init_varia.F`, `bling_light.F`, `bling_main.F`, `bling_min_val.F`, `bling_mixedlayer.F`, `bling_read_pickup.F`, `bling_readparms.F`, `bling_sgs.F`, `bling_tr_register.F`, `bling_write_pickup.F`

## Called from outside the package
- `BLING_FIELDS_LOAD` ← `pkg/gchem/gchem_fields_load.F:41`
- `BLING_MAIN` ← `pkg/gchem/gchem_forcing_sep.F:187,202,217`
- `BLING_INIT_FIXED` ← `pkg/gchem/gchem_init_fixed.F:49`
- `BLING_CARBONATE_INIT` ← `pkg/gchem/gchem_init_vari.F:81`
- `BLING_INIT_VARIA` ← `pkg/gchem/gchem_init_vari.F:79`
- `BLING_INI_FORCING` ← `pkg/gchem/gchem_init_vari.F:80`
- `BLING_READPARMS` ← `pkg/gchem/gchem_readparms.F:164`
- `BLING_TR_REGISTER` ← `pkg/gchem/gchem_tr_register.F:60`
- `BLING_WRITE_PICKUP` ← `pkg/gchem/gchem_write_pickup.F:50`

## Verification experiments compiling it (1)
`global_oce_biogeo_bling`
