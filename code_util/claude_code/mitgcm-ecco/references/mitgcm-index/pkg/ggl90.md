# pkg/ggl90

Gaspar-Gregoris-Lefevre (1990) TKE vertical mixing scheme (ECCO v4 default), with IDEMIX and Langmuir options.

**runtime switch:** `useGGL90`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.ggl90`
**manual:** `doc/phys_pkgs/ggl90.rst`, `doc/examples/examples.rst`
**adjoint support files:** ggl90_ad_check_lev1_dir.h, ggl90_ad_check_lev2_dir.h, ggl90_ad_check_lev3_dir.h, ggl90_ad_check_lev4_dir.h, ggl90_ad_diff.list

## Namelist parameters
### GGL90_PARM01
- `GGL90dumpFreq`
- `GGL90taveFreq`
- `GGL90diffTKEh`
- `GGL90mixingMaps`
- `GGL90writeState`
- `GGL90ck`
- `GGL90ceps`
- `GGL90alpha`
- `GGL90m2`
- `GGL90TKEmin`
- `GGL90TKEsurfMin`
- `GGL90TKEbottom`
- `GGL90mixingLengthMin`
- `mxlMaxFlag`
- `adMxlMaxFlag`
- `mxlSurfFlag`
- `GGL90viscMax`
- `GGL90diffMax`
- `GGL90TKEFile`
- `GGL90_dirichlet`
- `calcMeanVertShear` — calculate the mean (@ grid-cell center) of vertical shear compon. (instead of vert. shear of mean flow); also applies to surface stress (uStarSquare)
- `useIDEMIX`
- `useLANGMUIR`
### GGL90_PARM02
- `IDEMIX_tau_v` — time scale for vertical symmetrisation (s), def: 1d  _[ifdef ALLOW_GGL90_IDEMIX]_
- `IDEMIX_tau_h` — time scale for horizontal symmetrisation (s), def: 10d  _[ifdef ALLOW_GGL90_IDEMIX]_
- `IDEMIX_gamma` — scaling factor (see Olbers and Eden 2013, App) def: 1.57  _[ifdef ALLOW_GGL90_IDEMIX]_
- `IDEMIX_jstar` — spectral bandwidth in modes (default=10.)  _[ifdef ALLOW_GGL90_IDEMIX]_
- `IDEMIX_mu0` — dissipation parameter (wrong default=4/3, should be 1/3)  _[ifdef ALLOW_GGL90_IDEMIX]_
- `IDEMIX_tidal_file` — file containing tidal forcing  _[ifdef ALLOW_GGL90_IDEMIX]_
- `IDEMIX_wind_file` — file containing surface wind forcing  _[ifdef ALLOW_GGL90_IDEMIX]_
- `IDEMIX_mixing_efficiency` — used for diagnosing Osborn diff (def: 0.1666)  _[ifdef ALLOW_GGL90_IDEMIX]_
- `IDEMIX_diff_max` — maximum Osborn diffusivity (def: 1)  _[ifdef ALLOW_GGL90_IDEMIX]_
- `IDEMIX_diff_min` — minimum diffusivity, not used (def: 1e-9)  _[ifdef ALLOW_GGL90_IDEMIX]_
- `IDEMIX_frac_F_b` — scaling factor for  bottom forcing (def: 1)  _[ifdef ALLOW_GGL90_IDEMIX]_
- `IDEMIX_frac_F_s` — scaling factor for surface forcing (def: 0.2, because only 20% of the surface wind energy input is assumed to reach the interior below the mixed layer.)  _[ifdef ALLOW_GGL90_IDEMIX]_
- `IDEMIX_include_GM` — include eddy contribution parameterized by GMredi  _[ifdef ALLOW_GGL90_IDEMIX]_
- `IDEMIX_include_GM_bottom` — include eddy contribution only at the bottom  _[ifdef ALLOW_GGL90_IDEMIX]_
### GGL90_PARM03
- `LC_Gamma` — mixing-length Amplification factor from Langmuir Circ.  _[ifdef ALLOW_GGL90_LANGMUIR]_
- `LC_num` — Value for the Langmuir number (no unit)  _[ifdef ALLOW_GGL90_LANGMUIR]_
- `LC_lambda` — vertical scale for Stokes velocity profile ( m )  _[ifdef ALLOW_GGL90_LANGMUIR]_

## CPP options (defaults as shipped)
- `ALLOW_GGL90_HORIZDIFF` (undef, GGL90_OPTIONS.h) — Enable horizontal diffusion of TKE.
- `ALLOW_GGL90_SMOOTH` (undef, GGL90_OPTIONS.h) — Use horizontal averaging for viscosity and diffusivity as originally implemented in OPA.
- `ALLOW_GGL90_IDEMIX` (undef, GGL90_OPTIONS.h) — allow IDEMIX model
- `GGL90_IDEMIX_CVMIX_VERSION` (define, GGL90_OPTIONS.h) — The cvmix version of idemix uses different regularisations for the Coriolis parameter, buoyancy frequency etc, when used in the denominator
- `ALLOW_GGL90_LANGMUIR` (undef, GGL90_OPTIONS.h) — include Langmuir circulation parameterization
- `GGL90_REGULARIZE_MIXINGLENGTH` (undef, GGL90_OPTIONS.h) — Replace MAX(mxl,mxlMin) by SQRT(mxl**2+mxlMin**2) to help adjoint
- `GGL90_MISSING_HFAC_BUG` (undef, GGL90_OPTIONS.h) — recover old bug prior to Jun 2023

## Headers
- `GGL90.h` — BOP
- `GGL90_OPTIONS.h` — CPP options file for GGL90 package. Use this file for selecting options within the GGL90 package.
- `ggl90_ad_check_lev1_dir.h` — ADJ STORE GGL90viscArU       = comlev1, key=ikey_dynamics ADJ STORE GGL90viscArV       = comlev1, key=ikey_dynamics ADJ STORE GGL90diffKr        = com
- `ggl90_ad_check_lev2_dir.h` — ADJ STORE GGL90TKE           = tapelev2, key=ilev_2
- `ggl90_ad_check_lev3_dir.h` — ADJ STORE GGL90TKE           = tapelev3, key=ilev_3
- `ggl90_ad_check_lev4_dir.h` — ADJ STORE GGL90TKE           = tapelev4, key=ilev_4

## Routines (17)
`ggl90_add_stokesdrift.F`, `ggl90_calc.F`, `ggl90_calc_diff.F`, `ggl90_calc_visc.F`, `ggl90_check.F`, `ggl90_diagnostics_init.F`, `ggl90_exchanges.F`, `ggl90_idemix.F`, `ggl90_init_fixed.F`, `ggl90_init_varia.F`, `ggl90_mixinglength.F`, `ggl90_output.F`, `ggl90_read_pickup.F`, `ggl90_readparms.F`, `ggl90_write_pickup.F`

## Called from outside the package
- `GGL90_CALC_DIFF` ← `model/src/calc_3d_diffusivity.F:240`
- `GGL90_CALC_VISC` ← `model/src/calc_viscosity.F:114`
- `GGL90_CALC` ← `model/src/do_oceanic_phys.F:1010`
- `GGL90_EXCHANGES` ← `model/src/do_oceanic_phys.F:1109`
- `GGL90_OUTPUT` ← `model/src/do_the_model_io.F:175`
- `GGL90_CHECK` ← `model/src/packages_check.F:240`
- `GGL90_INIT_FIXED` ← `model/src/packages_init_fixed.F:311`
- `GGL90_INIT_VARIA` ← `model/src/packages_init_variables.F:245`
- `GGL90_READPARMS` ← `model/src/packages_readparms.F:201`
- `GGL90_WRITE_PICKUP` ← `model/src/packages_write_pickup.F:118`
- `GGL90_ADD_STOKESDRIFT` ← `pkg/mom_fluxform/mom_fluxform.F:1085`
- `GGL90_ADD_STOKESDRIFT` ← `pkg/mom_vecinv/mom_vecinv.F:692`

## Verification experiments compiling it (5)
`global_oce_latlon` `global_ocean.90x40x15` `global_ocean.cs32x15` `isomip` `vermix`
