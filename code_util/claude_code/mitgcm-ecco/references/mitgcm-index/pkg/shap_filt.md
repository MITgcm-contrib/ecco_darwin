# pkg/shap_filt

Shapiro filter for grid-scale noise.

**pkg_depend:** +mom_vecinv  (`+` requires, `-` excludes)
**in groups:** atmospheric
**runtime switch:** `useSHAP_FILT`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.shap`
**manual:** `doc/phys_pkgs/shap_filt.rst`, `doc/algorithm/algorithm.rst`, `doc/examples/held_suarez_cs/held_suarez_cs.rst`, `doc/outp_pkgs/outp_pkgs.rst`
**adjoint support files:** shap_filt_ad_diff.list

## Namelist parameters
### SHAP_PARM01
- `Shap_funct` — define which Shapiro Filter function is used = 1  (S1) : [1 - d_xx^n - d_yy^n] = 4  (S4) : [1 - d_xx^n][1- d_yy^n] = 2  (S2) : [1 - (d_xx+d_yy)^n]
- `shap_filt_uvStar` — filter applied to u*,v* (before SOLVE_FOR_P)
- `shap_filt_TrStagg` — if using a Stagger time-step, filter T,S before computing PhiHyd ; has no effect if syncr. time step is used
- `Shap_alwaysExchUV` — always call exch(U,V)    nShapUV times
- `Shap_alwaysExchTr` — always call exch(Tracer) nShapTr times
- `nShapT`
- `nShapS`
- `nShapTrPhys`
- `Shap_Trtau` — Time scale for tracer filter
- `Shap_TrLength` — Length scale for tracer filter
- `nShapUV`
- `nShapUVPhys`
- `Shap_uvtau` — Time scale for momentum filter
- `Shap_uvLength` — Length scale for momentum filter
- `Shap_noSlip` — No-slip parameter (=0 free sleep ; =1 No-slip)
- `Shap_diagFreq` — Frequency^-1 for diagnostic output (s)

## CPP options (defaults as shipped)
- `SEQUENTIAL_2D_SHAP` (define, SHAP_FILT_OPTIONS.h) — Use [1-d_yy^n)(1-d_xx^n] instead of [1-d_xx^n-d_yy^n] This changes the spectral response function dramatically. You need to do some analysis before changing this option. ;^)
- `USE_OLD_SHAPIRO_FILTERS` (undef, SHAP_FILT_OPTIONS.h) — overlap to be sufficiently wide and also does not work for arbitrarily arranged tiles (ie. as in cubed-sphere). *DO NOT USE THIS OPTION UNLESS YOU REALLY WANT TO*  :-(
- `NO_SLIP_SHAP` (undef, SHAP_FILT_OPTIONS.h) — Horizontal shear is calculated as if the boundaries are no-slip Note: option NO_SLIP_SHAP only used in OLD_SHAPIRO_FILTERS ; it is replaced by parameter "Shap_noSlip=1." in new S/R.
- `USE_SHAP_CALC_VORTICITY` (undef, SHAP_FILT_OPTIONS.h) — can be different from the masking applied for momentum advection; This option allows to use the local S/R: SHAP_FILT_RELVORT3 to compute vorticity, instead of pkg/mom_common S/R: MOM_CALC_RELVORT3

## Headers
- `SHAP_FILT.h` — -    Package flag and logical parameters : shap_filt_uvStar  :: filter applied to u*,v* (before SOLVE_FOR_P) shap_filt_TrStagg :: if using a Stagger t
- `SHAP_FILT_OPTIONS.h` — CPP options file for pkg SHAP_FILT

## Routines (17)
`shap_filt_apply_ts.F`, `shap_filt_apply_uv.F`, `shap_filt_computvort.F`, `shap_filt_diagnostics_init.F`, `shap_filt_init_fixed.F`, `shap_filt_readparms.F`, `shap_filt_relvort3.F`, `shap_filt_tracer_s1.F`, `shap_filt_tracer_s2.F`, `shap_filt_tracer_s4.F`, `shap_filt_tracerold.F`, `shap_filt_u.F`, `shap_filt_uv_s1.F`, `shap_filt_uv_s2.F`, `shap_filt_uv_s2c.F`, `shap_filt_uv_s4.F`, `shap_filt_v.F`

## Called from outside the package
- `SHAP_FILT_APPLY_UV` ← `model/src/forward_step.F:882`
- `SHAP_FILT_APPLY_UV` ← `model/src/momentum_correction_step.F:110`
- `SHAP_FILT_INIT_FIXED` ← `model/src/packages_init_fixed.F:233`
- `SHAP_FILT_READPARMS` ← `model/src/packages_readparms.F:171`
- `SHAP_FILT_APPLY_TS` ← `model/src/tracers_correction_step.F:73`

## Verification experiments compiling it (12)
`aim.5l_Equatorial_Channel` `aim.5l_LatLon` `aim.5l_cs` `atm_gray` `cpl_aim+ocn` `fizhi-cs-32x32x40` `fizhi-cs-aqualev20` `fizhi-gridalt-hs` `hs94.128x64x5` `hs94.1x64x5` `hs94.cs-32x32x5` `tutorial_held_suarez_cs`
