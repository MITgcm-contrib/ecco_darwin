# pkg/ebm

Energy-balance atmosphere model coupled to the ocean.

**runtime switch:** `useEBM`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.ebm`
**manual:** `doc/examples/examples.rst`
**adjoint support files:** ebm_ad_check_lev1_dir.h, ebm_ad_check_lev2_dir.h, ebm_ad_check_lev3_dir.h, ebm_ad_diff.list

## Namelist parameters
### EBM_PARM01
- `tauThetaZonRelax` — time-scale [s] for relaxation towards Zon.Aver. SST
- `scale_runoff`
- `RunoffFile` — Runoff input file name

## CPP options (defaults as shipped)
- `EBM_WIND_PERT` (undef, EBM_OPTIONS.h)
- `EBM_CLIMATE_CHANGE` (undef, EBM_OPTIONS.h)
- `EBM_VERSION_1BASIN` (undef, EBM_OPTIONS.h)

## Headers
- `EBM.h` — BOP
- `EBM_OPTIONS.h` — CPP options file for EBM package Use this file for selecting CPP options within the EBM package
- `ebm_ad_check_lev1_dir.h` — ADJ STORE zonalmeansst = comlev1, key = ikey_dynamics
- `ebm_ad_check_lev2_dir.h` — ADJ STORE zonalmeansst = tapelev2, key = ilev_2
- `ebm_ad_check_lev3_dir.h` — ADJ STORE zonalmeansst = tapelev3, key = ilev_3

## Routines (8)
`ebm_area_t.F`, `ebm_atmosphere.F`, `ebm_driver.F`, `ebm_ini_vars.F`, `ebm_load_climatology.F`, `ebm_readparms.F`, `ebm_wind_perturb.F`, `ebm_zonalmean.F`

## Called from outside the package
- `EBM_DRIVER` ← `model/src/forward_step.F:608`
- `EBM_INI_VARS` ← `model/src/packages_init_variables.F:310`
- `EBM_READPARMS` ← `model/src/packages_readparms.F:231`

## Verification experiments compiling it (1)
`global_oce_latlon`
