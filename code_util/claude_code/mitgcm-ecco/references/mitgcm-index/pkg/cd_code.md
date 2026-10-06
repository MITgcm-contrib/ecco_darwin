# pkg/cd_code

C-D scheme: carries D-grid velocities to stabilise Coriolis on the C-grid at coarse resolution (useCDscheme).

**pkg_depend:** +mom_common  (`+` requires, `-` excludes)
**runtime switch:** `useCD_CODE`-style flag in `data.pkg` (check exact name in packages_boot.F)
**manual:** `doc/algorithm/algorithm.rst`, `doc/getting_started/getting_started.rst`
**adjoint support files:** cd_code_ad_check_lev1_dir.h, cd_code_ad_check_lev2_dir.h, cd_code_ad_check_lev3_dir.h, cd_code_ad_check_lev4_dir.h, cd_code_ad_diff.list

## CPP options (defaults as shipped)
- `CD_CODE_NO_AB_MOMENTUM` (undef, CD_CODE_OPTIONS.h) — with Adams-Bashforth only on surface pressure term. Tests show that using AB on D-grid coriolis term improves stability (as expected from CD-scheme paper). The following 2 options allow to reproduce old results.
- `CD_CODE_NO_AB_CORIOLIS` (undef, CD_CODE_OPTIONS.h)

## Headers
- `CD_CODE_OPTIONS.h` — CPP options file for CD_CODE package Use this file for selecting CPP options within the cd_code package
- `CD_CODE_VARS.h` — uVelD  :: D grid zonal velocity vVelD  :: D grid meridional velocity
- `cd_code_ad_check_lev1_dir.h` — ADJ STORE uveld     = comlev1, key = ikey_dynamics, kind = isbyte ADJ STORE vveld     = comlev1, key = ikey_dynamics, kind = isbyte ADJ STORE etanm1  
- `cd_code_ad_check_lev2_dir.h` — ADJ STORE uveld      = tapelev2, key = ilev_2 ADJ STORE vveld     = tapelev2, key = ilev_2 ADJ STORE etanm1    = tapelev2, key = ilev_2 ADJ STORE unm1
- `cd_code_ad_check_lev3_dir.h` — ADJ STORE uveld      = tapelev3, key = ilev_3 ADJ STORE vveld     = tapelev3, key = ilev_3 ADJ STORE etanm1    = tapelev3, key = ilev_3 ADJ STORE unm1
- `cd_code_ad_check_lev4_dir.h` — ADJ STORE uveld     = tapelev4, key = ilev_4 ADJ STORE vveld     = tapelev4, key = ilev_4 ADJ STORE etanm1    = tapelev4, key = ilev_4 ADJ STORE unm1 

## Routines (5)
`cd_code_ini_vars.F`, `cd_code_init_fixed.F`, `cd_code_read_pickup.F`, `cd_code_scheme.F`, `cd_code_write_pickup.F`

## Called from outside the package
- `CD_CODE_INIT_FIXED` ← `model/src/packages_init_fixed.F:213`
- `CD_CODE_INI_VARS` ← `model/src/packages_init_variables.F:208`
- `CD_CODE_WRITE_PICKUP` ← `model/src/packages_write_pickup.F:104`
- `CD_CODE_SCHEME` ← `model/src/timestep.F:232`

## Verification experiments compiling it (14)
`bottom_ctrl_5x5` `cfc_example` `exp2` `global_oce_biogeo_bling` `global_oce_latlon` `global_ocean.90x40x15` `ideal_2D_oce` `isomip` `lab_sea` `so_box_biogeo` `tutorial_global_oce_biogeo` `tutorial_global_oce_latlon` `tutorial_global_oce_optim` `tutorial_tracer_adjsens`
