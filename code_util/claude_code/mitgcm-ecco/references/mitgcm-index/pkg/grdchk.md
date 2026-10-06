# pkg/grdchk

Gradient check: compares adjoint gradient against finite differences (data.grdchk).

**pkg_depend:** +autodiff +cost +ctrl  (`+` requires, `-` excludes)
**in groups:** adjoint
**runtime switch:** `useGRDCHK`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.grdchk`
**manual:** `doc/autodiff/autodiff.rst`

## Namelist parameters
### GRDCHK_NML
- `grdchk_eps`
- `nbeg`
- `nstep`
- `nend`
- `grdchkvarname`
- `grdchkvarindex`
- `useCentralDiff`
- `grdchkwhichproc`
- `iGloPos`
- `jGloPos`
- `kGloPos`
- `iGloTile`
- `jGloTile`
- `idep`
- `jdep`
- `obcsglo`
- `recglo`

## Headers
- `GRDCHK.h` — HEADER GRADIENT_CHECK Header for doing gradient checks with the ECCO ocean state estimation tool. started: Christian Eckert eckert@mit.edu  01-Mar-200
- `GRDCHK_OPTIONS.h` — BOP

## Routines (14)
`grdchk_check.F`, `grdchk_ctrl_fname.F`, `grdchk_get_mask.F`, `grdchk_get_obcs_mask.F`, `grdchk_get_position.F`, `grdchk_getadxx.F`, `grdchk_getxx.F`, `grdchk_init.F`, `grdchk_loc.F`, `grdchk_main.F`, `grdchk_print.F`, `grdchk_readparms.F`, `grdchk_setxx.F`, `grdchk_summary.F`

## Called from outside the package
- `GRDCHK_CHECK` ← `model/src/packages_check.F:445`
- `GRDCHK_READPARMS` ← `model/src/packages_readparms.F:333`
- `GRDCHK_MAIN` ← `model/src/the_model_main.F:735`
- `GRDCHK_MAIN` ← `pkg/openad/the_model_main.F:309`

## Verification experiments compiling it (16)
`1D_ocean_ice_column` `bottom_ctrl_5x5` `global_oce_biogeo_bling` `global_oce_latlon` `global_ocean.90x40x15` `global_ocean.cs32x15` `halfpipe_streamice` `hs94.1x64x5` `isomip` `lab_sea` `obcs_ctrl` `offline_exf_seaice` `tutorial_dic_adjoffline` `tutorial_global_oce_biogeo` `tutorial_global_oce_optim` `tutorial_tracer_adjsens`
