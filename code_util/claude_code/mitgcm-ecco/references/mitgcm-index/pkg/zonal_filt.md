# pkg/zonal_filt

Zonal (polar) Fourier filter for lat-lon grids.

**runtime switch:** `useZONAL_FILT`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.zonfilt`
**manual:** `doc/algorithm/algorithm.rst`

## Namelist parameters
### ZONFILT_PARM01
- `zonal_filt_uvStar` — filter applied to u*,v* (before SOLVE_FOR_P)
- `zonal_filt_TrStagg` — if using a Stager time-step, filter T,S before computing PhiHyd ; has no effect if syncr. time step is used
- `zonal_filt_lat` — Low latitude for FFT filtering of latitude circles
- `zonal_filt_cospow` — Latitude dependance of the damping function = ( cos Lat / cos zonal_filt_lat )**cospow
- `zonal_filt_sinpow` — zonal mode dependance of the damping function = 1 / ( sin pi.kx/Nx )**sinpow
- `zonal_filt_mode2dx` — to specify how to treat the 2.dx mode : = 0 : damped like other modes. = 1 : removed in regions where Zonal_filt apply = 2 : removed every where.

## Headers
- `FFTPACK.h` — Data structures/work-space for FFTPACK
- `ZONAL_FILT.h` — -    Package flag and logical parameters : zonal_filt_uvStar  :: filter applied to u*,v* (before SOLVE_FOR_P) zonal_filt_TrStagg :: if using a Stager 
- `ZONAL_FILT_OPTIONS.h` — Package-specific Options & Macros go here

## Routines (34)
`fftpack.F`, `zonal_filt_apply_ts.F`, `zonal_filt_apply_uv.F`, `zonal_filt_init.F`, `zonal_filt_nofill.F`, `zonal_filt_postsmooth.F`, `zonal_filt_presmooth.F`, `zonal_filt_readparms.F`, `zonal_filter.F`

## Called from outside the package
- `ZONAL_FILT_APPLY_UV` ← `model/src/forward_step.F:889`
- `ZONAL_FILT_APPLY_UV` ← `model/src/momentum_correction_step.F:119`
- `ZONAL_FILT_INIT` ← `model/src/packages_init_fixed.F:243`
- `ZONAL_FILT_READPARMS` ← `model/src/packages_readparms.F:176`
- `ZONAL_FILT_APPLY_TS` ← `model/src/tracers_correction_step.F:80`
- `ZONAL_FILTER` ← `pkg/ptracers/ptracers_zonal_filt_apply.F:44`

## Verification experiments compiling it (2)
`aim.5l_LatLon` `hs94.128x64x5`
