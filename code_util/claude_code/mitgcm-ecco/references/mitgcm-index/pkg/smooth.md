# pkg/smooth

Diffusion-operator smoothing for control/covariance (2-D/3-D correlation operators).

**runtime switch:** `useSMOOTH`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.smooth`
**manual:** `doc/ocean_state_est/ocean_state_est.rst`
**adjoint support files:** smooth_ad_diff.list

## Namelist parameters
### SMOOTH_NML
- `smooth2Dnbt`
- `smooth2Dtype`
- `smooth2Dsize`
- `smooth2D_Lx0`
- `smooth2D_Ly0`
- `smooth2Dfilter`
- `smooth2DmaskName`
- `smooth3Dnbt`
- `smooth3DtypeH`
- `smooth3DsizeH`
- `smooth3DtypeZ`
- `smooth3DsizeZ`
- `smooth3D_Lx0`
- `smooth3D_Ly0`
- `smooth3D_Lz0`
- `smooth3Dfilter`
- `smooth3DmaskName`
- `smoothDir`

## Headers
- `SMOOTH.h` — pkg/smooth constants
- `SMOOTH_OPTIONS.h` — CPP options file for SMOOTH package Use this file for selecting options within the SMOOTH package

## Routines (19)
`smooth2d.F`, `smooth3d.F`, `smooth_basic2d.F`, `smooth_check.F`, `smooth_correl2d.F`, `smooth_correl2dw.F`, `smooth_correl3d.F`, `smooth_diff2d.F`, `smooth_diff3d.F`, `smooth_filtervar2d.F`, `smooth_filtervar3d.F`, `smooth_hetero2d.F`, `smooth_impldiff.F`, `smooth_init2d.F`, `smooth_init3d.F`, `smooth_init_fixed.F`, `smooth_init_varia.F`, `smooth_readparms.F`, `smooth_rhs.F`

## Called from outside the package
- `SMOOTH_CHECK` ← `model/src/packages_check.F:360`
- `SMOOTH_INIT_FIXED` ← `model/src/packages_init_fixed.F:487`
- `SMOOTH_INIT_VARIA` ← `model/src/packages_init_variables.F:569`
- `SMOOTH_READPARMS` ← `model/src/packages_readparms.F:338`
- `SMOOTH2D` ← `pkg/ctrl/ctrl_get_gen.F:131,160`
- `SMOOTH2D` ← `pkg/ctrl/ctrl_map_genarr.F:139`
- `SMOOTH3D` ← `pkg/ctrl/ctrl_map_genarr.F:338`
- `SMOOTH_CORREL2D` ← `pkg/ctrl/ctrl_map_genarr.F:138`
- `SMOOTH_CORREL3D` ← `pkg/ctrl/ctrl_map_genarr.F:337`
- `SMOOTH2D` ← `pkg/ctrl/ctrl_map_ini_gentim2d.F:449`
- `SMOOTH_CORREL2D` ← `pkg/ctrl/ctrl_map_ini_gentim2d.F:448`
- `SMOOTH_BASIC2D` ← `pkg/ecco/cost_gencost_bpv4.F:315,333`
- `SMOOTH_HETERO2D` ← `pkg/ecco/cost_gencost_bpv4.F:311,329`
- `SMOOTH_HETERO2D` ← `pkg/ecco/cost_gencost_sshv4.F:580,947,950,1283`
- `SMOOTH_HETERO2D` ← `pkg/ecco/cost_gencost_sstv4.F:300,306`
- `SMOOTH_HETERO2D` ← `pkg/ecco/cost_generic.F:427`
