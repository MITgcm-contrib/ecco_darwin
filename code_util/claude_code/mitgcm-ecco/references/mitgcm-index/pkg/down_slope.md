# pkg/down_slope

Down-slope flow parameterisation (Campin & Goosse 1999) moving dense bottom water downslope (excludes bbl).

**runtime switch:** `useDOWN_SLOPE`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.down_slope`
**manual:** `doc/examples/examples.rst`
**adjoint support files:** dwnslp_ad_diff.list

## Namelist parameters
### DWNSLP_PARM01
- `DWNSLP_slope` — fixed slope (=0 => use the local slope)
- `DWNSLP_rec_mu` — reciprol friction parameter (unit = time scale [s]) used to compute the flow: U=dy*dz*(slope * g/mu * dRho / rho0)
- `DWNSLP_drFlow` — max. thickness [m] of the effective downsloping flow layer
- `temp_useDWNSLP` — true if Down-Sloping flow applies to temperature
- `salt_useDWNSLP` — true if Down-Sloping flow applies to salinity

## Headers
- `DWNSLP_OPTIONS.h` — CPP options file for Down-Slope package Use this file for selecting options within the Down-Slope package
- `DWNSLP_PARAMS.h` — -    Package flag and logical parameters : temp_useDWNSLP  :: true if Down-Sloping flow applies to temperature salt_useDWNSLP  :: true if Down-Sloping
- `DWNSLP_SIZE.h` — #ifdef ALLOW_DOWN_SLOPE
- `DWNSLP_VARS.h` — store the location of potential site where Down-Sloping Flow is applied DWNSLP_NbSite :: Number of bathymetry steps within each tile DWNSLP_ijDeep :: 

## Routines (7)
`dwnslp_apply.F`, `dwnslp_calc_flow.F`, `dwnslp_calc_rho.F`, `dwnslp_diagnostics_init.F`, `dwnslp_init_fixed.F`, `dwnslp_init_varia.F`, `dwnslp_readparms.F`

## Called from outside the package
- `DWNSLP_CALC_FLOW` ← `model/src/do_oceanic_phys.F:1052,1056`
- `DWNSLP_CALC_RHO` ← `model/src/do_oceanic_phys.F:736`
- `DWNSLP_INIT_FIXED` ← `model/src/packages_init_fixed.F:351`
- `DWNSLP_INIT_VARIA` ← `model/src/packages_init_variables.F:286`
- `DWNSLP_READPARMS` ← `model/src/packages_readparms.F:221`
- `DWNSLP_APPLY` ← `model/src/salt_integrate.F:449,456`
- `DWNSLP_APPLY` ← `model/src/temp_integrate.F:451,458`
- `DWNSLP_APPLY` ← `pkg/ptracers/ptracers_integrate.F:410,417`

## Verification experiments compiling it (2)
`global_ocean.90x40x15` `lab_sea`
