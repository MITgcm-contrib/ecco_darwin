# pkg/longstep

Longer time step for passive tracers than for dynamics (ptracers stepped every N steps with averaged velocities).

**runtime switch:** `useLONGSTEP`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.longstep`
**manual:** `doc/examples/examples.rst`

## README
```
Package longstep
================

This package allows the passive tracer time step to be longer than that for
dynamical fields: the ptracers are updated only every LS_nIter time step. 
Dynamical fields are averaged over LS_nIter time steps and are available as
fields LS_* (declared in LONGSTEP.h):

original fld.  averaged fld.
------------------------------               
UVEL           LS_uVel
VVEL           LS_vVel
WVEL           LS_wVel
THETA          LS_theta
SALT           LS_salt
IVDConvCount   LS_IVDConvCount
Qsw            LS_Qsw
               
Kwx            LS_Kwx
Kwy            LS_Kwy
Kwz            LS_Kwz
               
KPPdiffKzS     LS_KPPdiffKzS
KPPghat        LS_KPPghat

The T and S time step remains the same as that for u,v,...


Packages that use ptracers (like DIC) need to be adapted:

```

## Namelist parameters
### LONGSTEP_PARM01
- `LS_nIter` — number of dynamics time steps between ptracer steps
- `LS_whenToSample` — when to sample dynamical fields for the longstep average 0 - at beginning of timestep (reproduces offline results) 1 - after first THERMODYNAMICS but before DYNAMICS (use use old U,V,W for advection, but new T,S for GCHEM if

## Headers
- `LONGSTEP.h` — BOP
- `LONGSTEP_OPTIONS.h` — Package-specific options go here
- `LONGSTEP_PARAMS.h` — BOP

## Routines (17)
`longstep_average.F`, `longstep_average_3d.F`, `longstep_average_3d_fac.F`, `longstep_check.F`, `longstep_check_iters.F`, `longstep_correction_step.F`, `longstep_diagnostics_init.F`, `longstep_fill_3d.F`, `longstep_fill_3d_fac.F`, `longstep_fill_3d_rs.F`, `longstep_forcing_surf.F`, `longstep_init_fixed.F`, `longstep_init_varia.F`, `longstep_readparms.F`, `longstep_reset_3d.F`, `longstep_residual_flow.F`, `longstep_thermodynamics.F`

## Called from outside the package
- `LONGSTEP_AVERAGE` ← `model/src/forward_step.F:710,750,1040`
- `LONGSTEP_THERMODYNAMICS` ← `model/src/forward_step.F:719,759,1049`
- `LONGSTEP_CHECK` ← `model/src/packages_check.F:297`
- `LONGSTEP_INIT_FIXED` ← `model/src/packages_init_fixed.F:420`
- `LONGSTEP_INIT_VARIA` ← `model/src/packages_init_variables.F:339`
- `LONGSTEP_READPARMS` ← `model/src/packages_readparms.F:247`
- `LONGSTEP_CHECK_ITERS` ← `model/src/set_parms.F:333`

## Verification experiments compiling it (1)
`lab_sea`
