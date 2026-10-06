# pkg/flt

Lagrangian floats/particles advected online (2-D/3-D, profiling floats).

**pkg_depend:** +mdsio  (`+` requires, `-` excludes)
**runtime switch:** `useFLT`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.flt`
**manual:** `doc/outp_pkgs/flt.rst`, `doc/examples/examples.rst`, `doc/outp_pkgs/outp_pkgs.rst`

## README
```
c
c     ==============================
c     FLOAT Package for the MITgcmUV
c     ==============================
c
c
c     This package allows the advection of floats during a model run.
c     Although originally intended to simulate PALACE floats
c     (floats that drift in at depth and to the surface at a defined
c     time interval) it can also run ALACE floats (non-profiling) 
c     and surface drifters as well as sample moorings (simply a 
c     non-advective, profiling float).
c
c     The stepping of the float advection is done using a second
c     order Runga-Kutta scheme (Press et al., 1992, Numerical
c     Recipes), whereby velocities and positions are bilinear
c     interpolated between the grid points.
c
c     Current version: 1.0                              06-AUG-2001
c
c     Please report any bugs and suggestions to:
c            
c            Arne Biastoch (abiastoch@ucsd.edu)
c
c     
c     Implementation in MITgcmUV
c     --------------------------
c
c     The package has only few interfaces to the model. Despite a 
c     general introduction of the flag useFLT and an initialization in 
```

## Namelist parameters
### FLT_NML
- `flt_int_traj` — period between storing model state at float position, in s
- `flt_int_prof` — period between float vertical profiles, in s
- `flt_selectTrajOutp` — select which var. to output along trajectories
- `flt_selectProfOutp` — select which var. to output along profiles =0 : none ; =1 : position only ; =2 : +p,u,v,t,s
- `flt_noise` — range of noise added to the velocity component (randomly). The noise can be added or subtracted, the range is +/- flt_noise/2
- `flt_deltaT` — time-step to step forward floats (in flt_runga2.F) default is deltaTClock
- `FLT_Iter0` — timestep number when float are initialized
- `flt_file` — name of the file containing the initial positions. At initialization the program first looks for a global file flt_file.data. If that is not found it looks for tiled files flt_file.iG.jG.data.
- `mapIniPos2Index` — convert float initial position to (local) index map

## CPP options (defaults as shipped)
- `ALLOW_3D_FLT` (define, FLT_OPTIONS.h) — Include/Exclude part that allows 3-dimensional advection of floats
- `USE_FLT_ALT_NOISE` (define, FLT_OPTIONS.h) — Use the alternative method of adding random noise to float advection
- `ALLOW_FLT_3D_NOISE` (define, FLT_OPTIONS.h) — Add noise also to the vertical velocity of 3D floats
- `FLT_SECOND_ORDER_RUNGE_KUTTA` (undef, FLT_OPTIONS.h) — Define this to revert to old second-order Runge-Kutta integration
- `FLT_WITHOUT_X_PERIODICITY` (undef, FLT_OPTIONS.h) — Prevent floats to re-enter the opposite side of a periodic domain (stop instead)
- `FLT_WITHOUT_Y_PERIODICITY` (undef, FLT_OPTIONS.h)
- `DEVEL_FLT_EXCH2` (undef, FLT_OPTIONS.h) — Allow experimentation with pkg/flt + exch2 despite incomplete implementation

## Headers
- `FLT.h` — HEADER flt This header file contains variables that are used by the flt package. HEADER flt
- `FLT_BUFF.h` — BOP
- `FLT_OPTIONS.h` — CPP options file for FLT package
- `FLT_SIZE.h` — HEADER FLT_SIZE

## Routines (25)
`exch2_recv_get_vec.F`, `exch2_send_put_vec.F`, `exch_recv_get_vec.F`, `exch_send_put_vec.F`, `flt_down.F`, `flt_exch2.F`, `flt_exchg.F`, `flt_init_fixed.F`, `flt_init_varia.F`, `flt_interp_linear.F`, `flt_main.F`, `flt_mapping.F`, `flt_readparms.F`, `flt_runga2.F`, `flt_runga4.F`, `flt_traj.F`, `flt_up.F`, `flt_write_pickup.F`

## Called from outside the package
- `FLT_MAIN` ← `model/src/forward_step.F:1133`
- `FLT_INIT_FIXED` ← `model/src/packages_init_fixed.F:411`
- `FLT_INIT_VARIA` ← `model/src/packages_init_variables.F:323`
- `FLT_READPARMS` ← `model/src/packages_readparms.F:241`
- `FLT_WRITE_PICKUP` ← `model/src/packages_write_pickup.F:148`

## Verification experiments compiling it (2)
`MLAdjust` `exp4`
