# pkg/showflops

FLOP counting/timing (PAPI).

**pkg_depend:** +runclock  (`+` requires, `-` excludes)
**runtime switch:** `useSHOWFLOPS`-style flag in `data.pkg` (check exact name in packages_boot.F)

## CPP options (defaults as shipped)
- `USE_FLIPS` (undef, SHOWFLOPS_OPTIONS.h)

## Headers
- `SHOWFLOPS.h` — CE107 common block for per timestep timing !TIMING VARIABLES == Timing variables ==
- `SHOWFLOPS_INIT.h` — CE107 common block for per timestep timing !TIMING VARIABLES == Timing variables ==
- `SHOWFLOPS_OPTIONS.h` — CPP options file for SHOWFLOPS package Use this file for selecting options within the SHOWFLOPS package

## Routines (3)
`showflops_init.F`, `showflops_inloop.F`, `showflops_insolve.F`

## Called from outside the package
- `SHOWFLOPS_INLOOP` ← `model/src/forward_step.F:1213`
- `SHOWFLOPS_INSOLVE` ← `model/src/solve_for_pressure.F:463`
- `SHOWFLOPS_INIT` ← `model/src/the_main_loop.F:391`
