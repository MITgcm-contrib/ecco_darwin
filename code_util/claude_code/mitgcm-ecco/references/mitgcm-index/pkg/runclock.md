# pkg/runclock

Wall-clock based run termination (stop gracefully before walltime).

**runtime switch:** `useRUNCLOCK`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.runclock`

## Namelist parameters
### RUNCLOCK
- `RC_maxtime_hr`
- `RC_maxtime_mi`
- `RC_maxtime_sc`

## CPP options (defaults as shipped)
- `RUNCLOCK_USES_DATE_AND_TIME` (undef, RUNCLOCK_OPTIONS.h) — Define this macro if using an F90-compiler to compile and link the code

## Headers
- `RUNCLOCK.h` — Package flag
- `RUNCLOCK_OPTIONS.h` — CPP options file for RUNCLOCK package Use this file for selecting options within the RUNCLOCK package

## Routines (5)
`runclock_check.F`, `runclock_continue.F`, `runclock_gettime.F`, `runclock_init.F`, `runclock_readparms.F`

## Called from outside the package
- `RUNCLOCK_CHECK` ← `model/src/packages_check.F:496`
- `RUNCLOCK_INIT` ← `model/src/packages_init_fixed.F:152`
- `RUNCLOCK_READPARMS` ← `model/src/packages_readparms.F:418`
