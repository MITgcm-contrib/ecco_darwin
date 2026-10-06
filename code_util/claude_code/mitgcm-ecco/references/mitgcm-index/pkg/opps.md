# pkg/opps

OPPS convection (Paluszkiewicz & Romea penetrative plume scheme).

**runtime switch:** `useOPPS`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.opps`
**manual:** `doc/phys_pkgs/opps.rst`, `doc/examples/examples.rst`

## Namelist parameters
### OPPS_PARM01
- `MAX_ABE_ITERATIONS` — maximum for iteration on fractional size (default=1) In the present implementation, there is no iteration and max_abe_iterations should always be 1
- `OPPSdebugLevel` — sets internal debug level (default = 0) to produce some output for debugging
- `PlumeRadius` — default = 100 m
- `STABILITY_THRESHOLD` — threshold of vertical density difference, beyond which convection starts (default = -1.e-4 kg/m^3)
- `FRACTIONAL_AREA` — (initial) fractional area that plume(s) occupies (default = 0.1)
- `MAX_FRACTIONAL_AREA` — maximum of above (default = 0.8), not used
- `VERTICAL_VELOCITY` — initial (positive=downward) vertical velocity of plume (default=0.03m/s)
- `ENTRAINMENT_RATE` — default = - 0.05
- `useGCMwVel` — flag to replace VERTICAL_VELOCITY with actual vertical velocity of GCM, probably useless (default = .false.)

## CPP options (defaults as shipped)
- `ALLOW_OPPS_SNAPSHOT` (undef, OPPS_OPTIONS.h) — allow snap-shot OPPS output
- `ALLOW_OPPS_DEBUG` (define, OPPS_OPTIONS.h) — allow debugging OPPS_CALC

## Headers
- `OPPS.h` — BOP
- `OPPS_OPTIONS.h` — CPP options file for OPPS package. Use this file for selecting options within the OPPS package.

## Routines (7)
`opps_calc.F`, `opps_check.F`, `opps_init.F`, `opps_interface.F`, `opps_readparms.F`

## Called from outside the package
- `OPPS_CHECK` ← `model/src/packages_check.F:216`
- `OPPS_INIT` ← `model/src/packages_init_fixed.F:271`
- `OPPS_READPARMS` ← `model/src/packages_readparms.F:181`
- `OPPS_INTERFACE` ← `model/src/tracers_correction_step.F:109`

## Verification experiments compiling it (1)
`vermix`
