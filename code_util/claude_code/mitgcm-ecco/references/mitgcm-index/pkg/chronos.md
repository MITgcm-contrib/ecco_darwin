# pkg/chronos

Clock/date utilities (legacy).

**runtime switch:** `useCHRONOS`-style flag in `data.pkg` (check exact name in packages_boot.F)

## Headers
- `chronos.h` — *****                   Clock Variables                          *****

## Routines (16)
`chronos.F`

## Called from outside the package
- `ASTRO` ← `pkg/fizhi/do_fizhi.F:176,178`
- `SET_ALARM` ← `pkg/fizhi/fizhi_alarms.F:108,109,110,111`
- `GET_TIME` ← `pkg/fizhi/fizhi_clockstuff.F:154,182,209,239`
- `TICK` ← `pkg/fizhi/fizhi_clockstuff.F:242,1073,1074,1082`
- `SET_ALARM` ← `pkg/fizhi/fizhi_diagalarms.F:57,74`
- `GET_ALARM` ← `pkg/fizhi/fizhi_driver.F:125,126,127,128`
- `ASTRO` ← `pkg/fizhi/fizhi_swrad.F:134,138,140,144`
- `TICK` ← `pkg/fizhi/fizhi_swrad.F:137,143`
- `ASTRO` ← `pkg/fizhi/fizhi_turb.F:503,505`
- `GET_ALARM` ← `pkg/fizhi/fizhi_turb.F:282`
- `TICK` ← `pkg/fizhi/fizhi_update_time.F:26`
- `INTERP_TIME` ← `pkg/fizhi/fizhi_utils.F:854`
- `TIME_BOUND` ← `pkg/fizhi/fizhi_utils.F:853`
- `INTERP_TIME` ← `pkg/fizhi/update_chemistry_exports.F:79`
- `TIME_BOUND` ← `pkg/fizhi/update_chemistry_exports.F:78`
- `ASTRO` ← `pkg/fizhi/update_earth_exports.F:108,110`
- `INTERP_TIME` ← `pkg/fizhi/update_ocean_exports.F:346,608`
