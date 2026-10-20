# pkg/tapenade

Tapenade automatic-differentiation support.

**pkg_depend:** +autodiff  (`+` requires, `-` excludes)
**runtime switch:** `useTAPENADE`-style flag in `data.pkg` (check exact name in packages_boot.F)

## CPP options (defaults as shipped)
- `ALLOW_TAPENADE_ACTIVE_READ_XYZ` (define, TAPENADE_OPTIONS.h)
- `ALLOW_TAPENADE_ACTIVE_READ_XY` (define, TAPENADE_OPTIONS.h)
- `ALLOW_TAPENADE_ACTIVE_WRITE` (undef, TAPENADE_OPTIONS.h)

## Headers
- `COST_TAP_TLM.h` — HEADER COST_TAP_TLM
- `TAPENADE_OPTIONS.h` — BOP
- `adBinomial.h` — 
- `adComplex.h` — 
- `adFixedPoint.h` — 
- `adStack.h` — 

## Routines (35)
`active_read_tap.F`, `active_write_tap.F`, `dummy_tap.F`, `stubs_tap_adj.F`, `stubs_tap_tlm.F`

## Called from outside the package
- `ADEXCH_3D_RL` ← `pkg/autodiff/addummy_in_dynamics.F:88,89`
- `ADEXCH_3D_RL` ← `pkg/autodiff/addummy_in_stepping.F:130,131,132,151`
- `ADEXCH_UV_3D_RL` ← `pkg/autodiff/addummy_in_stepping.F:133`
- `ADEXCH_UV_XY_RS` ← `pkg/autodiff/addummy_in_stepping.F:108,143`
- `ADEXCH_XY_RL` ← `pkg/autodiff/addummy_in_stepping.F:124,164,165,167`
- `ADEXCH_XY_RS` ← `pkg/autodiff/addummy_in_stepping.F:109,110,144,145`
- `ADEXCH_UV_3D_RL` ← `pkg/autodiff/copy_ad_uv_outp.F:112,114`
- `ADEXCH_3D_RL` ← `pkg/autodiff/copy_advar_outp.F:97,99`
- `ADEXCH_3D_RL` ← `pkg/monitor/monitor_ad.F:116,118,119,120`
- `ADEXCH_UV_3D_RL` ← `pkg/monitor/monitor_ad.F:117`
- `ADEXCH_3D_RL` ← `pkg/ptracers/ptracers_ad_dump.F:80`
- `ADEXCH_UV_3D_RL` ← `pkg/seaice/seaice_ad_dump.F:99`
- `ADEXCH_XY_RL` ← `pkg/seaice/seaice_ad_dump.F:91,92,93,96`

## Verification experiments compiling it (8)
`global_oce_biogeo_bling` `global_oce_latlon` `global_ocean.cs32x15` `halfpipe_streamice` `isomip` `lab_sea` `tutorial_global_oce_biogeo` `tutorial_tracer_adjsens`
