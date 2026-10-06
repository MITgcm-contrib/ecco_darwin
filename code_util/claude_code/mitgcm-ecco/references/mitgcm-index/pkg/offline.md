# pkg/offline

Offline mode: reads precomputed velocities/diffusivities/forcing (e.g. ECCO archive) to drive passive tracers/BGC without dynamics.

**runtime switch:** `useOFFLINE`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.off`
**manual:** `doc/examples/examples.rst`
**adjoint support files:** offline_ad_check_lev1_dir.h, offline_ad_check_lev2_dir.h, offline_ad_check_lev3_dir.h, offline_ad_check_lev4_dir.h, offline_ad_diff.list

## Namelist parameters
### OFFLINE_PARM01
- `UvelFile`
- `VvelFile`
- `WvelFile`
- `ThetFile`
- `SaltFile`
- `GMwxFile`
- `GMwyFile`
- `GMwzFile`
- `ConvFile`
- `KPP_DiffSFile`
- `KPP_ghatKFile`
- `HFluxFile`
- `SFluxFile`
- `IceFile`
- `KPP_ghatFile`
### OFFLINE_PARM02
- `offlineIter0`
- `deltaToffline`
- `offlineTimeOffset`
- `offlineForcingPeriod`
- `offlineForcingCycle`
- `offlineLoadPrec`
- `offlineOffsetIter`

## CPP options (defaults as shipped)
- `NOT_MODEL_FILES` (undef, OFFLINE_OPTIONS.h)

## Headers
- `OFFLINE.h` — variable for forcing offline tracer
- `OFFLINE_OPTIONS.h` — BOP
- `OFFLINE_SWITCH.h` — variable for switching on/off some calculations
- `offline_ad_check_lev1_dir.h` — ADJ STORE save0 = comlev1, key = ikey_dynamics ADJ STORE save1 = comlev1, key = ikey_dynamics ADJ STORE tave0 = comlev1, key = ikey_dynamics ADJ STORE
- `offline_ad_check_lev2_dir.h` — ph( not sure exactly why these are needed. ADJ STORE phi0surf          = tapelev2, key = ilev_2 ADJ STORE rhoinsitu         = tapelev2, key = ilev_2 A
- `offline_ad_check_lev3_dir.h` — ph( not sure exactly why these are needed. ADJ STORE phi0surf          = tapelev3, key = ilev_3 ADJ STORE rhoinsitu         = tapelev3, key = ilev_3 A
- `offline_ad_check_lev4_dir.h` — ph( not sure exactly why these are needed. ADJ STORE phi0surf          = tapelev4, key = ilev_4 ADJ STORE rhoinsitu         = tapelev4, key = ilev_4 A

## Routines (6)
`offline_check.F`, `offline_fields_load.F`, `offline_get_diffus.F`, `offline_init_varia.F`, `offline_readparms.F`, `offline_reset_parms.F`

## Called from outside the package
- `OFFLINE_GET_DIFFUS` ← `model/src/do_oceanic_phys.F:1077`
- `OFFLINE_FIELDS_LOAD` ← `model/src/forward_step.F:824`
- `OFFLINE_CHECK` ← `model/src/packages_check.F:318`
- `OFFLINE_INIT_VARIA` ← `model/src/packages_init_variables.F:185`
- `OFFLINE_READPARMS` ← `model/src/packages_readparms.F:266`
- `OFFLINE_RESET_PARMS` ← `model/src/set_parms.F:47`
- `OFFLINE_FIELDS_LOAD` ← `pkg/bling/bling_carbonate_init.F:79`

## Verification experiments compiling it (2)
`tutorial_cfc_offline` `tutorial_dic_adjoffline`
