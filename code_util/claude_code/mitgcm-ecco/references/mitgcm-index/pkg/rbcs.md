# pkg/rbcs

Relaxation boundary conditions: 3-D relaxation of T, S, ptracers (and U/V) toward prescribed fields with masks/timescales.

**runtime switch:** `useRBCS`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.rbcs`
**manual:** `doc/phys_pkgs/rbcs.rst`, `doc/examples/examples.rst`, `doc/examples/reentrant_channel/reentrant_channel.rst`, `doc/getting_started/getting_started.rst`, `doc/phys_pkgs/obcs.rst`
**adjoint support files:** rbcs_ad_check_lev1_dir.h, rbcs_ad_check_lev2_dir.h, rbcs_ad_check_lev3_dir.h, rbcs_ad_check_lev4_dir.h, rbcs_ad_diff.list

## Namelist parameters
### RBCS_PARM01
- `tauRelaxU`
- `tauRelaxV`
- `tauRelaxT`
- `tauRelaxS`
- `relaxMaskUFile`
- `relaxMaskVFile`
- `relaxMaskFile`
- `relaxUFile`
- `relaxVFile`
- `relaxTFile`
- `relaxSFile`
- `useRBCuVel`
- `useRBCvVel`
- `useRBCtemp`
- `useRBCsalt`
- `useRBCptracers`
- `rbcsIniter`
- `rbcsForcingPeriod` — period of rbc data (in seconds)
- `rbcsForcingCycle` — cycle of rbc data (in seconds)
- `rbcsForcingOffset` — model time at beginning of first rbc period
- `rbcsVanishingTime` — when rbcsVanishingTime .NE. 0. the relaxation strength reduces
- `rbcsSingleTimeFiles` — if .TRUE., rbc fields are given 1 file per time
- `deltaTrbcs` — time step used to compute iteration numbers for singleTimeFiles
- `rbcsIter0` — singleTimeFile iteration number corresponding to rbcsForcingOffset
### RBCS_PARM02
- `useRBCpTrNum`  _[ifdef ALLOW_PTRACERS]_
- `tauRelaxPTR`  _[ifdef ALLOW_PTRACERS]_
- `relaxPtracerFile`  _[ifdef ALLOW_PTRACERS]_

## CPP options (defaults as shipped)
- `DISABLE_RBCS_MOM` (undef, RBCS_OPTIONS.h) — o disable relaxation conditions on momemtum

## Headers
- `RBCS_FIELDS.h` — BOP
- `RBCS_OPTIONS.h` — CPP options file for pkg RBCS Use this file for selecting options within package "RBCS"
- `RBCS_PARAMS.h` — BOP
- `RBCS_SIZE.h` — BOP
- `rbcs_ad_check_lev1_dir.h` — ADJ STORE rbct0 = comlev1, key = ikey_dynamics ADJ STORE rbct1 = comlev1, key = ikey_dynamics ADJ STORE rbcs0 = comlev1, key = ikey_dynamics ADJ STORE
- `rbcs_ad_check_lev2_dir.h` — ADJ STORE rbct0 = tapelev2, key = ilev_2 ADJ STORE rbct1 = tapelev2, key = ilev_2 ADJ STORE rbcs0 = tapelev2, key = ilev_2 ADJ STORE rbcs1 = tapelev2,
- `rbcs_ad_check_lev3_dir.h` — ADJ STORE rbct0 = tapelev3, key = ilev_3 ADJ STORE rbct1 = tapelev3, key = ilev_3 ADJ STORE rbcs0 = tapelev3, key = ilev_3 ADJ STORE rbcs1 = tapelev3,
- `rbcs_ad_check_lev4_dir.h` — ADJ STORE rbct0 = tapelev4, key = ilev_4 ADJ STORE rbct1 = tapelev4, key = ilev_4 ADJ STORE rbcs0 = tapelev4, key = ilev_4 ADJ STORE rbcs1 = tapelev4,

## Routines (5)
`rbcs_add_tendency.F`, `rbcs_fields_load.F`, `rbcs_init_fixed.F`, `rbcs_init_varia.F`, `rbcs_readparms.F`

## Called from outside the package
- `RBCS_ADD_TENDENCY` ← `model/src/apply_forcing.F:170,360,728,960`
- `RBCS_ADD_TENDENCY` ← `model/src/external_forcing.F:121,261,586,790`
- `RBCS_FIELDS_LOAD` ← `model/src/load_fields_driver.F:249`
- `RBCS_INIT_FIXED` ← `model/src/packages_init_fixed.F:447`
- `RBCS_INIT_VARIA` ← `model/src/packages_init_variables.F:354`
- `RBCS_READPARMS` ← `model/src/packages_readparms.F:261`
- `RBCS_ADD_TENDENCY` ← `pkg/ptracers/ptracers_apply_forcing.F:116`

## Verification experiments compiling it (2)
`exp4` `tutorial_reentrant_channel`
