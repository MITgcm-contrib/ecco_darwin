# pkg/salt_plume

Brine rejection from sea-ice growth distributed vertically (salt plume parameterisation, ECCO v4).

**runtime switch:** `useSALT_PLUME`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.salt_plume`
**manual:** `doc/examples/examples.rst`, `doc/ocean_state_est/ocean_state_est.rst`
**adjoint support files:** salt_plume_ad_check_lev1_dir.h, salt_plume_ad_check_lev2_dir.h, salt_plume_ad_check_lev3_dir.h, salt_plume_ad_check_lev4_dir.h, salt_plume_ad_diff.list

## Namelist parameters
### SALT_PLUME_PARM01
- `SaltPlumeSouthernOcean`
- `CriterionType`
- `PlumeMethod`
- `Npower`
- `SaltPlumeCriterion`
- `SPovershoot`
- `SPsalFRAC`
- `SPinflectionPoint`  _[ifdef SALT_PLUME_IN_LEADS]_
- `SaltPlumeSplitBasin`  _[ifdef SALT_PLUME_SPLIT_BASIN]_
- `SPbrineSconst` — salinity of brine pocket (g/kg)  _[ifdef SALT_PLUME_VOLUME]_
- `SPbrineSaltmax`  _[ifdef SALT_PLUME_VOLUME]_

## CPP options (defaults as shipped)
- `SALT_PLUME_IN_LEADS` (undef, SALT_PLUME_OPTIONS.h) — if seaice growth dh is from atmospheric cooling. if undefined: Activate pkg/salt_plume whenever seaice forms. This is the default of pkg/salt_plume.
- `SALT_PLUME_SPLIT_BASIN` (undef, SALT_PLUME_OPTIONS.h)
- `SALT_PLUME_VOLUME` (undef, SALT_PLUME_OPTIONS.h)

## Headers
- `SALT_PLUME.h` — --   SALT_PLUME parameters Find surface where the potential density (ref.lev=surface) is larger than surface density plus SaltPlumeCriterion.
- `SALT_PLUME_OPTIONS.h` — CPP options file for salt_plume package Use this file for selecting options within the salt_plume package
- `salt_plume_ad_check_lev1_dir.h` — ADJ STORE saltplumeflux   = comlev1, key = ikey_dynamics
- `salt_plume_ad_check_lev2_dir.h` — CADJ STORE saltplumeflux   = tapelev2, key = ilev_2
- `salt_plume_ad_check_lev3_dir.h` — CADJ STORE saltplumeflux   = tapelev3, key = ilev_3
- `salt_plume_ad_check_lev4_dir.h` — CADJ STORE saltplumeflux   = tapelev4, key = ilev_4

## Routines (15)
`salt_plume_apply.F`, `salt_plume_calc_depth.F`, `salt_plume_check.F`, `salt_plume_diagnostics_fill.F`, `salt_plume_diagnostics_init.F`, `salt_plume_do_exch.F`, `salt_plume_forcing_surf.F`, `salt_plume_frac.F`, `salt_plume_init_fixed.F`, `salt_plume_init_varia.F`, `salt_plume_mnc_init.F`, `salt_plume_readparms.F`, `salt_plume_tendency_apply_s.F`, `salt_plume_tendency_apply_t.F`, `salt_plume_volfrac.F`

## Called from outside the package
- `SALT_PLUME_TENDENCY_APPLY_S` ← `model/src/apply_forcing.F:952`
- `SALT_PLUME_TENDENCY_APPLY_T` ← `model/src/apply_forcing.F:720`
- `SALT_PLUME_APPLY` ← `model/src/do_oceanic_phys.F:911,915`
- `SALT_PLUME_CALC_DEPTH` ← `model/src/do_oceanic_phys.F:904`
- `SALT_PLUME_DIAGNOSTICS_FILL` ← `model/src/do_oceanic_phys.F:1127`
- `SALT_PLUME_DO_EXCH` ← `model/src/do_oceanic_phys.F:548`
- `SALT_PLUME_FORCING_SURF` ← `model/src/do_oceanic_phys.F:920`
- `SALT_PLUME_VOLFRAC` ← `model/src/do_oceanic_phys.F:908`
- `SALT_PLUME_TENDENCY_APPLY_S` ← `model/src/external_forcing.F:782`
- `SALT_PLUME_TENDENCY_APPLY_T` ← `model/src/external_forcing.F:578`
- `SALT_PLUME_FORCING_SURF` ← `model/src/external_forcing_surf.F:248`
- `SALT_PLUME_CHECK` ← `model/src/packages_check.F:339`
- `SALT_PLUME_INIT_FIXED` ← `model/src/packages_init_fixed.F:523`
- `SALT_PLUME_INIT_VARIA` ← `model/src/packages_init_variables.F:445`
- `SALT_PLUME_READPARMS` ← `model/src/packages_readparms.F:301`
- `SALT_PLUME_FRAC` ← `pkg/kpp/kpp_calc.F:691`
- `SALT_PLUME_FRAC` ← `pkg/kpp/kpp_routines.F:538,733,870`

## Verification experiments compiling it (3)
`cpl_aim+ocn` `lab_sea` `seaice_obcs`
