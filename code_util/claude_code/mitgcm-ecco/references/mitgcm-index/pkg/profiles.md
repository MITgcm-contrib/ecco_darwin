# pkg/profiles

Model-data comparison for in situ profiles (Argo, CTD, XBT...) for ECCO cost.

**pkg_depend:** +cal  (`+` requires, `-` excludes)
**runtime switch:** `usePROFILES`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.profiles`
**manual:** `doc/examples/examples.rst`, `doc/ocean_state_est/ocean_state_est.rst`, `doc/related_projects/related_projects.rst`
**adjoint support files:** profiles_ad_diff.list

## Namelist parameters
### PROFILES_NML
- `profilesDir`
- `profilesfiles`
- `mult_profiles`
- `mult_profiles_mean`
- `profiles_mean_indsamples`  _[ifdef ALLOW_PROFILES_SAMPLESPLIT_COST]_
- `prof_facmod`
- `prof_names`
- `prof_namesmod`
- `prof_namesclim`  _[ifdef ALLOW_PROFILES_CLIMMASK]_
- `prof_itracer`
- `profilesDoNcOutput`
- `profilesDoGenGrid`
- `prof_make_nc`
- `prof_dBugLevel` — control debug print to STDOUT or log file, higher -> more

## CPP options (defaults as shipped)
- `ALLOW_PROFILES_CLIMMASK` (undef, PROFILES_OPTIONS.h) — -- Undocumented Options:
- `ALLOW_PROFILES_EXCLUDE_CORNERS` (undef, PROFILES_OPTIONS.h)
- `ALLOW_PROFILES_SAMPLESPLIT_COST` (undef, PROFILES_OPTIONS.h)

## Headers
- `PROFILES_OPTIONS.h` — CPP options file for PROFILES package Use this file for selecting options within the PROFILES package
- `PROFILES_SIZE.h` — NOBSGLOB            :: maximum number of profiles per file and tile NFILESPROFMAX       :: maximum number of files NVARMAX             :: maximum numb
- `profiles.h` — --  PROF_PARAMS common block: prof_dBugLevel :: control debug print to STDOUT or log file, higher -> more

## Routines (25)
`active_file_control_profiles.F`, `active_file_profiles.F`, `active_file_profiles_ad.F`, `active_file_profiles_g.F`, `active_file_profiles_tap_adj.F`, `cost_profiles.F`, `profiles_cost.F`, `profiles_cost_final.F`, `profiles_findunit.F`, `profiles_ini_io.F`, `profiles_init_fixed.F`, `profiles_init_ncfile.F`, `profiles_init_varia.F`, `profiles_inloop.F`, `profiles_interp.F`, `profiles_make_ncfile.F`, `profiles_nc_utils.F`, `profiles_readparms.F`, `profiles_readvector.F`

## Called from outside the package
- `PROFILES_INI_IO` ← `model/src/ini_model_io.F:236`
- `PROFILES_INIT_FIXED` ← `model/src/packages_init_fixed.F:390`
- `PROFILES_INIT_VARIA` ← `model/src/packages_init_variables.F:576`
- `PROFILES_READPARMS` ← `model/src/packages_readparms.F:348`
- `PROFILES_COST` ← `model/src/the_main_loop.F:751`
- `PROFILES_INLOOP` ← `model/src/the_main_loop.F:682,748`
- `PROFILES_NC_CLOSE` ← `model/src/the_model_main.F:770`
- `PROFILES_COST_FINAL` ← `pkg/cost/cost_final.F:89`

## Verification experiments compiling it (2)
`global_oce_biogeo_bling` `global_oce_latlon`
