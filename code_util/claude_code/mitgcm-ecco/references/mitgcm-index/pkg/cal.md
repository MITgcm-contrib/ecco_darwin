# pkg/cal

Calendar tool (Gregorian/360-day/model calendar) used by exf, ecco, ctrl, diagnostics with calendarDumps.

**runtime switch:** `useCAL`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.cal`
**manual:** `doc/autodiff/autodiff.rst`, `doc/ocean_state_est/ocean_state_est.rst`, `doc/phys_pkgs/obcs.rst`
**adjoint support files:** cal_ad_diff.list

## Namelist parameters
### CAL_NML
- `TheCalendar` — type of calendar to use; available: 'model', 'gregorian' or 'noLeapYear'.
- `startDate_1`
- `startDate_2`
- `calendarDumps` — When set, approximate months (30-31 days) and years (360-372 days) for parameters chkPtFreq, pChkPtFreq and freq in pkg/diagnostics are converted to exact calendar months and years.

## Headers
- `CAL_OPTIONS.h` — BOP
- `cal.h` — HEADER calendar This header file contains variables that are used by the calendar tool. The calendar tool can be used in the ECCO SEALION release of t

## Routines (32)
`cal_addtime.F`, `cal_checkdate.F`, `cal_compdates.F`, `cal_convdate.F`, `cal_copydate.F`, `cal_daysformonth.F`, `cal_dayspermonth.F`, `cal_fulldate.F`, `cal_getdate.F`, `cal_getmonthsrec.F`, `cal_init_fixed.F`, `cal_intdays.F`, `cal_intmonths.F`, `cal_intyears.F`, `cal_isleap.F`, `cal_monthsforyear.F`, `cal_monthsperyear.F`, `cal_numints.F`, `cal_printdate.F`, `cal_printerror.F`, `cal_readparms.F`, `cal_set.F`, `cal_stepsforday.F`, `cal_stepsperday.F`, `cal_subdates.F`, `cal_summary.F`, `cal_time2dump.F`, `cal_timeinterval.F`, `cal_timepassed.F`, `cal_timestamp.F`, `cal_toseconds.F`, `cal_weekday.F`

## Called from outside the package
- `CAL_TIME2DUMP` ← `model/src/do_write_pickup.F:67,70`
- `CAL_GETDATE` ← `model/src/ini_mnc_vars.F:207`
- `CAL_INIT_FIXED` ← `model/src/packages_init_fixed.F:162`
- `CAL_READPARMS` ← `model/src/packages_readparms.F:156`
- `CAL_TIME2DUMP` ← `pkg/atm2d/atm2d_write_pickup.F:61,64`
- `CAL_GETDATE` ← `pkg/bling/bling_light.F:342`
- `CAL_ADDTIME` ← `pkg/ctrl/ctrl_get_gen_rec.F:155`
- `CAL_COPYDATE` ← `pkg/ctrl/ctrl_get_gen_rec.F:85`
- `CAL_GETDATE` ← `pkg/ctrl/ctrl_get_gen_rec.F:112,135`
- `CAL_GETMONTHSREC` ← `pkg/ctrl/ctrl_get_gen_rec.F:95`
- `CAL_TIMEINTERVAL` ← `pkg/ctrl/ctrl_get_gen_rec.F:154`
- `CAL_TIMEPASSED` ← `pkg/ctrl/ctrl_get_gen_rec.F:115,123,138,156`
- `CAL_TOSECONDS` ← `pkg/ctrl/ctrl_get_gen_rec.F:117,125,140,157`
- `CAL_FULLDATE` ← `pkg/ctrl/ctrl_init_rec.F:94,96`
- `CAL_TIMEPASSED` ← `pkg/ctrl/ctrl_init_rec.F:98`
- `CAL_TOSECONDS` ← `pkg/ctrl/ctrl_init_rec.F:100`
- `CAL_TIMEINTERVAL` ← `pkg/ctrl/ctrl_summary.F:256`
- `CAL_TIME2DUMP` ← `pkg/diagnostics/diagnostics_switch_onoff.F:113,203`
- `CAL_TIME2DUMP` ← `pkg/diagnostics/diagnostics_write.F:85,133`
- `CAL_TIME2DUMP` ← `pkg/diagnostics/diagnostics_write_adj.F:82`
- `CAL_ADDTIME` ← `pkg/ecco/cost_averagesflags.F:113`
- `CAL_COPYDATE` ← `pkg/ecco/cost_averagesflags.F:171,225,282`
- `CAL_FULLDATE` ← `pkg/ecco/cost_averagesflags.F:230,287`
- `CAL_GETDATE` ← `pkg/ecco/cost_averagesflags.F:109,110`
- `CAL_TIMEINTERVAL` ← `pkg/ecco/cost_averagesflags.F:112,186,242,292`
- `CAL_TIMEPASSED` ← `pkg/ecco/cost_averagesflags.F:178`
- `CAL_CONVDATE` ← `pkg/ecco/cost_gencal.F:84`
- `CAL_FULLDATE` ← `pkg/ecco/cost_gencal.F:91,93`
- `CAL_GETDATE` ← `pkg/ecco/cost_gencal.F:83`
- `CAL_TIMEPASSED` ← `pkg/ecco/cost_gencal.F:95`
- `CAL_TOSECONDS` ← `pkg/ecco/cost_gencal.F:96`
- `CAL_COPYDATE` ← `pkg/ecco/cost_gencost_sshv4.F:346,347,355,364`
- `CAL_CONVDATE` ← `pkg/ecco/cost_gencost_sstv4.F:205,360`
- `CAL_FULLDATE` ← `pkg/ecco/cost_gencost_sstv4.F:150,209,213,364`
- `CAL_GETDATE` ← `pkg/ecco/cost_gencost_sstv4.F:204,359`
- `CAL_PRINTDATE` ← `pkg/ecco/cost_gencost_sstv4.F:230,231,232,385`
- `CAL_TIMEPASSED` ← `pkg/ecco/cost_gencost_sstv4.F:215,370`
- `CAL_TOSECONDS` ← `pkg/ecco/cost_gencost_sstv4.F:216,371`
- `CAL_CONVDATE` ← `pkg/ecco/cost_sla_read.F:108`
- `CAL_FULLDATE` ← `pkg/ecco/cost_sla_read.F:112,115`
- `CAL_GETDATE` ← `pkg/ecco/cost_sla_read.F:107`
- `CAL_TIMEPASSED` ← `pkg/ecco/cost_sla_read.F:118`
- `CAL_TOSECONDS` ← `pkg/ecco/cost_sla_read.F:119`
- `CAL_COPYDATE` ← `pkg/ecco/ecco_cost_init_fixed.F:162,173`
- `CAL_FULLDATE` ← `pkg/ecco/ecco_cost_init_fixed.F:158,169`
- `CAL_PRINTDATE` ← `pkg/ecco/ecco_summary.F:78,79`
- `CAL_FULLDATE` ← `pkg/exf/exf_getffield_start.F:83`
- `CAL_GETDATE` ← `pkg/exf/exf_getffield_start.F:99`
- `CAL_TIMEPASSED` ← `pkg/exf/exf_getffield_start.F:90,100`
- `CAL_TOSECONDS` ← `pkg/exf/exf_getffield_start.F:92,102`
- `CAL_GETDATE` ← `pkg/exf/exf_getffieldrec.F:155`
- `CAL_TIMEPASSED` ← `pkg/exf/exf_getffieldrec.F:161`
- `CAL_TOSECONDS` ← `pkg/exf/exf_getffieldrec.F:162`
- `CAL_CONVDATE` ← `pkg/exf/exf_getmonthsrec.F:64`
- `CAL_GETDATE` ← `pkg/exf/exf_getmonthsrec.F:63`
- `CAL_GETMONTHSREC` ← `pkg/exf/exf_getmonthsrec.F:58`
- `CAL_GETMONTHSREC` ← `pkg/exf/exf_set_fld.F:137`
- `CAL_GETMONTHSREC` ← `pkg/exf/exf_set_uv.F:152`
- `CAL_GETDATE` ← `pkg/exf/exf_zenithangle.F:72`
- `CAL_TIMEPASSED` ← `pkg/exf/exf_zenithangle.F:82,91`
- `CAL_TOSECONDS` ← `pkg/exf/exf_zenithangle.F:83,92`
- `CAL_GETDATE` ← `pkg/mdsio/mdsio_write_meta.F:183`
- `CAL_GETMONTHSREC` ← `pkg/obcs/obcs_exf_load.F:288,380,587,679`
- `CAL_FULLDATE` ← `pkg/obsfit/obsfit_init_fixed.F:559`
- `CAL_TIMEPASSED` ← `pkg/obsfit/obsfit_init_fixed.F:561`
- `CAL_TOSECONDS` ← `pkg/obsfit/obsfit_init_fixed.F:563`
- `CAL_FULLDATE` ← `pkg/profiles/profiles_init_fixed.F:563`
- `CAL_TIMEPASSED` ← `pkg/profiles/profiles_init_fixed.F:565`
- `CAL_TOSECONDS` ← `pkg/profiles/profiles_init_fixed.F:567`
- `CAL_FULLDATE` ← `pkg/seaice/seaice_cost_init_fixed.F:44,56`
- `CAL_TIMEPASSED` ← `pkg/seaice/seaice_cost_init_fixed.F:46,58`
- `CAL_TOSECONDS` ← `pkg/seaice/seaice_cost_init_fixed.F:48,60`
- `CAL_TIME2DUMP` ← `pkg/seaice/seaice_jfnk.F:366`
- `CAL_TIME2DUMP` ← `pkg/seaice/seaice_krylov.F:430`
- `CAL_GETDATE` ← `pkg/seaice/seaice_readparms.F:619`

## Verification experiments compiling it (2)
`global_oce_biogeo_bling` `obcs_ctrl`
