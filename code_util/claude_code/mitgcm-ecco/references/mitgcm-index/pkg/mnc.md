# pkg/mnc

NetCDF I/O (per-tile files) for state, pickups, diagnostics when diag_mnc is on.

**runtime switch:** `useMNC`-style flag in `data.pkg` (check exact name in packages_boot.F)
**manual:** `doc/examples/baroclinic_gyre/baroclinic_gyre.rst`, `doc/examples/held_suarez_cs/held_suarez_cs.rst`, `doc/getting_started/getting_started.rst`, `doc/outp_pkgs/outp_pkgs.rst`, `doc/phys_pkgs/gmredi.rst`

## README
```
API Discussions:
================

As discussed in our group meeting of 2003-12-17 (AJA, CNH, JMC, AM,
PH, EH3), the NetCDF interface should resemble the following FORTRAN
subroutines:

  1) "stubs" of the form:  MNC_WV_[G|L]2D_R[S|L] ()

  2) MNC_INIT_VGRID('V_GRID_TYPE', nx, ny, nz, zc, zg)

  3) MNC_INIT_HGRID('H_GRID_TYPE', nx, ny, xc, yc, xg, yg)

  4) MNC_INIT_VAR('file', 'Vname', 'Vunits', 'H_GTYPE', 'V_GTYPE', PREC, FillVal)

  5) MNC_WRITE_VAR('file', 'Vname', var, bi, bj, myThid)

This is a reasonable start but its inflexible since it isn't easily
generalized to grids with dimensions other than [2,3,4] or grids with
non-horizontal orientations (eg. vertical slices).


Generalizing what we would like to write as "variables defined on 1-D
to n-D grids", one can imagine a small number of objects containing
all the relevant information:

  a dimension:   [ name, size, units ]
  a grid:        [ name, 1+ dim-ref ]
  a variable:    [ name, units, *1* grid-ref, data ]
  an attribute:  [ name, units, data ]
```

## Namelist parameters
### MNC_01
- `mnc_use_indir` — use "mnc_indir_str" as input filename prefix
- `mnc_use_outdir` — use "mnc_outdir_str" as output filename prefix
- `mnc_outdir_date` — use a date string within the output dir name
- `mnc_outdir_num` — use a seq. number within the output dir name
- `mnc_use_name_ni0` — use nIter0 in all the file names
- `mnc_echo_gvtypes` — echo type names (fails on many platforms)
- `pickup_write_mnc` — use mnc to write pickups
- `pickup_read_mnc` — use mnc to read  pickups
- `timeave_mnc`
- `snapshot_mnc`
- `monitor_mnc`
- `autodiff_mnc`
- `writegrid_mnc` — use mnc to write model-grid arrays to file
- `readgrid_mnc` — read INI_CURVILINEAR_GRID() info using mnc
- `mnc_outdir_str` — name of the output directory
- `mnc_indir_str` — name of the input directory
- `mnc_max_fsize` — maximum file size
- `mnc_filefreq`
- `mnc_read_bathy`
- `mnc_read_salt`
- `mnc_read_theta`

## CPP options (defaults as shipped)
- `MNC_DEF_FMNC` (define, MNC_OPTIONS.h) — These are the default minimum number of characters used for the per-file and per-tile file names
- `MNC_DEF_TMNC` (define, MNC_OPTIONS.h)

## Headers
- `MNC_BUFF.h` — The following is the size of the buffer used by MNC to read and write portions to/from NetCDF files.  The sizes of all reads and writes are checked an
- `MNC_COMMON.h` — MNC : an MITgcm wrapper package for NetCDF The following common block is the "state" for the MNC interface to NetCDF.  The intent is to keep track of 
- `MNC_OPTIONS.h` — EH3 package-specific options go here
- `MNC_PARAMS.h` — BOP
- `MNC_SIZE.h` — EH3 ;;; Local Variables: *** EH3 ;;; mode:fortran *** EH3 ;;; End: ***

## Routines (93)
`MNC_CW_READWRITE_I.F`, `MNC_CW_READWRITE_RL.F`, `MNC_CW_READWRITE_RS.F`, `mnc_cw_citer.F`, `mnc_cw_cvars.F`, `mnc_cw_init.F`, `mnc_cw_missingvals.F`, `mnc_cw_model_attr.F`, `mnc_cw_udim.F`, `mnc_cw_write_grid_info.F`, `mnc_cwrapper.F`, `mnc_dim.F`, `mnc_dump.F`, `mnc_file.F`, `mnc_grid.F`, `mnc_init.F`, `mnc_readparms.F`, `mnc_update_time.F`, `mnc_utils.F`, `mnc_var.F`

## Called from outside the package
- `MNC_UPDATE_TIME` ← `model/src/forward_step.F:813`
- `MNC_CW_SET_UDIM` ← `model/src/ini_cori.F:200`
- `MNC_CW_RS_R` ← `model/src/ini_curvilinear_grid.F:226,227,228,229`
- `MNC_CW_SET_CITER` ← `model/src/ini_curvilinear_grid.F:224`
- `MNC_CW_SET_UDIM` ← `model/src/ini_curvilinear_grid.F:223,225`
- `MNC_FILE_CLOSE_ALL_MATCHING` ← `model/src/ini_curvilinear_grid.F:222`
- `MNC_CW_ADD_VNAME` ← `model/src/ini_depths.F:107`
- `MNC_CW_DEL_VNAME` ← `model/src/ini_depths.F:114`
- `MNC_CW_RS_R` ← `model/src/ini_depths.F:112`
- `MNC_CW_SET_CITER` ← `model/src/ini_depths.F:110`
- `MNC_CW_SET_UDIM` ← `model/src/ini_depths.F:109,111`
- `MNC_FILE_CLOSE_ALL_MATCHING` ← `model/src/ini_depths.F:108,113`
- `MNC_CW_SET_UDIM` ← `model/src/ini_grid.F:201`
- `MNC_CW_ADD_VATTR_TEXT` ← `model/src/ini_mnc_vars.F:41,43,48,50`
- `MNC_CW_ADD_VNAME` ← `model/src/ini_mnc_vars.F:40,47,54,61`
- `MNC_CW_INIT` ← `model/src/ini_model_io.F:220`
- `MNC_INIT` ← `model/src/ini_model_io.F:219`
- `MNC_CW_RL_R` ← `model/src/ini_salt.F:72`
- `MNC_CW_SET_CITER` ← `model/src/ini_salt.F:70`
- `MNC_CW_SET_UDIM` ← `model/src/ini_salt.F:69,71`
- `MNC_FILE_CLOSE_ALL_MATCHING` ← `model/src/ini_salt.F:68,73`
- `MNC_CW_RL_R` ← `model/src/ini_theta.F:73`
- `MNC_CW_SET_CITER` ← `model/src/ini_theta.F:71`
- `MNC_CW_SET_UDIM` ← `model/src/ini_theta.F:70,72`
- `MNC_FILE_CLOSE_ALL_MATCHING` ← `model/src/ini_theta.F:69,74`
- `MNC_READPARMS` ← `model/src/packages_readparms.F:150`
- `MNC_CW_RL_R` ← `model/src/read_pickup.F:493,494,497,498`
- `MNC_CW_SET_CITER` ← `model/src/read_pickup.F:492`
- `MNC_CW_SET_UDIM` ← `model/src/read_pickup.F:491`
- `MNC_FILE_CLOSE_ALL_MATCHING` ← `model/src/read_pickup.F:490`
- `MNC_FILE_CLOSE_ALL` ← `model/src/the_model_main.F:762`
- `MNC_CW_RS_W` ← `model/src/write_grid.F:171,172,173,174`
- `MNC_CW_SET_CITER` ← `model/src/write_grid.F:169`
- `MNC_CW_SET_UDIM` ← `model/src/write_grid.F:168,170`
- `MNC_CW_WRITE_GRID_COORD` ← `model/src/write_grid.F:219,224,226,229`
- `MNC_CW_I_W_S` ← `model/src/write_pickup.F:406`
- `MNC_CW_RL_W` ← `model/src/write_pickup.F:407,408,411,412`
- `MNC_CW_RL_W_S` ← `model/src/write_pickup.F:405`
- `MNC_CW_SET_CITER` ← `model/src/write_pickup.F:399,401`
- `MNC_CW_SET_UDIM` ← `model/src/write_pickup.F:397,404`
- `MNC_FILE_CLOSE_ALL_MATCHING` ← `model/src/write_pickup.F:438`
- `MNC_CW_I_W_S` ← `model/src/write_state.F:187,199,207`
- `MNC_CW_RL_W` ← `model/src/write_state.F:189,190,191,192`
- `MNC_CW_RL_W_S` ← `model/src/write_state.F:185,197,205`
- `MNC_CW_SET_UDIM` ← `model/src/write_state.F:184,186,196,198`
- `MNC_CW_I_W_S` ← `pkg/aim_v23/aim_diagnostics.F:133,239`
- `MNC_CW_RL_W` ← `pkg/aim_v23/aim_diagnostics.F:136,137,138,139`
- `MNC_CW_RL_W_S` ← `pkg/aim_v23/aim_diagnostics.F:131,237`
- `MNC_CW_SET_UDIM` ← `pkg/aim_v23/aim_diagnostics.F:130,132,236,238`
- `MNC_CW_ADD_GNAME` ← `pkg/aim_v23/aim_mnc_init.F:62,73`
- `MNC_CW_ADD_VATTR_TEXT` ← `pkg/aim_v23/aim_mnc_init.F:76,79,82,85`
- `MNC_CW_ADD_VNAME` ← `pkg/aim_v23/aim_mnc_init.F:75,78,81,84`
- `MNC_CW_I_W_S` ← `pkg/autodiff/addummy_for_etan.F:119`
- `MNC_CW_RL_W` ← `pkg/autodiff/addummy_for_etan.F:122`
- `MNC_CW_RL_W_S` ← `pkg/autodiff/addummy_for_etan.F:117,120`
- `MNC_CW_SET_UDIM` ← `pkg/autodiff/addummy_for_etan.F:116,118`
- `MNC_CW_I_W_S` ← `pkg/autodiff/addummy_in_stepping.F:450`
- `MNC_CW_RL_W` ← `pkg/autodiff/addummy_in_stepping.F:457,458,460,462`
- `MNC_CW_RL_W_S` ← `pkg/autodiff/addummy_in_stepping.F:448,451`
- `MNC_CW_RS_W` ← `pkg/autodiff/addummy_in_stepping.F:504,505,506,507`
- `MNC_CW_SET_UDIM` ← `pkg/autodiff/addummy_in_stepping.F:447,449`
- `MNC_CW_ADD_VATTR_TEXT` ← `pkg/autodiff/autodiff_ini_model_io.F:124,125,127,131`
- `MNC_CW_ADD_VNAME` ← `pkg/autodiff/autodiff_ini_model_io.F:123,130,137,144`
- `MNC_CW_ADD_VATTR_TEXT` ← `pkg/cd_code/cd_code_init_fixed.F:26,28,30,35`
- `MNC_CW_ADD_VNAME` ← `pkg/cd_code/cd_code_init_fixed.F:24,33,42,51`
- `MNC_CW_RL_R` ← `pkg/cd_code/cd_code_read_pickup.F:57,58,59,60`
- `MNC_CW_SET_CITER` ← `pkg/cd_code/cd_code_read_pickup.F:55`
- `MNC_CW_SET_UDIM` ← `pkg/cd_code/cd_code_read_pickup.F:54,56`
- `MNC_FILE_CLOSE_ALL_MATCHING` ← `pkg/cd_code/cd_code_read_pickup.F:53`
- `MNC_CW_RL_W` ← `pkg/cd_code/cd_code_write_pickup.F:61,62,63,64`
- `MNC_CW_SET_CITER` ← `pkg/cd_code/cd_code_write_pickup.F:56,58`
- `MNC_CW_SET_UDIM` ← `pkg/cd_code/cd_code_write_pickup.F:54,60`
- `MNC_FILE_CLOSE_ALL_MATCHING` ← `pkg/cd_code/cd_code_write_pickup.F:66`
- `MNC_CW_ADD_GNAME` ← `pkg/diagnostics/diagnostics_mnc_out.F:115,150,345`
- `MNC_CW_ADD_VATTR_TEXT` ← `pkg/diagnostics/diagnostics_mnc_out.F:119,165,172,179`
- `MNC_CW_ADD_VNAME` ← `pkg/diagnostics/diagnostics_mnc_out.F:117,151,347`
- `MNC_CW_DEL_GNAME` ← `pkg/diagnostics/diagnostics_mnc_out.F:130,189,402`
- `MNC_CW_DEL_VNAME` ← `pkg/diagnostics/diagnostics_mnc_out.F:129,188,401`
- `MNC_CW_I_W_S` ← `pkg/diagnostics/diagnostics_mnc_out.F:103`
- `MNC_CW_RL_W` ← `pkg/diagnostics/diagnostics_mnc_out.F:126,393,397`
- `MNC_CW_RL_W_S` ← `pkg/diagnostics/diagnostics_mnc_out.F:101`
- `MNC_CW_RS_W` ← `pkg/diagnostics/diagnostics_mnc_out.F:187`
- `MNC_CW_SET_UDIM` ← `pkg/diagnostics/diagnostics_mnc_out.F:100,102`
- `MNC_CW_VATTR_MISSING` ← `pkg/diagnostics/diagnostics_mnc_out.F:123,185,377,386`
- `MNC_CW_ADD_GNAME` ← `pkg/diagnostics/diagnostics_read_pickup.F:93,116`
- `MNC_CW_ADD_VATTR_TEXT` ← `pkg/diagnostics/diagnostics_read_pickup.F:120`
- `MNC_CW_ADD_VNAME` ← `pkg/diagnostics/diagnostics_read_pickup.F:95,118`
- `MNC_CW_DEL_GNAME` ← `pkg/diagnostics/diagnostics_read_pickup.F:100,127`
- `MNC_CW_DEL_VNAME` ← `pkg/diagnostics/diagnostics_read_pickup.F:99,126`
- `MNC_CW_RL_R` ← `pkg/diagnostics/diagnostics_read_pickup.F:97`
- `MNC_CW_SET_UDIM` ← `pkg/diagnostics/diagnostics_read_pickup.F:69`
- `MNC_FILE_CLOSE_ALL_MATCHING` ← `pkg/diagnostics/diagnostics_read_pickup.F:68`
- `MNC_CW_ADD_GNAME` ← `pkg/diagnostics/diagnostics_write_pickup.F:124,152`
- `MNC_CW_ADD_VATTR_TEXT` ← `pkg/diagnostics/diagnostics_write_pickup.F:128,156`
- `MNC_CW_ADD_VNAME` ← `pkg/diagnostics/diagnostics_write_pickup.F:126,154`
- `MNC_CW_DEL_GNAME` ← `pkg/diagnostics/diagnostics_write_pickup.F:135,163`
- `MNC_CW_DEL_VNAME` ← `pkg/diagnostics/diagnostics_write_pickup.F:134,162`
- `MNC_CW_I_W` ← `pkg/diagnostics/diagnostics_write_pickup.F:159`
- `MNC_CW_I_W_S` ← `pkg/diagnostics/diagnostics_write_pickup.F:99`
- `MNC_CW_RL_W` ← `pkg/diagnostics/diagnostics_write_pickup.F:131`
- `MNC_CW_RL_W_S` ← `pkg/diagnostics/diagnostics_write_pickup.F:98`
- `MNC_CW_SET_CITER` ← `pkg/diagnostics/diagnostics_write_pickup.F:90,92`
- `MNC_CW_SET_UDIM` ← `pkg/diagnostics/diagnostics_write_pickup.F:88,95`
- `MNC_CW_ADD_GNAME` ← `pkg/diagnostics/diagstats_mnc_out.F:124,143,226,232`
- `MNC_CW_ADD_VATTR_TEXT` ← `pkg/diagnostics/diagstats_mnc_out.F:128,158,165,172`
- `MNC_CW_ADD_VNAME` ← `pkg/diagnostics/diagstats_mnc_out.F:126,144,247,289`
- `MNC_CW_DEL_GNAME` ← `pkg/diagnostics/diagstats_mnc_out.F:136,178,313,314`
- `MNC_CW_DEL_VNAME` ← `pkg/diagnostics/diagstats_mnc_out.F:135,177,279,307`
- `MNC_CW_I_W_S` ← `pkg/diagnostics/diagstats_mnc_out.F:106`
- `MNC_CW_RL_W` ← `pkg/diagnostics/diagstats_mnc_out.F:132,269,304`
- `MNC_CW_RL_W_S` ← `pkg/diagnostics/diagstats_mnc_out.F:104`
- `MNC_CW_RS_W` ← `pkg/diagnostics/diagstats_mnc_out.F:176`
- `MNC_CW_SET_UDIM` ← `pkg/diagnostics/diagstats_mnc_out.F:103,105`
- `MNC_CW_ADD_VATTR_TEXT` ← `pkg/dic/dic_mnc_init.F:32,34,39,41`
- `MNC_CW_ADD_VNAME` ← `pkg/dic/dic_mnc_init.F:30,37,44,51`
- `MNC_CW_I_W_S` ← `pkg/exf/exf_adjoint_snapshots_ad.F:164`
- `MNC_CW_RL_W` ← `pkg/exf/exf_adjoint_snapshots_ad.F:168,170,172,174`
- `MNC_CW_RL_W_S` ← `pkg/exf/exf_adjoint_snapshots_ad.F:162,165`
- `MNC_CW_SET_UDIM` ← `pkg/exf/exf_adjoint_snapshots_ad.F:161,163`
- `MNC_CW_APPEND_VNAME` ← `pkg/exf/exf_monitor.F:83`
- `MNC_CW_RL_W_S` ← `pkg/exf/exf_monitor.F:86`
- `MNC_CW_SET_UDIM` ← `pkg/exf/exf_monitor.F:85,88`
- `MNC_CW_APPEND_VNAME` ← `pkg/exf/exf_monitor_ad.F:85`
- `MNC_CW_RL_W_S` ← `pkg/exf/exf_monitor_ad.F:88`
- `MNC_CW_SET_UDIM` ← `pkg/exf/exf_monitor_ad.F:87,90`
- `MNC_CW_ADD_GNAME` ← `pkg/fizhi/fizhi_mnc_init.F:197,226,243,272`
- `MNC_CW_ADD_VNAME` ← `pkg/fizhi/fizhi_mnc_init.F:227,228,231,232`
- `MNC_CW_RL_R` ← `pkg/fizhi/fizhi_read_pickup.F:101,102,103,104`
- `MNC_CW_I_W` ← `pkg/fizhi/fizhi_readwrite_vegtiles.F:95,98`
- `MNC_CW_I_W_S` ← `pkg/fizhi/fizhi_readwrite_vegtiles.F:73`
- `MNC_CW_RL_R` ← `pkg/fizhi/fizhi_readwrite_vegtiles.F:339,340,341,342`
- `MNC_CW_RL_W` ← `pkg/fizhi/fizhi_readwrite_vegtiles.F:76,77,78,79`
- `MNC_CW_RL_W_S` ← `pkg/fizhi/fizhi_readwrite_vegtiles.F:72`
- `MNC_CW_SET_UDIM` ← `pkg/fizhi/fizhi_readwrite_vegtiles.F:71,336`
- `MNC_FILE_CLOSE_ALL_MATCHING` ← `pkg/fizhi/fizhi_readwrite_vegtiles.F:335`
- `MNC_CW_ADD_VATTR_TEXT` ← `pkg/gmredi/gmredi_mnc_init.F:30,31,38,39`
- `MNC_CW_ADD_VNAME` ← `pkg/gmredi/gmredi_mnc_init.F:29,37,43,49`
- `MNC_CW_I_W_S` ← `pkg/gmredi/gmredi_output.F:76`
- `MNC_CW_RL_W` ← `pkg/gmredi/gmredi_output.F:77,78,81,82`
- `MNC_CW_RL_W_S` ← `pkg/gmredi/gmredi_output.F:74`
- `MNC_CW_SET_UDIM` ← `pkg/gmredi/gmredi_output.F:73,75`
- `MNC_CW_ADD_VATTR_TEXT` ← `pkg/kpp/kpp_init_fixed.F:49,51,53,58`
- `MNC_CW_ADD_VNAME` ← `pkg/kpp/kpp_init_fixed.F:47,56,66,75`
- `MNC_CW_I_W_S` ← `pkg/kpp/kpp_output.F:171`
- `MNC_CW_RL_W` ← `pkg/kpp/kpp_output.F:172,174,176,178`
- `MNC_CW_RL_W_S` ← `pkg/kpp/kpp_output.F:169`
- `MNC_CW_SET_UDIM` ← `pkg/kpp/kpp_output.F:168,170`
- `MNC_CW_ADD_GNAME` ← `pkg/land/land_mnc_init.F:207`
- `MNC_CW_ADD_VATTR_TEXT` ← `pkg/land/land_mnc_init.F:223,227,231,234`
- `MNC_CW_ADD_VNAME` ← `pkg/land/land_mnc_init.F:221,225,229,233`
- `MNC_CW_APPEND_VNAME` ← `pkg/land/land_monitor.F:97`
- `MNC_CW_I_W_S` ← `pkg/land/land_monitor.F:100`
- `MNC_CW_SET_UDIM` ← `pkg/land/land_monitor.F:99,102`
- `MNC_CW_I_W_S` ← `pkg/land/land_output.F:120`
- `MNC_CW_RL_W` ← `pkg/land/land_output.F:122,124,126,129`
- `MNC_CW_RL_W_S` ← `pkg/land/land_output.F:118`
- `MNC_CW_SET_UDIM` ← `pkg/land/land_output.F:117,119`
- `MNC_CW_RL_R` ← `pkg/land/land_read_pickup.F:91,93,96,98`
- `MNC_CW_SET_CITER` ← `pkg/land/land_read_pickup.F:89`
- `MNC_CW_SET_UDIM` ← `pkg/land/land_read_pickup.F:88`
- `MNC_FILE_CLOSE_ALL_MATCHING` ← `pkg/land/land_read_pickup.F:87`
- `MNC_CW_I_W_S` ← `pkg/land/land_write_pickup.F:98`
- `MNC_CW_RL_W` ← `pkg/land/land_write_pickup.F:100,102,105,107`
- `MNC_CW_RL_W_S` ← `pkg/land/land_write_pickup.F:97`
- `MNC_CW_SET_CITER` ← `pkg/land/land_write_pickup.F:91,93`
- `MNC_CW_SET_UDIM` ← `pkg/land/land_write_pickup.F:89,95`
- `MNC_FILE_CLOSE_ALL_MATCHING` ← `pkg/land/land_write_pickup.F:88`
- `MNC_CW_ADD_VATTR_TEXT` ← `pkg/layers/layers_mnc_init.F:29,30,35,36`
- `MNC_CW_ADD_VNAME` ← `pkg/layers/layers_mnc_init.F:28,34`
- `MNC_CW_I_W_S` ← `pkg/mom_vecinv/mom_vecinv.F:177`
- `MNC_CW_RL_W_OFFSET` ← `pkg/mom_vecinv/mom_vecinv.F:715,717,807,809`
- `MNC_CW_RL_W_S` ← `pkg/mom_vecinv/mom_vecinv.F:175`
- `MNC_CW_SET_UDIM` ← `pkg/mom_vecinv/mom_vecinv.F:174,176`
- `MNC_CW_APPEND_VNAME` ← `pkg/monitor/mon_out.F:212`
- `MNC_CW_I_W` ← `pkg/monitor/mon_out.F:215`
- `MNC_CW_RL_W` ← `pkg/monitor/mon_out.F:218`
- `MNC_CW_APPEND_VNAME` ← `pkg/monitor/monitor.F:66`
- `MNC_CW_RL_W_S` ← `pkg/monitor/monitor.F:69`
- `MNC_CW_SET_UDIM` ← `pkg/monitor/monitor.F:68,71`
- `MNC_CW_APPEND_VNAME` ← `pkg/monitor/monitor_ad.F:81`
- `MNC_CW_RL_W_S` ← `pkg/monitor/monitor_ad.F:84`
- `MNC_CW_SET_UDIM` ← `pkg/monitor/monitor_ad.F:83,86`
- `MNC_CW_APPEND_VNAME` ← `pkg/monitor/monitor_g.F:77`
- `MNC_CW_RL_W_S` ← `pkg/monitor/monitor_g.F:80`
- `MNC_CW_SET_UDIM` ← `pkg/monitor/monitor_g.F:79,82`
- `MNC_CW_ADD_VATTR_TEXT` ← `pkg/mypackage/mypackage_mnc_init.F:31,33,38,40`
- `MNC_CW_ADD_VNAME` ← `pkg/mypackage/mypackage_mnc_init.F:29,36,43,50`
- `MNC_FILE_CLOSE_ALL` ← `pkg/openad/the_model_main.F:334`
- `MNC_CW_ADD_VATTR_TEXT` ← `pkg/ptracers/ptracers_mnc_init.F:60,67`
- `MNC_CW_ADD_VNAME` ← `pkg/ptracers/ptracers_mnc_init.F:52,55`
- `MNC_CW_APPEND_VNAME` ← `pkg/ptracers/ptracers_monitor.F:71`
- `MNC_CW_RL_W_S` ← `pkg/ptracers/ptracers_monitor.F:74`
- `MNC_CW_SET_UDIM` ← `pkg/ptracers/ptracers_monitor.F:73,76`
- `MNC_CW_APPEND_VNAME` ← `pkg/ptracers/ptracers_monitor_ad.F:79`
- `MNC_CW_RL_W_S` ← `pkg/ptracers/ptracers_monitor_ad.F:82`
- `MNC_CW_SET_UDIM` ← `pkg/ptracers/ptracers_monitor_ad.F:81,84`
- `MNC_CW_RL_R` ← `pkg/ptracers/ptracers_read_pickup.F:74,81`
- `MNC_CW_SET_CITER` ← `pkg/ptracers/ptracers_read_pickup.F:72`
- `MNC_CW_SET_UDIM` ← `pkg/ptracers/ptracers_read_pickup.F:71,79`
- `MNC_FILE_CLOSE_ALL_MATCHING` ← `pkg/ptracers/ptracers_read_pickup.F:70`
- `MNC_CW_I_W_S` ← `pkg/ptracers/ptracers_write_pickup.F:90,97`
- `MNC_CW_RL_W` ← `pkg/ptracers/ptracers_write_pickup.F:92,99`
- `MNC_CW_RL_W_S` ← `pkg/ptracers/ptracers_write_pickup.F:89,96`
- `MNC_CW_SET_CITER` ← `pkg/ptracers/ptracers_write_pickup.F:81,83`
- `MNC_CW_SET_UDIM` ← `pkg/ptracers/ptracers_write_pickup.F:79,86,95`
- `MNC_FILE_CLOSE_ALL_MATCHING` ← `pkg/ptracers/ptracers_write_pickup.F:77`
- `MNC_CW_I_W_S` ← `pkg/ptracers/ptracers_write_state.F:70`
- `MNC_CW_RL_W` ← `pkg/ptracers/ptracers_write_state.F:72`
- `MNC_CW_RL_W_S` ← `pkg/ptracers/ptracers_write_state.F:68`
- `MNC_CW_SET_UDIM` ← `pkg/ptracers/ptracers_write_state.F:67,69`
- `MNC_CW_ADD_VATTR_TEXT` ← `pkg/salt_plume/salt_plume_mnc_init.F:29,31`
- `MNC_CW_ADD_VNAME` ← `pkg/salt_plume/salt_plume_mnc_init.F:27`
- `MNC_CW_I_W_S` ← `pkg/seaice/seaice_ad_dump.F:148`
- `MNC_CW_RL_W` ← `pkg/seaice/seaice_ad_dump.F:154,157,160,166`
- `MNC_CW_RL_W_S` ← `pkg/seaice/seaice_ad_dump.F:146,149`
- `MNC_CW_SET_UDIM` ← `pkg/seaice/seaice_ad_dump.F:145,147`
- `MNC_CW_ADD_VATTR_TEXT` ← `pkg/seaice/seaice_mnc_init.F:33,34,37,38`
- `MNC_CW_ADD_VNAME` ← `pkg/seaice/seaice_mnc_init.F:32,36,40,47`
- `MNC_CW_APPEND_VNAME` ← `pkg/seaice/seaice_monitor.F:74`
- `MNC_CW_RL_W_S` ← `pkg/seaice/seaice_monitor.F:77`
- `MNC_CW_SET_UDIM` ← `pkg/seaice/seaice_monitor.F:76,79`
- `MNC_CW_APPEND_VNAME` ← `pkg/seaice/seaice_monitor_ad.F:80`
- `MNC_CW_RL_W_S` ← `pkg/seaice/seaice_monitor_ad.F:83`
- `MNC_CW_SET_UDIM` ← `pkg/seaice/seaice_monitor_ad.F:82,85`
- `MNC_CW_I_W_S` ← `pkg/seaice/seaice_output.F:74`
- `MNC_CW_RL_W` ← `pkg/seaice/seaice_output.F:79,81,83,87`
- `MNC_CW_RL_W_S` ← `pkg/seaice/seaice_output.F:72,75`
- `MNC_CW_RS_W` ← `pkg/seaice/seaice_output.F:101,102,103,104`
- `MNC_CW_SET_UDIM` ← `pkg/seaice/seaice_output.F:71,73`
- `MNC_CW_ADD_VATTR_TEXT` ← `pkg/shelfice/shelfice_mnc_init.F:28,29,33,34`
- `MNC_CW_ADD_VNAME` ← `pkg/shelfice/shelfice_mnc_init.F:27,32`
- `MNC_CW_I_W_S` ← `pkg/shelfice/shelfice_output.F:63`
- `MNC_CW_RL_W_S` ← `pkg/shelfice/shelfice_output.F:61,64`
- `MNC_CW_RS_W` ← `pkg/shelfice/shelfice_output.F:66,68`
- `MNC_CW_SET_UDIM` ← `pkg/shelfice/shelfice_output.F:60,62`
- `MNC_CW_ADD_VATTR_TEXT` ← `pkg/thsice/thsice_mnc_init.F:25,28,31,34`
- `MNC_CW_ADD_VNAME` ← `pkg/thsice/thsice_mnc_init.F:24,27,30,33`
- `MNC_CW_I_W_S` ← `pkg/thsice/thsice_monitor.F:88`
- `MNC_CW_RL_W_S` ← `pkg/thsice/thsice_monitor.F:91`
- `MNC_CW_SET_UDIM` ← `pkg/thsice/thsice_monitor.F:87,90`
- `MNC_CW_I_W_S` ← `pkg/thsice/thsice_output.F:117`
- `MNC_CW_RL_W` ← `pkg/thsice/thsice_output.F:120,121,122,123`
- `MNC_CW_RL_W_S` ← `pkg/thsice/thsice_output.F:119`
- `MNC_CW_SET_UDIM` ← `pkg/thsice/thsice_output.F:116,118`
- `MNC_CW_RL_R` ← `pkg/thsice/thsice_read_pickup.F:74,75,76,77`
- `MNC_CW_SET_CITER` ← `pkg/thsice/thsice_read_pickup.F:73`
- `MNC_CW_SET_UDIM` ← `pkg/thsice/thsice_read_pickup.F:72`
- `MNC_FILE_CLOSE_ALL_MATCHING` ← `pkg/thsice/thsice_read_pickup.F:71`
- `MNC_CW_I_W_S` ← `pkg/thsice/thsice_write_pickup.F:86`
- `MNC_CW_RL_W` ← `pkg/thsice/thsice_write_pickup.F:88,89,90,91`
- `MNC_CW_SET_CITER` ← `pkg/thsice/thsice_write_pickup.F:80,82`
- `MNC_CW_SET_UDIM` ← `pkg/thsice/thsice_write_pickup.F:78,85`

## Verification experiments compiling it (21)
`MLAdjust` `aim.5l_cs` `cheapAML_box` `fizhi-cs-aqualev20` `global_oce_biogeo_bling` `global_oce_latlon` `global_ocean.90x40x15` `global_ocean.cs32x15` `hs94.1x64x5` `internal_wave` `inverted_barometer` `isomip` `lab_sea` `tutorial_advection_in_gyre` `tutorial_baroclinic_gyre` `tutorial_dic_adjoffline` `tutorial_global_oce_biogeo` `tutorial_global_oce_latlon` `tutorial_held_suarez_cs` `tutorial_rotating_tank` `vermix`
