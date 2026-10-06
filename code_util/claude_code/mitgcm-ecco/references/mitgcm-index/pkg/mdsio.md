# pkg/mdsio

MDS (.data/.meta) binary I/O routines; always compiled.

**in groups:** gfd
**always-on utility package** (no data.pkg switch)
**manual:** `doc/examples/reentrant_channel/reentrant_channel.rst`, `doc/getting_started/getting_started.rst`, `doc/ocean_state_est/ocean_state_est.rst`, `doc/outp_pkgs/outp_pkgs.rst`, `doc/phys_pkgs/seaice.rst`

## CPP options (defaults as shipped)
- `SAFE_IO` (undef, MDSIO_OPTIONS.h) — -  Defining SAFE_IO stops the model from overwriting its own files
- `_NEW_STATUS` (define, MDSIO_OPTIONS.h)
- `_NEW_STATUS` (define, MDSIO_OPTIONS.h)
- `_OLD_STATUS` (define, MDSIO_OPTIONS.h)
- `_OLD_STATUS` (define, MDSIO_OPTIONS.h)
- `ALLOW_WHIO` (define, MDSIO_OPTIONS.h) — -  I/O that includes tile halos in the files Only used when pkg/autodiff is compiled:
- `ALLOW_WHIO_3D` (define, MDSIO_OPTIONS.h)
- `EXCLUDE_WHIO_GLOBUFF_2D` (undef, MDSIO_OPTIONS.h)
- `INCLUDE_WHIO_GLOBUFF_3D` (undef, MDSIO_OPTIONS.h)

## Headers
- `MDSIO_BUFF_3D.h` — BOP
- `MDSIO_BUFF_WH.h` — BOP
- `MDSIO_OPTIONS.h` — BOP

## Routines (56)
`mdsio_buffertorl.F`, `mdsio_buffertors.F`, `mdsio_check4file.F`, `mdsio_facef_read.F`, `mdsio_gl.F`, `mdsio_gl_slice.F`, `mdsio_pass_r4torl.F`, `mdsio_pass_r4tors.F`, `mdsio_pass_r8torl.F`, `mdsio_pass_r8tors.F`, `mdsio_rd_rec_rl.F`, `mdsio_rd_rec_rs.F`, `mdsio_read_field.F`, `mdsio_read_meta.F`, `mdsio_read_section.F`, `mdsio_read_tape.F`, `mdsio_read_whalos.F`, `mdsio_readvec_loc.F`, `mdsio_rw_field.F`, `mdsio_rw_slice.F`, `mdsio_seg4torl.F`, `mdsio_seg4tors.F`, `mdsio_seg8torl.F`, `mdsio_seg8tors.F`, `mdsio_segxtorx_2d.F`, `mdsio_wr_metafiles.F`, `mdsio_wr_rec_rl.F`, `mdsio_wr_rec_rs.F`, `mdsio_write_field.F`, `mdsio_write_meta.F`, `mdsio_write_section.F`, `mdsio_write_tape.F`, `mdsio_write_whalos.F`, `mdsio_writelocal.F`, `mdsio_writevec_loc.F`

## Called from outside the package
- `MDS_FACEF_READ_RS` ← `model/src/ini_cori.F:164`
- `MDS_FACEF_READ_RS` ← `model/src/ini_curvilinear_grid.F:294,297,300,303`
- `MDS_WR_METAFILES` ← `model/src/write_pickup.F:380`
- `MDS_READVEC_LOC` ← `pkg/aim_v23/aim_do_co2.F:112`
- `MDS_WRITEVEC_LOC` ← `pkg/aim_v23/aim_do_co2.F:170`
- `MDS_WR_METAFILES` ← `pkg/atm_compon_interf/cpl_write_pickup.F:166`
- `MDS_WR_METAFILES` ← `pkg/atm_phys/atm_phys_write_pickup.F:107`
- `MDS_READVEC_LOC` ← `pkg/autodiff/active_file_control.F:467,501,533,633`
- `MDS_READ_FIELD` ← `pkg/autodiff/active_file_control.F:108,145,191,285`
- `MDS_WRITEVEC_LOC` ← `pkg/autodiff/active_file_control.F:486,515,652,681`
- `MDS_WRITE_FIELD` ← `pkg/autodiff/active_file_control.F:132,167,309,344`
- `MDS_READ_SEC_XZ` ← `pkg/autodiff/active_file_control_slice.F:104,139,181,273`
- `MDS_READ_SEC_YZ` ← `pkg/autodiff/active_file_control_slice.F:441,476,518,610`
- `MDS_WRITE_SEC_XZ` ← `pkg/autodiff/active_file_control_slice.F:126,159,295,328`
- `MDS_WRITE_SEC_YZ` ← `pkg/autodiff/active_file_control_slice.F:463,496,632,665`
- `MDS_READ_TAPE` ← `pkg/autodiff/adread_adwrite.F:231,237`
- `MDS_READ_WHALOS` ← `pkg/autodiff/adread_adwrite.F:207`
- `MDS_WRITE_TAPE` ← `pkg/autodiff/adread_adwrite.F:446,452`
- `MDS_WRITE_WHALOS` ← `pkg/autodiff/adread_adwrite.F:424`
- `MDS_WR_METAFILES` ← `pkg/bbl/bbl_write_pickup.F:102`
- `MDS_CHECK4FILE` ← `pkg/bling/bling_read_pickup.F:58`
- `MDS_WR_METAFILES` ← `pkg/bling/bling_write_pickup.F:102`
- `MDS_WR_METAFILES` ← `pkg/cheapaml/cheapaml_write_pickup.F:136`
- `MDSREADFIELD_2D_GL` ← `pkg/ctrl/ctrl_set_pack_xy.F:215`
- `MDSREADFIELD_3D_GL` ← `pkg/ctrl/ctrl_set_pack_xy.F:141,148`
- `MDS_PASS_R8TORL` ← `pkg/ctrl/ctrl_set_pack_xy.F:363,370`
- `MDSREADFIELD_3D_GL` ← `pkg/ctrl/ctrl_set_pack_xyz.F:163,185,191`
- `MDS_PASS_R8TORL` ← `pkg/ctrl/ctrl_set_pack_xyz.F:344,351`
- `MDSREADFIELD_3D_GL` ← `pkg/ctrl/ctrl_set_pack_xz.F:181`
- `MDSREADFIELD_XZ_GL` ← `pkg/ctrl/ctrl_set_pack_xz.F:161,167,258`
- `MDSREADFIELD_3D_GL` ← `pkg/ctrl/ctrl_set_pack_yz.F:183`
- `MDSREADFIELD_YZ_GL` ← `pkg/ctrl/ctrl_set_pack_yz.F:161,167,262`
- `MDSREADFIELD_3D_GL` ← `pkg/ctrl/ctrl_set_unpack_xy.F:150`
- `MDSWRITEFIELD_2D_GL` ← `pkg/ctrl/ctrl_set_unpack_xy.F:337`
- `MDSWRITEFIELD_3D_GL` ← `pkg/ctrl/ctrl_set_unpack_xy.F:243`
- `MDS_PASS_R8TORL` ← `pkg/ctrl/ctrl_set_unpack_xy.F:399,446`
- `MDSREADFIELD_3D_GL` ← `pkg/ctrl/ctrl_set_unpack_xyz.F:177,199`
- `MDSWRITEFIELD_3D_GL` ← `pkg/ctrl/ctrl_set_unpack_xyz.F:285`
- `MDS_PASS_R8TORL` ← `pkg/ctrl/ctrl_set_unpack_xyz.F:356,402`
- `MDSREADFIELD_XZ_GL` ← `pkg/ctrl/ctrl_set_unpack_xz.F:173,179`
- `MDSWRITEFIELD_3D_GL` ← `pkg/ctrl/ctrl_set_unpack_xz.F:275`
- `MDSWRITEFIELD_XZ_GL` ← `pkg/ctrl/ctrl_set_unpack_xz.F:358`
- `MDSREADFIELD_YZ_GL` ← `pkg/ctrl/ctrl_set_unpack_yz.F:173,179`
- `MDSWRITEFIELD_3D_GL` ← `pkg/ctrl/ctrl_set_unpack_yz.F:275`
- `MDSWRITEFIELD_YZ_GL` ← `pkg/ctrl/ctrl_set_unpack_yz.F:352`
- `MDS_WR_METAFILES` ← `pkg/diagnostics/diagnostics_out.F:444`
- `MDS_CHECK4FILE` ← `pkg/dic/dic_read_co2_pickup.F:53`
- `MDS_CHECK4FILE` ← `pkg/dic/dic_read_pickup.F:115`
- `MDS_WRITEVEC_LOC` ← `pkg/dic/dic_write_pickup.F:68`
- `MDS_WR_METAFILES` ← `pkg/dic/dic_write_pickup.F:121`
- `MDS_READVEC_LOC` ← `pkg/ecco/cost_gencost_boxmean.F:126`
- `MDS_READVEC_LOC` ← `pkg/ecco/cost_gencost_moc.F:178`
- `MDS_READVEC_LOC` ← `pkg/ecco/ecco_read_pickup.F:78`
- `MDS_READVEC_LOC` ← `pkg/ecco/ecco_readparms.F:821`
- `MDS_WRITEVEC_LOC` ← `pkg/ecco/ecco_write_pickup.F:54`
- `MDS_WRITEVEC_LOC` ← `pkg/ecco/stergloh_output.F:89`
- `MDS_READVEC_LOC` ← `pkg/flt/flt_init_varia.F:122,137,184,197`
- `MDS_READVEC_LOC` ← `pkg/flt/flt_traj.F:176`
- `MDS_WRITEVEC_LOC` ← `pkg/flt/flt_traj.F:224,232`
- `MDS_READVEC_LOC` ← `pkg/flt/flt_up.F:174`
- `MDS_WRITEVEC_LOC` ← `pkg/flt/flt_up.F:222,230`
- `MDS_WRITEVEC_LOC` ← `pkg/flt/flt_write_pickup.F:72,91`
- `MDS_CHECK4FILE` ← `pkg/generic_advdiff/gad_read_pickup.F:85,136`
- `MDS_CHECK4FILE` ← `pkg/gmredi/gmredi_read_pickup.F:119`
- `MDS_WR_METAFILES` ← `pkg/gmredi/gmredi_write_pickup.F:207`
- `MDS_WR_METAFILES` ← `pkg/mypackage/mypackage_write_pickup.F:134`
- `MDS_READ_SEC_XZ` ← `pkg/obcs/obcs_init_fixed.F:461,465,469,473`
- `MDS_READ_SEC_YZ` ← `pkg/obcs/obcs_init_fixed.F:497,501,505,509`
- `MDS_CHECK4FILE` ← `pkg/ptracers/ptracers_read_pickup.F:303`
- `MDS_WR_METAFILES` ← `pkg/ptracers/ptracers_write_pickup.F:171`
- `MDS_READ_FIELD` ← `pkg/rw/read_fld_xy_rl.F:45`
- `MDS_READ_FIELD` ← `pkg/rw/read_fld_xy_rs.F:45`
- `MDS_READ_FIELD` ← `pkg/rw/read_fld_xyz_rl.F:45`
- `MDS_READ_FIELD` ← `pkg/rw/read_fld_xyz_rs.F:45`
- `MDS_READVEC_LOC` ← `pkg/rw/read_glvec_rl.F:61`
- `MDS_READVEC_LOC` ← `pkg/rw/read_glvec_rs.F:61`
- `MDS_READ_FIELD` ← `pkg/rw/read_mflds.F:348,477,606`
- `MDS_READ_META` ← `pkg/rw/read_mflds.F:151`
- `MDS_READ_FIELD` ← `pkg/rw/read_rec.F:65,124,183,242`
- `MDS_READ_SEC_XZ` ← `pkg/rw/read_rec.F:557,619`
- `MDS_READ_SEC_YZ` ← `pkg/rw/read_rec.F:681,743`
- `MDS_WRITE_FIELD` ← `pkg/rw/write_fld_3d_rl.F:50`
- `MDS_WRITE_FIELD` ← `pkg/rw/write_fld_3d_rs.F:50`
- `MDS_WRITE_FIELD` ← `pkg/rw/write_fld_xy_rl.F:50`
- `MDS_WRITE_FIELD` ← `pkg/rw/write_fld_xy_rs.F:50`
- `MDS_WRITE_FIELD` ← `pkg/rw/write_fld_xyz_rl.F:50`
- `MDS_WRITE_FIELD` ← `pkg/rw/write_fld_xyz_rs.F:50`
- `MDS_WRITEVEC_LOC` ← `pkg/rw/write_glvec_rl.F:60`
- `MDS_WRITEVEC_LOC` ← `pkg/rw/write_glvec_rs.F:60`
- `MDS_WRITELOCAL` ← `pkg/rw/write_local_rl.F:85`
- `MDS_WRITE_FIELD` ← `pkg/rw/write_local_rl.F:79`
- `MDS_WRITELOCAL` ← `pkg/rw/write_local_rs.F:85`
- `MDS_WRITE_FIELD` ← `pkg/rw/write_local_rs.F:79`
- `MDS_WRITE_FIELD` ← `pkg/rw/write_rec.F:132,195,258,321`
- `MDS_WRITE_SEC_XZ` ← `pkg/rw/write_rec.F:648,714`
- `MDS_WRITE_SEC_YZ` ← `pkg/rw/write_rec.F:780,846`
- `MDS_WRITEVEC_LOC` ← `pkg/sbo/sbo_output.F:96`
- `MDS_WR_METAFILES` ← `pkg/seaice/seaice_write_pickup.F:208`
- `MDS_WR_METAFILES` ← `pkg/shelfice/shelfice_write_pickup.F:120`
- `MDS_WR_METAFILES` ← `pkg/streamice/streamice_write_pickup.F:211`

## Verification experiments compiling it (58)
`1D_ocean_ice_column` `MLAdjust` `adjustment.cs-32x32x1` `advect_cs` `advect_xz` `aim.5l_Equatorial_Channel` `aim.5l_LatLon` `aim.5l_cs` `atm_gray` `bottom_ctrl_5x5` `cfc_example` `cheapAML_box` `cpl_aim+ocn` `deep_anelastic` `dome` `exp2` `exp4` `fizhi-cs-32x32x40` `fizhi-cs-aqualev20` `fizhi-gridalt-hs` `front_relax` `global_oce_biogeo_bling` `global_oce_latlon` `global_ocean.90x40x15` `global_ocean.cs32x15` `halfpipe_streamice` `hs94.128x64x5` `hs94.1x64x5` `hs94.cs-32x32x5` `ideal_2D_oce` `internal_wave` `inverted_barometer` `isomip` `lab_sea` `matrix_example` `obcs_ctrl` `offline_exf_seaice` `seaice_itd` `seaice_obcs` `shelfice_2d_remesh` `short_surf_wave` `so_box_biogeo` `solid-body.cs-32x32x1` `tutorial_advection_in_gyre` `tutorial_baroclinic_gyre` `tutorial_cfc_offline` `tutorial_deep_convection` `tutorial_dic_adjoffline` `tutorial_global_oce_biogeo` `tutorial_global_oce_in_p` `tutorial_global_oce_latlon` `tutorial_global_oce_optim` `tutorial_held_suarez_cs` `tutorial_plume_on_slope` `tutorial_reentrant_channel` `tutorial_rotating_tank` `tutorial_tracer_adjsens` `vermix`
