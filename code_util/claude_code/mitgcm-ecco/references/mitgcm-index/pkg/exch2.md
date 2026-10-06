# pkg/exch2

Generalised tile exchanges for cubed-sphere and LLC grids (facets, blank tiles via data.exch2, W2_mapIO).

**runtime switch:** `useEXCH2`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.exch2`
**manual:** `doc/phys_pkgs/exch2.rst`, `doc/examples/held_suarez_cs/held_suarez_cs.rst`, `doc/getting_started/getting_started.rst`, `doc/software_arch/software_arch.rst`
**adjoint support files:** exch2_ad_diff.list

## Namelist parameters
### W2_EXCH2_PARM01
- `preDefTopol` — pre-defined Topology selector:
- `dimsFacets`
- `facetEdgeLink`
- `blankList` — List of "Blank-Tiles" (non active)
- `W2_mapIO` — select option for global-IO mapping:
- `W2_printMsg` — select option for information messages printing
- `W2_useE2ioLayOut` — =T: use Exch2 global-IO Layout; =F: use model default

## CPP options (defaults as shipped)
- `W2_USE_E2_SAFEMODE` (define, W2_OPTIONS.h) — ... W2_USE_E2_SAFEMODE description ...
- `W2_E2_DEBUG_ON` (undef, W2_OPTIONS.h) — Debug mode option:
- `W2_USE_R1_ONLY` (undef, W2_OPTIONS.h) — Use only exch2_R1_cube (and avoid calling exch2_R2_cube)
- `W2_FILL_NULL_REGIONS` (undef, W2_OPTIONS.h) — Fill null regions (face-corner halo regions) with e2FillValue_RX (=0) notes: for testing (allow to check that results are not affected)
- `W2_CUMSUM_USE_MATRIX` (undef, W2_OPTIONS.h) — Process Global Cumulated-Sum using a Tile x Tile (x 2) Matrix notes: should be faster (vectorise) but storage of this matrix might become an issue on large set-up (with many tiles)

## Headers
- `W2_EXCH2_BUFFER.h` — BOP
- `W2_EXCH2_PARAMS.h` — BOP
- `W2_EXCH2_SIZE.h` — BOP
- `W2_EXCH2_TOPOLOGY.h` — BOP
- `W2_OPTIONS.h` — CPP options file for EXCH2 package

## Routines (127)
`exch2_3d_r4.F`, `exch2_3d_r8.F`, `exch2_3d_rl.F`, `exch2_3d_rs.F`, `exch2_ad_get_r41.F`, `exch2_ad_get_r42.F`, `exch2_ad_get_r81.F`, `exch2_ad_get_r82.F`, `exch2_ad_get_rl1.F`, `exch2_ad_get_rl2.F`, `exch2_ad_get_rs1.F`, `exch2_ad_get_rs2.F`, `exch2_ad_put_r41.F`, `exch2_ad_put_r42.F`, `exch2_ad_put_r81.F`, `exch2_ad_put_r82.F`, `exch2_ad_put_rl1.F`, `exch2_ad_put_rl2.F`, `exch2_ad_put_rs1.F`, `exch2_ad_put_rs2.F`, `exch2_check_depths.F`, `exch2_get_r41.F`, `exch2_get_r42.F`, `exch2_get_r81.F`, `exch2_get_r82.F`, `exch2_get_rl1.F`, `exch2_get_rl2.F`, `exch2_get_rs1.F`, `exch2_get_rs2.F`, `exch2_get_scal_bounds.F`, `exch2_get_uv_bounds.F`, `exch2_put_r41.F`, `exch2_put_r42.F`, `exch2_put_r81.F`, `exch2_put_r82.F`, `exch2_put_rl1.F`, `exch2_put_rl2.F`, `exch2_put_rs1.F`, `exch2_put_rs2.F`, `exch2_r41_cube.F`, `exch2_r41_cube_ad.F`, `exch2_r42_cube.F`, `exch2_r42_cube_ad.F`, `exch2_r81_cube.F`, `exch2_r81_cube_ad.F`, `exch2_r82_cube.F`, `exch2_r82_cube_ad.F`, `exch2_recv_r41.F`, `exch2_recv_r42.F`, `exch2_recv_r81.F`, `exch2_recv_r82.F`, `exch2_recv_rl1.F`, `exch2_recv_rl2.F`, `exch2_recv_rs1.F`, `exch2_recv_rs2.F`, `exch2_rl1_cube.F`, `exch2_rl1_cube_ad.F`, `exch2_rl1_cube_b.F`, `exch2_rl2_cube.F`, `exch2_rl2_cube_ad.F`, `exch2_rl2_cube_b.F`, `exch2_rs1_cube.F`, `exch2_rs1_cube_ad.F`, `exch2_rs1_cube_b.F`, `exch2_rs2_cube.F`, `exch2_rs2_cube_ad.F`, `exch2_rs2_cube_b.F`, `exch2_rs_rl_12_d.F`, `exch2_s3d_r4.F`, `exch2_s3d_r8.F`, `exch2_s3d_rl.F`, `exch2_s3d_rs.F`, `exch2_send_r41.F`, `exch2_send_r42.F`, `exch2_send_r81.F`, `exch2_send_r82.F`, `exch2_send_rl1.F`, `exch2_send_rl2.F`, `exch2_send_rs1.F`, `exch2_send_rs2.F`, `exch2_sm_3d_r4.F`, `exch2_sm_3d_r8.F`, `exch2_sm_3d_rl.F`, `exch2_sm_3d_rs.F`, `exch2_uv_3d_r4.F`, `exch2_uv_3d_r8.F`, `exch2_uv_3d_rl.F`, `exch2_uv_3d_rs.F`, `exch2_uv_agrid_3d_r4.F`, `exch2_uv_agrid_3d_r8.F`, `exch2_uv_agrid_3d_rl.F`, `exch2_uv_agrid_3d_rs.F`, `exch2_uv_bgrid_3d_r4.F`, `exch2_uv_bgrid_3d_r8.F`, `exch2_uv_bgrid_3d_rl.F`, `exch2_uv_bgrid_3d_rs.F`, `exch2_uv_cgrid_3d_r4.F`, `exch2_uv_cgrid_3d_r8.F`, `exch2_uv_cgrid_3d_rl.F`, `exch2_uv_cgrid_3d_rs.F`, `exch2_uv_dgrid_3d_r4.F`, `exch2_uv_dgrid_3d_r8.F`, `exch2_uv_dgrid_3d_rl.F`, `exch2_uv_dgrid_3d_rs.F`, `exch2_z_3d_r4.F`, `exch2_z_3d_r8.F`, `exch2_z_3d_rl.F`, `exch2_z_3d_rs.F`, `w2_cumulsum_z_tile.F`, `w2_e2setup.F`, `w2_eeboot.F`, `w2_map_procs.F`, `w2_print_comm_sequence.F`, `w2_print_e2setup.F`, `w2_readparms.F`, `w2_set_cs6_facets.F`, `w2_set_f2f_index.F`, `w2_set_gen_facets.F`, `w2_set_map_cumsum.F`, `w2_set_map_tiles.F`, `w2_set_myown_facets.F`, `w2_set_single_facet.F`, `w2_set_tile2tiles.F`

## Called from outside the package
- `W2_CUMULSUM_Z_TILE_RL` ← `eesupp/src/cumulsum_z_tile.F:71`
- `W2_EEBOOT` ← `eesupp/src/eeboot.F:154`
- `EXCH2_CHECK_DEPTHS` ← `model/src/ini_depths.F:294`

## Verification experiments compiling it (15)
`MLAdjust` `adjustment.cs-32x32x1` `advect_cs` `aim.5l_cs` `atm_gray` `cpl_aim+ocn` `fizhi-cs-32x32x40` `fizhi-cs-aqualev20` `fizhi-gridalt-hs` `global_oce_biogeo_bling` `global_ocean.90x40x15` `global_ocean.cs32x15` `hs94.cs-32x32x5` `solid-body.cs-32x32x1` `tutorial_held_suarez_cs`
