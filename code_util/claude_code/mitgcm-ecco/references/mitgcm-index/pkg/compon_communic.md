# pkg/compon_communic

MPI component communication layer for coupled (multi-executable) runs.

**runtime switch:** `useCOMPON_COMMUNIC`-style flag in `data.pkg` (check exact name in packages_boot.F)

## Headers
- `CPLR_SIG.h` — Special meanings/handles

## Routines (40)
`comprecv_i4vec.F`, `comprecv_r4.F`, `comprecv_r4tiles.F`, `comprecv_r8.F`, `comprecv_r8tiles.F`, `compsend_i4vec.F`, `compsend_r4.F`, `compsend_r4tiles.F`, `compsend_r8.F`, `compsend_r8tiles.F`, `couprecv_i4vec.F`, `couprecv_r4.F`, `couprecv_r4tiles.F`, `couprecv_r8.F`, `couprecv_r8tiles.F`, `coupsend_i4vec.F`, `coupsend_r4.F`, `coupsend_r4tiles.F`, `coupsend_r8.F`, `coupsend_r8tiles.F`, `generate_tag.F`, `mitcomponent_init.F`, `mitcomponent_register.F`, `mitcomponent_tile_register.F`, `mitcoupler_init.F`, `mitcoupler_register.F`, `mitcoupler_tile_register.F`, `mitcplr_all_check.F`, `mitcplr_char2dbl.F`, `mitcplr_char2int.F`, `mitcplr_char2real.F`, `mitcplr_dbl2char.F`, `mitcplr_init1.F`, `mitcplr_init2a.F`, `mitcplr_init2b.F`, `mitcplr_initcomp.F`, `mitcplr_int2char.F`, `mitcplr_match_comp.F`, `mitcplr_real2char.F`, `mitcplr_sortranks.F`

## Called from outside the package
- `COUPRECV_R8TILES` ← `pkg/atm2d/cpl_recv_ocn_fields.F:28,33,38,43`
- `COUPRECV_I4VEC` ← `pkg/atm2d/cpl_recv_ocn_ocnconfig.F:35`
- `COUPRECV_R8TILES` ← `pkg/atm2d/cpl_recv_ocn_ocnconfig.F:41`
- `MITCOUPLER_TILE_REGISTER` ← `pkg/atm2d/cpl_register_ocn.F:37`
- `COUPSEND_R8TILES` ← `pkg/atm2d/cpl_send_ocn_atmconfig.F:47`
- `COUPSEND_R8TILES` ← `pkg/atm2d/cpl_send_ocn_fields.F:27,31,35,39`
- `MITCPLR_ALL_CHECK` ← `pkg/atm2d/exch_component_configs.F:62`
- `MITCOUPLER_INIT` ← `pkg/atm2d/initialise.F:41,45`
- `COMPSEND_I4VEC` ← `pkg/atm_compon_interf/atm_export_atmconfig.F:56`
- `COMPSEND_R8TILES` ← `pkg/atm_compon_interf/atm_export_atmconfig.F:59`
- `COMPSEND_R8TILES` ← `pkg/atm_compon_interf/atm_export_fld.F:73`
- `COMPRECV_R8TILES` ← `pkg/atm_compon_interf/atm_import_fields.F:49,54,59,64`
- `COMPRECV_R8TILES` ← `pkg/atm_compon_interf/atm_import_ocnconfig.F:55`
- `MITCPLR_ALL_CHECK` ← `pkg/atm_compon_interf/cpl_exch_configs.F:78`
- `COMPRECV_I4VEC` ← `pkg/atm_compon_interf/cpl_import_cplparms.F:55`
- `MITCOMPONENT_INIT` ← `pkg/atm_compon_interf/cpl_init.F:38`
- `MITCOMPONENT_TILE_REGISTER` ← `pkg/atm_compon_interf/cpl_register.F:114`
- `COUPRECV_I4VEC` ← `pkg/atm_ocn_coupler/cpl_recv_atm_atmconfig.F:33`
- `COUPRECV_R8TILES` ← `pkg/atm_ocn_coupler/cpl_recv_atm_atmconfig.F:39`
- `COUPRECV_R8TILES` ← `pkg/atm_ocn_coupler/cpl_recv_atm_fields.F:36,41,46,51`
- `COUPRECV_R8TILES` ← `pkg/atm_ocn_coupler/cpl_recv_ocn_fields.F:37,42,47,52`
- `COUPRECV_I4VEC` ← `pkg/atm_ocn_coupler/cpl_recv_ocn_ocnconfig.F:34`
- `COUPRECV_R8TILES` ← `pkg/atm_ocn_coupler/cpl_recv_ocn_ocnconfig.F:40`
- `MITCOUPLER_TILE_REGISTER` ← `pkg/atm_ocn_coupler/cpl_register_atm.F:32`
- `MITCOUPLER_TILE_REGISTER` ← `pkg/atm_ocn_coupler/cpl_register_ocn.F:36`
- `COUPSEND_I4VEC` ← `pkg/atm_ocn_coupler/cpl_send_atm_cplparms.F:58`
- `COUPSEND_R8TILES` ← `pkg/atm_ocn_coupler/cpl_send_atm_fields.F:43,51,59,67`
- `COUPSEND_R8TILES` ← `pkg/atm_ocn_coupler/cpl_send_atm_ocnconfig.F:40`
- `COUPSEND_R8TILES` ← `pkg/atm_ocn_coupler/cpl_send_ocn_atmconfig.F:41`
- `COUPSEND_I4VEC` ← `pkg/atm_ocn_coupler/cpl_send_ocn_cplparms.F:58`
- `COUPSEND_R8TILES` ← `pkg/atm_ocn_coupler/cpl_send_ocn_fields.F:43,59,67,75`
- `MITCPLR_ALL_CHECK` ← `pkg/atm_ocn_coupler/exch_component_configs.F:61`
- `MITCOUPLER_INIT` ← `pkg/atm_ocn_coupler/initialise.F:49`
- `MITCPLR_ALL_CHECK` ← `pkg/ocn_compon_interf/cpl_exch_configs.F:69`
- `COMPRECV_I4VEC` ← `pkg/ocn_compon_interf/cpl_import_cplparms.F:55`
- `MITCOMPONENT_INIT` ← `pkg/ocn_compon_interf/cpl_init.F:38`
- `MITCOMPONENT_TILE_REGISTER` ← `pkg/ocn_compon_interf/cpl_register.F:114`
- `COMPSEND_R8TILES` ← `pkg/ocn_compon_interf/ocn_export_fields.F:49,53,57,61`
- `COMPSEND_I4VEC` ← `pkg/ocn_compon_interf/ocn_export_ocnconfig.F:54`
- `COMPSEND_R8TILES` ← `pkg/ocn_compon_interf/ocn_export_ocnconfig.F:57`
- `COMPRECV_R8TILES` ← `pkg/ocn_compon_interf/ocn_import_atmconfig.F:54`
- `COMPRECV_R8TILES` ← `pkg/ocn_compon_interf/ocn_import_fields.F:47,52,57,62`

## Verification experiments compiling it (1)
`cpl_aim+ocn`
