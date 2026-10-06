# pkg/mom_vecinv  (wad: ~/Documents/research/ECCO/wetting_drying/MITgcm)

S/R MOM_VECINV Form the right hand-side of the momentum equation. Terms are evaluated one layer at a time working from the bottom to the top. The vertically integrated barotropic flow tendency term is evluated by summing the

**vs its upstream base (merge-base, see README):** changed: mom_vecinv.F

**pkg_depend:** +mom_common  (`+` requires, `-` excludes)
**in groups:** gfd
**runtime switch:** `useMOM_VECINV`-style flag in `data.pkg` (check exact name in packages_boot.F)
**manual:** `doc/algorithm/algorithm.rst`, `doc/examples/held_suarez_cs/held_suarez_cs.rst`, `doc/getting_started/getting_started.rst`, `doc/outp_pkgs/outp_pkgs.rst`
**adjoint support files:** mom_vecinv_ad_diff.list

## CPP options (defaults as shipped)
- `MOM_VI_ORIGINAL_VISCA4` (undef, MOM_VECINV_OPTIONS.h) — use the original discretization (not recommended) for biharmonic viscosity that was in mom_vi_hdissip.F, version 1.1.2.1

## Headers
- `MOM_VECINV_OPTIONS.h` — CPP options file for mom_vecinv package Use this file for selecting CPP options within the mom_vecinv package

## Routines (12)
`mom_vecinv.F`, `mom_vi_coriolis.F`, `mom_vi_del2uv.F`, `mom_vi_hdissip.F`, `mom_vi_u_coriolis.F`, `mom_vi_u_coriolis_c4.F`, `mom_vi_u_grad_ke.F`, `mom_vi_u_vertshear.F`, `mom_vi_v_coriolis.F`, `mom_vi_v_coriolis_c4.F`, `mom_vi_v_grad_ke.F`, `mom_vi_v_vertshear.F`

## Called from outside the package
- `MOM_VECINV` ← `model/src/dynamics.F:527`
- `MOM_VI_U_CORIOLIS` ← `pkg/seaice/seaice_mom_advection.F:135`
- `MOM_VI_U_CORIOLIS_C4` ← `pkg/seaice/seaice_mom_advection.F:129`
- `MOM_VI_U_GRAD_KE` ← `pkg/seaice/seaice_mom_advection.F:171`
- `MOM_VI_V_CORIOLIS` ← `pkg/seaice/seaice_mom_advection.F:152`
- `MOM_VI_V_CORIOLIS_C4` ← `pkg/seaice/seaice_mom_advection.F:146`
- `MOM_VI_V_GRAD_KE` ← `pkg/seaice/seaice_mom_advection.F:177`
- `MOM_VI_DEL2UV` ← `pkg/shap_filt/shap_filt_uv_s2.F:178`

## Verification experiments compiling it (8)
`wad_balzano@wad` `wad_estuary_3d@wad` `wad_flat_xz@wad` `wad_iceground@wad` `wad_mangrove@wad` `wad_mudflat@wad` `wad_overflood@wad` `wad_thacker_1d@wad`
