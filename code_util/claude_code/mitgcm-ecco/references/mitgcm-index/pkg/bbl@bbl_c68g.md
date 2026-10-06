# pkg/bbl  (bbl_c68g: ~/Documents/research/ECCO/BBL/MITgcm_c68g)

bbl_wvel    :: default vertical entrainment velocity (m/s) bbl_hvvel   :: default horizontal velocity of BBL (m/s); upper bound of the speed when bbl_useNofSpeed=T bbl_initEta :: default initial thickness of BBL (m)

**vs its upstream base (merge-base, see README):** added: BBL_PTR.h, README.md, bbl_ad_check_lev1_dir.h, bbl_ad_check_lev2_dir.h, bbl_ad_check_lev3_dir.h, bbl_ad_check_lev4_dir.h; changed: BBL.h, BBL_OPTIONS.h, bbl_ad_diff.list, bbl_calc_rhs.F, bbl_check.F, bbl_diagnostics_init.F, bbl_init_varia.F, bbl_read_pickup.F, bbl_readparms.F, bbl_tendency_apply.F, bbl_write_pickup.F

**runtime switch:** `useBBL`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.bbl`
**adjoint support files:** bbl_ad_check_lev1_dir.h, bbl_ad_check_lev2_dir.h, bbl_ad_check_lev3_dir.h, bbl_ad_check_lev4_dir.h, bbl_ad_diff.list

## README
```
Package `BBL` is a simple bottom boundary layer scheme.  The initial
motivation is to allow dense water that forms on the continental shelf around
Antarctica (High Salinity Shelf Water) in the CS510 configuration to sink to
the bottom of the model domain and to become a source of Antarctic Bottom
Water.  The bbl package aims to address the following two limitations of
package down_slope:

(i) In `pkg/down_slope`, dense water cannot flow down-slope unless there is a
step, i.e., a change of vertical level in the bathymetry.  In pkg/bbl, dense
water can flow even on a slight incline or flat bottom.

(ii) In `pkg/down_slope`, dense water is diluted as it flows into grid cells
whose thickness depends on model configuration, typically much thicker than a
bottom boundary layer.  In pkg/bbl, dense water is contained in a thin
sub-layer and hence able to preserve its tracer properties.

Specifically, the bottommost wet grid cell of thickness

    thk = hFacC(kBot) * drF(kBot),

with properties `tracer`, and density `rho` is divided in two sub-levels:

1. A bottom boundary layer with T/S tracer properties `bbl_tracer`,
density `bbl_rho`, and thickness `bbl_eta`.

2. A residual thickness `resThk = thk - bbl_eta` with tracer properties

    resTracer = ( tracer * thk - bbl_tracer * bbl_eta ) / resThk

such that the volume integral of `bbl_tracer` and `resTracer` is consistent with
```

## Namelist parameters
### BBL_PARM01
- `bbl_wvel` — default vertical entrainment velocity (m/s)
- `bbl_hvel`
- `bbl_initEta` — default initial thickness of BBL (m)
- `bbl_useNofSpeed` — compute the lateral BBL speed from the density contrast and slope, u = bbl_nofCoeff*gPrime*slope/f
- `bbl_nofCoeff` — fraction of the Nof speed gPrime*slope/f that drains downslope (default 1; Campin & Goosse 1999 / MOM use 1/3, which gave too little export in tests)
- `bbl_fMin` — lower bound of |f| in the speed (1/s)
- `bbl_selectEntrain` — entrainment of ambient water into the BBL 0 = none; 1 = Turner (1986) E(Ri), w_E = E*U, with Ri = gPrime*bbl_eta/U^2 and U the BBL speed
- `bbl_selectDetrain` — detrainment of a dense BBL into its cell 0 = at bbl_wvel everywhere (original scheme); 1 = plume-like: at bbl_wvel*max(0,1-uOut/bbl_uDet), uOut the fastest outflow speed of the cell, so a descending plume keeps its water and loses it only
- `bbl_uDet` — outflow speed (m/s) above which a dense BBL no longer detrains (bbl_selectDetrain = 1)
- `bbl_gpTaper` — reduced-gravity scale (m/s^2) over which the BBL fades out as it approaches neutral density: for 0 < gPrime < bbl_gpTaper entrainment is scaled by gPrime/bbl_gpTaper and the BBL dissolves into its cell at the complementary rate. 0 (default) = the
- `bbl_tauMix` — timescale (s) on which a neutral or light BBL dissolves into its cell; 0 (default) = within one time step
- `bbl_coldStart` — on a restart, do not read pickup_bbl but start the BBL as at iteration 0 (BBL T/S = bottom cell, zero thickness, or bbl_*File if set)
- `bbl_doPtracers` — the BBL also carries the passive tracers of pkg/ptracers (default: .TRUE. when usePTRACERS)
- `bbl_thetaFile`
- `bbl_saltFile`
- `bbl_etaFile`

## CPP options (defaults as shipped)
- `BBL_CHECK_BUDGET` (undef, BBL_OPTIONS.h) — Print global sums of BBL T/S tendency x volume x tracer time step (should be zero for a conservative lateral exchange)

## Headers
- `BBL.h` — bbl_wvel    :: default vertical entrainment velocity (m/s) bbl_hvvel   :: default horizontal velocity of BBL (m/s); upper bound of the speed when bbl_
- `BBL_OPTIONS.h` — CPP options file for BBL Use this file for selecting options within package "BBL"
- `BBL_PTR.h` — Passive tracers (pkg/ptracers) carried by the bottom boundary layer. Requires PTRACERS_SIZE.h to be included before this file.
- `bbl_ad_check_lev1_dir.h` — BBL state carried between time steps (pkg/bbl) ADJ STORE bbl_eta   = comlev1, key = ikey_dynamics, kind = isbyte ADJ STORE bbl_theta = comlev1, key = 
- `bbl_ad_check_lev2_dir.h` — BBL state carried between time steps (pkg/bbl) ADJ STORE bbl_eta   = tapelev2, key = ilev_2 ADJ STORE bbl_theta = tapelev2, key = ilev_2 ADJ STORE bbl
- `bbl_ad_check_lev3_dir.h` — BBL state carried between time steps (pkg/bbl) ADJ STORE bbl_eta   = tapelev3, key = ilev_3 ADJ STORE bbl_theta = tapelev3, key = ilev_3 ADJ STORE bbl
- `bbl_ad_check_lev4_dir.h` — BBL state carried between time steps (pkg/bbl) ADJ STORE bbl_eta   = tapelev4, key = ilev_4 ADJ STORE bbl_theta = tapelev4, key = ilev_4 ADJ STORE bbl

## Routines (17)
`bbl_calc_rho.F`, `bbl_calc_rhs.F`, `bbl_check.F`, `bbl_diagnostics_init.F`, `bbl_diagnostics_state.F`, `bbl_init_fixed.F`, `bbl_init_varia.F`, `bbl_read_pickup.F`, `bbl_readparms.F`, `bbl_tendency_apply.F`, `bbl_write_pickup.F`

## Called from outside the package
- `BBL_TENDENCY_APPLY_S` ← `model/src/apply_forcing.F:998`
- `BBL_TENDENCY_APPLY_T` ← `model/src/apply_forcing.F:766`
- `BBL_CALC_RHO` ← `model/src/do_oceanic_phys.F:697`
- `BBL_CALC_RHS` ← `model/src/do_oceanic_phys.F:1057`
- `BBL_DIAGNOSTICS_STATE` ← `model/src/do_statevars_diags.F:90`
- `BBL_TENDENCY_APPLY_S` ← `model/src/external_forcing.F:807`
- `BBL_TENDENCY_APPLY_T` ← `model/src/external_forcing.F:603`
- `BBL_CHECK` ← `model/src/packages_check.F:260`
- `BBL_INIT_FIXED` ← `model/src/packages_init_fixed.F:336`
- `BBL_INIT_VARIA` ← `model/src/packages_init_variables.F:522`
- `BBL_READPARMS` ← `model/src/packages_readparms.F:212`
- `BBL_WRITE_PICKUP` ← `model/src/packages_write_pickup.F:123`
- `BBL_TENDENCY_APPLY_PTR` ← `pkg/ptracers/ptracers_apply_forcing.F:116`
