# pkg/bbl: an embedded bottom boundary layer for MITgcm and ECCO-Darwin

A thin sub-layer of the bottom grid cell lets dense shelf water flow downslope while
keeping its properties, so it can become Antarctic Bottom Water. The bottom boundary
layer (BBL) carries T, S and every passive tracer (Darwin fields included), conserves
all of them to round-off, gives bit-identical results on any tiling, MPI/OpenMP
decomposition and restart, and has a verified TAF adjoint. Optional plume physics:
Nof (1983) speed, Turner (1986) entrainment, a fade-out near neutral density, and
plume-like detrainment.

This package replaces the earlier, undocumented MITgcm pkg/bbl of the same name; copy
`src/` over `pkg/bbl` in your MITgcm or darwin3 tree.

Author: Dustin Carroll (SJSU / JPL), 2026. Method, parameters, diagnostics and the
full list of tests: `src/README_pkg_bbl.md`.

## Contents

```
bbl/
  README.md                this file
  src/                     the package: 11 .F files, BBL.h, BBL_PTR.h, BBL_OPTIONS.h,
                           adjoint files (bbl_ad_*), README_pkg_bbl.md
  darwin3_changes/         the 2 darwin3 routines that call the package (forward runs)
    bbl_darwin3_hooks.patch    both changes as one patch
    model/src/packages_init_variables.F   complete, at darwin3 0078141 (2025-03-06)
    pkg/ptracers/ptracers_apply_forcing.F complete, at darwin3 0078141
    older_darwin3/         the same two routines for the older darwin3 of the
                           ECCO v4r5 LLC90 (v05) set-up
  mitgcm_changes/          the same hooks plus the TAF adjoint hooks for MITgcm
    bbl_hooks.patch        against MITgcm master 38a538975 (2026-10-13)
    model/ pkg/            the 8 changed files, complete, at that commit
    checkpoint68g/bbl_hooks_c68g.patch   the same against checkpoint68g
                           (ECCO v4r5–r7 physics)
  verification/global_oce_latlon/        test input.bblphys (Nof speed + entrainment)
                           and its reference output, for MITgcm testreport
```

## What the hooks do

1. `pkg/ptracers/ptracers_apply_forcing.F`: call `BBL_TENDENCY_APPLY_PTR`, so the BBL
   moves passive tracers (otherwise dense water moved by the BBL would leave its DIC,
   O2, nutrients and plankton behind).
2. `model/src/packages_init_variables.F`: call `BBL_INIT_VARIA` after
   `PTRACERS_INIT_VARIA` and `CTRL_INIT_VARIABLES`, so the BBL starts from the
   bottom-cell tracers, including initial-condition controls in an optimization.
3. Adjoint only (`mitgcm_changes/`): `BBL.h` and `BBL_PTR.h` in `the_main_loop.F` and
   `forward_step.F`, and the BBL state in the `pkg/autodiff` checkpoint directives.

## Install in darwin3 (ECCO-Darwin)

From the top of a darwin3 checkout:

```sh
rm -r pkg/bbl && cp -r <path>/bbl/src pkg/bbl
git apply <path>/bbl/darwin3_changes/bbl_darwin3_hooks.patch
```

The patch applies to darwin3 `0078141` (2025-03-06) and `28ac947a3` (2024-09-24). For
the older darwin3 of the v05 LLC90 set-up use the files in
`darwin3_changes/older_darwin3/`, or make the two edits by hand (see the patch). You
can also put the package files and the two routines in your experiment's `code/`
directory instead of editing darwin3. Rebuild from scratch (`make CLEAN`) after
changing headers or adding files to `code/`.

## Install in MITgcm

From the top of an MITgcm checkout:

```sh
rm -r pkg/bbl && cp -r <path>/bbl/src pkg/bbl
git apply <path>/bbl/mitgcm_changes/bbl_hooks.patch          # master
# or: git apply <path>/bbl/mitgcm_changes/checkpoint68g/bbl_hooks_c68g.patch
cp -r <path>/bbl/verification/global_oce_latlon/* verification/global_oce_latlon/
cd verification && ./testreport -t global_oce_latlon
```

`global_oce_latlon` compiles pkg/bbl; `input.bblphys` runs it with the Nof speed and
entrainment.

## Using it

Add `bbl` to `code/packages.conf`, `useBBL = .TRUE.` to `data.pkg`, and a `data.bbl`.
Recommended settings (used in the ECCO v4r5 LLC90 + darwin3 runs):

```
 &BBL_PARM01
 bbl_useNofSpeed   = .TRUE.,
 bbl_nofCoeff      = 1.,
 bbl_selectEntrain = 1,
 bbl_hvel          = 1.,
 bbl_wvel          = 1.E-5,
 bbl_initEta       = 20.,
 bbl_coldStart     = .TRUE.,
 &
```

`bbl_coldStart = .TRUE.` is needed when starting from pickups without `pickup_bbl`
(e.g. switching the BBL on in a spun-up ECCO-Darwin run). Passive tracers are carried
by default when `usePTRACERS = .TRUE.` (`bbl_doPtracers`), and their BBL values are
saved in `pickup_bbl`. For adjoint runs also set the fade-out (`bbl_gpTaper = 1.E-5`,
`bbl_tauMix = 86400.`). Diagnostics: `BBLeta`, `BBLtheta`, `BBLsalt`, `BBLuTrp`,
`BBLvTrp`, `BBLwEnt`, ... `#define BBL_CHECK_BUDGET` in `BBL_OPTIONS.h` prints a
conservation check for T, S and every tracer each step (testing only).

## Tests

- MITgcm: bit-identical across tilings (1x1, 3x2, 6x4), MPI, OpenMP and restarts;
  heat, salt and passive tracers conserved to ≤ 2.5e-16 per step with all options on;
  TAF adjoint vs finite differences 1e-9 to 1e-7 (4° `global_oce_latlon`), adjoint =
  tangent linear to 12 digits; ECCO v4r6 LLC90 adjoint within 0.3% with GM/Redi off
  (3–11% with GM/Redi on, from ECCO's stabilized GM/Redi adjoint).
- darwin3: ECCO-Darwin v06 3° (36 tracers), 24 steps, darwin3 `0078141`: BBL budget per
  step T 1.9e-16, S 8.7e-17, all 36 tracers ≤ 2.0e-16; with the fade-out and plume-like
  detrainment also on, T 1.8e-16, S 1.1e-16, tracers ≤ 2.0e-16. Every tracer changes
  with the BBL on; the BBL adds ~3% to the run time.
- ECCO v4r5 LLC90 + darwin3 (`older_darwin3/` files), 1992–2001: budgets ≤ 5e-16 per
  step over 10 years. With denser shelf water (shelf GM/Redi off + salt-plume brine to
  the bottom) the BBL moves 5.0 Sv south of 60°S (1997–2001) and makes the deep
  Southern Ocean (below 2000 m, south of 60°S) denser by 3.3e-4 kg m⁻³ after 10 years,
  with DIC −0.11 and O2 +0.08 mmol m⁻³.

## Notes

- `mitgcm_changes/bbl_hooks.patch` was regenerated for MITgcm master 38a538975; the
  package was developed and tested on master 63cdc0b96 and on checkpoint68g. pkg/bbl
  itself has no upstream changes in between, but the master build has not been rerun.
- The darwin3 adjoint with the BBL has not been tested.
- The shelf changes used in the LLC90 runs (GM/Redi off on the shelf, brine to the
  bottom) are experiments, not part of this package.
