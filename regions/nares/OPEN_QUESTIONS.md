# Open questions before this tree goes to GitHub (2026-11-16)

Decisions for Dustin. Nothing below is decided in the staged files.

## 1. The physics is not frozen (DRAFT everywhere)

- Stop rule adopted 2026-11-16: at most 2 more full Jan-Aug 2015 runs, one
  change per run, then freeze the run with the lowest total error.
- Run 1 of 2 is `run_arch_D10` (PBS 25399446, started 2026-11-16): ice-ocean
  drag on the top-10-m mean ocean velocity (`code_drag10`, 6 files). The staged
  `code/` INCLUDES these 6 files because `job_build_drag10.sh` is the newest
  physics build. If D10 fails the stop rule, remove them from `code/` (the
  previous binary, `build_swfix`, = `code/` minus the 6 `code_drag10` files).
- `input/` = `run_arch_Y69m_sw` namelists (the last scored full run). D10 uses
  the same namelists plus one daily diagnostics list (`uv10`, levels 1-7).
- Values that depend on this: everything in `data.seaice`, `data.ggl90`, the
  sea-ice source overrides, the `code_drag10` files, and the readme numbers
  for ice-arch timing.

## 2. pkg/iceplume -- DECIDED 2026-11-16 (Dustin: add it to the tree)

- `v05/pkg/iceplume/` holds the 15 package files as compiled (byte-identical
  to the pfe tree), plus Mike Wood's BSD-3 `LICENSE` and a `README.md`.
  Upstream: https://github.com/mhwood/iceplume at e4d4ac9.
- The SEAICE_SWFrac backport (2 pkg/seaice files) is shipped as `code/`
  overrides. With these, public darwin3 at tag `darwin_ckpt69k` builds the
  configuration. Not verified by a build yet.
- Still open: the ODEPACK files (`opkd*.F`) carry no license header.
- Still open: `useICEPLUME = .FALSE.` in every run so far. Decide whether the
  shipped configuration should say "iceplume on" or "iceplume off".

## 3. BGC production tiling and period

- The BGC chain Jan 1 to Sep 1 2015 ran in 3 segments on 192 ranks
  (`run_bgc15_s1..s3`, 2.337 s per step, 1.07 model days per wall hour).
- A scaling test for the 2016-2025 production is in progress (20x20 tiles,
  678 ranks; 40x20 tiles, 365 ranks). `build_bgc_s20x20`, `build_bgc_s40x20`,
  `code_size20x20`, `code_size40x20`, `prep/scaling/`, `make_build_scaling.py`,
  `stage_scaling_1d.sh` and `job_run_scaling.sh` are NOT staged, per the
  instruction to ignore them. If a new tiling is adopted, `code/SIZE.h` and
  `input/data.exch2` change together (and the build pre-flight checks that).
- The production period (2016-2025 proposed) follows the measured rate.

## 4. The staged BGC build has never been built as staged

- `build_bgc` (2026-10-29) was `-mo "code_capload code_bgc code"`. It has NO
  SW-heating fix (`code_swfix`) and NO drag10 (`code_drag10`). Both came
  later. For BGC runs (`useGCHEM = .TRUE.`) the swfix is a no-op by
  construction (it only guards the CDOMdecay path with `IF (useGCHEM)`), so
  the staged `code_darwin + code` build should give the same BGC physics as
  `build_bgc` plus drag10. This is not verified by a build or run.
- `input_darwin/` comes from `run_bgc15_s1`, whose physics namelists are the
  older `run_arch_Y16` set (P* 8.25e4, presPow 1/3, no parent GGL90 set, no
  SImelt diagnostics). The staged `input/` is the newer Y69m_sw set.
  `input_darwin/data.diagnostics` = Y16 list + BGC streams; it lacks the
  `SImelt` list that Y69m_sw added. Decide: regenerate `input_darwin/
  data.diagnostics` from `input/data.diagnostics` + the BGC streams, or keep
  it as run.

## 5. Version label -- DECIDED 2026-11-16 (Dustin: v05)

- Renamed `v06/` to `v05/`, after the parent ECCO-Darwin version, as
  kerguelen does. (ecco_darwin now has a top-level `v06/` for ECCO-Darwin
  v06, which this configuration is not built on.)

## 6. Files I could not source (not in the tree)

- The script that wrote `Nares_grid18.face001.bin` (16 mitgrid fields +
  AngleCS/AngleSN from `Nares_ncgrid.nc`, 2026-10-15). DECISIONS.md records
  the rule but not the script name.
- The dv-mask generation for Nares (`parent/inputs/{east,west,north,south}_BC_mask.bin`).
  DECISIONS.md says `gen_dvmasks` (the toolkit's `utils/gen_dvmasks.py`) on
  the Nares ncgrid. No Nares wrapper script was found under `prep/`.
- `python/nares_strait/gen_exch2_blanklist.py` and `gen_ob_index.py` are
  staged from the local Mac copy (not from pfe). `gen_obcs_tides.py` and
  `make_ic_thin_mouth.py` likewise.
- The thin-mouth sea-ice IC (`pickup_seaice_siHEFF_ICT`) script
  `make_ic_thin_mouth.py` is staged, but its inputs (the original IC and the
  mouth mask) are documented only in DECISIONS.md 2026-10-20.
- The nutrient-fix factor file `nutfix_factor_v3.npz` and `glodap_north50.npz`
  come from the pan_AO project (not in this tree).

## 7. Things to clean before copying to GitHub

- `scripts/` and `scripts/pipeline/` are shipped AS RUN with absolute
  Pleiades paths and PBS resource lines. They are records, not portable tools.
- `code/` contains the 30 Darwin files of the OLD (2026-10-07) `code_darwin`
  set, because the physics builds compile them (packages.conf lists darwin,
  gchem, radtrans, oasim). 17 of them are stale (17-argument
  `darwin_plankton.F`, old RADIv2 sediment) and are REPLACED by `code_darwin/`
  in a BGC build. They are inactive in physics runs (`useGCHEM = .FALSE.`).
  Option: drop them from `code/` and make `code_darwin/` the only Darwin set.
  That changes which Darwin sources the physics binary compiles, so it needs
  one physics build to check (nothing runtime should change).
- `code/seaice_growth.F` and `code/apply_forcing.F` are the `code_swfix`
  versions (6-m ice-load cap + `IF (useGCHEM)` guard around the CDOMdecay
  shortwave path). The BGC chain's winners for these two names
  (`code_capload/seaice_growth.F`, `code_bgc/apply_forcing.F`) are the same
  code WITHOUT the guard. With `useGCHEM = .TRUE.` the guard is a no-op, so
  they are NOT staged in `code_darwin/` (staging them would undo the fix for
  physics-only runs of a BGC binary). See PROVENANCE.md.
