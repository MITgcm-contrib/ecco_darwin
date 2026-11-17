# Open questions before this tree goes to GitHub (2026-11-16)

Decisions for Dustin. Nothing below is decided in the staged files.

## 1. The physics is not frozen (DRAFT everywhere)

- Stop rule adopted 2026-11-16: at most 2 more full Jan-Aug 2015 runs, one
  change per run, then freeze the run with the lowest total error.
- Run 1 of 2 is `run_arch_D10` (PBS 25399446, submitted 2026-11-16; still
  QUEUED at 18:00 on 2026-11-16, "no available resources"): ice-ocean
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
- ODEPACK (`opkd*.F`): no license header, but the LLNL ODEPACK page states
  the solvers are in the public domain (read 2026-11-16). Resolved.
- `useICEPLUME = .FALSE.` in every run so far. DECIDED 2026-11-16 (Dustin):
  keep it off in `input/data.pkg`, as run; `readme_pleiades.txt` section 5
  says how to switch it on.

## 3. BGC production tiling and period

- The BGC chain Jan 1 to Sep 1 2015 ran in 3 segments on 192 ranks
  (`run_bgc15_s1..s3`, 2.337 s per step, 1.07 model days per wall hour).
- Scaling test for the 2016-2025 production (1-day BGC runs, Jan 1 cold start):
  - 40x40 tiles, 192 ranks (run_bgc_test2d, 2026-10-30): 2.337 s per step.
  - 20x20 tiles, 678 ranks (run_bgc_sc_20x20, PBS 25399447, 2026-11-16,
    ended normally, 1440 steps): 1385.7 s wall = 0.962 s per step, so 2.43x
    faster for 3.53x the ranks (69 % parallel efficiency). BLOCKING_EXCHANGES
    134 s, THERMODYNAMICS 429 s, GCHEM_FORCING_SEP 353 s of the 1386 s.
    One model day = 23 min; 10 years (2016-2025) = about 1400 wall hours =
    12 long-queue segments of 120 h. (My arithmetic from STDOUT.0000; the
    BGC session owns the scaling analysis.)
  - 40x20 tiles, 365 ranks (run_bgc_sc_40x20, PBS 25399711): queued. `build_bgc_s20x20`, `build_bgc_s40x20`,
  `code_size20x20`, `code_size40x20`, `prep/scaling/`, `make_build_scaling.py`,
  `stage_scaling_1d.sh` and `job_run_scaling.sh` are NOT staged, per the
  instruction to ignore them. If a new tiling is adopted, `code/SIZE.h` and
  `input/data.exch2` change together (and the build pre-flight checks that).
- The production period (2016-2025 proposed) follows the measured rate.

## 4. The staged BGC build -- test build DECIDED 2026-11-16 (Dustin: stage and qsub)

- `build_bgc` (2026-10-29) was `-mo "code_capload code_bgc code"`, before
  `code_swfix` and `code_drag10`. The staged `code_darwin + code` = build_bgc +
  swfix + drag10. For BGC runs the swfix guard is a no-op by construction.
- Test build: `regional_run/v05test/` (text copy of the staged `v05/`) +
  `regional_run/job_build_v05test.sh` -> `build_v05test/`. RESULT (PBS
  25399785, 2026-11-16 18:48): BUILD OK, 0 error lines, mitgcmuv 291 MB,
  all 8 expected objects, packages darwin gchem ptracers radtrans oasim
  iceplume shelfice obcs exch2 seaice ggl90 configured. Every override came
  from the staged tree (12 physics files incl. the SWFrac backport, 5 Darwin
  files checked; darwin_plankton md5 6cc6dcc5). Compiled sources differ from
  build_bgc ONLY in the 2 swfix + 4 drag10 files (+ the 2 drag10 files new to
  it), and from build_drag10 ONLY in the 18 code_darwin files. Not run.
  (The job's own assert on packages.conf failed because genmake2 does not
  link that file; the asserts were rerun by hand and the script fixed.)
- `input_darwin/data.diagnostics` regenerated from `input/data.diagnostics` +
  the run_bgc15_s1 BGC streams (DECIDED 2026-11-16; see PROVENANCE.md 4).
  `input_darwin` otherwise = run_bgc15_s1, whose physics was the older Y16 set.

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

## 7. Clean-ups -- DECIDED 2026-11-16 (Dustin)

- `scripts/` and `scripts/pipeline/` stay AS RUN, with absolute Pleiades paths
  and PBS resource lines. They are records, not portable tools.
- `code/` keeps the 30 Darwin files of the old (2026-10-07) `code_darwin` set,
  so that `code/` reproduces `build_drag10` file for file. 17 of them are stale
  and are replaced by `code_darwin/` in a BGC build; they never execute in a
  physics run (`useGCHEM = .FALSE.`). Listed in PROVENANCE.md section 2.
- `code/seaice_growth.F` and `code/apply_forcing.F` are the `code_swfix`
  versions. The BGC chain's winners for these two names are the same code
  without the `IF (useGCHEM)` guard, so they are not staged in `code_darwin/`.
