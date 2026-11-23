# Open questions before this tree goes to GitHub (2026-11-16, updated 2026-11-22)

Decisions for Dustin, and what was decided.

## 1. Frozen physics -- DECIDED 2026-11-20 (Dustin: "Freeze D10")

- Stop rule (adopted 2026-11-16): at most 2 more full Jan-Aug 2015 runs, one
  change per run, then freeze the run with the lowest total error.
- Run 1 = `run_arch_D10` (PBS 25399446): ice-ocean drag on the top-10-m mean
  ocean velocity (`code_drag10`, 6 files). Normal end, 28:11 wall. Score:
  closure Feb 8 (obs Jan 30, +9 d), polynya 90/94/121 % of MODIS (Feb/Mar/Apr),
  summer date errors 26 d in sum. It does not meet the stop rule.
- Run 2 = `run_arch_K10` (PBS 25409065): D10 + GGL90 diffusivity x0.1 in
  Kennedy Channel. Normal end, 29.1 h wall. Its summer error was worse (31 d
  vs 26 d); it gained only the April polynya point (117 vs 121 %). Not kept.
- So `code/` = the `build_drag10` chain (it already holds the 6 `code_drag10`
  files), and the physics namelists in `input/` are the D10 namelists
  byte for byte. One exception: `input/data.diagnostics` leaves out D10's
  evaluation-only list 12 (daily `uv10`, UVEL/VVEL levels 1-7). It has no
  effect on the solution.

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

## 3. BGC production tiling -- DECIDED 2026-11-20 (Dustin: 20x20)

1-day BGC runs, Jan 1 2015 cold start, 1440 steps each:

| tiling | ranks | s per step | speed-up | run |
|---|---|---|---|---|
| 40x40 | 192 | 2.337 | 1.0x | run_bgc_test2d (2026-10-30) |
| 40x20 | 365 | 1.425 | 1.64x | run_bgc_sc_40x20 (PBS 25399711) |
| 20x20 | 678 | 0.962 | 2.43x | run_bgc_sc_20x20 (PBS 25399447) |

- Estimate for 2016-2025 (3653 model days) on 20x20: about 1405 wall hours
  = 12 long-queue segments of 120 h on 18 sky_ele nodes, about 25,300
  node-hours (my arithmetic from the timers; queue waits not included).
- Done in the freeze commit (2026-11-22): `code/SIZE.h` (sNx = sNy = 20,
  nPx = 678) and `input/data.exch2` (56 x 18 = 1008 tiles, 330 blank) are
  the `code_size20x20/` and `prep/scaling/data.exch2_20x20` files. New
  `scripts/job_build_v05.sh` checks 678 ranks. `scripts/job_run_scaling.sh`
  runs with NP=678. The older build and run scripts in `scripts/` are the
  as-run record of the 40x40 tuning runs and stop unless they see 192 ranks.

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
- Freeze test build (2026-11-22): this tree with the 20x20 `SIZE.h` and
  `data.exch2`, built by `scripts/job_build_v05.sh` -> `build_v05freeze`
  (PBS 25423792, 1 cas_ait node; the Electra nodes took no jobs that day).
  RESULT: `V05_BUILD_OK`, exit 0, 0 compiler errors, mitgcmuv 137 MB (the
  20x20 `build_bgc_s20x20` binary is also 137 MB). Pre-flight: 33/33 tracers,
  678 ranks = 1008 - 330 blank. Sources that differ from `build_drag10`: the
  17 changed `code_darwin` files + `SIZE.h` (+ `BUILD_INFO.h`). From
  `build_bgc_s20x20`: the 2 swfix + 4 drag10 files (+ `BUILD_INFO.h`, and
  `SIZE.h` in its comment lines only). Not run.
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
