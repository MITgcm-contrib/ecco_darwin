# Regional downscaling (ECCO-Darwin parent → regional MITgcm/darwin3)

Dustin nests km-scale regional MITgcm/darwin3 configurations inside ECCO-Darwin (LLC270 v05,
also LLC4320/AO set-ups) with the `ecco_darwin/regions/downscaling` toolkit. This file is the
map of where the work lives, the pipeline, the working rules and the failure modes that
actually happened. Status lines are a snapshot (2026-10-06): **the region's `DECISIONS.md` is
the live record, so read it before touching a region.**

## Where things live

| What | Where |
|---|---|
| Region working dirs (Pleiades) | `/nobackup/<nas_user>/downscaling/<region>/`; large `run`/`forcings` dirs are often symlinks to `/nobackupp19/<nas_user>/downscaling/...` |
| Per-region working log | `<region>/DECISIONS.md`: dated `### ` entries saying what was run, what failed and why; append, never rewrite |
| Toolkit docs (generic STEP1–4) | `ecco_darwin/regions/downscaling/{README,STEP1..4}.md`; local clone `~/Documents/research/downscaling/method/regions/downscaling/` |
| Worked, documented example | `/nobackup/<nas_user>/downscaling/repo_staging/kerguelen/v05/` (`README.md`, `PIPELINE.md`, `readme_pleiades.txt`, `scripts/pipeline/`) |
| Production-chain supervisor | `/nobackup/<nas_user>/downscaling/tools/chain_watchdog.sh` (+ `chain_watchdog_selftest.sh`) |
| Analysis python (pfe) | `/nobackup/<nas_user>/downscaling/micromamba_root/envs/downscaling/bin/python` (numpy, scipy, netCDF4, xarray, matplotlib, cartopy). **Never `source setup_env.sh`**: it runs `micromamba install` through the proxy and hangs the ssh call |
| Bathymetry source | `/nobackup/<nas_user>/downscaling/gebco/GEBCO_2026_CF.nc` (2024 also present) |
| Scratch (pfe) | `/nobackup/<nas_user>/claude_scratch/`; move finished work into `<region>/analysis/` |
| Mac side | `~/Documents/research/downscaling/{figures,movies,python,model_setup,grid,m_files,mat,notes,papers,method}/<region>/` (organised by type, then region) |
| OIF / Kerguelen ecosystem work | `~/Documents/research/OIF/model_setup/kerguelen_eco/` (Mac master), mirrored at `/nobackup/<nas_user>/downscaling/kerguelen_4km/kerguelen_eco/` (newer for some scripts, so diff before copying either way) |

## Regions (snapshot 2026-10-06)

| Region | Grid / ranks | Parent | Status / notes |
|---|---|---|---|
| `kerguelen` | 1170×1008×90, ~1 km, dt 90 s, 720 ranks (30×24) | LLC270 v05, 31 tracers | 3-yr production chain (2020–2022) under the watchdog; validation in `kerguelen/validation/` (obs/, model_v05/, scripts/validate_phaseA.py) |
| `kerguelen_4km` | 292×252×90, 84 ranks (4×21 of 73×12) on 1 rom_ait node | coarsened from the 1 km inputs (`kerguelen_eco/fourkm/make_4km.py`) | kerguelen_eco test bed: darwin3 `3f0529872` (`darwin3_2021/`), radtrans/OASIM backport, pkg/profiles, build in `build_profiles/`; inputs `inputs_v2/` (eco = seeded) |
| `fiji` | 2240-wide grid, production on 2688 ranks (`job_production_chain_r2688.sh`) | LLC270 | production under the watchdog; OBCS regenerated on the extended grid (`OBCS_1700_preext` kept) |
| `pan_AO` | 1680×1680 polar disc, 60×60 tiles, 413 of 784 active (`grid/blanklist_sNx60.txt`), dt 300 | LLC270 | eta blow-up fixed (wet pockets masked → `grid_nopocket/`, per-arc OBCS correction, observed-2020 Bering target); Bering sits ~0.6 m above parent but stable; analysis in `pan_AO/analysis/` |
| `nares`, `west_AO` | new (2026-10-06) | | scaffold, grid, bathymetry; west_AO has pkg/wad + sediment ported |
| `GoM`, `GoM_jra55do_nutrients`, `mac_delta` | darwin3 + ecco_darwin trees | | older set-ups (build notes in `~/Documents/research/downscaling/notes/{GOM,mac_delta}.txt`) |
| Mac-only | `disko_bay`, `nares_strait`, `CCS_surf`, `pan_greenland` | | local grids/figures/m_files; Mike Wood's Greenland L0–L2 runs: `notes/mike_wood_nobackup_path.txt` |

## The pipeline (STEP1–4) and the one rule that matters

```
grid -> dv masks -> parent extraction -> OBCS + pickups -> (transport correction) -> run
```
**A change upstream silently invalidates everything downstream.** Shapes don't change, only
coordinates, so size checks pass and the model starts. Rebuild everything after the grid changes.

1. **STEP1 grid:** mitgrid (`run_step1_mitgrid.sh`), bathymetry from GEBCO
   (`job_step1_bathy.sh`), tiles → `*_ncgrid.nc` (`job_step1_stitch.sh`), dv masks selecting the
   parent cells that feed each boundary (`job_step1_dvmasks*.sh`).
2. **STEP2 parent extraction:** a `pkg/diagnostics_vec` LLC270 run restarted from the
   period's pickup (2020-01-01 = iter 736344; 2023-01-01 = 815256) writing hourly values at the
   mask points. `DIAGNOSTICS_VEC_SIZE.h` `nVEC_mask` must equal the number of `nml_vecFiles`
   entries in `data.diagnostics_vec` (a matched pair); `VEC_points` is per process per mask. A restart
   writes a new file set tagged `<nIter0+1>` with records renumbered, so stitch segments and drop the overlap.
3. **STEP3 BCs/ICs:** `gen_obcs_fast.py` (griddata horizontally, linear vertically, nearest
   outside the hull) → `{VAR}_{north,south,east,west}.bin`; `gen_pickups.py` → cold-start
   `pickup_*.<iter>.data` read via `hydrogThetaFile` etc. at `nIter0=0`; derived tracers
   (`gen_derived_tracer_obcs.py`/`_pickup.py`, e.g. rDOC, CDOM); sea-ice ICs `gen_seaice_ic.py`.
   Then the **barotropic transport correction** (`gen_obcs_transport_correction.py`):
   interpolation doesn't conserve net transport, which shows up as secular ETAN drift. Weight each parent point by
   its position along the true boundary coordinate, not native DXG/DYG: on rotated LLC faces those
   overcount the boundary length 2–3×.
4. **STEP4 run:** stage (`stage_regional_run.sh` / `setup_run.sh` / `stage_run_4km.sh`), short
   test, then `job_production_chain.sh` (self-resubmitting segments under the `long` queue's
   120 h limit, resolving restarts with `find_restart_iter.py`).

## Working rules (on top of the SKILL.md ground rules)

- **Read `DECISIONS.md` first; append a dated entry for every change**, giving what, why, the backup
  name and the job ID. User decisions are labelled `USER DECISION <date>`.
- Work in copies (`grid_nopocket/`, `input_nopocket/`, `OBCS_arcfix/`); keep the live
  `grid/`, `input/`, `forcings/OBCS` untouched until a test passes. Backups are
  `<file>.bak_pre_<change>_<yyyymmdd>`. Generators with hardcoded output paths have
  overwritten live files before (`make_pan_AO_grid_nc.py` GRID_DIR), so check output paths before running.
- Heavy python on the pfe front end gets OOM-killed. Use a small PBS job, or chunk with memmap.
- The ssh ControlMaster wedges under concurrent traffic: serialise ssh/scp, don't run parallel agents
  against pfe, and if `ssh -O check pfe` fails ask the user to reopen it (`! ssh -fN pfe`). Never
  handle passcodes.
- Ask before every qsub/qdel/qalter, including resubmitting a fixed test.

## Failure modes seen (symptom → cause → fix)

**Build / staging**
- Build silently ignores an override → genmake2 keeps an existing symlink; `job_build.sh` now deletes stale links and `*.o *.f`.
- Wrong OBCS/diagnostics options → `-mo "code_darwin code"` order: the first match wins; `code_darwin` must come first.
- Physics-only run that "works" → copy `input/` then `input_darwin/` (the latter overwrites `data.pkg`, `data.obcs`, `data.diagnostics`).
- Traps at start on AMD (`rom_ait`) → remove `-axCORE-AVX2 -xSSE4.2`, add `-mcmodel=medium` (+ `-Wl,--no-relax`). The LLC270 parent binary is vendor-locked: run it on `bro_ele`.
- Test runs for days → `endtime` and `nTimeSteps` are exclusive: delete `endtime` (also lowercase `endtime`) when patching `nTimeSteps`.
- `forrtl severe (29)` on the first diagnostics write → the `diags/<stream>/` dirs don't exist; staging creates every dir named in `data.diagnostics`.
- `X is not a Diagnostic` → names differ between darwin3 versions (the 2021 tree has Rirr/Ed/Es/Eu/PARF/surfPAR but no aplk/bbplk/a/bb). Check every field against the names registered in `build/*diagnostics_init.f` before submitting.
- `numDiags=... but needs at least ...` → raise `numDiags` in the governing `DIAGNOSTICS_SIZE.h` (`code_darwin`), then clean rebuild.
- EOF reading `data.traits` on some ranks → with radtrans, `data.traits` needs an (empty) `&DARWIN_RADTRANS_TRAITS` group.
- SIGFPE in `oasim_calc_solz` → `oasim_dTsolz` > `deltaT` gives a zero stride.
- `mpiexec` returns 0 on an aborted model: judge by `Execution ended Normally` / `ABNORMAL END` in STDOUT.0000. Don't grep bare `STOP`: namelist comments are echoed into STDOUT.

**Grid / bathymetry / inputs**
- Antimeridian lookup errors → GEBCO lon is −180..180, a 0–360 mitgrid isn't (fiji GOTCHA 12).
- Lakes after a bathymetry regrid → fill closed basins not face-connected to the open ocean (`scipy.ndimage.label`, `fill_lakes()` in `make_4km.py`).
- Wet pockets with OB cells but no connection to the interior → blow-up/eta hot spots; mask them (pan_AO: 507 cells).
- Coarsened ICs/OBCS with 8 psu salt and w blow-up → block means included cells the fine model had made land via hFacMin; weight with the fine model's hFacC/W/S.
- Wrong dtype is never caught by size checks: OBCS via exf are usually `>f4`; pickups follow `readBinaryPrec` (1 km `>f8`, 4 km `>f4`). Check bytes = nx·ny·nz·itemsize.

**Boundaries / sea level**
- Secular ETAN drift → transport correction; `useOBCSbalance` with `OBCSbalanceSurf=.TRUE.`.
- A local head that won't go away (pan_AO Bering) → arc-scale imbalance in the OBCS velocities; correct per arc and set targets from observations (Bering A3 2020 monthly). A residual 0.6 m offset with the right transport points at strait resistance, not the boundary.
- Arctic ETAN 1–2 m below the subpolar seas → not drift: with sea ice + real freshwater, ETAN is the surface under the ice. Interpret `ETAN + (910·SIheff + 330·SIhsnow)/1027` (corr 0.95 with the load in pan_AO).

**Ecosystem inputs**
- Plankton that never appear in a regional run → mapped from parent types that are extinct (v05 types 21–24 and microzooplankton 25 are ≤ 0 almost everywhere at Kerguelen). Clip parent plankton at 0 and seed from diatom C (kerguelen_eco `eco_tracers.SEED_FRAC`, `ZOO_SEED_FRAC`), and give mortality a refuge (`XMIN`; darwin3 applies mortality to X − XMIN). Box tests miss this because they seed every type themselves.
- Remapping an ecosystem (e.g. v05 → AO v6 types): see `notes/clement_ecosystem.txt`.

**PBS / chains**
- Chain silently dead → a system hold (`Hold_Types=s`) after repeated failed starts; it can't be released by the user (`qrls`, `qalter` → Unauthorized), so qdel + resubmit. That is what `chain_watchdog.sh` does (opt-in per region: `--arm <region>`, `--status`, `--once --dry-run`). It self-chains with `qsub -a` (no cron here) and needs `/PBS/bin` on PATH. Chains re-arm it every ~30 min.
- `qstat` elapsed time lags; for progress use `time_tsnumber` in STDOUT.0000.

## Analysis recipes

- **Sparse global 2-D files** (blank tiles at the end are never written): read with
  `read_etan()` in `pan_AO/diag_eta_tendency_map.py`. It pads with NaN only at the tail and fails if a
  padded cell is wet.
- **Boundary vs parent:** `pan_AO/analysis/scripts/panao30d_rim_vs_parent.py` (model 2–4 cells
  inside vs the OBCS ETAN, by sector and day).
- **Movies:** render frames on pfe with matplotlib (Agg). pfe has no ffmpeg, so copy the frames and encode on the
  Mac: `ffmpeg -framerate N -pattern_type glob -i 'frame_*.png' -vf "pad=ceil(iw/16)*16:ceil(ih/16)*16" -c:v libx264 -pix_fmt yuv420p`.
  Save to `~/Documents/research/downscaling/movies/<region>/`. Example: `pan_AO/analysis/scripts/panao_movie_frames.py`.
- **Ecosystem validation:** `kerguelen/validation/` (OC-CCI, CMEMS PFT, MODIS PIC, BGC-Argo,
  SWINGS); Nature-style figures and deck builders in `kerguelen_eco/docs/`; pkg/profiles input via
  `kerguelen_eco/scripts/make_profiles.py` (prof names ≤ 8 characters).
