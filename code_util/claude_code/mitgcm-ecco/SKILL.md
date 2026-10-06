---
name: mitgcm-ecco
description: Expert knowledge of the whole MITgcm code base (every package, namelist parameter, CPP option, the time-step call tree and all verification experiments, via a bundled source index) and of MITgcm, ECCO (v4r4–v4r6, LLC90/LLC270) and ECCO-Darwin / darwin3 as Dustin Carroll uses them — building with genmake2 and optfiles, running on NASA Pleiades (pfe, PBS) and the local Mac, checking code_*/run_* configurations for cross-file consistency, reading MDS/NetCDF output and LLC tiles, regional downscaling from ECCO-Darwin parents (grids, dv masks, diagnostics_vec extraction, OBCS/pickup generation, transport correction, production chains and the chain watchdog for kerguelen, kerguelen_4km/kerguelen_eco, fiji, pan_AO, nares, west_AO, GoM, mac_delta and others under /nobackup/<nas_user>/downscaling and ~/Documents/research/downscaling), writing and extending MITgcm and Darwin packages (pkg/bbl, pkg/wad, OBCS, offline, radtrans, adjoint with TAF/Tapenade), verification with testreport, and upstreaming to MITgcm or ecco_darwin. Use this skill for ANY task that touches MITgcm, ECCO, Darwin/darwin3, ECCO-Darwin, a SIZE.h / data.* namelist / packages.conf / optfile, mitgcmuv, a verification experiment, or model output in .data/.meta form — even when the user doesn't say "MITgcm", and in any project.
---

# MITgcm + ECCO + Darwin

You are working as a co-developer on Dustin's ocean-model code. He builds and extends
MITgcm packages himself (pkg/bbl, pkg/wad, pkg/overflood, OBCS sea-ice sponges), runs
ECCO v4r6 and ECCO-Darwin on NASA Pleiades, and publishes configurations to
MITgcm-contrib/ecco_darwin. Treat the model source as the ground truth: MITgcm's own
code is the documentation, so read the routine before explaining or changing it.

**If you are running on Pleiades** (hostname `pfe*`, home `/home1/<nas_user>`): the Mac paths in
these references don't exist and the index source clones aren't there. Use the index as a map only
and read source in the run's own MITgcm/darwin3 tree under `/nobackup/<nas_user>/`. PBS commands run
directly (still ask before `qsub`/`qdel`), and the Mac-only rules don't apply.

## Ground rules (these come from real incidents)

1. **Runs are precious.** Never `pkill`/`killall mitgcmuv` or kill by pattern — other
   sessions run models on this Mac at the same time, and a pattern kill once destroyed an
   LR17 run. Kill only a PID you launched (saved to a pid file) after checking its cwd.
   Never delete or overwrite a run directory without first grepping the figure/render
   scripts that may read it.
2. **Ask before every `qsub`/`qdel`** on Pleiades, and before anything that spends
   allocation. Treat the user's existing run directories as read-only and work in copies.
3. **Short test before long run.** Run a few steps (a temporary `nTimeSteps`/`nEndIter`
   that you then restore) and check for `NORMAL END` in STDOUT before launching years.
4. **`make CLEAN` after header changes** or after adding a new file to a `code/` dir.
   Incremental builds with stale common blocks have produced fake-good runs, and genmake2
   never replaces an existing symlink, so a new override silently isn't compiled
   (check `ls -la build/<file>.F`).
5. **Keep separate MITgcm clones separate.** `ECCO/BBL/MITgcm`, `ECCO/wetting_drying/MITgcm`,
   `ECCO/sea_ice_BCs/MITgcm`, `ECCO/MITgcm_wad_checkin` etc. each hold uncommitted or
   branch-specific work. Don't edit one to fix another's problem; don't push — Dustin pushes.
6. **Verify claims against source.** When a diagnosis is relayed (from a colleague, a
   log, another session), confirm it in the code before acting. Before handing over a
   result, audit it adversarially and retract anything that fails.
7. **Give shell commands for the user as plain copy-paste text** in your reply, not via
   `echo` in a tool call.

## How to be the MITgcm expert in the room

You have a generated index of the **whole** upstream MITgcm tree (every package, every
verification experiment, every core namelist parameter, the timestep call tree, the manual)
plus darwin3's `pkg/darwin` and `pkg/radtrans`, in `references/mitgcm-index/`. Use it as a
map, then read the actual source before answering:

1. **Locate.** `grep -ril <term> references/mitgcm-index/` finds which package owns a
   parameter, CPP flag, routine or diagnostic; `packages.md` / `verification.md` are the
   catalogues; `pkg/<name>.md` gives a package's namelist params with descriptions, CPP
   options and defaults, headers, and every call site in the core (file:line);
   `core-calltree.md` shows where in the time step something happens; `core-params.md`
   covers `data`/`data.pkg`/`eedata`.
2. **Read the source.** Line numbers drift; open the real file in the relevant clone
   (upstream: `git -C ~/Documents/research/ECCO/BBL/MITgcm show origin/master:<path>`;
   the user's branches: the matching clone in configs.md; darwin3: `~/Documents/GitHub/darwin3`).
   Cite `path:line` in answers.
3. **Learn from a working example.** For "how do I configure X", find verification
   experiments that compile X (listed at the bottom of `pkg/X.md`) and read their
   `code/` and `input*/` files. That's the canonical usage, and testreport-checked.
4. **Check the manual** for equations and intent: `docs.md` lists every manual section with its line number and the parameters it cites, so `grep -n <param> references/mitgcm-index/docs.md` gives file + line to read. The manual can lag the code (e.g. it still mentions the removed `EXACT_CONSERV`).
   `docs-ext.md` does the same for docs outside MITgcm: ECCO v4 configurations (r1–r7 READMEs, flux-forced,
   Docker), the ECCO-v4-Python tutorial, Gael's ECCOv4 docs, every ecco_darwin readme (with the
   checkpoint/optfile/grid it names), and the darwin3 manual's darwin/radtrans pages.
5. When the user's branch differs from upstream, say which you're describing. Branch pages are
   named `<pkg>@<label>` / `<exp>@<label>` (labels → clones in the index `README.md`); they cover
   pkg/wad, overflood, mangrove, sediment, the bbl and obcs/seaice changes, and darwin3 `backport_ckpt68y`.
6. **Freshness:** the index `README.md` gives the indexed commit and date. If that's more than
   ~2 months old, or the user mentions new upstream changes or branch work, offer to run
   `scripts/refresh_index.sh` (~2 min, read-only apart from `git fetch` in `ECCO/BBL/MITgcm`).
   After changing the skill on the Mac, push it to Pleiades with `scripts/sync_to_pfe.sh`.
7. **Community history:** `references/troubleshooting/` is distilled from the mitgcm-support
   list (2003+) and MITgcm GitHub issues/PRs, as symptom → cause → fix with source links and era.
   It records what core developers advised at the time, so confirm parameter names and behaviour in
   current source before relying on an entry, and check the "Era" line (fixed upstream? old checkpoint?).
   Its `README.md` says how far it has been distilled; `scripts/refresh_troubleshooting.sh` pulls newer threads.

## Route to the right reference

Read only the reference(s) the task needs:

| Task | Read |
|---|---|
| What does package/parameter/flag/routine X do; where is it called; which experiment uses it | `references/mitgcm-index/` (start at `README.md`, `packages.md`, `verification.md`) |
| How an ECCO release or ecco_darwin config is set up (checkpoint, optfile, namelists, forcing, ctrl/cost); darwin3 ecosystem equations | `references/mitgcm-index/docs-ext.md`, then the file it points to |
| Which paper to cite / read for a scheme, parameter or package; ECCO, ECCO-Darwin or Dustin's papers | `references/mitgcm-index/bibliography.md` (manual bib + where each paper is cited, by parameter), `references/literature.md` (verified ECCO / ECCO-Darwin / fjord papers). Never cite from memory: use these or resolve the DOI |
| Compile, optfiles, genmake2 errors, Mac/arm64 builds, adjoint builds | `references/build.md` |
| Pleiades access, PBS scripts, modules, remote paths, local long runs | `references/hpc.md` |
| Checking or editing a `code_*`/`run_*` pair, namelists, calendars, offline forcing, diagnostics set-up | `references/config-consistency.md` |
| Reading/plotting output, LLC tiles, MDS format, budgets, figures/movies | `references/analysis.md` |
| Writing a new package, adding a feature/hook to MITgcm or darwin3, diagnostics, pickups, adjoint-readiness | `references/package-development.md` |
| testreport, verification experiments, restart/tiling tests, check-in/PRs, repos | `references/verification-upstream.md` |
| Which config/run lives where (ECCO v4r6, offline Darwin, BBL, WAD, sea-ice BCs) | `references/configs.md` |
| Regional downscaling: any region under `/nobackup/<nas_user>/downscaling` or `~/Documents/research/downscaling` (pipeline STEP1–4, OBCS/pickups, transport correction, eco tracer seeding, production chains, watchdog, region status, known failure modes) | `references/downscaling.md`, then that region's `DECISIONS.md` |
| An error message, blow-up/NaN, hang, wrong-looking output, or "has anyone hit this before?" | `references/troubleshooting/`: grep the exact error text or parameter first (`grep -ri <text> references/troubleshooting/`); `README.md` lists topics |

Project-level `CLAUDE.md` files and the per-project memory directories hold the latest
project specifics. When any of them disagrees with these references or with the source,
the **source code wins**: check it, then say which note is wrong so it can be fixed.

## Fast facts worth always having in mind

- The binary is `mitgcmuv` (`mitgcmuv_ad` for adjoint, `mitgcmuv_ftl` for tangent linear).
- `nPx*nPy` in `SIZE.h` must equal the MPI rank count, so changing ranks means
  recompiling. On LLC grids with exch2, blank (all-land) tiles listed in `data.exch2`
  are dropped: ECCO v4 LLC90 = 117 tiles of 30×30 minus 4 blanks = 113 ranks, Nr=50.
- MDS output is big-endian, usually float32 (`readBinaryPrec`/`writeBinaryPrec`);
  diagnostics in `.data/.meta`, unless `diag_mnc=.TRUE.`.
- Namelists: comment with `!`, never a trailing `#`; keep lines under ~200 chars.
- Diagnostic names change between darwin3 versions — grep `available_diagnostics.log`
  in the run dir rather than trusting memory; compare tracers by name, never by index.
- Check `STDERR.0000` as well as `STDOUT.0000`: bad diagnostic names, parameter
  conflicts and `*_CHECK` failures land there.
- Local darwin3 source: `~/Documents/GitHub/darwin3` (darwinproject/darwin3); its
  `pkg/darwin/darwin_check.F` is generated by `cog` (see package-development.md).
