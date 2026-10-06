# MITgcm troubleshooting history

Symptom → cause → fix entries distilled from the mitgcm-support mailing list and MITgcm GitHub
issues/PRs. 589 entries in total. **Grep first**: `grep -ri '<exact error text or parameter>' references/troubleshooting/`.

Distilled through: **2026-10-05** (mail archive 2003-08 to 2026-09, GitHub issues/PRs to 2026-10-05).
Names checked against MITgcm origin/master ae4c03af7 (2026-10-03) and darwin3 28ac947a3.

| File | Entries | Covers |
|---|---|---|
| `build.md` | 50 | genmake2, optfiles, compilers, MPI/NetCDF linking, start-up segfaults, TAF/Tapenade builds |
| `parallel.md` | 51 | SIZE.h, nPx/nPy, tiling, exch2/blank tiles, OpenMP, hangs, scaling, reproducibility across tilings |
| `instability.md` | 53 | blow-ups, NaNs, CFL, time-step choices, plus a "How to diagnose a blow-up" checklist |
| `restart.md` | 25 | pickups, nIter0, changing deltaT/tiling on restart, automated restarts |
| `grid.md` | 51 | cubed sphere/LLC/curvilinear grid files, bathymetry, hFac, vertical coordinates, r* |
| `obcs.md` | 43 | open boundaries, Orlanski/Stevens, sponge, tides, OB balance, OB input files |
| `seaice.md` | 39 | pkg/seaice, thsice, LSR/JFNK/EVP, very thick ice |
| `forcing.md` | 43 | pkg/exf, cal, interpolation, runoff, restoring, flux forcing with sea ice |
| `diagnostics.md` | 59 | data.diagnostics, output timing, MDS/MNC, **budget closure** (linear FS, r*/nonlinear FS) |
| `adjoint.md` | 68 | TAF errors, tapes/store directives, ctrl/cost, grdchk, optim, Tapenade, ECCO |
| `tracers.md` | 53 | ptracers, gchem/DIC/BLING/darwin, GMRedi, KPP, GGL90, viscosity |
| `physics.md` | 54 | free surface/cg2d, EOS/TEOS-10, shelfice, RBCS, drag, advection schemes, nonhydrostatic |

Entry format: `### <symptom / exact error>`, then `Cause`, `Fix`, `Era` (years seen, checkpoint, whether it was
fixed upstream or a name was renamed or removed) and `Src` (GitHub URL, or `mitgcm-support YYYY-Month 'Subject'`:
the month's thread list is http://mailman.mitgcm.org/pipermail/mitgcm-support/YYYY-Month/thread.html).

## How far to trust an entry

- It records what was advised at the time, mostly by core developers. Read the **Era** line, and confirm
  names and behaviour in current source (or the user's checkpoint, e.g. checkpoint68g for ECCO v4r5) before acting.
- Cited parameter/CPP/routine names were grepped against master; names that don't exist are flagged
  in the entry (renamed, removed, never merged). Linker symbols, error text and placeholders are left as quoted.
- Coverage is uneven. build, tracers, diagnostics, instability, restart and physics read every thread.
  seaice/forcing read only core-answered threads (long posts truncated). parallel/adjoint/grid/obcs read
  the topic-relevant threads after skimming all titles. Absence of an entry doesn't mean nobody hit it:
  grep the raw digests too (below).
- The same problem can appear in more than one file (e.g. CALC_R_STAR in instability, seaice and adjoint).

## Raw material and updating

- Raw cache (not in the skill, not synced to pfe): `~/.cache/mitgcm-ecco-sources/` (`mail/*.txt.gz`,
  `github/*.jsonl`, `digest/<topic>.md`). Grep the digests for full threads when an entry is thin.
- `scripts/fetch_support_sources.sh`: download/refresh the archive and GitHub issues/PRs (~5 min).
- `scripts/digest_support_sources.py [cache] --since YEAR`: answered threads → topic digests.
- `scripts/refresh_troubleshooting.sh [YEAR]`: both of the above. Then read the new digest items (newer than
  "Distilled through"), merge them into the matching file in the same format, and update the date and counts above.
