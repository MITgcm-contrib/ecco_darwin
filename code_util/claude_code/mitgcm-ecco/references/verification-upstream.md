# Verification, testreport and upstreaming

## testreport

```
cd <MITgcm>/verification
./testreport -of ../../build_options/darwin_x86_gfortran_local -t global_oce_latlon -devel
./testreport -of <optfile> -j 8 -t 'exp1 exp2'      # several, parallel make
./testreport -of <optfile> -MPI 4 -t <exp>           # MPI with 4 procs (-mpi alone = 2)
./testreport -clean                                   # tidy afterwards
```

- Compares `%MON` stats and cg2d output against `verification/<exp>/results/output.txt`
  and reports matching digits; default pass threshold is `MATCH_CRIT=10` digits (`-match N` to change). Expect ~13+ with -devel/-ieee, fewer when optimised.
- `-devel` adds bounds checking and IEEE settings. Use it for development.
- On this Mac, `internal_wave` and `global_ocean.cs32x15` fail even in stock MITgcm; that
  is not your change.
- Each experiment can have secondary tests, `input.<suffix>/` → `results/output.<suffix>.txt`.
- A full sweep takes a while (`global_oce_latlon` alone takes ~11 min). Run the experiments
  touched by the change first.
- Regenerating reference output for the user's new experiments: `wetting_drying/scripts/make_wad_refs.sh`;
  checking: `check_wad_refs.sh` / `check_wad_refs_x86.sh`.

## A new verification experiment

```
verification/<exp>/
  code/      packages.conf, SIZE.h, *_OPTIONS.h, overrides (small; serial tiling)
  input/     data, data.pkg, eedata, data.<pkg>, small binaries (or prepare_run script)
  results/   output.txt from a -devel build
  README / description in the experiment's docs
```

Keep it small (seconds to a couple of minutes), deterministic, and designed so the new
feature changes `%MON` output. Otherwise testreport can't detect a regression.

## Correctness checks the user expects

- `NORMAL END` and the expected record counts; no `ABNORMAL END`, no NaN in `%MON`.
- Restart: run 2N steps vs N + N from pickup → bit-identical fields.
- Tiling: compare fields (not `%MON`, whose global-sum order differs) across tilings/MPI.
- OpenMP vs serial bit-identity when threaded.
- Conservation/budget closure flags; `conscheck` for darwin3.
- Keep an adversarial audit note (as in `wetting_drying/docs/audit_*.md`) and retract
  claims that fail it.

## Repositories

| Repo | Use |
|---|---|
| https://github.com/MITgcm/MITgcm | upstream model; tags `checkpoint68g` (ECCO v4r5–r7 base) etc. |
| https://github.com/darwinproject/darwin3 | Darwin; branches `backport_ckpt68y`, `backport_ckpt68g` |
| https://github.com/MITgcm-contrib/ecco_darwin | ECCO-Darwin v04/v05/v06 configs, regional set-ups, `offline/V4r6_darwin_offline` |
| https://github.com/ECCO-GROUP/ECCO-v4-Configurations | official ECCOv4 code/namelists (Release 6); NAS reproduction doc |

## Upstreaming decisions so far

- **pkg/wad**: check-in branch `wad` in `ECCO/MITgcm_wad_checkin` (based on `d861cd501`,
  sea-ice coupling removed for the first submission at a colleague's suggestion). There's
  also a drop-in `ECCO/wad_for_ecco_darwin/pkg/wad` kept in sync by
  `ECCO/wetting_drying/scripts/sync_wad_checkin.sh`. Change the check-in tree only on purpose, and re-sync after.
- **pkg/bbl**: no MITgcm PR. It goes to MITgcm-contrib/ecco_darwin. Backport for
  checkpoint68g lives in `ECCO/BBL/MITgcm_c68g` (branch `bbl-c68g`).
- The user pushes to their own fork and opens PRs. Claude prepares commits, the PR text and
  the `doc/tag-index` entry, but doesn't push. Don't suggest contacting individual MITgcm
  maintainers (e.g. Jean-Michel Campin) directly.

## MITgcm PR conventions

- One logical change per PR; new behaviour off by default so existing results don't change.
  If reference outputs must change, update the affected `results/output*.txt` in the same
  PR and say why.
- Add a `doc/tag-index` entry (`o pkg/<name>:` + indented description) and
  `doc/phys_pkgs/<name>.rst` docs for a new package.
- GitHub Actions (`.github/workflows/build_testing.yml`) runs testreport on a set of
  experiments. Run the relevant ones locally first.
- Commit messages: `pkg/<name>: <what changed>` (the style used in the user's branches).
