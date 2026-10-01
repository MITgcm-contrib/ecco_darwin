# Pipeline scripts (STEP1–STEP3)

The scripts that generated this configuration's grid, boundary conditions and
initial conditions. They are the thin Kerguelen-specific wrappers; the actual
work is done by `ecco_darwin/regions/downscaling/utils/*.py`. See
`../../PIPELINE.md` for what each step produces and why.

**These are shipped as run, not generalised.** Every one of them carries
absolute Pleiades paths (`/nobackup/dcarrol2/...`, `/nobackupp19/dcarrol2/...`),
PBS `-l select=` lines sized for this domain, and in some cases assertions
tuned to this configuration's exact counts. Read and edit before reuse —
running them unmodified will either fail on a missing path or write into
someone else's directory.

| script | step | produces |
|---|---|---|
| `run_step1_mitgrid.sh` | 1 | `kerguelen.mitgrid` |
| `run_step1_mitgrid_repad.sh` | 1 | re-padded mitgrid (wet-point fix) |
| `job_step1_bathy.sh` | 1 | `kerguelen_bathymetry.bin` |
| `job_step1_stitch.sh` | 1 | `kerguelen_ncgrid.nc` |
| `job_step1_dvmasks.sh` | 1 | `parent/inputs/*_BC_mask.bin` |
| `job_step2_parent_dv.sh` | 2 | parent `diagnostics_vec` extraction |
| `job_gen_obcs_fast.sh` | 3 | `forcings/OBCS/*.bin` |
| `job_gen_obcs_transport_correction.sh` | 3 | rewrites `UVEL_*`/`VVEL_*` in place |
| `job_gen_pickups.sh` | 3 | `forcings/pickups/*` |

`job_step2_parent_dv.sh` is the parent-side run script (originally
`parent_run/job_ECCO_darwin`). It is a **csh** script, unlike the rest.

Order is strictly sequential and a change upstream silently invalidates
everything downstream — the arrays keep their shapes and the model still
starts. See the regeneration-order note at the end of `../../PIPELINE.md`.
