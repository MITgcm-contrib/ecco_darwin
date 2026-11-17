# Pipeline scripts (grid to BGC inputs)

The scripts that generated this configuration's grid, bathymetry, parent
extraction, boundary conditions, initial conditions, tides, iceplume forcing
and BGC inputs. Most are thin Nares-specific wrappers; the heavy work is in
`ecco_darwin/regions/downscaling/utils/*.py` (gen_pickups.py, gen_dvmasks.py,
stitch_ncgrid.py). See `../../PIPELINE.md` for what each step produces and why.

**These are shipped as run, not generalized.** Every one of them carries
absolute Pleiades paths (`/nobackup/dcarrol2/...`, `/nobackupp19/dcarrol2/...`),
PBS `-l select=` lines sized for this domain, and assertions tuned to this
configuration's exact counts. Read and edit before reuse.

| folder | script | step | produces |
|---|---|---|---|
| grid | `gen_nares_mitgrid.py` | 1 | `Nares.mitgrid` (rotated lat-lon, analytic metrics, one facet) |
| grid | `gen_bathy_nares.py` | 1 | `Nares_bathymetry.bin`, `Nares_shelficeTopo.bin` (GEBCO 2026 + BedMachine v6 Petermann cavity) |
| grid | `wall_fill_bridges_nares.py` | 1 | `Nares_bathymetry_walled.bin` (571 min-depth fill cells walled off) = the bathymetry the runs use |
| grid | `deepen_min_depth_nares.py`, `job_min20m_nares.sh` | 1 | 20 m min-depth variant (tested, NOT adopted) |
| grid | `gen_exch2_blanklist.py` | 1 | `input/data.exch2` (252 tiles, 60 blank) |
| grid | `gen_ob_index.py` | 1 | the OB_I*/OB_J* arrays in `input/data.obcs` |
| grid | `gen_humboldt_plumes.py` | 1 | Humboldt plume mask/length files and `Humboldt_Qsg_daily.npz` |
| grid | `plot_nares_grid.py` | 1 | review figure |
| gen_ncgrid | `code/`, `input/`, `job_build_ncgrid`, `job_run_ncgrid` | 1 | `Nares_ncgrid.nc` (1-step grid dump on all 252 tiles) |
| parent_run | `code_darwin/`, `input/`, `job_build.sh`, `job_ECCO_darwin`, `chain_next.sh`, `check_dv.py` | 2 | shared west_AO + Nares LLC270 V5r1 `diagnostics_vec` extraction (40 vectors, hourly, 2015-2025) |
| obcs_ic | `stage_nares_config.sh` | 3 | symlink layout (parent grid, parent pickups, dv masks) |
| obcs_ic | `job_nares_arch_inputs.sh` | 3 | first physics inputs: snapshot_dv, gen_pickups, sea-ice IC, tRef/sRef, gen_obcs_fast, throughflow offset, shelfice load |
| obcs_ic | `snapshot_dv.py` | 3 | fixed-length copies of the growing extraction |
| obcs_ic | `gen_obcs_fast.py`, `fast_extrapol.py` | 3 | OBCS files (fast replacement for gen_obcs.py) |
| obcs_ic | `rotate_obcs_velocity.py`, `job_rotate_obcs.sh` | 3 | OBCS u/v rotated from east/north to the grid axes |
| obcs_ic | `nares_throughflow_offset.py`, `job_offset_rot.sh`, `obcs_boundary_budget.py`, `job_obcs_budget.sh` | 3 | +0.331 Sv southward throughflow offset and its budget check |
| obcs_ic | `job_obcs_v2.sh`, `job_obcs_2015JanAug.sh` | 3 | the full chain as run: `OBCS2_rot_tf` (Jan-Jun), `OBCS2_rot_tf_JanAug` (Jan 1 to Sep 1, 5856 records) |
| obcs_ic | `nares_seaice_ic.py`, `make_ic_thin_mouth.py` | 3 | sea-ice IC files; thinned-mouth HEFF (`_ICT`) |
| obcs_ic | `nares_tref_sref.py`, `nares_shelfice_load.py`, `jmd95.py` | 3 | tRef/sRef for `data`; `Nares_shelficeLoadAnomaly.bin` |
| tides | `nares_ob_points.py`, `gen_obcs_tides.py`, `check_tide_files.py` | 3 | `OBCS_tides_2015Jan_c2` (8 constituents, FES2014 currents, FES2022 prism scaling) |
| iceplume | `gen_runoffQsg_files.py`, `job_runoffQsg.sh` | 3 | yearly `runoffQsg_YYYY` (Humboldt daily discharge) |
| bgc | `job_bgc_inputs.sh` and the 6 tools it calls | 3 | BGC pickups and OBCS (Jan 1 to Sep 1 2015): derived rDOC/CDOM, clipped negatives, seeded plankton, nutrient fix v3 |

`job_ECCO_darwin` is a **csh** script, unlike the rest. The parent extraction
is shared with the west_AO configuration (`parent_V5r1_westAO_nares`): its
`data.diagnostics_vec` lists 8 masks (4 west_AO + 4 Nares) x 5 groups.

Order is strictly sequential, and a change upstream silently invalidates
everything downstream. See the regeneration-order note at the end of
`../../PIPELINE.md`.
