# Generate regional set-up, initial and boundary conditions

---

## Preliminary information

**Requirement**:
> Before proceeding to the following instructions, you will need to complete steps in the README, STEP 1, and STEP 2.
> The anaconda python environment will be needed to run the code.

At the end of this step, you will have generated the initial and boundary conditions required to run the regional set-up.

---

## I. Generate the initial conditions

**Requirement**:
> First, you need to copy the parent model files with initial and boundary conditions to your local machine (see code below).

```
cd /path/where/region/files/are/stored/.../parent/output
mkdir pickups OBCS
===========================
from super computer to local machine:
cp pickup*.XXXXXXXXXX.* output/pickups/.
cp -r dv output/.
===========================
mv dv OBCS
```

Note: replace "XXXXXXXXXX" by the iteration of the pickup file you want your regional model to start with  
**Example:** In v05 ECCO-Darwin: pickup*.0000683784* corresponds to the begining of the year: 1992 + (683784 x timestep) / (365.25 x 86400) = 2018) with timestep = 1200.

---

By running `gen_pickups.py` in the `utils` folder, you will generate the initial conditions for your regional set-up (`pickup` files)  
Below are the instructions to run the code in a terminal with the anaconda envrironment:

```
conda activate downscaling
cd ecco_darwin/regions/downscaling/utils/
python3 gen_pickups.py -d /path/to/regional/files -n name_of_the_region\
                       -i pickup_iteration_number -sg Sigma_gaussian_filter\
                       -bgc -v -nc  
```

**Example:**

```
python3 gen_pickups.py -d ecco_darwin/regions/downscaling/ -n NorthSlope \
                      -i 683784 -bgc -v -nc
```

To get more information about the options required for this code run `gen_pickup.py -h`. Here are additional details about the options:
> - -d: The directory where mitgrid is stored.
> - -n: The name of your region.
> - -i: iteration number of the pickup file chosen (see requirements).
> - -sg: sigma for the Gaussian Filter applied to the interpolation. This parameter is **optional**, if not prescribed, sigmaG will take the value of the average size of the grid cell of the parent grid in the selected region.
> - -bgc: generate Darwin pickup files if running the ECCO-Darwin model.
> - -v: verbose.
> - -nc: generate a NetCDF file containing all forcing matrices (**optional**). This is for control purposes only, it is not required to run the simulation.
---

## II. Generate the boundary conditions

Running `gen_obcs.py` in the `utils` folder, you will generate the boundary conditions for your regional set-up (`OBCS` files)  
Below are the instructions to run the code in a terminal with the anaconda envrironment:

```
mkdir path/to/save/the/grid/forcings/OBCS
conda activate downscaling
cd ecco_darwin/regions/downscaling/utils/
python3 gen_obcs.py -d /path/to/regional/files -n name_of_the_region\
                    -bnd EWNS -i iter_list -seaice -bgc -v  
```

**Example:**

```
python3 gen_obcs.py -d ecco_darwin/regions/downscaling/ -n NorthSlope \
                    -bnd EWN -i 683785 710065 -seaice -bgc -v
```

To get more information about the options required for this code run `gen_pickup.py -h`. Here are additional details about the options:
> - -d: The directory where mitgrid is stored.
> - -n: The name of your region.
> - -bnd: Open boundaries of your domain where you need to generate boundary conditions >> 'E' (East), 'W' (West), 'N' (North), 'S' (South)
> - -i: iteration of the obcs files extracted from diagnostics_vec.
> - -seaice: generate sea-ice obcs files.
> - -bgc: generate Darwin obcs files (if running the ECCO-Darwin model).
> - -v: verbose.
---

## III. Verify the boundary conditions (recommended)

Before launching the regional simulation, it is strongly recommended to verify that the
generated boundary-normal velocities are consistent with the parent model. A downscaled model
whose open-boundary volume transport does not match the parent will typically be unstable or
require `OBCSbalance` as a workaround. Two diagnostic scripts in the `utils` folder help catch
problems (missing/masked boundary velocities, transport mismatches, parent/child coverage
gaps) **before** you spend compute on a run.

Both scripts are read-only with respect to the OBCS files that `gen_obcs.py` produced (the
transport check writes only if you explicitly ask it to correct), and both re-use
`gen_obcs.py`'s own grid- and velocity-reading helpers so their results stay consistent with
what was generated. Run them after `gen_obcs.py` has written the `OBCS` files. Pass the **same
iteration(s)** to `-i` that you gave `gen_obcs.py`.

### III.1 Boundary velocity sections — `plot_boundary_velocity_profile.py`

Plots the time-mean, depth-resolved boundary-normal velocity at each open boundary as a
depth × along-boundary "section" — the **parent (native LLC270)** panel beside the
**downscaled (child)** panel, on a shared axis. This exposes the vertical (baroclinic)
structure directly, so you can see exactly where — both along the boundary and at depth — the
child departs from the parent. It only saves PNGs; it never writes grid or OBCS files.

```
conda activate downscaling
cd ecco_darwin/regions/downscaling/utils/
python3 plot_boundary_velocity_profile.py -d /path/to/regional/files -n name_of_the_region\
                                          -bnd ENS -i iter_list -o /path/to/plots -v
```

**Example:**

```
python3 plot_boundary_velocity_profile.py -d ecco_darwin/regions/downscaling/ -n NorthSlope \
                                          -bnd ENS -i 683785 -v
```

> - -d: Directory where `gen_obcs.py`'s inputs/outputs live (the same `-d` you gave it).
> - -n: The name of your region.
> - -bnd: Open boundaries to plot, e.g. `ENS` ('E' East, 'W' West, 'N' North, 'S' South).
> - -i: Same iteration(s) used when running `gen_obcs.py`/`gen_obcs2.py` (**required**).
> - -o: Where to save plots (**optional**, default `<config_dir>`); plots are written under `<output_dir>/diagnostics/`.
> - -legacy-bug-mode: Replicate the original `gen_obcs.py`'s multi-file merge bug when reading the native files. Only needed if the OBCS files you are comparing against were generated by the **unfixed** `gen_obcs.py` (see note in III.2). **Optional.**
> - -v: verbose.

**What to look for:** the child panel should reproduce the parent panel — same sign, same
vertical and along-boundary structure. If a child panel is blank (or covers a much narrower
strip) while the parent clearly shows flow, that boundary is not carrying its velocity;
inspect the boundary mask/geometry before running (see the note at the end of this section).

### III.2 Boundary transport check — `check_obcs_transport.py`

Integrates the volume transport across each boundary and compares the **parent** against the
**downscaled** OBCS, reporting the mismatch per boundary. With `-plot`, it also saves a
time-mean transect (transport per unit length vs. latitude/longitude) plus a
boundary-integrated net-transport time series (Sv) for each boundary. By default it only
**reads** the OBCS files and reports/plots — it writes nothing unless you pass `-correct`.

```
conda activate downscaling
cd ecco_darwin/regions/downscaling/utils/
python3 check_obcs_transport.py -d /path/to/regional/files -n name_of_the_region\
                                -bnd ENS -i iter_list -plot -v
```

**Example (report + plots only):**

```
python3 check_obcs_transport.py -d ecco_darwin/regions/downscaling/ -n NorthSlope \
                                -bnd ENS -i 683785 -plot -v
```

**Example (also write a transport-corrected OBCS set):**

```
python3 check_obcs_transport.py -d ecco_darwin/regions/downscaling/ -n NorthSlope \
                                -bnd ENS -i 683785 -correct -plot -v
```

> - -d: Directory where `gen_obcs.py`'s inputs/outputs live (the same `-d` you gave it).
> - -n: The name of your region.
> - -bnd: Open boundaries to check, e.g. `ENS`. Only boundaries with a normal-velocity component (E/W → `UVEL`, N/S → `VVEL`) are meaningful here.
> - -i: Same iteration(s) used when running `gen_obcs.py` (**required**).
> - -correct: Write transport-corrected `UVEL`/`VVEL` files (plus unchanged copies of every other boundary file) to `--output_dir`. Without this flag nothing is written — you get only the mismatch report (and plots, with `-plot`). The correction is a single per-boundary **barotropic** velocity offset, uniform across depth and along the whole boundary, so the baroclinic shape and along-boundary pattern that `gen_obcs.py` produced are left untouched. **Optional.**
> - -o: Where to write the corrected files (**optional**, default `<config_dir>/forcings/OBCS_corrected/`). It never overwrites `gen_obcs.py`'s original output, so you can point `data.obcs` at whichever set you trust.
> - -plot: Save the transect + time-series diagnostic plot per boundary (requires matplotlib). **Optional.**
> - -maxcorr: Safety cap (m/s) on the barotropic correction added to any single column (**optional**, default 2.0).
> - -veltol: How far beyond the parent's own observed velocity range (as a fraction of that range, **optional**, default 0.2 = 20%) the correction may push any single depth level. Keeps the corrected min/max velocities close to the parent's, not just the net transport.
> - -legacy-bug-mode: Replicate a known bug in the **original** `gen_obcs.py` — when 3+ native `diagnostics_vec` files are present, only the first and last survive its merge loop. Use this **only** when checking OBCS that were generated by that unfixed `gen_obcs.py`; leave it off (default) for files from `gen_obcs2.py` or any fixed generator, which correctly use every available file. **Optional.**
> - -v: verbose.

**What to look for:** the downscaled net transport should track the parent net transport at
each boundary. A large, persistent shortfall (the downscaled boundary carrying much less
transport than the parent) points to a boundary that is masked to too few cells, or a
parent/child coverage gap, and is best resolved at its source rather than masked by
`OBCSbalance`. If the residual mismatch is genuinely physical (e.g. the coarse parent pushes
flow over cells your finer child resolves as land), `-correct` writes a barotropic-matched
set to `forcings/OBCS_corrected/` that you can point `data.obcs` at.

> **Note (boundary masks / C-grid staggering):** both scripts read the child grid from the
> regional `_ncgrid.nc` file. On a C-grid the meridional face field `HFacS` has one extra row
> (`ny+1`) relative to the cell-centered fields, so the open southern/northern boundary lives
> on the **interior** face, not the outer domain edge. The scripts therefore trim `HFacS` as
> `[:, 1:, :]` (interior face) — not `[:, :-1, :]` (closed outer edge). If a boundary transport
> or section looks starved even though the grid and OBCS files are correct, confirm that any
> script reloading the grid uses the interior-face trim; a stale copy of `[:, :-1, :]` in a
> diagnostic will mask a correct boundary down to a thin strip along one end.

---

**CONGRATULATIONS!!** You have now generated (and verified) the pickup and boundary files required to run your regional set-up.\

Considering the vast capabilities of ECCO/MITgcm and the number of parameters that can be tuned, we let the users generate and customize their own regional set-up.  
You can however start from the `gen_mod` folder used in STEP1 as a starting point.  
You can also check the `gen_mod/reg_example_files`.  
We placed important files for setting up the regional model (such as `data.obcs`) where you set the OBCS files to read.