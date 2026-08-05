# Extract regional model boundary information

---
## Preliminary information
**Requirement**: 
> Before proceeding to the following instructions, you will need to complete the steps in README and STEP1.\
> Supercomputing capacities will be needed to run any ECCO model.

At the end of this step, you will have run ECCO global configuration using ``diagnostics_vec`` package and extracted vectors along the boundaries of your regional model.\
<u>Note:</u> ``diagnostics_vec`` package was designed to output model diagnostics from ECCO LLC global model in a subset of the model domain, e.g., along a vector (or "vec").

---
## I. Prepare the simulation for ``diagnostics_vec`` extraction

### a. Turn on ``diagnostics_vec`` package

> - Enable the package in ``packages.conf``
```
cd regions/configs/parent_run/code_darwin
vim packages.conf
##### add the following line to package.conf #####
diagnostics_vec
```
> - Turn on the package in ``data.pkg``
```
cd ../input
vim data.pkg
##### add the following line to data.pkg #####
 useDiagnostics_vec=.TRUE.,
```

### b. Generate a ``data.diagnostics_vec`` parameter file
> - copy the ``data.diagnostics_vec`` file from this repository's ``utils`` folder into the ``input`` directory (paths below are relative to ``downscaling/``)
```
cd regions/configs/parent_run/input/
cp ../../../../ecco_darwin/regions/downscaling/utils/data.diagnostics_vec .
```
> - Modify the ``data.diagnostics_vec`` file according to the specificities of your domain and your requirements (see details inside the files).

<u>Warning:</u> the mask filenames you list in ``nml_vecFiles`` determine the names of the output
files, and ``gen_obcs.py`` in STEP3 looks for files named ``<boundary>_BC_mask_<FIELD>.<iter>.bin``.
Set ``nml_vecFiles`` to the **exact filenames written by ``gen_dvmasks.py``** in STEP1 — that is
``dv/east_BC_mask.bin``, ``dv/west_BC_mask.bin``, ``dv/north_BC_mask.bin``, ``dv/south_BC_mask.bin``
(lowercase, ``_BC_mask``). Renaming them here will produce output that STEP3 cannot find.

<u>Note:</u> the number of entries in ``nml_vec_iters_per_file`` and ``nml_vec_avg_periods`` must
cover every mask you declare in ``nml_vecFiles``. If you list 20 masks, both arrays must be
``(1:20)``.

<u>Note:</u>  You can modify the ``data.diagnostics`` file with only the diagnostics you want to save.\
This won't affect ``diagnostics_vec``, and the fewer diagnostics saved the faster the simulation can be integrated.

### c. Set the compile time ``DIAGNOSTICS_VEC_SIZE.h`` file
> - copy the ``DIAGNOSTICS_VEC_SIZE.h`` file into ``code_darwin``, so it is picked up at compile time

<u>Note:</u> this is a compile-time header, so it belongs in a ``code``/``code_darwin`` directory,
**not** in ``input``. You can take it either from the merged package in ``darwin3`` (shown below)
or from this repository's ``utils`` folder.
```
cd regions/configs/parent_run/code_darwin/
cp ../../../../darwin3/pkg/diagnostics_vec/DIAGNOSTICS_VEC_SIZE.h .
```
> - Modify the ``DIAGNOSTICS_VEC_SIZE.h`` file as follows:
>    - VEC_points: the maximum number of points that can be stored in a ``diagnostics_vec`` mask. Set it to at least the largest of the per-boundary point counts printed by ``gen_dvmasks.py -v`` in STEP1, Section V.c.
>    - nVEC_mask: number of lateral boundary masks used in the ``data.diagnostics_vec`` file — i.e. the number of ``nml_vecFiles`` entries (the example file declares 20).
>    - nSURF_mask: number of surface boundary masks used in the ``data.diagnostics_vec`` file — i.e. the number of ``nml_sfFiles`` entries (the example file declares none, so 0).

---
## II. Compile and run the simulation
Before proceeding, copy the mask files generated in STEP1 (``<config_dir>/parent/inputs/*_BC_mask.bin``)
to the supercomputer, into the folder "downscaling/regions/configs/parent_run/dv".\
This must match the ``dv/`` prefix used in ``nml_vecFiles`` in ``data.diagnostics_vec``.\
Below is an example to extract vectors from v05 ECCO-Darwin.\
On the supercomputer (example on Pleiades) run:

<u>Note:</u> remove the existing ``parent_run`` build/run directory first if it is left over from STEP1, Section V.
```
cd darwin3/regions
rm -r parent_run
mkdir parent_run
cd parent_run
mkdir build run
cd build
```
> - Compile the code
```
for Pleiades users only:
module purge
module load comp-intel mpi-hpe/mpt hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt python3/3.9.12
../../../tools/genmake2 -of ../../../../regions/configs/parent_run/code/linux_amd64_ifort+mpi_ice_nas \
-mpi -mo '../../../../regions/configs/parent_run/code_darwin ../../../../regions/configs/parent_run/code'
make depend
make -j 16
```
> - Run the simulation
```
cd ../run
ln -sf ../build/mitgcmuv .
ln -sf /nobackup/hzhang1/pub/llc270_FWD/nbp19_dmenemen_public_llc270/* .
ln -sf /nobackup/dcarrol2/forcing/apCO2/NOAA_MBL/* .
ln -sf /nobackup/hzhang1/forcing/era_xx .
ln -sf /nobackup/hzhang1/pub/llc270_FWD/input/19920101/to2023/xx*42.data .
cp -r ../../../../regions/configs/parent_run/dv .
cp ../../../../regions/configs/parent_run/input/* .
mkdir diags
# modify data as needed to run the model only on the years required
# modify job_ECCO_darwin as needed
qsub job_ECCO_darwin
```

<u>Note:</u> At the end of the simulation, you will get a binary file for every parameter set in ``data.diagnostics_vec``, each containing the number of iterations chosen for the parameter.

The files are written to the ``dv/`` output directory and named after the mask they came from,
following the pattern ``<boundary>_BC_mask_<FIELD>.<iteration>.bin`` — for example
``east_BC_mask_THETA.0000683784.bin``. These are the files STEP3 consumes: copy the whole set
into ``<config_dir>/parent/outputs/OBCS/`` on the machine running your anaconda environment
(see STEP3, Section I).

<u>Note:</u> if a boundary or field you expected is missing, check that its mask is listed in
``nml_vecFiles`` and that ``nVEC_mask`` in ``DIAGNOSTICS_VEC_SIZE.h`` is large enough to hold
every declared mask — a value that is too small silently truncates the list.

