# Extract regional model boundary information

---
## Preliminary information
**Requirement**: 
> Before proceding to the following intructions, you will need to complete steps in README and STEP1.\
> Supercomputing capacities will be needed to run any ECCO model.

At the end of this step, you will have ran ECCO global configuration using ``diagnostic_vec`` package and extracted vectors along the boundaries of your regional model.\
<u>Note:</u> ``diagnostic_vec`` package was designed to output model diagnostics from ECCO llc global model in a subset of the model domain e.g. along a vector (or "vec").

---
## I. Prepare the simulation for ``diagnostic_vec`` extraction

### a. Turn on ``diagnostic_vec`` package

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
> - copy ``data.diagnostics_vec`` file from the ``utils`` folder into the ``input`` directory
```
cd regions/configs/parent_run/input/
cp ../../ecco_darwin/regions/downscaling/utils/data.diagnostics_vec .
```
> - Modify the ``data.diagnostics_vec`` file according to the specificities of your domain and your requirements (see details inside the files).

<u>Note:</u>  You can modify ``data.diagnostic`` file with only the diagnostics you want to save.\
This won't affect ``diagnostic_vec`` and the fewer diagnostics saved the faster the simulation can be integrated.

### c. Set the compile time ``DIAGNOSTICS_VEC_SIZE.h`` file
> - copy ``DIAGNOSTICS_VEC_SIZE.h`` file from the ``utils`` folder into the ``input`` directory
```
cd regions/configs/parent_run/code_darwin/
cp ../../darwin3/pkg/diagnostics_vec/DIAGNOSTICS_VEC_SIZE.h .
```
> - Modify the ``DIAGNOSTICS_VEC_SIZE.h`` file as follows:
    - VEC_points: correspond to the max number of points that can be stored in a diagnostics_vec mask. See maximum number given in STEP 1 runing ``gen_dvmasks``.
    - nVEC_mask: number of lateral boundary mask used in the ``data.diagnostics_vec`` file (in the example file nVEC_mask=20)
    - nSURF_mask: number of surface boundary mask used in the ``data.diagnostics_vec`` file (in the example file nSURF_mask=1)

---
## II. Compile and run the simulation
Before proceding copy the mask files on the supercomputer capability in the folder: "regions/configs/parent_run/dv".\
Below is an exemple to extract vectors from v05 ECCO-Darwin.\
On the supercomputer (example on Pleiades) run:
```
cd darwin3/regions
rm -r parent_run (if extsting from STEP1)
mkdir parent_run
cd parent_run
mkdir build run
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

<u>Note:</u> At the end of the simulation you will get a binary file for every parameter set in ``data.diagnostics_vec``, each containing the number of iterations chosen for the parameter.

