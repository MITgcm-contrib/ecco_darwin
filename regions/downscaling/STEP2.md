# Step II: Extracting regional model boundary information

## I. Preliminary information
**Requirement**: Before proceding to the following intructions, you will need to complete steps in the README and in STEP1.

This step explains how to exctract vectors from any llc global configuration along the boundaries of the required regional model using the ``diagnostic_vec`` package. This package was designed to output model diagnostics from llc global model in a subset of the model domain e.g. along a vector (or "vec"). This package is not included in the official MITgcm realease but can be easly merged to it. More information on ``diagnostic_vec`` package at https://github.com/mhwood/diagnostics_vec (*Credit*: Mike Wood). 

## II. Prepare the simulation for ``diagnostic_vec`` extraction

### a. Turn on ``diagnostic_vec`` package

> - Enable the package in ``packages.conf``
```
cd config/code_darwin
vim packages.conf
##### add the following line to package.conf #####
diagnostic_vec
```
> - Turn on the package in ``data.pkg``
```
cd ../input
vim data.pkg
##### add the following line to data.pkg #####
 useDiagnostics_vec=.TRUE.,
```

### b. Generate a ``data.diagnostics_vec`` parameter file
> - copy ``data.diagnostics_vec`` file into the ``input`` directory
```
cd config/input
cp ../../ecco_darwin/regions/downscaling/utils/data.diagnostics_vec .
```
> - Modify the ``data.diagnostics_vec`` file according to the specificities of your domain and your requirements.

**Note:** You can modify data.diagnostic file with only the diagnostics you want to save. This won't affect diagnostic_vec and the fewer diagnostics saved the faster the simulation

### c. Set the compile time ``DIAGNOSTICS_VEC_SIZE.h`` file
> - copy ``DIAGNOSTICS_VEC_SIZE.h`` file into the ``code_darwin`` directory
```
cd ../code_darwin
cp ../../darwin3/pkg/diagnostics_vec/DIAGNOSTICS_VEC_SIZE.h .
```
> - Modify the ``DIAGNOSTICS_VEC_SIZE.h`` file as follows:
    - VEC_points: correspond to the max number of points that can be stored in a diagnostics_vec mask (the shorter this list, the less time is spent on I/O).
    - nVEC_mask: number of lateral boundary mask used in the ``data.diagnostics_vec`` file (in the example file nVEC_mask=20)
    - nSURF_mask: number of surface boundary mask used in the ``data.diagnostics_vec`` file (in the example file nSURF_mask=1)

## IV. Compile and run the simulation
```
cd darwin3
mkdir build run
```
> - Compile the code
```
cd build
module purge
module load comp-intel mpi-hpe/mpt hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt python3/3.9.12
../tools/genmake2 -of ../../config/code/linux_amd64_ifort+mpi_ice_nas \
  -mpi -mo '../../config/code_darwin/ ../../config/code/'
make depend
make -j 16
```
> - Run the simulation
```
cd ../run
ln -sf ../build/mitgcmuv .
ln -sf /nobackupp19/dmenemen/public/llc_270/iter42/input/* .
ln -sf /nobackupp19/dmenemen/public/llc_270/ecco_darwin_v5/input/darwin_initial_conditions/* .
ln -sf /nobackupp19/dmenemen/public/llc_270/ecco_darwin_v5/input/darwin_forcing/* .
ln -sf /nobackup/dcarrol2/forcing/apCO2/NOAA_MBL/* .
ln -sf /nobackup/hzhang1/forcing/era_xx .
ln -sf /nobackup/hzhang1/pub/llc270_FWD/input/19920101/to2023/xx*42.data .
cp ../../config/input/* .
cp -r ../../config/dv .
mkdir diags
# modify job_ECCO_darwin as needed
qsub job_ECCO_darwin
```

**Note:** At the end of the simulation you will get a binary file for every parameter set in ``data.diagnostics_vec`` file each containing the number of iterations chosen for the parameter.

