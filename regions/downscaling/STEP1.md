# Step I: Extracting regional model boundary information

## I. Preliminary information
**Requirement**: Before proceding to the following intructions, you will need to complete steps in README.

Following this instructions you will be able to exctract vectors from any llc global configuration along the boundaries of the required regional model using the ``diagnostic_vec`` package. This package was designed to output model diagnostics from llc global model in a subset of the model domain e.g. along a vector (or "vec"). This package is not include in the official realease of MITgcm but can be easly merged to it. More information on diagnostic_vec package at https://github.com/mhwood/diagnostics_vec.

## II. Merging Instructions
To perform the following steps you need to connect on a capability able to run a llc global simulation *(Below, we detail the instructions to run ECCO-Darwin model on NASA Pleiades supercomputer)*
> - clone darwin checkpointv67x github
```
mkdir downscalling 
cd downscalling
git clone https://github.com/darwinproject/darwin3
cd darwin3
git checkout 24885b71
```
> - clone diagnostic_vec github
```
cd ..
git clone https://github.com/mhwood/diagnostics_vec.git
```
> - Merge diagnostic_vec package to darwin3
```
cd diagnostics_vec/utils/
python3 copy_doc_files_to_MITgcm.py -m ../../darwin3/
python3 copy_pkg_files_to_MITgcm.py -m ../../darwin3/
python3 copy_verification_files_to_MITgcm.py -m ../../darwin3/
```

## III. Set diagnostic_vec parameterization files
Below are the instructions to run the ECCO-Darwin v5 [llc270] simulation with the ``diagnostics_vec`` package and extract a choosen region.

### a. Set the darwin5 configuration to run
```
cd downscalling
git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
mkdir config
cp -r ecco_darwin/v05/llc270/code config/.
cp -r ecco_darwin/v05/llc270/code_darwin config/.
cp -r ecco_darwin/v05/llc270/input config/.
```

### b. Turn on ``diagnostic_vec`` package
> - Enable the package in ``packages.conf``
```
cd config/code_darwin
vim package.conf
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

### c. Generate "masks" files to isolate the region
```

```

### d. Generate a ``data.diagnostics_vec`` parameter file
```

```
**Note:** You can modify data.diagnostic file with only the diagnostics you want to save. This won't affect diagnostic_vec and the fewer diagnostics saved the faster the simulation

### e. Set the compile time ``DIAGNOSTICS_VEC_SIZE.h`` file
```
cp ../../darwin3/pkg/diagnostics_vec/DIAGNOSTICS_VEC_SIZE.h .
```

### f. Compile and run the simulation
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
mkdir diags
# modify job_ECCO_darwin as needed
qsub job_ECCO_darwin
```

