# Build and run your downscaled regional set-up

---
## Preliminary information
**Requirement**: 
> Before proceeding to the following instructions, you will need to complete steps in the README, STEP 1, STEP 2 and STEP 3.

---
## I. Prepare code directories

First, you need to build the code and code_darwin directories containing specific routines of the model in the ecco_darwin repo you cloned.\
You can follow or copy the architecture used in existing regional set-ups:

```
mkdir ecco_darwin/regions/YOURSETUP
cp -r ecco_darwin/regions/GoM/GoM_1km/code ecco_darwin/regions/YOURSETUP/.
cp -r ecco_darwin/regions/GoM/GoM_1km/code_darwin ecco_darwin/regions/YOURSETUP/.
```

Copy the SIZE.h file you configured in STEP 1, Section IV.a (and compiled with in Section IV.b). It must match your downscaled regional set-up dimensions.

```
cp downscaling/darwin3/regions/gen_ncgrid/build/SIZE.h ecco_darwin/regions/YOURSETUP/code/.
```

Make sure the packages.conf file in ecco_darwin/regions/YOURSETUP/code/ contains the following packages:
```
cal
diagnostics
exf
gchem
generic_advdiff
ggl90
mdsio
mom_vecinv
monitor
rw
salt_plume
obcs
```
These packages will be necessary to run the downscaled model. If you want to run your model with the Darwin package, add the following packages to packages.conf:
```
ptracers
darwin
```

Make sure the version of Darwin you are using (based on your parent run) matches routines in ecco_darwin/regions/YOURSETUP/code_darwin/. See more details on compiling options for the Darwin package: https://darwin3.readthedocs.io/en/latest/phys_pkgs/darwin.html#compiling.

---
## II. Prepare namelists

Then, you need to prepare data namelists files as inputs to copy when preparing the model for run. Once again, you can follow and copy examples from existing regional set-ups:

```
cp -r ecco_darwin/regions/GoM/GoM_1km/input ecco_darwin/regions/YOURSETUP/.
cp -r ecco_darwin/regions/GoM/GoM_1km/input_darwin ecco_darwin/regions/YOURSETUP/.
```

### a. data namelist (ecco_darwin/regions/YOURSETUP/input/data)

Like for SIZE.h, you can re-use the data namelist file you configured in STEP 1, Section IV.a.

Make sure that hFacMin and hFacMinDr in the data file match the values you provided in STEP 1, Section II when generating the bathymetry file. 

Set your iteration number to zero, remove the pickupSuff line and set deltaT to a small value to respect CFL conditions in # Time stepping parameters (&PARM03):

```
 nIter0=0,
 deltaT=30.,
```

Specify the duration of the model integration with endtime or nTimeSteps. See documentation for more details (https://mitgcm.readthedocs.io/en/latest/getting_started/getting_started.html#run-start-and-duration).

Make sure tRef, sRef and delR are matching your vertical grid and add the following lines in # Gridding parameters (&PARM04):

```
 xgOrigin =   -97.9849,
 ygOrigin =  16.7545,
 delX     = 2130*0.009381140787627662
 delYfile = 'delYFile',
```

with xgOrigin and ygOrigin, the minimum longitude and latitude of your domain, delX, the grid spacing in x-direction (use grid dimension in x-direction and the grid spacing value as STEP 1, Section I). delYFile is a simple binary file containing the grid spacing in the y-direction in vector format. You can generate your own delYfile by following these instructions:
```
conda activate downscaling
cd ecco_darwin/regions/downscaling/utils/
python3 gen_delYFile.py -d /path/to/config_dir -n name_of_the_region -x n_rows -y n_cols
```
**Example:**
```
python3 gen_delYFile.py -d ecco_darwin/regions/downscaling/ -n NorthSlope -x 840 -y 510
```
> - -d: Your `config_dir` — the directory where `<name_of_the_region>.mitgrid` is stored; `delYFile` is written back into the same directory.
> - -n: The name of your region (must match the `.mitgrid` filename).
> - -x: Grid size in the x-direction (`n_rows` from STEP 1, Section I).
> - -y: Grid size in the y-direction (`n_cols` from STEP 1, Section I).

<u>Note:</u> `-x` and `-y` are both required; the script reads the raw `.mitgrid` binary and needs the domain shape to reshape it.

The delYFile will be saved in /path/to/config_dir.

Finally, add the following lines in # Input datasets (&PARM05) to restart the model from a specific iteration in your pickup files generated in STEP 3, Section I:

```
 bathyFile      ='GoM_1km_bathymetry.bin',
 hydrogThetaFile = 'pickup_THETA.0000026352.data',
 checkIniTemp=.FALSE.,
 hydrogSaltFile  = 'pickup_SALT.0000026352.data',
 checkIniSalt=.FALSE.,
 uVelInitFile    = 'pickup_U.0000026352.data',
 vVelInitFile    = 'pickup_V.0000026352.data',
 pSurfInitFile   = 'pickup_ETAN.0000026352.data',
```

The files that `readBinaryPrec` governs — the bathymetry (`gen_bathy.py`), `delYFile` (`gen_delYFile.py`), and the pickup/initial-condition files (`gen_pickups.py`) — are all written in 64-bit precision, so set:
```
readBinaryPrec=64,
```

<u>Warning:</u> this does **not** apply to the OBCS files. `gen_obcs.py` writes them in **32-bit** precision, and they are read by the EXF package, which is governed by `exf_iprec_obcs` in `data.exf` (default `32`). Leave that default alone — setting `exf_iprec_obcs=64` to "match" `readBinaryPrec` will make the model misread the boundary conditions.
### b. data.cal namelist (ecco_darwin/regions/YOURSETUP/input/data.cal)

Make sure to set startDate_1 in data.cal to the corresponding iteration you would like to restart the model from (suffix in the pickup files set in data).\
**Example:** In v05 ECCO-Darwin: pickup*.0000683784* corresponds to the beginning of the year: 1992 + (683784 x timestep) / (365.25 x 86400) = 2018) with timestep of 1200 seconds. Here, the timestep is from the parent run.

### c. data.exf namelist (ecco_darwin/regions/YOURSETUP/input/data.exf)

This is where you configure the forcing fields to be read by the model. Make sure to either provide global forcings along with their interpolation parameters or to provide forcing fields adapted to your grid in &EXF_NML_02 and &EXF_NML_04 (see https://mitgcm.readthedocs.io/en/latest/phys_pkgs/exf.html). It is also in data.exf that you configure OBCS date, time and period parameters in &EXF_NML_OBCS. Be sure that your OBCS either starts before or at the same time as your startDate_1 date in data.cal. obcsXXstartdate1 is given by the iteration number used in STEP 3, Section II, gen_obcs.py.
Make sure to have lines for each of your OBCS (here only 3 for East, North and South):

```
 &EXF_NML_OBCS
  obcsNstartdate1   = 19921231,
  obcsNstartdate2   = 002000,
  obcsNperiod       = 2629800.0,
#
  obcsSstartdate1   = 19921231,
  obcsSstartdate2   = 002000,
  obcsSperiod       = 2629800.0,
#
  obcsEstartdate1   = 19921231,
  obcsEstartdate2   = 002000,
  obcsEperiod       = 2629800.0,
```

### d. data.obcs namelist (ecco_darwin/regions/YOURSETUP/input/data.obcs)

Here, you provide the size and direction of your OBCS and corresponding files for each variable (see https://mitgcm.readthedocs.io/en/latest/phys_pkgs/obcs.html).

### e. data.ggl90 namelist (ecco_darwin/regions/YOURSETUP/input/data.ggl90)

Finally, provide the name of the pickup_ggl90 file generated in STEP 3, Section I:

```
GGL90TKEFile = 'pickup_ggl90.0000026352.data',
```

### f. data.diagnostics namelist (ecco_darwin/regions/YOURSETUP/input/data.diagnostics)

You can specify your desired outputs that the model will produce during the integration in data.diagnostics. Follow documentation for more details on available diagnostics:  https://mitgcm.readthedocs.io/en/latest/outp_pkgs/outp_pkgs.html#usage-notes.

### g. data.pkg namelist (ecco_darwin/regions/YOURSETUP/input/data.pkg)

This namelist defines packages used for the model integration. Make sure it looks like this:

```
# Packages
 &PACKAGES
 useCAL         = .TRUE.,
 useEXF         = .TRUE.,
 useOBCS        = .TRUE.,
 useDiagnostics = .TRUE.,
 useGGL90 = .TRUE.,
 /
```
usePTRACERS and useGCHEM can be undefined if running WITHOUT Darwin.

### h. Additional namelists for Darwin run (ecco_darwin/regions/YOURSETUP/input_darwin/)

If you want to run your downscaled regional set-up with the Darwin package, make sure you generated initial and boundary conditions for all the ptracers and Darwin state variables in STEP 3.

Make sure the GCHEM and Ptracers packages are enabled in data.pkg:
```
 useGCHEM       = .TRUE.,
 usePTRACERS    = .TRUE.,
```

Then, modify data.ptracers according to the version of Darwin used in the parent model when proceeding to the downscaling (number of tracers, names, units).

In STEP 3, Section I, you generated initial condition pickup files for each Ptracers. In data.ptracers, add the following lines in &PTRACERS_PARM01 to restart from specific initial conditions for each Ptracers:
```
 PTRACERS_initialFile( 1)= 'pickup_pTr01.0000026352.data',
 PTRACERS_initialFile( 2)= 'pickup_pTr02.0000026352.data',
 ...
```
Make sure PTRACERS_Iter0 is absent or commented.

In data.darwin, make sure darwin_pickupSuff is absent or commented and set darwin_chlIter0 to 0. In data.darwin, you also specify forcings for iron dust deposition and atmospheric surface pCO2. These can be either already adapted to your downscaled regional set-up grid or, if provided on a different grid, interpolated during the model integration if you provide associated interpolation parameters (&DARWIN_INTERP_PARAMS, similar to data.exf https://mitgcm.readthedocs.io/en/latest/phys_pkgs/exf.html).

Modify data.diagnostics according to the Ptracers or Darwin state variables you want as outputs. See more details:
https://darwin3.readthedocs.io/en/latest/outp_pkgs/outp_pkgs.html#using-available-diagnostics.

---
## III. Prepare forcing files

### a. Forcing files for physics

To constrain your downscaled regional set-up during its integration, you need physical forcings for atmospheric surface temperature, specific humidity, precipitation, wind horizontal velocity components, freshwater runoff and short- and long-wave downwelling radiation. Provide files either already adapted to your downscaled regional set-up grid or on a different grid. If provided on a different grid, make sure to fill out the interpolation parameters (&EXF_NML_04) so they can be interpolated during the model integration. You can find a set of forcing files for the global LLC270 solution here: https://ecco.jpl.nasa.gov/drive/files/Version5/Alpha/era_xx.

Follow EXF package documentation for more details: https://mitgcm.readthedocs.io/en/latest/phys_pkgs/exf.html.

### b. Additional forcings if running with Darwin

If running the model with the Darwin package, you have to provide additional forcings for iron dust deposition and atmospheric surface pCO2. Make sure you specify these files in data.darwin in &DARWIN_FORCING_PARAMS:
```
 ironfile = 'exf/L1_exf_IRONDUST',
 ironstartdate1 = 19920101,
 ironstartdate2 = 030000,
 ironperiod = 21600.0,
#
 pCO2File = 'apCO2',
 pCO2startdate1 = 19920101,
 pCO2startdate2 = 000000,
 pCO2period     = 86400.0,
```

If darwin_useQsw, darwin_useSEAICE and darwin_useEXFwind are set to TRUE, you don't need to provide files for Photosynthetically active radiation (PAR), fraction of surface covered by ice and wind speed, as the model will derive these forcings from the EXF package. See running parameters for the Darwin package: https://darwin3.readthedocs.io/en/latest/phys_pkgs/darwin.html#running.

Once again, depending on the grid of your forcings, you might need to provide interpolation parameters in &DARWIN_INTERP_PARAMS. If provided on your downscaled regional set-up grid, you can leave iron_interpMethod and pCO2_interpMethod to zero.

Follow Darwin and EXF packages documentation for more details:
https://darwin3.readthedocs.io/en/latest/phys_pkgs/darwin.html
https://mitgcm.readthedocs.io/en/latest/phys_pkgs/exf.html

---
## IV. Compile and run

The following example assumes you intend to run the model on a high-performance computer in parallel mode. To run the model on your local machine, I recommend MITgcm_on_Mac, MITgcm_on_Ubuntu and MITgcm_on_Windows documentation in https://github.com/MITgcm-contrib/ecco_darwin/tree/master/doc.

### a. Compile

To compile the model, first, create new directories for compiling and running your set-up in your downscaling repo (or create a new one):

```
mkdir downscaling/darwin3/regions/downscaled
cd downscaling/darwin3/regions/downscaled
mkdir build run
cd build
```

Load packages, build the makefile and compile your executable with:

```
module purge
module load comp-intel/2020.4.304 mpi-hpe/mpt hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt python3/3.9.12
../../../tools/genmake2 -of ../../../../ecco_darwin/v05/llc270/code/linux_amd64_ifort+mpi_ice_nas -mo '../../../../ecco_darwin/regions/YOURSETUP/code_darwin ../../../../ecco_darwin/regions/YOURSETUP/code' -mpi
make depend
make -j 16
```

A mitgcmuv executable should be produced in your build directory.


### b. Run

Now, go into the run directory:
```
cd ../run
```

Create a symbolic link pointing to your newly created mitgcmuv executable:
```
ln -sf ../build/mitgcmuv .
```

Create a symbolic link pointing to your forcing files:
```
ln -sf /PATH2FORCINGS/* . 
```
Here, you can either create the link directly into your run directory or create a forcings directory in your run directory. Make sure you provide the correct path in data.exf and data.darwin.

Then, create a symbolic link pointing to your initial, open boundary conditions, bathymetry and delYFile files for the downscaled regional set-up. These are the files produced in STEP 1 and STEP 3, so the paths below all point into the `config_dir` you have been passing to the `utils` scripts as `-d` (see "Directory layout" in the README):
```
ln -sf /path/to/config_dir/forcings/pickups/*0000026352* .
ln -sf /path/to/config_dir/forcings/OBCS/* .
ln -sf /path/to/config_dir/YOURSETUP_bathymetry.bin .
ln -sf /path/to/config_dir/delYFile .
```

<u>Note:</u> if you ran `check_obcs_transport.py -correct` in STEP 3, Section III and want to use the corrected set, link `forcings/OBCS_corrected/*` instead of `forcings/OBCS/*`.

Copy your namelist files:
```
cp ecco_darwin/regions/YOURSETUP/input/* .
cp ecco_darwin/regions/YOURSETUP/input_darwin/* .
```

If using the diagnostics package, make sure to create a diagnostic directory that matches the paths provided in data.diagnostics in your run directory.

Finally, prepare a batch file that follows this example:
https://github.com/MITgcm-contrib/ecco_darwin/blob/master/v05/llc270/input/job_ECCO_darwin

When allocating resources in line, make sure you provide enough cores for the parallelization defined by nPx*nPy from SIZE.h file in your build directory:
```
select=<nodes>:ncpus=<cores_per_node>
```
You can reserve slightly more cores than required.

At the last line, make sure to specify the exact number of cores you want to use:
```
mpiexec -np XXX ./mitgcmuv
```
with XXX = nPx*nPy from your SIZE.h provided in your build directory.

To start the integration and run your model, submit your batch file with the qsub command:
```
qsub batch_file
```

---

Before committing supercomputer time to a long production integration, run the short
verification experiment in **STEP5**. It confirms that the model integrated stably (no NaN,
CFL respected, no runaway drift) and — just as important — that it actually ingested the
initial and boundary conditions you generated in STEP3. Both failure modes are silent: a run
can complete cleanly while ignoring your pickups, or while carrying no flow at all through
one of the open boundaries.

