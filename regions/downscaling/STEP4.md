# Build and run your downscaled regional setup (in progress)

---
## Preliminary information
**Requirement**: 
> Before proceding to the following intructions, you will need to complete steps in the README, STEP 1, STEP 2 and STEP 3.

---
## I. Prepare code directories

First you need to build code and ``code_darwin`` directories containing specific routines of the model in the ecco_darwin repo you cloned.\
You can follow or copy the architecture used in existing regional setups:

```
mkdir ecco_darwin/regions/YOURSETUP
cp -r ecco_darwin/regions/GoM_1km/code ecco_darwin/regions/YOURSETUP/.
cp -r ecco_darwin/regions/GoM_1km/code_darwin ecco_darwin/regions/YOURSETUP/.
```

Copy your previously configured SIZE.h file from STEP 1, Section IV, b. It must match your downscaled setup dimensions.

```
cp downscalling/darwin3/regions/gen_ncgrid/build/SIZE.h ecco_darwin/regions/YOURSETUP/.
```

---
## II. Prepare namelists

Then, you need to prepare data namelists files as inputs to copy when preparing the model for run. Once again, you can follow and copy examples from existing regional setups:

```
cp -r ecco_darwin/regions/GoM_1km/inputs ecco_darwin/regions/YOURSETUP/.
cp -r ecco_darwin/regions/GoM_1km/inputs_darwin ecco_darwin/regions/YOURSETUP/.
```

### a. data namelist (ecco_darwin/regions/YOURSETUP/code/data)

Like for SIZE.h, you can re-use the data namelist file configured for STEP 1, Section IV, b.

Set your iteration number to zero, remove pickupSuff line and set deltaT to samll value to respect CFL conditions\ 
in # Time stepping parameters (&PARM03):

```
 nIter0=0,
 deltaT=30.,
```

Make sure tRef, sRef and delR are matching your vertical grid and add the following lines in # Gridding parameters (&PARM04):

```
 xgOrigin =   -97.9849,
 ygOrigin =  16.7545,
 delX     = 2130*0.009381140787627662
 delYfile = 'delYFile',
```

with xgOrigin and ygOrigin, the minimum longitude and latitude of your domain, delX, the grid spacing in x-direction (use same value as STEP 1, Section I).\ delYFile is a simple binary file containing the grid spacing in the y-direction in a vector format(see Python function gen_delYFile.py).

Finally, add the following lines in # Input datasets (&PARM05) to restart the model from a specific iteration in your pickup files generated in STEP 3, Section I:

```
 bathyFile      ='GoM_1km_bathymetry_8.bin',
 hydrogThetaFile = 'pickup_THETA.0000026352.data',
 checkIniTemp=.FALSE.,
 hydrogSaltFile  = 'pickup_SALT.0000026352.data',
 checkIniSalt=.FALSE.,
 uVelInitFile    = 'pickup_U.0000026352.data',
 vVelInitFile    = 'pickup_V.0000026352.data',
 pSurfInitFile   = 'pickup_ETAN.0000026352.data',
````

### b. data.cal namelist (ecco_darwin/regions/YOURSETUP/inputs/data.cal)

Make sure to set startDate_1 in data.cal to the corresponding iteration you would like to restart the model from (suffix in the pickup files set in data).\
**Example:** In v05 ECCO-Darwin: pickup*.0000683784* corresponds to the begining of the year: 1992 + (683784 x timestep) / (365.25 x 86400) = 2018) with timestep = 1200. Here, the timestep is from the parent run.

### c. data.exf namelist (ecco_darwin/regions/YOURSETUP/inputs/data.exf)

This is where you configure the forcing to be read by the model. Make sure to either provide global forcing with their interpolation parameters or to provide forcings adapted to your grid in &EXF_NML_02 and &EXF_NML_04 (see https://mitgcm.readthedocs.io/en/latest/phys_pkgs/exf.html). It is also in data.exf that you configure OBCs date, time and period parameters in &EXF_NML_OBCS. Be sure that your OBCs either start before or at the same time than your startDate_1 in data.cal. obcsXXstartdate1 is given by iteration number used in STEP 3, Section II, gen_obcs.py.
Make sure to have lines for each of your OBCs (here only 3 for East, North and South):

```
 &EXF_NML_OBCS
  obcsNstartdate1   = 19921231,
  obcsNstartdate2   = 002000,
  obcsNperiod       = 2629800.0,
\#
  obcsSstartdate1   = 19921231,
  obcsSstartdate2   = 002000,
  obcsSperiod       = 2629800.0,
\#
  obcsEstartdate1   = 19921231,
  obcsEstartdate2   = 002000,
  obcsEperiod       = 2629800.0,
```

### d. data.obcs namelist (ecco_darwin/regions/YOURSETUP/inputs/data.obcs)

Here, you provide the size and direction of your OBCS and corresponding files for each variables (see https://mitgcm.readthedocs.io/en/latest/phys_pkgs/obcs.html).

### e. data.ggl90 namelist (ecco_darwin/regions/YOURSETUP/inputs/data.ggl90)

Finally, provide name of pickup_ggl90 file generated in STEP 3, Section I:

```
GGL90TKEFile = 'pickup_ggl90.0000026352.data',
```

### f. data.pkg namelist (ecco_darwin/regions/YOURSETUP/inputs/data.pkg)

This namelist defines packages used for the model integration. Make sure it looks like this:
```
\# Packages
 &PACKAGES
 useCAL         = .TRUE.,
 useEXF         = .TRUE.,
 useOBCS        = .TRUE.,
 useDiagnostics = .TRUE.,
 useGCHEM       = .TRUE.,
 usePTRACERS    = .TRUE.,
 useGGL90 = .TRUE.,
 /
```

### g. Additional namelists for Darwin run (ecco_darwin/regions/YOURSETUP/inputs_darwin/)

---
## III. Prepare forcing files


---
## IV. Compile and run


