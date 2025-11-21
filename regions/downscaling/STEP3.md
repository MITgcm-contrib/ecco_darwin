# Generate the regional setup

---
## Preliminary information
**Requirement**: 
> Before proceding to the following intructions, you will need to complete steps in the README, STEP 1 and STEP 2.\
> The anaconda python environment will be needed to run the codes.

At the end of this step, you will have configured your regional setup and generated the initial and boundary conditions required to run it.

---
## I. Generate the Initial conditions

**Requirement**:
> First you need to copy on your local machine the parent model files with initial and boundary conditions information (see code below) 

```
cd /path/where/region/files/are/strored/.../parent/ouput
mkdir pickups OBCS
===========================
from super computer to local machine:
cp pickup*.XXXXXXXXXX.* ouput/pickups/.
cp -r dv ouput/.
===========================
mv dv OBCS
```

<u>Note:</u>: replace "XXXXXXXXXX" by the iteration of the pickup file you want your regional model to start\
**Example:** In ECCO-Darwinv5: pickup*.0000683784* corresponds to the begining of year: 1992+(683784xtimestep)/(365.25x86400) = 2018) with timestep=1200

---

Runing ``gen_pickups.py`` in the ``utils`` folder, you will generate the initial conditions for your regional setup (``pickup`` files)\
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
To get more information about the options required for this code run ``gen_pickup.py -h``. Here are additional details about the options:
> - -d: The directory where to store the mitgrid
> - -n: The name of your region
> - -i: iteration number of the pickup file chosen. (see requirements)
> - -sg: Sigma for the Gaussian Filter applied to the interpolation. This parameter is **optional**, if not prescribed sigmaG will take the value of the average size of the grid cell of the parent grid in the selected region.
> - -bgc: generate darwin pickup files if runing the ECCO-Darwin model.
> - -v: verbose
> - -nc: generate a netcdf file with all the forcing matrices inside (**optional**). This is for a control purpose only. It is not required to run the simulation. 

---
## II. Generate the Boundary conditions

**Work in Progress**


