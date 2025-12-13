# Generate the regional set-up forcing conditions

---
## Preliminary information
**Requirement**: 
> Before proceding to the following intructions, you will need to complete steps in the README, STEP 1, and STEP 2.\
> The anaconda python environment will be needed to run the code.

At the end of this step, you will have generated the initial and boundary conditions required to run the regional set-up.

---
## I. Generate the initial conditions

**Requirement**:
> First you need to copy the parent model files with initial and boundary conditions to your local machine (see code below). 

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

<u>Note:</u>: replace "XXXXXXXXXX" by the iteration of the pickup file you want your regional model to start with\
**Example:** In v05 ECCO-Darwin: pickup*.0000683784* corresponds to the begining of year: 1992 + (683784 x timestep) / (365.25 x 86400) = 2018) with timestep = 1200

---

Runing ``gen_pickups.py`` in the ``utils`` folder, you will generate the initial conditions for your regional set-up (``pickup`` files)\
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
> - -d: The directory where mitgrid is stored.
> - -n: The name of your region.
> - -i: iteration number of the pickup file chosen (see requirements).
> - -sg: sigma for the Gaussian Filter applied to the interpolation. This parameter is **optional**, if not prescribed, sigmaG will take the value of the average size of the grid cell of the parent grid in the selected region.
> - -bgc: generate Darwin pickup files if runing the ECCO-Darwin model.
> - -v: verbose.
> - -nc: generate a NetCDF file containing all forcing matrices (**optional**). This is for a control purpose only, it is not required to run the simulation. 

---
## II. Generate the Boundary conditions

Runing ``gen_obcs.py`` in the ``utils`` folder, you will generate the booundary conditions for your regional set-up (``OBCS`` files)\
Below are the instructions to run the code in a terminal with the anaconda envrironment:

```
conda activate downscaling
cd ecco_darwin/regions/downscaling/utils/
python3 gen_obcs.py -d /path/to/regional/files -n name_of_the_region\
                    -bnd EWNS -i iter_list  -bgc -v          
```
**Example:**
```
 python3 gen_obcs.py -d ecco_darwin/regions/downscaling/ -n NorthSlope \
                     -bnd EWN -i 683785 710065 -bgc -v
```

To get more information about the options required for this code run ``gen_pickup.py -h``. Here are additional details about the options:
> - -d: The directory where mitgrid is stored.
> - -n: The name of your region.
> - -bnd: Open boundaries of your domain where you need to generate boundary conditions >> 'E' (East), 'W' (West), 'N' (North), 'S' (South)
> - -i: iteration of the obcs files extracted from diagnostics_vec.
> - -bgc: generate Darwin pickup files if runing the ECCO-Darwin model.
> - -v: verbose.

---

**CONGRATULATIONS!!** You have now generated the pickup and boundary files required to run your regional setup.\

Considering the vast capabilities of ECCO/MITgcm and the amount of parameters that can be tuned, we let the users generate their own regional setup.\
You can however start from the ``gen_mod`` folder used in STEP1 as a starting point.\
You can also check the ``gen_mod/reg_example_files``.\
We put important files to setup the regional model such as ``data.obcs`` where you set the OBCS files to read.

