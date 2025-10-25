# Step III: Generate the Regional setup

## I. Preliminary information
**Requirement**: Before proceding to the following intructions, you will need to complete steps in the README, STEP 1 and STEP 2.

Following these instructions you will convert the output vector from STEP 2 into initial and boundary conditions of your cut out.

## II. Get output from parent grid

Following these instructions you will get the output files from the parent grid simulation (STEP 2) that are going to be used as input file in the following python routines.

Go in the run directory where you ran the global model simulation (here ECCO-Darwinv5 LLC270):

```
mkdir dv/outputs
cd dv/outputs
mkdir grid pickups OBCS
cd ..
```

> - Get the grid files (**Warning:** These files are available only if you set debugLevel=1 in data – see STEP 2)

```
cp XC* YC* AngleCS* AngleSN* hFacC* DRF* dv/outputs/grid/.
```
> - Get the initial and boundary (dv) conditions files
```
cp pickup*.XXXXXXXXXX.* dv/outputs/pickups/.
mv dv OBCS
mv OBCS dv/outputs/.
```

**Note**: replace "XXXXXXXXXX" by the iteration of the pickup file you want to start the regional model (Example: In ECCO-Darwinv5: pickup*.0000683784* corresponds to the begining of year: 1992+(683784xtimestep)/(365.25x86400) = 2018) with timestep=1200

Once all the files moved in the folder, make the output_dv folder accessible by the python3 codes available in ``utils``. **Warning:** The folder can be very heavy depending on the resolution of the regional model.

### III. Generate the Initial conditions

Following these instructions you will generate the pickup files necessary to run your regional setup later. Below are the instructions to run the code in a terminal:

```
conda activate downscaling
cd ecco_darwin/regions/downscaling/utils/
python3 gen_pickups.py -d /path/to/save/the/grid -n name_of_the_region\
                      -i teration_nubmer -sg Sigma_gaussian_filter\
                      -v -nc             
```
**Example:**
```
 python3 gen_pickups.py -d ecco_darwin/regions/downscaling/ -n NorthSlope \
                       -i 683784 -v -nc
```
To get more information about the options required for this code run ``gen_pickup.py -h``. Here are additional details about the options:
> - -d: The directory where to store the mitgrid
> - -n: The name of your region
> - -i: iteration number of the pickup file chosen. (see section II.)
> - -sg: Sigma for the Gaussian Filter applied to the interpolation. This parameter is **optional**, if not prescribed sigmaG will take the value of the average size of the grid cell of the parent grid in the selected region.
> - -v: verbose
> - -nc: generate a netcdf file with all the forcing matrices inside (**optional**). This is for a control purpose only. It is not required to run the simulation. 

## III. Generate the Boundary conditions




