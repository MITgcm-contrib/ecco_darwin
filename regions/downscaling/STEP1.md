# Generate regional cut out input files

---
## Preliminary information
**Requirement**: 
> Before proceding to the following intructions, you will need to complete steps in README.\
> The anaconda python environment will be needed to run the codes.

At the end of this step, you will have generated several model input files needed to run ECCO model with ``diagnostics_vec``.

---
## I. Generate the ``.mitgrid`` file 

``.mitgrid`` file contains all the grid information of the regional domain.\
Runing ``gen_mitgrid.py`` in the ``utils`` folder, you will generate the ``.mitgrid`` file\
Below are the instructions to run the code in a terminal with the envrironment previously created:

```
conda activate downscaling
cd ecco_darwin/regions/downscaling/utils/
python3 gen_mitgrid.py -d /path/to/save/the/grid -n name_of_the_region\ 
                       -c left_lon right_lon left_lat right_lat -r dlon dlat\
                       -e espg -v
```  
**Example:**
```
python3 gen_mitgrid.py -d ecco_darwin/regions/downscaling/ -n NorthSlope\ 
                       -c -155.071750 -143.017750 69.940900 72.256850\
                       -r 0.01435 0.00455 -e 4326 -v
```
To get more information about the options required for this code run ``gen_mitgrid.py -h``.\
Below are additional details about the option:
> - -d: The directory where to store the mitgrid
> - -n: The name of your region. It will be the name of the mitgrid
> - -c: The coordinates of the corners of your region (for the center of the cell). It can be in any ESPG.
> - -r: The resolution of the grid in the format of the chosen ESPG
> - -e: The ESPG of the coordinates input
> - -v: Verbose

<u>Note</u>: The code will return the output shape information (n_rows, n_cols). Note this information it will be needed on the next step.

---
## II. Generate the bathymetry file 

**Requirement**: 

> The code generating the bathymetry file uses GEBCO data to get the bathymetry information.\
> This step then require to first download the [netcdf GEBCO  dataset](https://www.gebco.net/data_and_products/gridded_bathymetry_data/#global).\
> <u>Note:</u> you can use any bathymetry file as long as depth variable is named “elevation”.

Runing ``gen_bathy.py`` in the ``utils`` folder, you will generate the bathymetry file.\
Below are the instructions to run the code in a terminal:

```
conda activate downscaling
cd ecco_darwin/regions/downscaling/utils/
python3 gen_bathy.py -d /path/to/save/the/grid -g /path/to/gebco/file/ 
                     -n name_of_the_region -s n_rows n_cols -cs min_Scell_size max_Scell_size\
                     -cw first_surface_cell_height second_surface_cell_height\
                     -wp row_wet_central_cell col_wet_central_cell -v 
```
**Example:**
```
python3 gen_bathy.py -d ecco_darwin/regions/downscaling/ -g /data/gebco/gebco.nc\
                     -n NorthSlope -s 840 510 -cs 0.3 1.0 -cw 1.0 1.14\
                     -wp 421 254 -v
``` 

To get more information about the options required for this code run ``gen_bathy.py -h``.\
Here are additional details about the option:
> - -d: The directory where the mitgrid is stored and where the bathymetry file will be stored
> - -g: Path to the GEBCO netcdf file.
> - -n: The name of your region. It will be the name of the bathymetryfile.
> - -s: Size of the downscaled domain. It correspond to n_rows and n_cols information from mitgrid code.
> - -cs: Minimum fraction size of a cell (hFacMin) and minimum dimension size of a cell (hFacMinDr). This is two paramters you will set later for in the ``data`` file of your downscaled model. More information on these parameters [here](https://darwin3.readthedocs.io/en/latest/algorithm/vert-grid.html#topography-partially-filled-cells).
> - -cw: Height of the 2 first surface layers of the downscaled model. This is the 2 first numbers you will set later for the deltaR parameter in the ``data`` file of your model.  
> - -wp: Indexes (row and column) of the central cells of the downscaled domain.  **Warning**: This should be a wet cell, so select the wet cell closest to the center.
> - -v: Verbose

---
## III. Generate tiles for mutliprocessing

ECCO simulations run on multiprocessors by divinding the regional domain in tiles for which on processor will be dedicated.\
Runing ``mitgrid2tiles.py`` in the ``utils`` folder, you will split the mitgrid onto the processors depending on your choice of sNx & sNy [(see MITgcm doc)](https://darwin3.readthedocs.io/en/latest/getting_started/getting_started.html#customizing-the-model-configuration-code-parameters-and-compilation-options) and generate the tile files.\
Below are the instructions to run the code in a terminal:

```
conda activate downscaling
cd ecco_darwin/regions/downscaling/utils/
python3 mitgrid2tiles.py -d /path/to/save/the/grid -n name_of_the_region\
                         -s n_rows n_cols -p sNx sNy
```
**Example:**
```
 python3 mitgrid2tiles.py -d ecco_darwin/regions/downscaling/ -n NorthSlope \
                          -s 840 510 -p 30 30
```
To get more information about the options required for this code run ``mitgrid2tiles.py -h``. Here are additional details about the option:
> - -d: The directory where to store the mitgrid
> - -n: The name of your region
> - -s: Size of the downscaled domain. It correspond to n_rows and n_cols information from mitgrid code.
> - -p: Number of processors you want to slpit the region on the width (sNx) and the height (sNy)

---
## IV. Generate grid file (netcdf)

The following steps will permit to generate a netcdf file containing the different grid information as they will be interpreted by the model for you regional cut out.\
<u>Note:</u> For this step you will need to run the ECCO-Darwin model on the supercomputer for 1 time step. 

### a. Set the cut-out configuration parameters

Here you will to start setting up the configuration parameters for your downscalled model.\
On your supercomputer:
```
cd downscalling/regions/configs/
cp -r ecco_darwin/regions/downscaling/gen_ncgrid .
cd gen_ncgrid/code/
vim SIZE.h
    ---> Modify sNx to Nr parameters according to your choices in III.
cd ../input/
vim data
    ---> Modify parameters according to your downscalled setup (see instructions below)
    ---> Add the name of you bathyfile: "bathyFile="
cd ..

==============================
for Pleiades users only:
vim job_downsc_ivy 
    ---> Modify the job file according to you setup 
```

<u>Note:</u> See how to modify the job file [here](https://www.nas.nasa.gov/hecc/support/kb/running-jobs-with-pbs-121/)\
In the data file change the parameters according to the downsacled model configuration you want to create.\
Below are the details of importante parameters to set:
> - tRef: reference vertical profile for potential temperature (size = Nr the number of vertical levels)
> - sRef: reference vertical profile for salinity/specific humidity (size = Nr the number of vertical levels)
> - hFacMinDr: Minimum dimension size of a cell (same as bathy file option -cs)
> - hFacMin: Minimum fraction size of a cell (same as bathy file option -cs)
> - delR: Vertical grid spacing (size = Nr the number of vertical levels)
> - bathyFile: name of the bathymetry file generated previously

### b. Compile and run the regional model on 1 step

**Requirement**:
> The bathymetry and tiles files generated in eralier should be accessible for the simulation.

Uplpoad ``regnm_bathymetry.bin``, ``regnm.mitgrid`` files and ``tiles`` folder on your supercomputer capability in the following folder "downscalling/regions/configs/".\
Then follow these instructions (example on Pleiades):

```
cd ../../darwin3/regions/
mkdir gen_ncgrid
cd gen_ncgrid
mkdir run build
cd build
```
> - Compile the code
```
module purge
module load comp-intel mpi-hpe/mpt hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt python3/3.9.12
../../../../tools/genmake2 -of ../../../../tools/build_options/linux_amd64_ifort+mpi_ice_nas\
-mpi -mo '../../../../regions/configs/gen_ncgrid/code/'
make depend
make -j 16
```
> - Run the model
```
cd ../run
ln -sf ../build/mitgcmuv .
ln -sf ../../../../regions/regnm_bathymetry.bin .
ln -sf ../../../../regions/tiles/* .
cp ../../../../regions/configs/gen_ncgrid/input/* .
cp ../../../../regions/configs/gen_ncgrid/job_downsc_ivy .
qsub job_downsc_ivy
```

### c. Stitch the grid tiles in a netcdf file

**Requirement**:
> In part b. you have generated mnc_XXX folders. Download them back to your local machine (so it is accessible by anaconda environment).

Runing ``stitch_ncgrid.py`` in the ``utils`` folder, you will generate the netcdf grid file\
Below are the instructions to run the code in a terminal:

```
conda activate downscaling
cd ecco_darwin/regions/downscaling/utils/
python3 stitch_ncgrid.py -d /path/to/save/the/grid -n name_of_the_region\
                         -z Nr -s n_rows n_cols -p sNx sNy
```
**Example:**
```
 python3 stitch_ncgrid.py -d ecco_darwin/regions/downscaling/ -n NorthSlope \
                          -z 81 -s 840 510 -p 30 30
```
To get more information about the options required for this code run ``stitch_ncgrid.py -h``. Here are additional details about the options:
> - -d: The directory where to store the mitgrid
> - -n: The name of your region
> - -z: Number of vertical levels on the regional domain.
> - -s: Size of the downscaled domain. This corresponds to n_rows and n_cols information from mitgrid code.
> - -p: Number of processors you want to split the region using width (sNx) and height (sNy)

---
## V. Generate mask files for ``diagnostics_vec``

### a. Compile and run ECCO global state estimate on 1 step
The following steps will permit to generate the grid files from the ECCO global model ("parent_run") required to generate the masks files.\
<u>Note:</u> For this step you will need to run the ECCO-Darwin model on the supercomputer for 1 time step. 

Below is an exemple with ECCO-Darwin v5. On the supercomputer (example on Pleiades) run:
```
cd downscalling/darwin3/regions/parent_run
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
cp ../../../../regions/configs/parent_run/input/* .
vim data # set debugLevel = 1 (Important here to extract the grid files)
         # set endtime to 24000
# modify job_ECCO_darwin as needed
qsub job_ECCO_darwin
```
<u>Note:</u> You can modify ``job_ECCO_darwin`` and use the debug loop to go faster here.

Once the run finished copy following grid files generated on your local machine (so it is accessible by anaconda environment).
on your local machine:
```
cd /path/where/region/files/are/strored/
mkdir parents/outputs/grid/
cp XC* YC* DXC* DYC* AngleCS* AngleSN* hFacC* DRF* parents/outputs/grid/.
cp bathy270_filled_noCaspian_r4 parents/outputs/grid/.
```

### c. Generate the masks files
Runing ``gen_dvmasks.py`` in the ``utils`` folder, you will generate the dv mask files required for Step 2.\
Below are the instructions to run the code in a terminal:

```
conda activate downscaling
cd ecco_darwin/regions/downscaling/utils/
python3 gen_dvmasks.py -d /path/to/save/the/masks -n name_of_the_region\
                         -bfl bathy_file_name -bnd EWNS -r reso -v
```
**Example:**
```
python3 gen_dvmasks.py -d ecco_darwin/regions/downscaling/ -n NorthSlope\
                       -bfl bathy270_filled_noCaspian_r4 -bnd ENW -r 18 -v
```
To get more information about the options required for this code run ``gen_dvmasks.py -h``. Here are additional details about the options:
> - -d: The directory where "parent" folder is stored or should be stored (parent folder nust have the grid files in it)
> - -n: The name of your region
> - -bfl: name of the ECCO global model bathymetry file to use
> - -bnd: Open boundaries of your domain where you want to extract boundary conditions. Can be either 'E' (East), 'W' (West), 'N' (North), 'S' (South)
> - -r: Horizontal grid resolution in meters of the global model in km (**Warning**: This must be an integer). This information is used the searching radius for matching global model coordinates with regional setup boundaries coordinates.

