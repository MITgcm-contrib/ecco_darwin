# Step I: Generate regional model input files

## I. Preliminary information
**Requirement**: Before proceding to the following intructions, you will need to complete steps in README. The anaconda python environment will be needed to run the codes.

Following these instructions you will be able to generate several model input files needed to run the model with ``diagnostics_vec`` in the following steps.

## II. Generate the ``.mitgrid`` file 

The ``.mitgrid`` files contain all the grid information of your regional domain.
The python code ``gen_mitgrid.py`` in the ``utils`` folder will help you generate this file. Below are the instructions to run the code in a terminal with the anaconda prviously created:

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
To get more information about the options required for this code run ``gen_mitgrid.py -h``. Here are additional details about the option:
> - -d: The directory where to store the mitgrid
> - -n: The name of your region. It will be the name of the mitgrid
> - -c: The coordinates of the corners of your region (for the center of the cell). It can be in any ESPG.
> - -r: The resolution of the grid in the format of the chosen ESPG
> - -e: The ESPG of the coordinates input
> - -v: Verbose

**Note**: The code will return the output shape information (n_rows, n_cols). Note this information it will be needed on the next step.

## II. Generate the bathymetry file 

### a. Download GEBCO bathymetry files.

The code generating the bathymetry file uses GEBCO data to get the bathymetry information. This step then require to first download the [netcdf GEBCO  dataset](https://www.gebco.net/data_and_products/gridded_bathymetry_data/#global).

### b. Generate the bathymetry file

The python code ``gen_bathy.py`` in the ``utils`` folder will generate the bathymetry file. Below are the instructions to run the code in a terminal:

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

To get more information about the options required for this code run ``gen_bathy.py -h``. Here are additional details about the option:
> - -d: The directory where the mitgrid is stored and where the bathymetry file will be stored
> - -g: Path to the GEBCO netcdf file.
> - -n: The name of your region. It will be the name of the bathymetryfile.
> - -s: Size of the downscaled domain. It correspond to n_rows and n_cols information from mitgrid code.
> - -cs: Minimum fraction size of a cell (hFacMin) and minimum dimension size of a cell (hFacMinDr). This is two paramters you will set later for in the ``data`` file of your downscaled model. More information on these parameters [here](https://darwin3.readthedocs.io/en/latest/algorithm/vert-grid.html#topography-partially-filled-cells).
> - -cw: Height of the 2 first surface layers of the downscaled model. This is the 2 first numbers you will set later for the deltaR parameter in the ``data`` file of your model.  
> - -wp: Indexes (row and column) of the central cells of the downscaled domain.  **Warning**: This should be a wet cell, so select the wet cell closest to the center.
> - -v: Verbose

## III. Generate the tile files for mutliprocessing

**Note:** If you plan to run the model on a single processor you can jump to the following step.

ECCO-Darwin model allow multiprocessing by divinding the regional domain in tiles for which on processor will be dedicated. The python code ``mitgrid2tiles.py`` in the ``utils`` folder will help you split the mitgrid onto the processors depending on your choice of sNx & sNy [(see MITgcm doc)](https://darwin3.readthedocs.io/en/latest/getting_started/getting_started.html#customizing-the-model-configuration-code-parameters-and-compilation-options) and generate the tilefiles. Below are the instructions to run the code in a terminal:

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

## IV. Generate a netcdf grid file with the model grid fields

Here, we will generate a netcdf file containing the different grid information as they will be interpreted by the model for you regional cut out.\
**Note:** For this step you will need to run the ECCO-Darwin model on the supercomputer for 1 time step. 

### a. Set the ecco-darwin v5 configurationparameters

Here you will to start setting up the configuration parameters for your downscalled model.\
**Requirement**: The bathymetry and tiles files generated in II. should be accessible for the simulation.

```
cd downscalling
cp -r ecco_darwin/regions/downscaling/gen_ncgrid .
cd gen_ncgrid/code/
vim SIZE.h
    ---> Modify sNx to Nr parameters according to your choices in III.
cd ../input/
vim data
    ---> Modify parameters according to your downscalled setup (see instructions below)
    ---> Add the name of you bathyfile: "bathyFile="
cd ..
vim job_downsc_ivy (for Pleiades users)
    ---> Modify the job file ccording to you setup 
```

**Note:** See how to modify the job file [here](https://www.nas.nasa.gov/hecc/support/kb/running-jobs-with-pbs-121/)\
In the data file change the parameters according to the downsacled model configuration you want to create. Below are the details of importante parameters to set:
> - tRef: reference vertical profile for potential temperature (size = Nr the number of vertical levels)
> - sRef: reference vertical profile for salinity/specific humidity (size = Nr the number of vertical levels)
> - hFacMinDr: Minimum dimension size of a cell (same as bathy file option -cs)
> - hFacMin: Minimum fraction size of a cell (same as bathy file option -cs)
> - delR: Vertical grid spacing (size = Nr the number of vertical levels)
> - bathyFile: name of the bathymetry file generated previously

### b. Compile and run the model on 1 step

```
cd ../../darwin3/
mkdir reg_mod/gen_nc
cd reg_mod/gen_nc
mkdir run build
cd build
```
> - Compile the code
```
module purge
module load comp-intel mpi-hpe/mpt hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt python3/3.9.12
../../../tools/genmake2 -of ../../../tools/build_options/linux_amd64_ifort+mpi_ice_nas\
-mpi -mo '../../../../config/gen_ncgrid/code/
make depend
make -j 16
```
> - Run the model
```
cd ../run
ln -sf ../build/mitgcmuv .
ln -sf ../../../../config/bathymetry_file.bin .
ln -sf ../../../../config/tiles/* .
cp ../../../../config/gen_ncgrid/input/* .
cp ../../../../config/job_downsc_ivy .
qsub job_downsc_ivy
```

### d. Stitch the grid tiles in a netcdf file
**Requirement**: You should have the mnc folder generated by the model accessible for the python environment.

Below are the instructions to run the code in a terminal to generate the netcdf grid file:

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

### e. Generate mask files for ``diagnostics_vec`` package
**Requirement**: you need to create the netcdf grid file (see above).

Below are the instructions to run the code in a terminal to generate the dv mask files:

```
conda activate downscaling
cd ecco_darwin/regions/downscaling/utils/
python3 gen_dvmasks.py -d /path/to/save/the/grid -n name_of_the_region\
                         -b EWNS -r reso
```
**Example:**
```
 python3 gen_dvmasks.py -d ecco_darwin/regions/downscaling/ -n NorthSlope \
                          -b ENW -r 500 -v
```
To get more information about the options required for this code run ``gen_dvmasks.py -h``. Here are additional details about the options:
> - -d: The directory where to store the mitgrid
> - -n: The name of your region
> - -b: Open boundaries of your domain where you want to extract boundary conditions. Can be either 'E' (East), 'W' (West), 'N' (North), 'S' (South)
> - -r: Horizontal grid resolution in meters. **Warning**: This must be an integer.

**CONGRATULATIONS!!** You have generated all the files necessary to proceed to the 2nd step
