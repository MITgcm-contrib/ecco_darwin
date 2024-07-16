# Generating mask file for ``diagnostics_vec``

## I. Preliminary information
**Requirement**: Before proceding to the following intructions, you will need to complete steps in README. The anaconda python environment will be needed to run the codes.

Following these instructions you will be able to generate a mask file of your regional
domain that ``diagnostics_vec`` will use to save the boundary conditions of the domain. Several files will also be generated (mitgrid, bathymetry), which will be required for STEP2. 

## II. Generate the ``.mitgrid`` file 

You will need to generate a file with all the grid information of your regional domain. The python code ``gen_mitgrid.py`` in the ``utils`` folder will generate this file. Below are the instructions to run the code in a terminal:

```
conda activate downscaling
cd ecco_darwin/regions/downscaling/utils/
python3 gen_mitgrid.py -d /path/to/save/the/grid -n name_of_the_region\ 
                       -c left_lon right_lon left_lat right_lat -r dlon dlat\
                       -e espg -v
```  
**Example:**
```
python3 gen_mitgrid.py -d ecco_darwin/regions/downscaling/ -n Mackenzie_Delta\ 
                       -c -153.1435 -118.0147 68.4727 75.2386 -r 0.0287 0.0091\
                       -e 4326 -v
```
To get more information about the options required for this code run ``gen_mitgrid.py -h``. Here are additional details about the option:
> - -d: The directory where to store the mitgrid
> - -n: The name of your region. It will be the name of the mitgrid
> - -c: The coordinates of the corners of your region. It can be in any ESPG.
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
                     -n Mackenzie_Delta -s 1224 744 -cs 0.3 1.0 -cw 1.0 1.14\
                     -wp 372 612 -v
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

## III. Generate an netcdf grid file with the model grid fields 

On this step we generate a netcdf file with the model grid fields as they will be interpreted by the model. **Note:** For this step you will need to run the ECCO-Darwin model for on time step. 

### a. Set the ecco-darwin v5 configuration to generate the netcdf

You will need here to start setting up the configuration you planed for your downscalled model, including the number of processors you to run the simulation on and how they will be split. 

```
cd downscalling
cp -r ecco_darwin/regions/downscaling/gen_ncgrid .
cd gen_ncgrid/code/
vim SIZE.h
##### Modify sNx to Nr parameters according to your downscalled setup #####
cd ../namelist/
vim data
##### Modify parameters according to your downscalled setup (see instructions below) #####
```

In the data file change the parameters according to the downsacled model configuration you want to create. Below are the details of importante parameters to set:
> - tRef: reference vertical profile for potential temperature (size = Nr the number of vertical levels)
> - sRef: reference vertical profile for salinity/specific humidity (size = Nr the number of vertical levels)
> - hFacMinDr: Minimum dimension size of a cell (same as bathy file option -cs)
> - hFacMin: Minimum fraction size of a cell (same as bathy file option -cs)
> - delR: Vertical grid spacing (size = Nr the number of vertical levels)
> - bathyfile: name of the bathymetry file generated previously

### b. Split mitgrid on the processors

Then, we need to split the mitgrid onto the processors depending on sNx and sNy choices. The python code ``mitgrid2tiles.py`` in the ``utils`` folder will generate these files in a ``mitgrids`` folder. Below are the instructions to run the code in a terminal:

```
conda activate downscaling
cd ecco_darwin/regions/downscaling/utils/
python3 mitgrid2tiles.py -d /path/to/save/the/grid -n name_of_the_region\
                         -s n_rows n_cols -p sNx sNy
```
**Example:**
```
 python3 mitgrid2tiles.py -d ecco_darwin/regions/downscaling/ -n SBS \
                          -s 1224 744 -p 34 24
```

### c. Run the model on 1 step


### d. Stictch the netcdf file

















