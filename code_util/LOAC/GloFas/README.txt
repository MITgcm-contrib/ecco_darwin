Create daily freshwater runoff forcing files for ECCO from Global Flood Awareness System (GloFas), Copernicus, 3 min or 0.05x0.05 degree, 
historical run available at: https://ewds.climate.copernicus.eu/datasets/cems-glofas-historical?tab=overview

Global Flood Awareness System (GloFas) runoff daily outputs come in m3/s. Runoff has to be divided by area (m2).

Auxiliary files are available on: https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/CEMS-GLOFAS/LISFLOOD_static_and_parameter_maps_for_GloFAS/
Pixel area and Water mask can be found at:
/nobackup/rsavelli/GloFas/GLOFAS_pixarea_Global_03min.nc    ::  pixel area m2
/nobackup/rsavelli/GloFas/GLOFAS_fracwater_Global_03min.nc  ::  water mask [0-1]

################## make_coastal_GLOFAS.py ##################

make_coastal_GLOFAS.py is a Python script that keeps only the coastal grid cells, converts daily runoff into m/s, and saves it into a binary 
file. The script also removes duplicated rivers that are accounted for several times along the coastline. Run the following command in the 
terminal:

python3 make_coastal_GLOFAS.py glofas_input_NCDFfile True

Inputs:
  glofas_input_NCDFfile : input netcdf file. m3/sec
  arg_plot : save map of annual discharge
Outputs:
  Binary file GloFas_<year> compatible with load_GLOFAS_to_ECCO.py

################## load_GLOFAS_to_ECCO.py ##################

load_GLOFAS_to_ECCO.py is a Python script that loads GloFas daily freshwater runoff from previously created binary files GloFas_<year> to wet 
coastal grid cells of ECCO and writes binary files GloFas_runoff_<year> that can be used as forcings in the EXF package. It computes matching 
indexes and area weights between GloFas and ECCO grids and aggregates runoff to ECCO grid cells. Set directories and ECCO grid dimensions in 
the script's header. Script saves grid mapping into a pickle file. For example, run the following command in the terminal:

python3 load_GLOFAS_to_ECCO.py --compute-grid true --coast-mask ECCO_V4r5_coastMask_orig.mat --mapping-file glofas_ECCO_V4r5_grid.pkl --out-prefix GloFAS_runoff_ECCO_V4r5_

Inputs:
  Automatically scans and iterates through GloFas_<year> files
  --compute-grid Whether to compute the GloFAS to ECCO mapping (true/false) if already exists (uses existing pickle mapping file)
  --coast-mask Name of ECCO coastal mask file (e.g. ECCO_V4r5_coastMask_orig.mat)
  --mapping-file Mapping pickle file name (e.g. glofas_ECCO_V4r5_grid_orig.pkl)
  --out-prefix Output binary file prefix (e.g. GloFAS_runoff_ECCO_V4r5_)
Outputs:
  Binary file GloFas_runoff_<grid>_<year> compatible with EXF package.

################## patch_jra55_do_antarctica_to_glofas.py ##################

patch_jra55_do_antarctica_to_glofas.py is a Python script that patches regridded JRA55-do freshwater runoff the ECCO grid into regridded GloFas 
on the same ECCO grid for latitudes south of 60°S. JRA55-do at these latitudes is a climatology repeating over the years. This step necessitates 
that you have already loaded JRA55-DO onto your ECCO grid (for example, see load_jra55_do_ECCO_V4r5_NAS.m at 
https://github.com/MITgcm-contrib/ecco_darwin/blob/master/code_util/LOAC/GlobalNews/load_jra55_do_ECCO_V4r5_NAS.m). Set directories and ECCO grid
dimensions in the script's header. Run the following command in the terminal:

python3 patch_jra55_do_antarctica_to_glofas.py --glofas_prefix GloFAS_runoff_ECCO_V4r5_

Inputs:
  Automatically scans and iterates through glofas_prefix_<year>
  --glofas_prefix Prefix of previously regridded GloFAS to ECCO runoff binary files
Outputs:
  Binary file glofas_prefix_<year> compatible with EXF package. Overwrites existing glofas_prefix_<year> files.

