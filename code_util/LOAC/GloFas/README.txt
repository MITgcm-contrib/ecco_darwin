Create freshwater runoff forcing files from Global Flood Awareness System (GloFas), Copernicus, 3 min or 0.05x0.05 degree
https://global-flood.emergency.copernicus.eu

Global Flood Awareness System (GloFas) runoff outputs come in m3/s.

make_GLOFAS.py is a Python script that keeps only the coastal grid cells, converts runoff into m/s, and saves it into a binary file. Call it in the terminal with:
python3 make_GLOFAS.py glofas_input_NCDFfile True
Inputs:
  glofas_input_NCDFfile : input netcdf file. Must be one full year. m3/sec
  arg_plot : save map of annual discharge

Auxiliary files:
  GLOFAS_pixarea_Global_03min.nc: pixel area m2
  GLOFAS_fracwater_Global_03min.nc: water mask [0-1]

Outputs:
  Binary file GloFas_<year> compatible with load_glofas_ECCO_V4r5.m

load_glofas_ECCO_V4r5.m is a Matlab script that loads GloFas freshwater runoff to wet coastal grid cells of ECCO_V4r5 and writes binary files GloFas_ECCO_V4r5_<year> that can be used as forcings in data.exf.
