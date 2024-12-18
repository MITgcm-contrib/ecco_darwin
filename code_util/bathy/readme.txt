This folder contains MATLAB code written by Dustin Carroll for accurately re-gridding bathymetric products onto various LLC and lon-lat grids

Note: you will need to adjust the file paths at the top of the various MATLAB files for your local machine

Processing workflow:
1. Generate model cell corners for your grid of choice under cell_corners/
2. Convert bathymetric/ice base products from NetCDF to MATLAB files under raw_bathy/
3. Generate product->model grid indices under indices/
4. Reconstruct global bathymetry/ice base from indices under reconstruct_bathy/ or reconstruct_ice/
5. If needed, merge bathymetry/ice base under postprocess/
