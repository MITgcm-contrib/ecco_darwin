This folder contains MATLAB code written by Dustin Carroll for accurately re-gridding bathymetric products onto various LLC and lon-lat grids

Note: you will need to adjust the file paths at the top of the various MATLAB scripts for your local machine

MATLAB packages required: cmocean for plotting and gcmfaces + MATLAB Image Processing Toolbox for step 5 below

cmocean: https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
gcmfaces: https://github.com/MITgcm/gcmfaces

Example bathymetry/ice base products used are:
GEBCO: https://www.gebco.net/
Schaffer et al. (2016): https://doi.pangaea.de/10.1594/PANGAEA.856844?format=html#download
Schaffer et al. (2019): https://doi.pangaea.de/10.1594/PANGAEA.905295

BedMachine Greenland: https://sites.ps.uci.edu/morlighem/dataproducts/bedmachine-greenland/
BedMachine Antarctica: https://sites.ps.uci.edu/morlighem/dataproducts/bedmachine-antarctica/

Processing workflow:
1. Generate model cell corners for your grid of choice under cell_corners/
2. Convert bathymetric/ice base products from NetCDF to MATLAB files under raw_bathy/
3. Generate product->model grid indices under indices/
4. Reconstruct global bathymetry/ice base from indices under reconstruct_bathy/ or reconstruct_ice/
5. If needed, merge bathymetry/ice base under postprocess/
