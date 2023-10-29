This file is available in


============
Mediteranean and Black Sea cutout of the global ECCO-Darwin version 5
https://onlinelibrary.wiley.com/doi/10.1029/2019MS001888
https://onlinelibrary.wiley.com/doi/10.1029/2021GB007162

Extracted for Louisa Giannoudi and Aleka Pavlidou, HCMR,
on October 25, 2023, using
https://github.com/MITgcm-contrib/ecco_darwin/tree/master/regions/Med

Extracted region is available on pfe:~dmenemen/ecco_darwin/Med
and on https://nasa-ext.box.com/s/7z3jsoymr39c7jov81cmzafqqfx6922j

============
Grid information is in ../Med/grid

RC_47           :: vertical coordinate of center of cell (m)
DRF_47          :: Cell face separation along Z axis (m)

XC_144x72       :: longitude East of center of grid cell
XG_144x72       :: longitude East of southwest corner of grid cell
YC_144x72       :: latitude North of center of grid cell
YG_144x72       :: latitude North of southwest corner of grid cell

DXC_144x72      :: Cell center separation in X across western cell wall (m)
DXG_144x72      :: Cell face separation in X along southern cell wall (m)
DYC_144x72      :: Cell center separation in Y across southern cell wall (m)
DYG_144x72      :: Cell face separation in Y along western cell wall (m)

Depth_144x72    :: Model bathymetry (m)

RAC_144x72      :: vertical face area of tracer cell (m^2)
RAZ_144x72      :: vertical face area of vorticity points (m^2)

hFacC_144x72x47 :: mask of tracer cell (0 is land, >0 is wet)
hFacS_144x72x47 :: mask of v cell (0 is land, >0 is wet)
hFacW_144x72x47 :: mask of u cell (0 is land, >0 is wet)

More details about MITgcm grid are available in:
https://github.com/MITgcm/MITgcm/blob/master/model/inc/GRID.h

============
Available monthly-mean output is:
Eta.*.data      sea surface height (m)
Salt.*.data     salinity (g/kg)
Theta.*.data    potential temperature (deg C)
U.*.data        zonal (relative to grid) velocity, >0 from West to East (m/s)
                specified on Southwest C-grid U point
V.*.data        merid. (relative to grid) velocity, >0 from South to North (m/s)
                specified on Southwest C-grid V point
