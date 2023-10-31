https://github.com/MITgcm-contrib/ecco_darwin/tree/master/regions/GoA

Gulf of Alaska cutout of the global ECCO-Darwin version 5
https://onlinelibrary.wiley.com/doi/10.1029/2019MS001888
https://onlinelibrary.wiley.com/doi/10.1029/2021GB007162
https://ecco-group.org/storymaps.htm?id=45
https://ecco-group.org/storymaps.htm?id=71
https://darwin3.readthedocs.io/en/latest/phys_pkgs/darwin.html

Extracted for Takamitsu Ito on October 31, 2023

Extracted region is available on pfe:~dmenemen/ecco_darwin/GoA
and on https://nasa-ext.box.com/s/7z3jsoymr39c7jov81cmzafqqfx6922j

============
Grid information is in ../GoA/grid

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
Extracted monthly-mean output for January 1992 to June 2023
ETAN  sea surface height (m)
SALT  salinity (g/kg)
THETA potential temperature (deg C)
U     zonal velocity, >0 from West to East (m/s) specified on Southwest C-grid U point
V     meridional velocity, >0 from South to North (m/s) specified on Southwest C-grid V point

DIC   concentration of dissolved inorganic carbon       (mmol C   m^-3)
NO3   concentration of nitrate                          (mmol N   m^-3)
NO2   concentration of nitrite                          (mmol N   m^-3)
NH4   concentration of ammonia                          (mmol N   m^-3)
PO4   concentration of phosphate                        (mmol P   m^-3)
FeT   concentration of total dissolved iron             (mmol Fe  m^-3)
SiO2  concentration of inorganic silica                 (mmol Si  m^-3)
DOC   concentration of dissolved organic carbon         (mmol C   m^-3)
DON   concentration of dissolved organic nitrogen       (mmol N   m^-3)
DOP   concentration of dissolved organic phosphorus     (mmol P   m^-3)
DOFe  concentration of dissolved organic iron           (mmol Fe  m^-3)
POC   concentration of particulate organic carbon       (mmol C   m^-3)
PON   concentration of particulate organic nitrogen     (mmol N   m^-3)
POP   concentration of particulate organic phosphorus   (mmol P   m^-3)
POFe  concentration of particulate organic iron         (mmol Fe  m^-3)
POSi  concentration of particulate organic silica       (mmol Si  m^-3)
PIC   concentration of particulate inorganic carbon     (mmol C   m^-3)
ALK   alkalinity                                        (meq      m^-3)
O2    concentration of oxygen                           (mmol O2  m^-3)
c1    concentration of carbon in plankton type 1        (mmol C   m^-3)
c2    concentration of carbon in plankton type 2        (mmol C   m^-3)
c3    concentration of carbon in plankton type 3        (mmol C   m^-3)
c4    concentration of carbon in plankton type 4        (mmol C   m^-3)
c5    concentration of carbon in plankton type 5        (mmol C   m^-3)
c6    concentration of carbon in plankton type 6        (mmol C   m^-3)
c7    concentration of carbon in plankton type 7        (mmol C   m^-3)
Chl1  concentration of Chlorophyll-a in plankton type 1 (mg Chl a m^-3)
Chl2  concentration of Chlorophyll-a in plankton type 2 (mg Chl a m^-3)
Chl3  concentration of Chlorophyll-a in plankton type 3 (mg Chl a m^-3)
Chl4  concentration of Chlorophyll-a in plankton type 4 (mg Chl a m^-3)
CHl5  concentration of Chlorophyll-a in plankton type 5 (mg Chl a m^-3)

apCO2    atmospheric pCO2               (atm)
pH       potential of hydrogen          (pH)
fugCO2   fugacity of CO2                (atm)
CO2_flux flux of CO2 - air-sea exchange (mmol C/m^2/s)

Additional model variables that can easily be diagnosed
are defined in the file available_diagnostics.log

============
Model output naming convention and format

All files are plain binary, real*4, IEEE big-endian,
with dimensions 47 (with file name termination _47),
144x72 (with file name termination _144x72), and
144x72x47 (with file name termination _144x72x47).

Monthly means are indicated with end of averaging period,
for example, ".19920201T000000" indicates a monthly mean
field for January 1992.
