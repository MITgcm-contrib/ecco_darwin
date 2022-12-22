GlobalNEWS biogeochemical river exports (DIN, DON, DIP, DOP, DOC, PN, PP, POC, DSi)

GlobalNEWS2__RH2000Dataset-version1.0.xls
was obtained by Tom Van der Stocken from Emilio Mayorga
(mayorga@marine.rutgers.edu) on February 11, 2019.

The GlobalNEWS analysis is described here:
https://marine.rutgers.edu/globalnews/datasets.htm
Mayorga et al. (2010): https://www.sciencedirect.com/science/article/pii/S1364815210000186
Beusen et al. (2009): https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2008GB003281
Seitzinger et al. (2010): https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2009GB003587

DIC runoff is based on forcings from GlobalNEWS (dominant lithology, natural runoff and discharge)
and is computed from the consumption of CO2 by rock weathering with the formula of Suchet et al 2003
that is used as a predictor of DIC fluxes along with discharge in the formula of Li et al 2017
Suchet et al 2003: https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2002GB001891
Li et al 2017: https://www.sciencedirect.com/science/article/pii/S1470160X17302352
Details are in BDIC_GlobalNEWS.xlsx and final DIC fluxes in DIC_final_globalnews.xlsx

DIN is distributed in NO3, NO2 and NH4 according to a ratio calculated on the Glorich database. 
DIC to Alkalinity ratio is also computed from this database.
https://www.geo.uni-hamburg.de/en/geologie/forschung/aquatische-geochemie/glorich.html

Fe runoff is based on a constant P : Fe ratio from Lacroix et al. 2020:
https://doi.org/10.5194/bg-17-55-2020

Files:

mk_jra55_2000.m
Compute jra55_do time-mean year-2000 runoff

mk_SnapGlobalNEWS_indices.m 
Compute indices for snapping river point sources from GlobalNEWS to JRA-55 point sources

Snap_Examples.m
Example of using GlobalNews_to_JRA55.mat to add GlobalNEWS2 to JRA55 locations

SnapGlobalNEWS.m
Snap river point sources from GlobalNEWS to JRA-55 point sources

load_jra55_do_LLC_270_Nutrients.m
Snap-sum Global NEWS nutrient runoff from JRA-55 to LLC_270
Initially copied on April 2, 2022 from
/nobackup/dcarrol2/LOAC/m_files/jra55_do/LLC_270
outputs directory: /nobackup/rsavelli/LOAC/write_bin/jra55_do/v1.4.0/LLC_270/
