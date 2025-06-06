Python Idealized Estuary C-GEM model for the North Slope (Alaska) with 
carbonate chemistry and air-sea CO2 flux in ./code_python_FCO2

This version reads external forcings as 365-day climatology or continuous 
with repeatYear parameter in config.py. It also accounts for ice cover where
biogeochenical reactions are stopped if cumulated water temperature over 
user-specified window (nbday_ice) is negative.

Original C version from Volta et al., 2014 www.geosci-model-dev.net/7/1271/2014/
Converted to Python in ./code_python

FCO2 computed from:
Laruelle et al 2017 and Regnier et al 2002 for parametrizing the air-sea fluxes
Regnier et al 2013 to parametrize DIC/ALK
Follows et al 2006 for computing alkalinity and aqueous CO2

Cite following paper for base version:
Volta et al., 2014 www.geosci-model-dev.net/7/1271/2014/

Cite following papers for FCO2 version:
Volta et al., 2014 - www.geosci-model-dev.net/7/1271/2014/
Volta et al., 2016 - www.hydrol-earth-syst-sci.net/20/991/2016/
Regnier et al., 2002 - https://doi.org/10.1016/S0307-904X(02)00047-1
Follows et al., 2006 - https://doi.org/10.1016/j.ocemod.2005.05.004
Regnier et al., 2013 - https://link.springer.com/article/10.1007/s10498-013-9218-3
Laruelle et al., 2017 - www.biogeosciences.net/14/2441/2017/

See full readme file for more details on the model:
https://github.com/MITgcm-contrib/ecco_darwin/blob/master/code_util/LOAC/C_GEM/readme.txt

