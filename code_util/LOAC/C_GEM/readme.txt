Python Idealized Estuary C-GEM model
Original C version from Volta et al., 2014 www.geosci-model-dev.net/7/1271/2014/
Converted to Python in ./code_python

code_python_v2 includes an updated version of the code. It is faster and includes a built-in function for activating CO2 fluxes in config.py in the main code (USE_CO2_FLUX). This new version uses more vectorization, removes excess use of loops, and uses Numba when possible. It also relaxes the tolerance criterion for flux convergence with possibility to relax when CFL condition allows. Thresholds can be set in config.py.

Addition of carbonate chemistry and air-sea CO2 flux in ./code_python_FCO2
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


To set-up model, follow Volta et al., 2014:

config.py
- Specification of estuarine geometry (protocol step 1) and physical parameters, such as 
pure water density and gravity acceleration
  The parameter distance delineates the number of grid point in the saline area. It is
  used as a fixed threshold to determine what forcings to use: upstream or 
  downstream.
- Specification of hydrodynamic and sediment parameters (protocol steps 2.1 and 4.1)
- Specification of biogeochemical parameters (protocol step 5.2)
- Specification of external forcing (protocol step 5.4)
- Specification of mathematical constants (Euler’s constant, number Pi), setup of model
resolution (time, space, warm up period for biogeochemistry) and specification of
libraries

init_module.py
- Flag to use a variable or a constant estuarine depth. The user can decide to use a 
variable depth (include_constant_depth=0) as in the simulations presented in our study or a
constant depth (include_constant_depth=1).
If:
· include_constant_depth=0. The tidally-averaged water depth decreases linearly in the 
upstream direction and the Chézy coefficient, as well as the critical shear for erosion 
and deposition are constant until a set distance. Further upstream, they increase 
linearly towards the inland boundary. In this case, the user needs to specify in config.py 
the number of grid points (parameter named distance) corresponding to this distance.
· include_constant_depth=1. The tidally-averaged water depth is constant along the 
estuarine domain and two different constant values, one for the saline estuary and 
one for the tidal river, are assigned to the Chézy coefficient and the critical shear 
stress for erosion and deposition.
- Formulation of dispersion coefficient D (protocol step 3.1)
- Formulation of estuarine convergence length (protocol step 1)
- Initializing numerical arrays for hydrodynamics, transport and biogeochemical reaction
network
- Specification of boundary conditions for salinity, chemical species and SPM 
(protocol steps 3.2, 4.2, 5.3)
- Saving salinity and chemical species and SPM concentrations in outputs files. All 
variables considered in our manuscript are written in output files

variables.py
- List and description of all model variables

hyd_module.py
- Main hydrodynamic routine 

transport_module.py
- Main transport routine 

fun_module.py
- Formulation of piston velocity for O exchange across air-water interface 2
- Formulation of tidal elevation, absolute temperature and light intensity (protocol
 step 5.4)
- Formulation of temperature-dependence for heterotrophic reactions and nitrification 
- Compute total/carbonate alkalinity, coefficients of dissociation, pH for FCO2

schemes_module.py
- Formulations of advection and dispersion numerical schemes 

tridag_module.py
- Solution of the tridiagonal matrix 

uphyd_module.py
- Update hydrodynamic variables 

biogeo_module.py
- Formulation of underwater light field and nutrient dependence for primary production
- Formulation of biogeochemical reaction network (protocol step 5.1)
- Update of biogeochemical state variables. Note that phytoplankton concentration is
referred as diatom concentration (DIA)
- Saving biogeochemical process rates in outputs files. All biogeochemical processes 
considered in our manuscript are written in output files

sed_module.py
- Formulation of longitudinal variation in SPM parameters (protocol step 4.1)
- Formulation of SPM erosion and deposition rates
- Update of SPM concentrations
        
file_module.py
- Saving geometrical variables in output files
- Format specification for concentration and biogeochemical reaction rate output files

To run model, run main.py

Output Files
C-GEM output data are printed in *.dat type files. Output files, produced by default by 
C-GEM, with specification of the corresponding modeled data are:
- depth.dat: water depth data
- width.dat: estuarine width data
_ U.dat: horizontal velocity (-: flow towards the ocean; +: flow towards river source)
- S.dat: salinity data
- SPM.dat: suspended particulate matter concentrations
- DIA.dat: diatoms concentrations
- NH4.dat: ammonium concentrations
- NO3.dat: nitrate concentrations
- O2.dat: oxygen concentrations
- PO4.dat: phosphate concentrations
- DSi.dat: dissolved silica concentrations
- TOC.dat: total organic carbon concentrations
- NPP.dat: Net Primary Production rates
- aer_deg.dat: aerobic degradation rates
- denit.dat: denitrification rates
- nit.dat: nitrification rates
- O2_ex.dat: air-water O exchange rates 
- NEM.dat: Net Ecosystem Metabolism rates

if FCO2 version, additional outputs:
- DIC.dat: dissolved inorganic carbon concentration
- ALK.dat: alkalinity
- pH.dat: pH
- FCO2.dat: air-water CO2 exchange rates (negative raw output means CO2 outgassing to the atmosphere)

The first column of the output file corresponds to the time expressed in second and other
columns correspond to the selected data along the estuarine length from the mouth of the 
estuary to the
upstream cell. 
