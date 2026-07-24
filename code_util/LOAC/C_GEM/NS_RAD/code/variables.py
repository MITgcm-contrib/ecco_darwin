"""
Global variables and structures (translated from variables.h)
"""

import numpy as np

from config import M, MAXV

# Hydrodynamic Variables
# ----------------------
# These are numpy float64 arrays rather than Python lists so the hot solver kernels
# in tridag_module / uphyd_module can be @njit-compiled -- the hydrodynamic loop is
# ~60% of runtime and runs ~44 iterations per timestep under the shallow
# observation-based geometry. Item assignment (`arr[i] = x`) behaves identically to
# a list, so the surrounding code is unchanged; only the mutate-in-place contract
# matters, and that is preserved. Do NOT rebind these to lists.
_z = lambda: np.zeros(M + 1, dtype=np.float64)

Y = _z()                      # Convergence test
E = _z()                      # Convergence test
ZZ = _z()                     # Cross-section at reference level [m^2]
D = _z()                      # Total cross-section [m^2]
Dold = _z()                   # Previous total cross-section [m^2]
Dold2 = _z()                  # Additional previous cross-section [m^2]
DEPTH = _z()                  # Water depth [m]
Chezy = _z()                  # Chezy coefficient [m^-1/2 s^-1] (units per config.py)
H = _z()                      # Free cross-section [m^2]
U = _z()                      # Flow velocity [m/s]
B = _z()                      # Width [m] -- TOTAL conveyance/surface width, summed over
                              #   parallel threads where the channel is braided or deltaic
B_thread = _z()               # PER-THREAD width [m] = B / n_chan. Equals B for a
                              #   single-thread channel (and whenever MULTICHANNEL is off).
                              #   Used only by the Seo & Cheong shear-dispersion closure,
                              #   which is a WITHIN-channel process -- see config.MULTICHANNEL.
n_chan = _z()                 # Parallel thread count at each grid point [-] (>= 1)
TH = _z()                     # Temporary free cross-section [m]
TU = _z()                     # Temporary velocity [m/s]
LC = 0.0                      # Width convergence length [m]

# Transport Variables
C = np.zeros((M + 1, 5), dtype=np.float64)  # Coefficients for tridiagonal matrix
Z = _z()                                    # Tridiagonal matrix coefficients
fl = _z()                   # Advective flux [mmol/s]
aold = _z()                 # Dispersion scheme variable
cold = _z()                 # Dispersion scheme variable
rold = _z()                 # Dispersion scheme variable
ccold = _z()                # Dispersion scheme variable
dispersion = _z()           # Dispersion coefficient [m^2/s]

# Piston Velocity Variables
kwind = _z()                # Wind component for piston velocity [m/s]
kflow = _z()                # Current component for piston velocity [m/s]
vp = _z()                   # Piston velocity [m/s]

# Ice Variables (item 2 -- river ice model)
ice_thickness = _z()        # ice thickness [m], prognostic per cell
ice_frac = _z()             # areal ice cover fraction [0-1], diagnosed from thickness

# Sediment Variables
Mero = _z()                 # Erosion coefficient [mg/m^2/s]
tau_ero = _z()              # Critical shear stress for erosion [N/m^2]
tau_dep = _z()              # Critical shear stress for deposition [N/m^2]
tau_b = _z()                # Bottom shear stress [N/m^2]
erosion = _z()              # Sediment erosion rate [mg/m^2/s]
deposition = _z()           # Sediment deposition rate [mg/m^2/s]

# Biogeochemical Variables
NPP_NO3 = _z()              # NPP rate using NO3 [mmol C/m^3/s]
NPP_NH4 = _z()              # NPP rate using NH4 [mmol C/m^3/s]
GPP = _z()                  # Gross Primary Production rate [mmol C/m^3/s]
phy_death = _z()            # Phytoplankton death rate [mmol C/m^3/s]
NPP = _z()                  # Net Primary Production rate [mmol C/m^3/s]
Si_consumption = _z()       # Silica consumption rate [mmol Si/m^3/s]
aer_deg = _z()              # Aerobic degradation rate [mmol C/m^3/s]
denit = _z()                # Denitrification rate [mmol C/m^3/s]
nit = _z()                  # Nitrification rate [mmol N/m^3/s]
O2_ex = _z()                # O2 flux across air-water interface [mmol O2/m^3/s]
NEM = _z()                  # Net Ecosystem Metabolism [mmol C/m^3/s]
FCO2 = _z()                  # air-sea CO2 fluxes [mmol C/m^3/s] || CO2 outgassing (FCO2 > 0)
Hplus = _z()                  # Hydrogen ions concentration [mol kg-1]

# Arctic biogeochemistry extension (docs/arctic_biogeochemistry.md) -- process-rate
# diagnostics for the added terms, same convention as the rates above [mmol m^-3 s^-1].
rdoc_ox = _z()              # refractory-DOC aerobic oxidation
photo = _z()                # CDOM photomineralisation (RDOC -> DIC + labile DOC)
ch4_ox = _z()               # methanotrophy (CH4 -> DIC, consumes O2)
ch4_ex = _z()               # air-water CH4 flux  (>0 outgassing)
n2o_prod = _z()             # N2O production from nitrification + denitrification
n2o_ex = _z()               # air-water N2O flux  (>0 outgassing)
sod = _z()                  # sediment oxygen demand -> benthic DIC efflux
# Distributed lateral inflow per cell [m^3 s^-1], built in init from config.LATERAL_INFLOW
# (zero unless the active site enables it). Injected as a mixing source in main.py.
q_lat = _z()
#mol kg-1
#mol kg-1 x 1000 = mmol kg-1
#mmol kg-1 = mmol L-1
#mmol L-1 / 1000 = mmol m-3
#1 mol kg-1 = 1 mmol m-3

# Chemical Species and SPM
# 'T' (water temperature, degC) is a TRANSPORTED scalar, advected and dispersed by
# exactly the same TVD + dispersion machinery as salinity (env=1). It replaces the
# single scalar `water_temp` that used to be applied uniformly to all M grid points.
# Its boundary values are time-varying and are refreshed every timestep in main.py:
# the downstream end takes SEA temperature (PRDA2), the upstream end river
# temperature -- so watertemp.csv, which was wrong as a river forcing, is correct
# here as the marine boundary condition.
#
# T is transported like salinity -- step (a) of docs/ice_model_plan.md 4b. Steps (b) and
# (c) are ALSO implemented now and run every timestep in main.py: the surface heat budget
# (heat_module.py) heats/cools the channel interior, and the prognostic ice model
# (ice_module.py) grows/melts a cover from it. So T is no longer a pure mixing tracer
# between the two boundary end-members -- the interior is atmospherically forced.
# TOC is the LABILE / bioavailable DOC pool (its original role). RDOC is the added
# refractory + chromophoric pool (slow oxidation + photomineralisation); CH4 and N2O
# are the added dissolved gases. All are env=1 (transported) like the rest.
names = ['DIA', 'dSi', 'NO3', 'NH4', 'PO4', 'O2', 'TOC', 'RDOC', 'CH4', 'N2O',
         'S', 'SPM', 'DIC', 'ALK', 'pH', 'T']

# The per-species arrays are numpy for the same reason as the hydrodynamic state
# above: `biogeo` indexes them ~60-70 times per grid point per timestep, and its loop
# is the largest remaining hotspot. As with the hydrodynamic arrays, NEVER rebind
# these -- assign in place (`arr[:] = 0.0`), or every module holding a reference
# silently detaches. transport_module used to rebind `c` to a fresh list when the
# ice gate closed; that is now an in-place fill.
v = {name: {
        "name": name,
        "env": 0,
        "c": np.zeros(M + 1, dtype=np.float64),
        "clb": 0.0,
        "cub": 0.0,
        "avg": np.zeros(M + 1, dtype=np.float64),
        "concflux": np.zeros(M + 1, dtype=np.float64),
        "advflux": np.zeros(M + 1, dtype=np.float64),
        "disflux": np.zeros(M + 1, dtype=np.float64),
    } for name in names}

# Global Flag
include_constantDEPTH = 0  # Include constant depth formulation

