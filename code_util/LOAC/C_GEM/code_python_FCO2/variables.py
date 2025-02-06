"""
Global variables and structures (translated from variables.h)
"""

from config import M, MAXV

# Hydrodynamic Variables
Y = [0.0] * (M + 1)           # Convergence test
E = [0.0] * (M + 1)           # Convergence test
ZZ = [0.0] * (M + 1)          # Cross-section at reference level [m^2]
D = [0.0] * (M + 1)           # Total cross-section [m^2]
Dold = [0.0] * (M + 1)        # Previous total cross-section [m^2]
Dold2 = [0.0] * (M + 1)       # Additional previous cross-section [m^2]
DEPTH = [0.0] * (M + 1)       # Water depth [m]
Chezy = [0.0] * (M + 1)       # Chezy coefficient [m s-2]
H = [0.0] * (M + 1)           # Free cross-section [m^2]
U = [0.0] * (M + 1)           # Flow velocity [m/s]
B = [0.0] * (M + 1)           # Width [m]
TH = [0.0] * (M + 1)          # Temporary free cross-section [m]
TU = [0.0] * (M + 1)          # Temporary velocity [m/s]
LC = 0.0                      # Width convergence length [m]

# Transport Variables
C = [[0.0] * 5 for _ in range(M + 1)]  # Coefficients for tridiagonal matrix
Z = [0.0] * (M + 1)                    # Tridiagonal matrix coefficients
fl = [0.0] * (M + 1)                   # Advective flux [mmol/s]
aold = [0.0] * (M + 1)                 # Dispersion scheme variable
cold = [0.0] * (M + 1)                 # Dispersion scheme variable
rold = [0.0] * (M + 1)                 # Dispersion scheme variable
ccold = [0.0] * (M + 1)                # Dispersion scheme variable
dispersion = [0.0] * (M + 1)           # Dispersion coefficient [m^2/s]

# Piston Velocity Variables
kwind = [0.0] * (M + 1)                # Wind component for piston velocity [m/s]
kflow = [0.0] * (M + 1)                # Current component for piston velocity [m/s]
vp = [0.0] * (M + 1)                   # Piston velocity [m/s]

# Sediment Variables
Mero = [0.0] * (M + 1)                 # Erosion coefficient [mg/m^2/s]
tau_ero = [0.0] * (M + 1)              # Critical shear stress for erosion [N/m^2]
tau_dep = [0.0] * (M + 1)              # Critical shear stress for deposition [N/m^2]
tau_b = [0.0] * (M + 1)                # Bottom shear stress [N/m^2]
erosion = [0.0] * (M + 1)              # Sediment erosion rate [mg/m^2/s]
deposition = [0.0] * (M + 1)           # Sediment deposition rate [mg/m^2/s]

# Biogeochemical Variables
NPP_NO3 = [0.0] * (M + 1)              # NPP rate using NO3 [mmol C/m^3/s]
NPP_NH4 = [0.0] * (M + 1)              # NPP rate using NH4 [mmol C/m^3/s]
GPP = [0.0] * (M + 1)                  # Gross Primary Production rate [mmol C/m^3/s]
phy_death = [0.0] * (M + 1)            # Phytoplankton death rate [mmol C/m^3/s]
NPP = [0.0] * (M + 1)                  # Net Primary Production rate [mmol C/m^3/s]
Si_consumption = [0.0] * (M + 1)       # Silica consumption rate [mmol Si/m^3/s]
aer_deg = [0.0] * (M + 1)              # Aerobic degradation rate [mmol C/m^3/s]
denit = [0.0] * (M + 1)                # Denitrification rate [mmol C/m^3/s]
nit = [0.0] * (M + 1)                  # Nitrification rate [mmol N/m^3/s]
O2_ex = [0.0] * (M + 1)                # O2 flux across air-water interface [mmol O2/m^3/s]
NEM = [0.0] * (M + 1)                  # Net Ecosystem Metabolism [mmol C/m^3/s]
FCO2 = [0.0] * (M + 1)                  # air-sea CO2 fluxes [mmol C/m^3/s] || CO2 outgassing (FCO2 < 0)
Hplus = [0.0] * (M + 1)                  # Hydrogen ions concentration [mol kg-1]
#mol kg-1
#mol kg-1 x 1000 = mmol kg-1
#mmol kg-1 = mmol L-1
#mmol L-1 / 1000 = mmol m-3
#1 mol kg-1 = 1 mmol m-3

# Chemical Species and SPM
names = ['DIA', 'dSi', 'NO3', 'NH4', 'PO4', 'O2', 'TOC', 'S', 'SPM', 'DIC', 'ALK', 'pH']
v = {name: {
        "name": name,
        "env": 0,
        "c": [0.0] * (M + 1),
        "clb": 0.0,
        "cub": 0.0,
        "avg": [0.0] * (M + 1),
        "concflux": [0.0] * (M + 1),
        "advflux": [0.0] * (M + 1),
        "disflux": [0.0] * (M + 1),
    } for name in names}

# Global Flag
include_constantDEPTH = 0  # Include constant depth formulation