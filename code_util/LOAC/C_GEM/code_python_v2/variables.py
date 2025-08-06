"""
Global variables and structures (translated from variables.h)
"""
import numpy as np
from config import M, MAXV

# Create a helper function for consistent zero-based indexing with ignored 0 index
def zeros_array():
    return np.zeros(M + 1, dtype=np.float64)

# Hydrodynamic Variables
Y = zeros_array()           # Convergence test
E = zeros_array()           # Convergence test
ZZ = zeros_array()          # Cross-section at reference level [m^2]
D = zeros_array()           # Total cross-section [m^2]
Dold = zeros_array()        # Previous total cross-section [m^2]
Dold2 = zeros_array()       # Additional previous cross-section [m^2]
DEPTH = zeros_array()       # Water depth [m]
Chezy = zeros_array()       # Chezy coefficient [m s-2]
H = zeros_array()           # Free cross-section [m^2]
U = zeros_array()           # Flow velocity [m/s]
B = zeros_array()           # Width [m]
TH = zeros_array()          # Temporary free cross-section [m]
TU = zeros_array()          # Temporary velocity [m/s]
LC = 0.0                   # Width convergence length [m]

# Transport Variables
# C was 2D list [M+1][5], now NumPy 2D array
C = np.zeros((M + 1, 5), dtype=np.float64)  # Coefficients for tridiagonal matrix
Z = zeros_array()                           # Tridiagonal matrix coefficients
fl = zeros_array()                          # Advective flux [mmol/s]
aold = zeros_array()                        # Dispersion scheme variable
cold = zeros_array()                        # Dispersion scheme variable
rold = zeros_array()                        # Dispersion scheme variable
ccold = zeros_array()                       # Dispersion scheme variable
dispersion = zeros_array()                  # Dispersion coefficient [m^2/s]

# Piston Velocity Variables
kwind = zeros_array()                       # Wind component for piston velocity [m/s]
kflow = zeros_array()                       # Current component for piston velocity [m/s]
vp = zeros_array()                          # Piston velocity [m/s]

# Sediment Variables
Mero = zeros_array()                        # Erosion coefficient [mg/m^2/s]
tau_ero = zeros_array()                     # Critical shear stress for erosion [N/m^2]
tau_dep = zeros_array()                     # Critical shear stress for deposition [N/m^2]
tau_b = zeros_array()                       # Bottom shear stress [N/m^2]
erosion = zeros_array()                     # Sediment erosion rate [mg/m^2/s]
deposition = zeros_array()                  # Sediment deposition rate [mg/m^2/s]

# Biogeochemical Variables
NPP_NO3 = zeros_array()                     # NPP rate using NO3 [mmol C/m^3/s]
NPP_NH4 = zeros_array()                     # NPP rate using NH4 [mmol C/m^3/s]
GPP = zeros_array()                         # Gross Primary Production rate [mmol C/m^3/s]
phy_death = zeros_array()                   # Phytoplankton death rate [mmol C/m^3/s]
NPP = zeros_array()                         # Net Primary Production rate [mmol C/m^3/s]
Si_consumption = zeros_array()              # Silica consumption rate [mmol Si/m^3/s]
aer_deg = zeros_array()                     # Aerobic degradation rate [mmol C/m^3/s]
denit = zeros_array()                       # Denitrification rate [mmol C/m^3/s]
nit = zeros_array()                         # Nitrification rate [mmol N/m^3/s]
O2_ex = zeros_array()                       # O2 flux across air-water interface [mmol O2/m^3/s]
NEM = zeros_array()                         # Net Ecosystem Metabolism [mmol C/m^3/s]


# Chemical Species and SPM
names = ['DIA', 'dSi', 'NO3', 'NH4', 'PO4', 'O2', 'TOC', 'S', 'SPM']
v = {
    name: {
        "name": name,
        "env": 0,
        "c": np.zeros(M + 1, dtype=float),
        "clb": 0.0,
        "cub": 0.0,
        "avg": np.zeros(M + 1, dtype=float),
        "concflux": np.zeros(M + 1, dtype=float),
        "advflux": np.zeros(M + 1, dtype=float),
        "disflux": np.zeros(M + 1, dtype=float),
    }
    for name in names
}

# Global Flag
include_constantDEPTH = 0  # Include constant depth formulation