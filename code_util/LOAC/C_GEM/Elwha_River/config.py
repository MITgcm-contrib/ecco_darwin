import math
from datetime import datetime


# GEOMETRICAL AND PHYSICAL PARAMETERS
EL = 72000  # Estuarine length [m]
DEPTH_lb = 30.0  # Depth at downstream boundary [m]
DEPTH_ub = 3.2  # Depth at upstream boundary [m]
B_lb = 30.0  # Width at downstream boundary [m]
B_ub = 3.2  # Width at upstream boundary [m]
RS = 1.0  # Storage width ratio
rho_w = 1000.0  # Density of pure water [kg/m^3]
G = 9.81  # Gravity acceleration [m/s^2]
distance = 5  # Grid points in saline zone  (for experimentation set up to be half the distance of river in km (36 km))

# HYDRODYNAMIC AND SEDIMENT PARAMETERS
Chezy_lb = 60  # Chezy coefficient (downstream) [m^-1/2 s^-1]
Chezy_ub = 40  # Chezy coefficient (upstream) [m^-1/2 s^-1]
ws = 1.0e-3  # Settling velocity for sediment [m/s]
tau_ero_lb = 0.4  # Erosion shear stress (downstream) [N/m^2]
tau_dep_lb = 0.4  # Deposition shear stress (downstream) [N/m^2]
tau_ero_ub = 1.0  # Erosion shear stress (upstream) [N/m^2]
tau_dep_ub = 1.0  # Deposition shear stress (upstream) [N/m^2]
Mero_lb = 3.5e-6  # Erosion coefficient (downstream) [mg/m^2 s]
Mero_ub = 6.0e-8  # Erosion coefficient (upstream) [mg/m^2 s]

# BIOGEOCHEMICAL PARAMETERS
Pbmax = 2.58e-5  # Max photosynthetic rate [s^-1]
alpha = 4.11e-7  # Photosynthetic efficiency [m^2 s/muE]
KdSi = 1.07  # Silica constant [muM Si]
KPO4 = 0.2  # Phosphate constant [muM P]
KNH4 = 228.9  # Ammonium constant [muM N]
KNO3 = 26.07  # Nitrate constant [muM N]
KTOC = 186.25  # Organic matter constant [muM C]
KO2_ox = 31.0  # Oxygen constant [muM O2]
KO2_nit = 51.25  # Oxygen constant [muM O2]
KN = 1.13  # Dissolved nitrogen constant [muM N]
KinO2 = 33.0  # Denitrification inhibition term [muM O2]
redsi = 16.0 / 80.0  # Redfield ratio for silica
redn = 16.0 / 106.0  # Redfield ratio for nitrogen
redp = 1.0 / 106.0  # Redfield ratio for phosphorus
kmaint = 4.6e-7  # Maintenance rate constant [s^-1]
kmort = 1.56e-6  # Phytoplankton mortality [s^-1]
kexcr = 5e-2  # Excretion constant [-]
kgrowth = 2.9e-1  # Growth constant [-]
KD1 = 1.3  # Background extinction coefficient [m^-1]
KD2 = 0.06  # Suspended matter attenuation [mg^-1 m^-1]
kox = 6.08e-4  # Aerobic degradation rate [muM C s^-1]
kdenit = 5.05e-4  # Denitrification rate [muM C s^-1]
knit = 2.73e-5  # Nitrification rate [muM C s^-1]

# EXTERNAL FORCINGS
Qr = -47.7  # River discharge [m^3/s]  ################# (negative == downstream) (input value as negated)
AMPL = 3.5  # Tidal amplitude at mouth [m]
pfun = 0.0830  # Tidal frequency [cycle/hr]
Uw_sal = 4.39  # Wind speed in saline estuary [m/s]
Uw_tid = 2.195  # Wind speed in tidal river [m/s]
water_temp = 10  # Water temperature [C]
#pCO2 = 380 * 1e-6  # CO2 partial pressure in the atmosphere [atm]  ################## (change to time series) (test first with sine
use_real_pCO2 = True # toggle wheather to use the real world pCO2 time series or dummy sine wave function 

# OTHER PARAMETERS
Euler = 0.5772156649  # Euler's constant
PI = math.pi  # Pi value
pH_ite = 50  # number of iterations to converge to pH following Follows et al., 2006
mass_mol_B = 10.8110  # molar mass of Boron g/mol

# NUMERICAL INTEGRATION
SIM_START_DATETIME = datetime(2011, 1, 1)  # start date of the pCO2 for sine wavesimulation
MAXT = (365*3) * 24 * 60 * 60  # Max time [s]  (roughly 94 million sec)
WARMUP = (365*2) * 24 * 60 * 60  # Warmup period [s]
DELTI = 150  # Delta t [s]
TS = 12  # Save every TS timesteps
DELXI = 1000  # Delta x [m]
TOL = 1e-10  # Convergence criterion
M = int(EL / DELXI) + 1  # Max even grid points
if M % 2 == 0:
    pass
else:
    M = M - 1
M1 = M - 1  # Max odd grid points
M2 = M - 2  # Last even grid point
M3 = M - 3  # Last odd grid point
MAXV = 12  # Max species in chemical array
