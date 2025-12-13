import math

# GEOMETRICAL AND PHYSICAL PARAMETERS
EL = 855000  # Estuarine length [m]
DEPTH_lb = 15.0  # Depth at downstream boundary [m]
DEPTH_ub = 30.0  # Depth at upstream boundary [m]
B_lb = 12000.0  # Width at downstream boundary [m]
B_ub = 2500.0  # Width at upstream boundary [m]
RS = 1.0  # Storage width ratio
rho_w = 1000.0  # Density of pure water [kg/m^3]
G = 9.81  # Gravity acceleration [m/s^2]
distance = 57  # Grid points in saline zone 
# The parameter distance delineates the number of
# grid point in the saline area. It is
# used as a fixed threshold to determine what
# forcings to use: upstream or downstream.

# HYDRODYNAMIC AND SEDIMENT PARAMETERS
Chezy_lb = 60  # Chezy coefficient (downstream) [m^-1/2 s^-1] 
Chezy_ub = 15  # Chezy coefficient (upstream) [m^-1/2 s^-1]  
ws = 1.0e-4  # Settling velocity for sediment [m/s]  #
tau_ero_lb = 0.6  # Erosion shear stress (downstream) [N/m^2] 
tau_dep_lb = 0.6  # Deposition shear stress (downstream) [N/m^2] 
tau_ero_ub = 0.25  # Erosion shear stress (upstream) [N/m^2] 
tau_dep_ub = 0.25  # Deposition shear stress (upstream) [N/m^2] 
Mero_lb = 1e-6  # Erosion coefficient (downstream) [mg/m^2 s] 
Mero_ub = 6.0e-6  # Erosion coefficient (upstream) [mg/m^2 s]

# BIOGEOCHEMICAL PARAMETERS
Pbmax = 5.58e-5  # Max photosynthetic rate [s^-1] 
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
kmort = 37e-8  # Phytoplankton mortality [s^-1] 
kexcr = 5e-2  # Excretion constant [-]
kgrowth = 0.3  # Growth constant [-] 
KD1 = 1.3  # Background extinction coefficient [m^-1]
KD2 = 0.06  # Suspended matter attenuation [mg^-1 m^-1]
kox = 1.44e-4  # Aerobic degradation rate [muM C s^-1] 
kdenit = 5.0e-4  # Denitrification rate [muM C s^-1] 
knit = 4.62e-4  # Nitrification rate [muM C s^-1] 

# --- physics/constants and forcings used by CO2 chemistry ---
USE_CO2_FLUX = True          # ⇦ turn ON to enable air–sea CO2 exchange
MOLAR_MASS_B = 10.811          # g mol^-1 (atomic B)
PH_ITERS = 50                   # pH fixed-point iterations (fast + stable)
CO2_PISTON_FROM_O2 = 0.915
pCO2 = 400 * 1e-6  # CO2 partial pressure in the atmosphere [atm] 

# EXTERNAL FORCINGS
Qr = -203782  # River discharge [m^3/s] 
AMPL = 5  # Tidal amplitude at mouth [m] 
pfun = 0.083  # Tidal frequency [cycle/hr] 
Uw_sal = 1  # Wind speed in saline estuary [m/s] 
Uw_tid = 1.78  # Wind speed in tidal river [m/s] 
water_temp = 30  # Water temperature [C] 
I = 367.6 # PAR [muE m^-2 s^-1] (0.4 x total downward shortwave radiation (200W.m-2)

# OTHER PARAMETERS
Euler = 0.5772156649  # Euler's constant
PI = math.pi  # Pi value

# NUMERICAL INTEGRATION
MAXT = (3*365) * 24 * 60 * 60  # Max time [s]
WARMUP = (2*365) * 24 * 60 * 60  # Warmup period [s]
DELTI = 150  # Delta t [s]
TS = 12  # Save every TS timesteps
DELXI = 5000  # Delta x [m]
# Convergence criterion
TOL = 1e-2  # Convergence criterion
ITE = 1000
TH_ABS_FLOOR = TOL   # m
TU_ABS_FLOOR = TOL   # m/s
TH_REL       = TOL   # relative per |TH|~1 m → 1e-10
TU_REL       = TOL   # relative per |TU|~1 m/s → 1e-10
LOOSE_CAP    = TOL*10    # allow up to 10× looser when CFL/shallow
M = int(EL / DELXI) + 1  # Max even grid points
if M % 2 == 0:
    pass
else:
    M = M - 1
M1 = M - 1  # Max odd grid points
M2 = M - 2  # Last even grid point
M3 = M - 3  # Last odd grid point
MAXV = 12  # Max species in chemical array
