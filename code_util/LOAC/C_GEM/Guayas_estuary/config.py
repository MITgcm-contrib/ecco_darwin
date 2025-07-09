import math

# GEOMETRICAL AND PHYSICAL PARAMETERS
EL = 74000  # Estuarine length [m] x Twilley et al
DEPTH_lb = 10.0  # Depth at downstream boundary [m] x Twilley et al
DEPTH_ub = 8.0  # Depth at upstream boundary [m] x Twilley et al
B_lb = 3000.0  # Width at downstream boundary [m] x riverlines database
B_ub = 600.0  # Width at upstream boundary [m] x riverlines database (74 km from puna island shore)
RS = 1.0  # Storage width ratio
rho_w = 1000.0  # Density of pure water [kg/m^3]
G = 9.81  # Gravity acceleration [m/s^2]
distance = 22  # Grid points in saline zone (45 km from BElliard et al between 60 and 30 km in dry and rainy season)
# The parameter distance delineates the number of 
# grid point in the saline area. It is
# used as a fixed threshold to determine what 
# forcings to use: upstream or downstream.

# HYDRODYNAMIC AND SEDIMENT PARAMETERS
Chezy_lb = 60  # Chezy coefficient (downstream) [m^-1/2 s^-1] # Nguyen et al
Chezy_ub = 15  # Chezy coefficient (upstream) [m^-1/2 s^-1]  # Nguyen et al
ws = 1.0e-4  # Settling velocity for sediment [m/s]  # Nguyen et al
tau_ero_lb = 0.6  # Erosion shear stress (downstream) [N/m^2] # Nguyen et al
tau_dep_lb = 0.6  # Deposition shear stress (downstream) [N/m^2] # Nguyen et al
tau_ero_ub = 0.25  # Erosion shear stress (upstream) [N/m^2] # Nguyen et al
tau_dep_ub = 0.25  # Deposition shear stress (upstream) [N/m^2] # Nguyen et al
Mero_lb = 1e-6  # Erosion coefficient (downstream) [mg/m^2 s] # Nguyen et al
Mero_ub = 6.0e-6  # Erosion coefficient (upstream) [mg/m^2 s] # Nguyen et al

# BIOGEOCHEMICAL PARAMETERS
Pbmax = 5.58e-5  # Max photosynthetic rate [s^-1] # Nguyen et al
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
kmort = 37e-8  # Phytoplankton mortality [s^-1] # Nguyen et al
kexcr = 5e-2  # Excretion constant [-]
kgrowth = 0.3  # Growth constant [-] # Nguyen et al
KD1 = 1.3  # Background extinction coefficient [m^-1]
KD2 = 0.06  # Suspended matter attenuation [mg^-1 m^-1]
kox = 1.44e-4  # Aerobic degradation rate [muM C s^-1] # Nguyen et al
kdenit = 5.0e-4  # Denitrification rate [muM C s^-1] # Nguyen et al
knit = 4.62e-4  # Nitrification rate [muM C s^-1] # Nguyen et al

# EXTERNAL FORCINGS 
Qr = -1143  # River discharge [m^3/s] x Twilley et al
AMPL = 1.8  # Tidal amplitude at mouth [m] x Twilley et al
pfun = 0.0805  # Tidal frequency [cycle/hr] x Twilley et al
Uw_sal = 1.5  # Wind speed in saline estuary [m/s] x Twilley et al
Uw_tid = 1.5  # Wind speed in tidal river [m/s] x Twilley et al
water_temp = 27  # Water temperature [C] x Twilley et al
pCO2 = 390 * 1e-6  # CO2 partial pressure in the atmosphere [atm] NOAA MBL 2000_2020 lat -3 to -2

# OTHER PARAMETERS
Euler = 0.5772156649  # Euler's constant
PI = math.pi  # Pi value
pH_ite = 50  # number of iterations to converge to pH following Follows et al., 2006
mass_mol_B = 10.8110  # molar mass of Boron g/mol

# NUMERICAL INTEGRATION
MAXT = (365*3) * 24 * 60 * 60  # Max time [s]
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
