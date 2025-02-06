import math

# GEOMETRICAL AND PHYSICAL PARAMETERS
EL = 158000  # Estuarine length [m]
DEPTH_lb = 11.5  # Depth at downstream boundary [m]
DEPTH_ub = 1.0  # Depth at upstream boundary [m]
B_lb = 6952.0  # Width at downstream boundary [m]
B_ub = 30.0  # Width at upstream boundary [m]
RS = 1.0  # Storage width ratio
rho_w = 1000.0  # Density of pure water [kg/m^3]
G = 9.81  # Gravity acceleration [m/s^2]
distance = 50  # Grid points in saline zone

# HYDRODYNAMIC AND SEDIMENT PARAMETERS
Chezy_lb = 70  # Chezy coefficient (downstream) [m^-1/2 s^-1]
Chezy_ub = 40  # Chezy coefficient (upstream) [m^-1/2 s^-1]
ws = 1.0e-3  # Settling velocity for sediment [m/s]
tau_ero_lb = 0.4  # Erosion shear stress (downstream) [N/m^2]
tau_dep_lb = 0.4  # Deposition shear stress (downstream) [N/m^2]
tau_ero_ub = 1.0  # Erosion shear stress (upstream) [N/m^2]
tau_dep_ub = 1.0  # Deposition shear stress (upstream) [N/m^2]
Mero_lb = 3.5e-6  # Erosion coefficient (downstream) [mg/m^2 s]
Mero_ub = 6.0e-8  # Erosion coefficient (upstream) [mg/m^2 s]

# BIOGEOCHEMICAL PARAMETERS
Pbmax = 1.157e-4  # Max photosynthetic rate [s^-1]
alpha = 5.8e-7  # Photosynthetic efficiency [m^2 s/muE]
KdSi = 20.0  # Silica constant [muM Si]
KPO4 = 0.5  # Phosphate constant [muM P]
KNH4 = 100.0  # Ammonium constant [muM N]
KNO3 = 45.0  # Nitrate constant [muM N]
KTOC = 60.0  # Organic matter constant [muM C]
KO2 = 15.0  # Oxygen constant [muM O2]
KN = 5.0  # Dissolved nitrogen constant [muM N]
KinO2 = 50.0  # Denitrification inhibition term [muM O2]
redsi = 16.0 / 80.0  # Redfield ratio for silica
redn = 16.0 / 106.0  # Redfield ratio for nitrogen
redp = 1.0 / 106.0  # Redfield ratio for phosphorus
kmaint = 9.2593e-7  # Maintenance rate constant [s^-1]
kmort = 7.1e-7  # Phytoplankton mortality [s^-1]
kexcr = 0.03  # Excretion constant [-]
kgrowth = 0.3  # Growth constant [-]
KD1 = 1.3  # Background extinction coefficient [m^-1]
KD2 = 0.06  # Suspended matter attenuation [mg^-1 m^-1]
kox = 2.0e-4  # Aerobic degradation rate [muM C s^-1]
kdenit = 1.0e-4  # Denitrification rate [muM C s^-1]
knit = 1.5e-4  # Nitrification rate [muM C s^-1]

# EXTERNAL FORCINGS
Qr = -32.0  # River discharge [m^3/s]
AMPL = 3.7  # Tidal amplitude at mouth [m]
pfun = 0.080536912751677847  # Tidal frequency [cycle/hr]
Uw_sal = 4.0  # Wind speed in saline estuary [m/s]
Uw_tid = 1.0  # Wind speed in tidal river [m/s]
water_temp = 17  # Water temperature [C]

# OTHER PARAMETERS
Euler = 0.5772156649  # Euler's constant
PI = math.pi  # Pi value

# NUMERICAL INTEGRATION
MAXT = 90 * 24 * 60 * 60  # Max time [s]
WARMUP = 3 * 24 * 60 * 60  # Warmup period [s]
DELTI = 150  # Delta t [s]
TS = 12  # Save every TS timesteps
DELXI = 2000  # Delta x [m]
TOL = 1e-10  # Convergence criterion
M = int(EL / DELXI) + 1  # Max even grid points
M1 = M - 1  # Max odd grid points
M2 = M - 2  # Last even grid point
M3 = M - 3  # Last odd grid point
MAXV = 10  # Max species in chemical array