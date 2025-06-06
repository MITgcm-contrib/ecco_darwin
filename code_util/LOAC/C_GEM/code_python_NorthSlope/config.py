import math

# GEOMETRICAL AND PHYSICAL PARAMETERS
EL = 27175  # Estuarine length [m] - delta segment data from SWORD
DEPTH_lb = 15  # Depth at downstream boundary [m] - Sagavanirktok R NR Pump Sta 3 AK - 15908000 gauge data
DEPTH_ub = 15  # Depth at upstream boundary [m] - Sagavanirktok R NR Pump Sta 3 AK - 15908000 gauge data
B_lb = 1215  # Width at downstream boundary [m] - delta segment data from SWORD
B_ub = 859  # Width at upstream boundary [m] - delta segment data from SWORD
RS = 1.0  # Storage width ratio
rho_w = 1000.0  # Density of pure water [kg/m^3]
G = 9.81  # Gravity acceleration [m/s^2]
distance = 1  # Grid points in saline zone - almost no saline intrusion

# HYDRODYNAMIC AND SEDIMENT PARAMETERS
Chezy_lb = 60  # Chezy coefficient (downstream) [m^-1/2 s^-1] Panchenko and Alabyan, 2022 https://doi.org/10.1016/j.mex.2022.101669
Chezy_ub = 40  # Chezy coefficient (upstream) [m^-1/2 s^-1] Panchenko and Alabyan, 2022 https://doi.org/10.1016/j.mex.2022.101669
#ws = 1.0e-3  # Settling velocity for sediment [m/s]
kISS = 51.0e-3  # Inorganic sinking velocity half saturation concentration [g L^-1] - Clark et al 2022 https://doi.org/10.1029/2022JG007139
wMAX = 2.0 / 86400  #	Maximum inorganic suspended sediment sinking velocity [m s^−1] - Clark et al 2022 https://doi.org/10.1029/2022JG007139
tau_ero_lb = 0.005  # Erosion shear stress (downstream) [N/m^2] - Clark et al 2022 https://doi.org/10.1029/2022JG007139
tau_dep_lb = 0.4  # Deposition shear stress (downstream) [N/m^2]
tau_ero_ub = 0.005  # Erosion shear stress (upstream) [N/m^2] - Clark et al 2022 https://doi.org/10.1029/2022JG007139
tau_dep_ub = 1.0  # Deposition shear stress (upstream) [N/m^2]
Mero_lb = 5.787037037037037e-05  # Erosion coefficient (downstream) [mg/m^2 s] - Clark et al 2020 https://doi.org/10.1029/2019JG005442
Mero_ub = 5.787037037037037e-05  # Erosion coefficient (upstream) [mg/m^2 s] - Clark et al 2020 https://doi.org/10.1029/2019JG005442

# BIOGEOCHEMICAL PARAMETERS
Pbmax = 1.3888888888888888e-05  # Max photosynthetic rate [s^-1] - Le Fouest 2013 www.biogeosciences.net/10/4785/2013/
Chla2CMIN = 0.0125  # Minimum Chla : C ratio for phytoplankton [g/g] - Le Fouest 2013 www.biogeosciences.net/10/4785/2013/
KE = 8.0  # Photoacclimation parameter [E/m^2/d^1] - Le Fouest 2013 www.biogeosciences.net/10/4785/2013/
alpha = 2.0  # Photosynthetic efficiency [mg C/(mg Chl)/(E/m^2/d^1)] - Le Fouest 2013 www.biogeosciences.net/10/4785/2013/
KdSi = 1.07  # Silica constant [muM Si]
KPO4 = 0.2  # Phosphate constant [muM P]
KNH4 = 0.5  # Ammonium constant [muM N] - Le Fouest 2013 www.biogeosciences.net/10/4785/2013/
KNO3 = 1  # Nitrate constant [muM N] - Le Fouest 2013 www.biogeosciences.net/10/4785/2013/
KTOC = 186.25  # Organic matter constant [muM C]
KO2_ox = 31.0  # Oxygen constant [muM O2]
KO2_nit = 51.25  # Oxygen constant [muM O2]
KN = 1.13  # Dissolved nitrogen constant [muM N]
KinO2 = 33.0  # Denitrification inhibition term [muM O2]
redsi = 16.0 / 80.0  # Redfield ratio for silica
redn = 16.0 / 106.0  # Redfield ratio for nitrogen
redp = 1.0 / 106.0  # Redfield ratio for phosphorus
kmaint = 4.6e-7  # Maintenance rate constant [s^-1]
kmort = 1.1574074074074074e-07  # Phytoplankton mortality [s^-1] - Le Fouest 2013 www.biogeosciences.net/10/4785/2013/
kexcr = 5e-2  # Excretion constant [-]
kgrowth = 2.9e-1  # Growth constant [-]
KD1 = 1.3  # Background extinction coefficient [m^-1]
KD2 = 0.06  # Suspended matter attenuation [mg^-1 m^-1]
kox = 6.08e-4  # Aerobic degradation rate [muM C s^-1]
kdenit = 5.05e-4  # Denitrification rate [muM C s^-1]
knit = 3.472222222222222e-06  # Nitrification rate [muM C s^-1] - Le Fouest 2013 www.biogeosciences.net/10/4785/2013/

# EXTERNAL FORCINGS
# Qr = -177.0  # River discharge [m^3/s]
AMPL = 0.2  # Tidal amplitude at mouth [m]
pfun = 0.0787  # Tidal frequency [cycle/hr]
#Uw_sal = 8  # Wind speed in saline estuary [m/s]
#Uw_tid = 8  # Wind speed in tidal river [m/s]
#water_temp = 12  # Water temperature [C]
#pCO2 = 380 * 1e-6  # CO2 partial pressure in the atmosphere [atm]
repeatYear = 1  # equal to 1 if 365-day climatology (repeat forcings)
nbday_ice = 10  # number of day of negative water temperature that will stop all reactions (transport, biogeo, sed)


# OTHER PARAMETERS
Euler = 0.5772156649  # Euler's constant
PI = math.pi  # Pi value
pH_ite = 50  # number of iterations to converge to pH following Follows et al., 2006
mass_mol_B = 10.8110  # molar mass of Boron g/mol

# NUMERICAL INTEGRATION
MAXT = 365*2 * 24 * 60 * 60  # Max time [s]
WARMUP = 365 * 24 * 60 * 60  # Warmup period [s]
DELTI = 75  # Delta t [s]
TS = 12  # Save every TS timesteps
DELXI = 200  # Delta x [m]
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