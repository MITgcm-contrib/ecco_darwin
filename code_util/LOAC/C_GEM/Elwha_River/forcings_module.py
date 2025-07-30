# forcing_module.py
import numpy as np
import pandas as pd
from datetime import timedelta, datetime
import csv
import matplotlib.pyplot as plt
from config import Qr, Uw_sal, pCO2, WARMUP, DELTI, DEBUG_PLOT_FOR_INTERP_ARRAYS, series_info, water_temp


# Global variable to hold the interpolated discharge array
interpolated_discharge = None
interpolated_wind_speed = None
interpolated_pCO2 = None
interpolated_water_temp = None
interpolated_sediment = None



###################################################################################### pCO2
'''
def get_dummy_pCO2(t_seconds, sim_start_dt_dt):
    """
    Returns a synthetic, time-varying pCO2 [atm] using a sinusoidal pattern.
    
    Parameters:
    
        t_seconds : float
            Simulation time in seconds since t = 0
        sim_start_dt_dt : datetime
            Start datetime of the simulation

    Returns:
        float : Atmospheric pCO2 [atm]
    """
    base_ppm = 400     # central value (ppm)
    amplitude = 80 # 20 # +/- swing (ppm)
    period_days = 365  # seasonal fluctuation
    omega = 2 * np.pi / (period_days * 86400)  # convert to radians/sec

    ppm = base_ppm + amplitude * np.sin(omega * t_seconds)
    atm = ppm / 1_000_000  # convert to atm

    return atm
'''

'''
# Load real-world pCO2 time series once
csv_path = "usgs_elwha_pCO2_timeseries_2011_2016.csv"
try:
    pCO2_df = pd.read_csv(csv_path)
    pCO2_df["datetime"] = pd.to_datetime(pCO2_df[["Year", "Month", "Day"]])
    pCO2_df.set_index("datetime", inplace=True)
    pCO2_df.sort_index(inplace=True)
except Exception as e:
    print(f"Warning: Failed to load real pCO2 data: {e}")
    pCO2_df = None
'''


def get_real_pCO2(t, sim_start_dt=None):
    global interpolated_pCO2
    if t <= WARMUP:
        return pCO2  # use constant from config (already in atm)
    index = int((t - WARMUP) // DELTI) # model start time at sim_start_date
    if index < 0:
        index = 0
    elif index >= len(interpolated_pCO2):
        index = len(interpolated_pCO2) - 1
    return interpolated_pCO2[index] / 1_000_000  # convert ppm to atm if needed


###################################################################################### discharge

def get_discharge(t, sim_start_dt=None):
    global interpolated_discharge
    if t <= WARMUP:
        return abs(Qr)  # use constant from config and ensure positive
    index = int((t - WARMUP) // DELTI) # model start time at sim_start_date
    if index < 0:
        index = 0
    elif index >= len(interpolated_discharge):
        index = len(interpolated_discharge) - 1
    return interpolated_discharge[index]
    

###################################################################################### wind speed

def get_wind_speed(t, sim_start_dt=None):
    global interpolated_wind_speed
    if t <= WARMUP:
        return Uw_sal  # use constant from config
    index = int((t - WARMUP) // DELTI) # model start time at sim_start_date
    if index < 0:
        index = 0
    elif index >= len(interpolated_wind_speed):
        index = len(interpolated_wind_speed) - 1
    return interpolated_wind_speed[index]


###################################################################################### water temprature


def get_water_temp(t, sim_start_dt=None):
    global interpolated_water_temp
    if t <= WARMUP:
        return water_temp  # use constant from config
    index = int((t - WARMUP) // DELTI) # model start time at sim_start_date
    if index < 0:
        index = 0
    elif index >= len(interpolated_water_temp):
        index = len(interpolated_water_temp) - 1
    return interpolated_water_temp[index]

###################################################################################### suspended sediment

def get_sediment(t, sim_start_dt=None):
    global interpolated_sediment
    if t <= WARMUP:
        return 0.01  # use constant from init_module
    index = int((t - WARMUP) // DELTI)
    if index < 0:
        index = 0
    elif index >= len(interpolated_sediment):
        index = len(interpolated_sediment) - 1
    return interpolated_sediment[index] / 1000.0  # Convert mg/L to g/L
    #sediment_value = interpolated_sediment[index] / 1000.0  # Convert mg/L to g/L
    # Cap at reasonable maximum (1 g/L = 1000 mg/L)
    #return min(sediment_value, 1.0)



################################################# Loading and linear interpolation #####################################################
def load_and_interpolate_timeseries(series_info, sim_start_dt, delti, maxt):
    model_times = np.arange(0, maxt+1, delti)
    results = {}

    for varname, (csv_path, date_cols, value_col) in series_info.items():
        data_times = []
        data_values = []
        with open(csv_path, 'r', newline='') as f:
            reader = csv.DictReader(f)
            for row in reader:
                try:
                    if isinstance(date_cols, str):
                        date = row[date_cols]
                        try:
                            dt = datetime.fromisoformat(date)
                        except ValueError:
                            dt = datetime.strptime(date, "%Y-%m-%d %H:%M:%S")
                    else:
                        year = int(row[date_cols[0]])
                        month = int(row[date_cols[1]])
                        day = int(row[date_cols[2]])
                        dt = datetime(year, month, day)
                    seconds_since_start = int((dt - sim_start_dt).total_seconds())
                    value = float(row[value_col])
                    if value == -999:
                        continue
                    data_times.append(seconds_since_start)
                    data_values.append(value)
                except Exception as e:
                    print(f"Skipping row in {csv_path} due to error: {e}")
        if not data_times:
            raise ValueError(f"No valid data found for {varname} in {csv_path}")
        data_times = np.array(data_times)
        data_values = np.array(data_values)
        interp_values = np.interp(model_times, data_times, data_values)
        results[varname] = interp_values

    return results

def initialize_forcings(sim_start, delti, maxt):
    global interpolated_discharge, interpolated_wind_speed, interpolated_pCO2, interpolated_water_temp, interpolated_sediment
    results = load_and_interpolate_timeseries(series_info, sim_start, delti, maxt)
    interpolated_discharge = results['discharge']
    interpolated_wind_speed = results['wind_speed']
    interpolated_pCO2 = results['pCO2']
    interpolated_water_temp = results['water_temp']
    interpolated_sediment = results['sediment']


    
    if DEBUG_PLOT_FOR_INTERP_ARRAYS:
        # print lines for array verification 
        YELLOW = "\033[93m"
        RESET = "\033[0m"
        # print lines for array verification 
        print(f"\n{YELLOW}Discharge array: len = {len(interpolated_discharge)}, min = {np.min(interpolated_discharge)}, max = {np.max(interpolated_discharge)}{RESET}")
        print(f"{YELLOW}Wind speed array: len = {len(interpolated_wind_speed)}, min = {np.min(interpolated_wind_speed)}, max = {np.max(interpolated_wind_speed)}{RESET}")
        print(f"{YELLOW}pCO2 array: len = {len(interpolated_pCO2)}, min = {np.min(interpolated_pCO2)}, max = {np.max(interpolated_pCO2)}{RESET}")
        print(f"{YELLOW}water_temp array: len = {len(interpolated_water_temp)}, min = {np.min(interpolated_water_temp)}, max = {np.max(interpolated_water_temp)}{RESET}")
        print(f"{YELLOW}sediment array: len = {len(interpolated_sediment)}, min = {np.min(interpolated_sediment)}, max = {np.max(interpolated_sediment)}{RESET}\n")


        
        # plots to varify arrays

        fig, axs = plt.subplots(5, 1, figsize=(12, 10), sharex=True)

        axs[0].plot(interpolated_discharge)
        axs[0].set_title('Discharge')
        axs[0].set_ylabel('Discharge (m^3/s)')

        axs[1].plot(interpolated_wind_speed)
        axs[1].set_title('Wind Speed')
        axs[1].set_ylabel('Wind Speed (m/s)')

        axs[2].plot(interpolated_pCO2)
        axs[2].set_title('pCO2')
        axs[2].set_ylabel('pCO2 (atm)')

        axs[3].plot(interpolated_water_temp)
        axs[3].set_title('water_temp')
        axs[3].set_ylabel('temp (C)')

        axs[4].plot(interpolated_sediment)
        axs[4].set_title('Suspended Sediment')
        axs[4].set_ylabel('SSC (mg/L)')

        plt.xlabel("Model Timestep")
        plt.tight_layout()
        plt.show()
