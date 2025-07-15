# forcing_module.py
import numpy as np
import pandas as pd
from datetime import timedelta
from config import WARMUP, USE_WARMUP_DISCHARGE
###################################################################################### pCO2

def get_dummy_pCO2(t_seconds, sim_start_dt):
    """
    Returns a synthetic, time-varying pCO2 [atm] using a sinusoidal pattern.
    
    Parameters:
    
        t_seconds : float
            Simulation time in seconds since t = 0
        sim_start_dt : datetime
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



def get_real_pCO2(t_seconds, sim_start_dt):
    """
    Return interpolated atmospheric pCO2 [atm] from real data.
    Falls back to dummy value if data unavailable.
    """
    if pCO2_df is None or pCO2_df.empty:
        return get_dummy_pCO2(t_seconds, sim_start_dt)

    current_dt = sim_start_dt + timedelta(seconds=t_seconds)

    if current_dt < pCO2_df.index[0]:
        ppm = pCO2_df["CO2_Value"].iloc[0]
    elif current_dt > pCO2_df.index[-1]:
        ppm = pCO2_df["CO2_Value"].iloc[-1]
    else:
        ppm = pCO2_df.loc[:current_dt, "CO2_Value"].iloc[-1]

    return ppm / 1_000_000

###################################################################################### pCO2

###################################################################################### discharge

# Load real-world discharge time series once
discharge_df = None
try:
    discharge_df = pd.read_csv("Elwha_Cleaned_Discharge.csv")
    discharge_df["datetime"] = pd.to_datetime(discharge_df["datetime"])
    discharge_df.set_index("datetime", inplace=True)
    discharge_df.sort_index(inplace=True)
except Exception as e:
    print(f"[WARN] Failed to load discharge data: {e}")

def get_discharge(t_seconds, sim_start_dt):
    """
    Return constant discharge during warmup, dynamic discharge after.
    If USE_WARMUP is False, always use dynamic discharge (CSV or fallback).
    """
    if USE_WARMUP_DISCHARGE and t_seconds <= WARMUP:
        return 47.7  # constant during warmup (positive)
    else:
        # Use dynamic discharge after warmup or always if USE_WARMUP is False
        if discharge_df is None or discharge_df.empty:
            return 47.7  # fallback (positive)

        current_dt = sim_start_dt + timedelta(seconds=t_seconds)
        
        if current_dt < discharge_df.index[0]:
            return abs(discharge_df["discharge_cms"].iloc[0])  # ensure positive
        elif current_dt > discharge_df.index[-1]:
            return abs(discharge_df["discharge_cms"].iloc[-1])  # ensure positive

        subset = discharge_df.loc[:current_dt]
        if subset.empty:
            return 47.7

        val = subset["discharge_cms"].iloc[-1]
        if pd.isna(val):
            return 47.7

        return abs(val)  # ensure positive
    
###################################################################################### discharge

###################################################################################### wind speed

# Load wind speed CSV 
try:
    wind_df = pd.read_csv("POWER_Point_Daily_CLEANED.csv", parse_dates=['datetime'])
    wind_df.set_index('datetime', inplace=True)
except Exception as e:
    print(f"[WARN] Failed to load wind speed data: {e}")
    wind_df = None

def get_wind_speed(t_seconds, sim_start_dt):
    """
    Return wind speed (WS2M) for the current simulation time.
    Looks up the most recent value at or before the current datetime.
    Falls back to a default if data is missing.
    """
    DEFAULT_WIND = 4.0  # fallback value (adjust as needed)
    if wind_df is None or wind_df.empty:
        return DEFAULT_WIND

    current_dt = sim_start_dt + timedelta(seconds=t_seconds)

    # If before/after data range, use first/last value
    if current_dt < wind_df.index[0]:
        return wind_df["WS2M"].iloc[0]
    elif current_dt > wind_df.index[-1]:
        return wind_df["WS2M"].iloc[-1]

    # find the most recent value at or before current_dt
    subset = wind_df.loc[:current_dt]
    if subset.empty:
        return DEFAULT_WIND

    val = subset["WS2M"].iloc[-1]
    if pd.isna(val):
        return DEFAULT_WIND

    return val

###################################################################################### wind speed
