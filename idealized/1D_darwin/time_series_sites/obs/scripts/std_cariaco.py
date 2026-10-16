"""Reader for the CARIACO Ocean Time-Series (Cariaco Basin, Venezuela; station 10.50N 64.67W).

load(site, site_dir) -> tidy DataFrame (see standardize_obs.COLUMNS). Only applies to site "CARIACO".
Raw files (BCO-DMO, CC BY 4.0) under <site_dir>/raw/<source>/:
    bcodmo_bottle  dataset 3093 v1 (Niskin bottles, Nov 1995 - Jan 2017), 3093_v1_niskin.csv
    bcodmo_flux    dataset 3094 (moored sediment traps 150-1200 m, 1995-2010), 3094_v1_Sediment_Trap.csv
    bcodmo_hplc    datasets 3293 (BBRS), 3292 (Mote), 3235 (HPL), 777689 (NASA GSFC) HPLC pigments

Flag convention (3093 q_* columns, strings):
    0 = data correct; 1 = outlier (> 3 SD from the long-term mean at that depth); suffix 'v' = evidence of
    ventilation (0v / 1v).
  KEPT: flag in {'0', '0v'}; dropped '1' and '1v'. 'nd' (not measured) and 'BM_master' dropped.
  Bacterial abundance/biomass have no flags (qc '').
Traps (3094) and HPLC: no flags; 'nd' and -9999 dropped (HPLC -9999 is undocumented; treated as missing).

Column-label error in 3093 (verified 2026-10-09 from the data themselves): the five particulate C/N columns
are permuted relative to their headers. Actual content:
    POC_ug_kg -> POC ug/L       PON_ug_kg -> PON ug/L       POC_ug_L -> molar C:N
    PN_ug_L   -> POC ug/kg      C_N_particulate -> PON ug/kg
Evidence: PN_ug_L / POC_ug_kg = 0.97590 (= 1/1.0247, i.e. per-kg/per-L) for all 3967 samples,
C_N_particulate / PON_ug_kg = 0.97590 for all 3920, and (PN_ug_L/12)/(C_N_particulate/14) == POC_ug_L
to 1e-8. We output POC/PON in ug/kg from the PN_ug_L and C_N_particulate columns (flags travel with them).

Times: hydrocast start UTC (ISO_DateTime_start_hc_UTC) for hydrographic/chemical variables; biocast start UTC
(ISO_DateTime_start_bc_UTC) for PP, Chl and phaeopigments; fallback Year/Month/Day 00:00. Traps: midpoint of
the collection interval from date_open (local date, 00:00) + duration_d/2. HPLC: Date_time_UTC.
HPLC positions in dataset 3235 (and part of 777689) are degrees + decimal minutes (e.g. 10.302 = 10 deg 30.2',
-64.4018 = 64 deg 40.18' W); converted to decimal degrees when the station would otherwise be ~37 km off.
HPLC: the two CAR-144 surface red-tide ("Marea Roja") samples are dropped; Mote (3292) Tot/MV/DV chl a
dropped where the lab comment says chl a precipitated (other pigments kept, per PI comment). Negative
pigment values (a few MV/DV chl a) are kept as reported.
"""
import glob
import os

import numpy as np
import pandas as pd

SITE_POS = (10.50, -64.67)
KEEP_FLAGS = {"0", "0v"}


def _iso(t):
    t = pd.to_datetime(t, utc=True, errors="coerce", format="mixed")
    return t.dt.strftime("%Y-%m-%dT%H:%M:%SZ")


def _num(s):
    s = pd.to_numeric(s, errors="coerce")
    return s.where(s > -998)


def _melt(df, specs, source, time_col="time"):
    """specs: list of (column, qf_column_or_None, variable, units[, time_col])."""
    out = []
    for spec in specs:
        col, qf, var, units = spec[:4]
        tcol = spec[4] if len(spec) > 4 else time_col
        if col not in df.columns:
            continue
        v = _num(df[col])
        ok = v.notna() & df[tcol].notna() & df["depth_m"].notna()
        if qf is not None and qf in df.columns:
            q = df[qf].astype(str).str.strip()
            ok &= q.isin(KEEP_FLAGS)
            qcs = q
        else:
            qcs = pd.Series("", index=df.index)
        out.append(pd.DataFrame({
            "time": df.loc[ok, tcol], "depth_m": df.loc[ok, "depth_m"],
            "variable": var, "value": v[ok], "units": units, "source": source,
            "qc": qcs[ok], "lat": df.loc[ok, "lat"], "lon": df.loc[ok, "lon"],
        }))
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def _bottle(site_dir):
    fs = sorted(glob.glob(os.path.join(site_dir, "raw", "bcodmo_bottle", "*niskin*.csv")))
    if not fs:
        return pd.DataFrame()
    df = pd.read_csv(fs[-1], dtype=str, keep_default_na=False)
    df = df.replace({"nd": np.nan, "BM_master": np.nan, "": np.nan})
    ymd = pd.to_datetime(dict(year=_num(df["Year"]), month=_num(df["Month"]), day=_num(df["Day"])),
                         errors="coerce", utc=True).dt.strftime("%Y-%m-%dT%H:%M:%SZ")
    t_hc = _iso(df["ISO_DateTime_start_hc_UTC"])
    t_bc = _iso(df["ISO_DateTime_start_bc_UTC"])
    df["time"] = t_hc.fillna(t_bc).fillna(ymd)
    df["time_bio"] = t_bc.fillna(t_hc).fillna(ymd)
    df["depth_m"] = _num(df["Depth_real"]).fillna(_num(df["Depth_target"]))
    df["lat"] = _num(df["Latitude"])
    df["lon"] = _num(df["Longitude"])
    # salinity: salinometer where kept, else CTD salinity at the bottle (as for BATS)
    qsb = df["q_Salinity_bottles"].astype(str).str.strip()
    qsc = df["q_Salinity_CTD"].astype(str).str.strip()
    sal = _num(df["Salinity_bottles"]).where(qsb.isin(KEEP_FLAGS))
    ctd = _num(df["Salinity_CTD"]).where(qsc.isin(KEEP_FLAGS))
    df["salt_best"] = sal.fillna(ctd)
    df["q_salt_best"] = np.where(sal.notna(), qsb, qsc)
    u = " (CARIACO; per-L = mmol/m^3)"
    specs = [
        ("Temperature", "q_Temperature", "temp", "degC (ITS-90, in-situ, CTD at bottle)"),
        ("salt_best", "q_salt_best", "salt", "PSU (bottle salinometer; CTD where bottle missing)"),
        ("O2_umol_kg", "q_O2_umol_kg", "O2", "umol/kg (Winkler; 0 in anoxic water below ~250 m)"),
        ("NO3_NO2_USF", "q_NO3_NO3_USF", "NO3", "uM (NO3+NO2, USF)" + u),
        ("NO2_USF", "q_NO2_USF", "NO2", "uM (USF)" + u),
        ("NH4_USF", "q_NH4_USF", "NH4", "uM (USF)" + u),
        ("PO4_USF", "q_PO4_USF", "PO4", "uM (USF)" + u),
        ("SiO4_USF", "q_SiO4_USF", "SiO2", "uM (USF; frozen samples pre-Aug 2001 low below 300 m)" + u),
        ("NO3_UDO", "q_NO3_UDO", "NO3_UDO", "uM (nitrate, UDO/EDIMAR lab, 1995-2003)" + u),
        ("PO4_UDO", "q_PO4_UDO", "PO4_UDO", "uM (UDO/EDIMAR lab, 1995-2003)" + u),
        ("SiO4_UDO", "q_SiO4_UDO", "SiO2_UDO", "uM (UDO/EDIMAR lab, 1995-2003)" + u),
        ("pH_corrected", "q_pH_corrected", "pH",
         "total scale at 25 degC (spectrophotometric, corrected for dye impurities)"),
        ("Alkalinity_umol_kg", "q_Alkalinity_umol_kg", "ALK", "umol/kg (spectrophotometric)"),
        ("TCO2_corrected", "q_TCO2_corrected", "DIC",
         "umol/kg (TCO2_corrected: calculated from dye-corrected pH + ALK, not measured)"),
        ("fCO2_corrected", "q_fCO2_corrected", "fCO2",
         "uatm (fCO2_corrected: calculated from dye-corrected pH + ALK; reference temperature not stated)"),
        ("PN_ug_L", "q_PN_ug_L", "POC", "ug/kg (from column labeled PN_ug_L; see reader docstring)"),
        ("C_N_particulate", "q_C_N_particulate", "PON",
         "ug/kg (from column labeled C_N_particulate; see reader docstring)"),
        ("TPP", "q_TPP", "TPP", "nM (total particulate P)"),
        ("PIP", "q_PIP", "PIP", "nM (particulate inorganic P)"),
        ("DOC", "q_DOC", "DOC", "uM (1995-2001 'irregular', suspended; restarted Mar 2005)"),
        ("TOC", "q_TOC", "TOC", "uM (total organic C, ~DOC)"),
        ("DON", "q_DON", "DON", "uM (by difference from TDN, USF)"),
        ("DOP", "q_DOP", "DOP", "uM (by difference from TDP, USF)"),
        ("PrimaryProductivity", "q_PrimaryProductivity", "PP",
         "mgC/m^3/hr (14C in-situ, 'daily mean hourly average'; x12 = daily per methods)", "time_bio"),
        ("Chlorophyll", "q_Chlorophyll", "Chl", "mg/m^3 (fluorometric extracted chl a)", "time_bio"),
        ("Phaeopigments", "q_Phaeopigments", "Phaeo", "mg/m^3 (fluorometric phaeopigments)", "time_bio"),
        ("Total_Prokaryotes", None, "BactAbund", "1e8 cells/L (total prokaryotes)"),
        ("Bact_Biomass_mgC_m3", None, "BactBiomass", "mgC/m^3 (bacterial biomass)"),
    ]
    return _melt(df, specs, "bcodmo_bottle")


def _flux(site_dir):
    fs = sorted(glob.glob(os.path.join(site_dir, "raw", "bcodmo_flux", "*Sediment_Trap*.csv")))
    if not fs:
        return pd.DataFrame()
    df = pd.read_csv(fs[-1], na_values=["nd"], dtype={"date_open": str}, low_memory=False)
    t0 = pd.to_datetime(df["date_open"].astype(str).str[:8], format="%Y%m%d", errors="coerce", utc=True)
    df["time"] = (t0 + pd.to_timedelta(_num(df["duration_d"]) / 2.0, unit="D")).dt.strftime(
        "%Y-%m-%dT%H:%M:%SZ")
    df["depth_m"] = _num(df["depth_trap"])
    df["lat"] = _num(df["lat"])
    df["lon"] = _num(df["lon"])
    u = " (moored trap, ~14-day cups; time = cup midpoint, local date)"
    specs = [
        ("MF_Total", None, "mass_flux", "g/m^2/day" + u),
        ("MF_Corg", None, "POC_flux", "gC/m^2/day (organic C)" + u),
        ("MF_N", None, "PON_flux", "gN/m^2/day (total N)" + u),
        ("MF_CaCO3", None, "PIC_flux", "g CaCO3/m^2/day" + u),
        ("MF_Sibio", None, "bSi_flux", "g opal/m^2/day (biogenic silica)" + u),
        ("MF_Terr", None, "lith_flux", "g/m^2/day (terrigenous)" + u),
    ]
    return _melt(df, specs, "bcodmo_flux")


HPLC_PIG = {
    "Tot_Chl_b": "chl_b", "Tot_Chl_c": "chl_c", "MV_Chl_a": "mv_chl_a", "DV_Chl_a": "dv_chl_a",
    "Chlide_a": "chlide_a", "Chl_c3": "chl_c3", "Chl_c1c2": "chl_c1c2", "Perid": "perid",
    "But_fuco": "but_fucox", "Fuco": "fucox", "Hex_fuco": "hex_fucox", "Pras": "prasino",
    "Diadino": "diadino", "Allo": "allo", "Diato": "diato", "Zea": "zea", "Lut": "lut",
    "alpha_beta_Car": "ab_carot", "alpha_Car": "a_carot", "beta_Car": "b_carot",
    "Neo": "neo", "Viola": "viola",
}


def _dist_km(lat, lon):
    la0, lo0 = SITE_POS
    return 111.2 * np.hypot(lat - la0, (lon - lo0) * np.cos(np.radians(la0)))


def _hplc(site_dir):
    fs = sorted(glob.glob(os.path.join(site_dir, "raw", "bcodmo_hplc", "*hplc_*.csv")))
    if not fs:
        return pd.DataFrame()
    parts = []
    for f in fs:
        d = pd.read_csv(f, na_values=["nd"], low_memory=False)
        d["_file"] = os.path.basename(f)
        parts.append(d)
    df = pd.concat(parts, ignore_index=True)
    df["time"] = _iso(df["Date_time_UTC"])
    df["depth_m"] = _num(df["depth"])
    lat, lon = _num(df["lat"]), _num(df["lon"])
    # degrees + decimal minutes recorded as decimal degrees (e.g. 10.302 / -64.4018)
    dm = (_dist_km(lat, lon) > 10)
    lat_dm = np.trunc(lat) + (lat - np.trunc(lat)) * 100 / 60
    lon_dm = np.trunc(lon) + (lon - np.trunc(lon)) * 100 / 60
    fix = dm & (_dist_km(lat_dm, lon_dm) < 5)
    df["lat"] = lat.where(~fix, lat_dm)
    df["lon"] = lon.where(~fix, lon_dm)
    # Mote (3292): filters not dry -> chl a precipitated; PI says other pigments reliable
    bad_chla = df["Comments"].astype(str).str.contains("Chl a to pricipitate", na=False)
    for c in ["Tot_Chl_a", "MV_Chl_a", "DV_Chl_a"]:
        df.loc[bad_chla, c] = np.nan
    # CAR-144 (2008-04-08) surface "Marea Roja" (red tide) grab samples: not station profile samples
    red = df["Comments"].astype(str).str.contains("red tide", case=False, na=False)
    red |= df["sample_number"].astype(str).str.contains("Marea Roja", case=False, na=False)
    df = df[~red].copy()
    specs = [("Tot_Chl_a", None, "Chl_HPLC", "mg/m^3 (HPLC Tot_Chl_a)")]
    specs += [(k, None, "HPLC_" + v, "mg/m^3") for k, v in HPLC_PIG.items()]
    return _melt(df, specs, "bcodmo_hplc")


def load(site, site_dir):
    if site != "CARIACO":
        return pd.DataFrame()
    parts = [r(site_dir) for r in (_bottle, _flux, _hplc)]
    parts = [p for p in parts if len(p)]
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
