"""Readers for BIOS sites: BATS and Hydrostation S (raw files from BCO-DMO / NCEI).

load(site, site_dir) -> tidy DataFrame (see standardize_obs.COLUMNS).

Flag convention (BATS / Hydro S BCO-DMO files, per-parameter QF_* columns):
    1 = unverified, 2 = verified/acceptable, 3 = questionable, 4 = bad, 9 = no data.
    Bottle flag (QF_bottle): -3 = suspect bottle, 1 = unverified, 2 = verified.
  KEPT: parameter flag in {1, 2} AND bottle flag != -3 AND depth flag not in {3, 4}.
  Missing values (-999 / NA / blank) are dropped.
BTM mooring (MAPCO2): xCO2 QF 2 = good (WOCE), 3 = questionable, 4 = bad -> keep 2 only.
PP: QF10_pp kept in {1, 2}. Pigments: QF_p* kept in {1, 2}.
Scripps CO2 surface DIC/ALK: pre-screened by provider, no flags.
Traps (BATS PITS, OFP): no flags; values 'nd' / blank dropped.
"""
import glob
import os

import numpy as np
import pandas as pd

GOOD = {1, 2}
SITE_POS = {"BATS": (31.67, -64.17), "HydroS": (32.17, -64.5)}
MAX_KM = {"BATS": 100.0, "HydroS": 50.0}


def _iso(t):
    t = pd.to_datetime(t, utc=True, errors="coerce", format="mixed")
    return t.dt.strftime("%Y-%m-%dT%H:%M:%SZ")


def _dist_km(lat, lon, site):
    la0, lo0 = SITE_POS[site]
    return 111.2 * np.hypot(lat - la0, (lon - lo0) * np.cos(np.radians(la0)))


def _num(s):
    s = pd.to_numeric(s, errors="coerce")
    return s.where(s > -998)


def _melt(df, specs, source, time_col="time", depth_col="depth_m"):
    """specs: list of (column, qf_column_or_None, variable, units)."""
    out = []
    for col, qf, var, units in specs:
        if col not in df.columns:
            continue
        v = _num(df[col])
        ok = v.notna()
        if qf is not None and qf in df.columns:
            q = pd.to_numeric(df[qf], errors="coerce")
            ok &= q.isin(GOOD)
            qcs = q.astype("Int64").astype(str)
        else:
            qcs = pd.Series("", index=df.index)
        if "_mask" in df.columns:
            ok &= df["_mask"]
        d = pd.DataFrame({
            "time": df.loc[ok, time_col], "depth_m": df.loc[ok, depth_col],
            "variable": var, "value": v[ok], "units": units, "source": source,
            "qc": qcs[ok],
            "lat": df.loc[ok, "lat"] if "lat" in df.columns else np.nan,
            "lon": df.loc[ok, "lon"] if "lon" in df.columns else np.nan,
        })
        out.append(d)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def _bottle_common(df, site):
    df["time"] = _iso(df["ISO_DateTime_UTC"])
    df["depth_m"] = _num(df["Depth"])
    df["lat"] = _num(df["Latitude"])
    df["lon"] = _num(df["Longitude"])
    qb = df["QF_bottle"] if "QF_bottle" in df.columns else df["QF_Bottle"]
    mask = (pd.to_numeric(qb, errors="coerce") != -3)
    mask &= ~pd.to_numeric(df["QF_Depth"], errors="coerce").isin([3, 4])
    dist = _dist_km(df["lat"], df["lon"], site)
    mask &= ~(dist > MAX_KM[site])  # NaN position -> kept (assumed on station)
    df["_mask"] = mask & df["depth_m"].notna()
    return df


def _bats_bottle(site_dir, site):
    fs = sorted(glob.glob(os.path.join(site_dir, "raw", "bcodmo_bottle", "*bats_bottle*.csv")))
    if not fs:
        return pd.DataFrame()
    df = _bottle_common(pd.read_csv(fs[-1], low_memory=False), site)
    # Salinity: bottle (salinometer) where valid, else CTD salinity at bottle
    sal = _num(df["Salinity"]).where(pd.to_numeric(df["QF_Salinity"], errors="coerce").isin(GOOD))
    ctd = _num(df["CTD_Salinity"]).where(pd.to_numeric(df["QF_CTD_Sal"], errors="coerce").isin(GOOD))
    df["salt_best"] = sal.fillna(ctd)
    df["QF_salt_best"] = np.where(sal.notna(), df["QF_Salinity"], df["QF_CTD_Sal"])
    # TN column: DON reported before BATS cruise 121, TN (total dissolved N) after
    cn = pd.to_numeric(df["Cruise_num"], errors="coerce")
    early = (cn >= 10000) & (cn < 10121)  # core cruise numbers 10001..; bloom cruises are 2xxxx
    df["DON_early"] = _num(df["TN"]).where(early)
    df["TN_late"] = _num(df["TN"]).where(~early)
    specs = [
        ("Temperature", "QF_Temp", "temp", "degC (ITS-90, in-situ, CTD at bottle)"),
        ("salt_best", "QF_salt_best", "salt", "PSS-78 (bottle; CTD where bottle missing)"),
        ("Oxygen_1", "QF_Oxygen", "O2", "umol/kg"),
        ("CO2", "QF_DIC", "DIC", "umol/kg"),
        ("Alkalinity", "QF_Alk", "ALK", "umol/kg"),
        ("NO3_plus_NO2", "QF_NO3_NO2", "NO3", "umol/kg (NO3+NO2)"),
        ("NO2", "QF_NO2", "NO2", "umol/kg"),
        ("PO4", "QF_PO4", "PO4", "umol/kg"),
        ("Silicate", "QF_Silicate", "SiO2", "umol/kg"),
        ("POC", "QF_POC", "POC", "ug/kg"),
        ("PON", "QF_PON", "PON", "ug/kg"),
        ("POP", "QF_POP", "POP", "umol/kg"),
        ("TOC", "QF_TOC", "TOC", "umol/kg (total organic C, ~DOC)"),
        ("DON_early", "QF_TN", "DON", "umol/kg (DON, cruises <121)"),
        ("TN_late", "QF_TN", "TDN", "umol/kg (total N incl. NO3, cruises >=121)"),
        ("TDP", "QF_TDP", "TDP", "nmol/kg"),
        ("SRP", "QF_LLP_SRP", "SRP_LL", "nmol/kg (low-level SRP)"),
        ("Bio_Si", "QF_bio_Si", "bSi", "umol/kg"),
        ("Litho_Si", "QF_litho_Si", "lSi", "umol/kg (lithogenic Si)"),
        ("Bact_Enum", "QF_Bact_enum", "BactAbund", "1e8 cells/kg"),
        ("Prochlorococcus", "QF_Prochloro", "Prochlorococcus", "cells/mL"),
        ("Synechococcus", "QF_Synecho", "Synechococcus", "cells/mL"),
        ("Picoeukaryotes", "QF_Picoeuk", "Picoeuk", "cells/mL"),
        ("Nanoeukaryotes", "QF_Nanoeuk", "Nanoeuk", "cells/mL"),
    ]
    return _melt(df, specs, "bcodmo_bottle")


def _hydro_bottle(site_dir, site):
    fs = sorted(glob.glob(os.path.join(site_dir, "raw", "bcodmo_bottle", "*hydrostation_s_bottle*.csv")))
    if not fs:
        return pd.DataFrame()
    df = _bottle_common(pd.read_csv(fs[-1], low_memory=False), site)
    sal = _num(df["Salinity_1"]).where(pd.to_numeric(df["QF_Salinity"], errors="coerce").isin(GOOD))
    ctd = _num(df["CTD_Salinity"]).where(pd.to_numeric(df["QF_CTD_Sal"], errors="coerce").isin(GOOD))
    df["salt_best"] = sal.fillna(ctd)
    df["QF_salt_best"] = np.where(sal.notna(), df["QF_Salinity"], df["QF_CTD_Sal"])
    specs = [
        ("Temperature", "QF_Temp", "temp", "degC (in-situ; reversing thermometer pre-CTD era)"),
        ("salt_best", "QF_salt_best", "salt", "PSS-78 (bottle; CTD where bottle missing)"),
        ("Oxygen", "QF_Oxygen", "O2", "umol/kg"),
    ]
    return _melt(df, specs, "bcodmo_bottle")


def _bats_pp(site_dir, site):
    fs = sorted(glob.glob(os.path.join(site_dir, "raw", "bcodmo_pp", "*primary_production*.csv")))
    if not fs:
        return pd.DataFrame()
    df = pd.read_csv(fs[-1], low_memory=False)
    # incubation start: in-situ array deployment time if present, else CTD-in time, else date
    t = df["ISO_DateTime_UTC_Array_in"].fillna(df["ISO_DateTime_UTC_CTD_in"]).fillna(df["Date"])
    df["time"] = _iso(t)
    df["depth_m"] = _num(df["Depth"])
    df["lat"] = _num(df["Latitude_Array_in"]).fillna(_num(df["Latitude_CTD_in"]))
    # some post-2013 array longitudes are reported positive (sign error); site is at ~64 W
    df["lon"] = -_num(df["Longitude_Array_in"]).fillna(_num(df["Longitude_CTD_in"])).abs()
    df["_mask"] = ~(_dist_km(df["lat"], df["lon"], site) > MAX_KM[site]) & df["depth_m"].notna()
    df["_mask"] &= pd.to_numeric(df["QF_Niskin_GoFlo"], errors="coerce") != -3
    specs = [("pp", "QF10_pp", "PP", "mgC/m^3/day (14C, in-situ array, mean light - dark)")]
    return _melt(df, specs, "bcodmo_pp")


PIG = {
    "p1": "chl_c3", "p2": "chlide_a", "p3": "chl_c1c2", "p4": "perid", "p5": "but_fucox",
    "p6": "fucox", "p7": "hex_fucox", "p8": "prasino", "p9": "diadino", "p10": "allo",
    "p11": "diato", "p12": "zea_lut", "p13": "chl_b", "p15": "ab_carot", "p18": "lut",
    "p19": "zea", "p20": "a_carot", "p21": "b_carot",
}


def _bats_pigments(site_dir, site):
    fs = sorted(glob.glob(os.path.join(site_dir, "raw", "bcodmo_pigments", "*pigments*.csv")))
    if not fs:
        return pd.DataFrame()
    df = pd.read_csv(fs[-1], low_memory=False)
    df["time"] = _iso(df["ISO_DateTime_UTC"])
    df["depth_m"] = _num(df["Depth"])
    df["lat"] = _num(df["Latitude"])
    df["lon"] = _num(df["Longitude"])
    m = ~(_dist_km(df["lat"], df["lon"], site) > MAX_KM[site]) & df["depth_m"].notna()
    m &= pd.to_numeric(df["QF_Niskin_GoFlo"], errors="coerce") != -3
    m &= ~pd.to_numeric(df["QF_depth"], errors="coerce").isin([3, 4])
    df["_mask"] = m
    specs = [
        ("p16_Chl", "QF_p16_Chl", "Chl", "ug/kg (fluorometric chl a)"),
        ("p17_Phae", "QF_p17_Phae", "Phaeo", "ug/kg (fluorometric phaeopigments)"),
        ("p14", "QF_p14", "Chl_HPLC", "ng/kg (HPLC total chl a, MV+DV)"),
    ] + [(k, "QF_" + k, "HPLC_" + v, "ng/kg") for k, v in PIG.items()]
    return _melt(df, specs, "bcodmo_pigments")


def _bats_flux(site_dir, site):
    fs = sorted(glob.glob(os.path.join(site_dir, "raw", "bcodmo_flux", "*particle_flux*.csv")))
    if not fs:
        return pd.DataFrame()
    df = pd.read_csv(fs[-1], low_memory=False)
    t0 = pd.to_datetime(df["Date_deployed"], utc=True, errors="coerce")
    t1 = pd.to_datetime(df["Date_recovered"], utc=True, errors="coerce")
    df["time"] = (t0 + (t1 - t0) / 2).dt.strftime("%Y-%m-%dT%H:%M:%SZ")  # deployment midpoint
    df["depth_m"] = _num(df["Depth"])
    df["lat"] = _num(df["Latitude_deployed"])
    df["lon"] = _num(df["Longitude_deployed"])
    specs = [
        ("M_avg", None, "mass_flux", "mg/m^2/day (PITS, ~3-day deployment)"),
        ("C_avg", None, "POC_flux", "mgC/m^2/day (PITS, not blank-corrected)"),
        ("N_avg", None, "PON_flux", "mgN/m^2/day (PITS, not blank-corrected)"),
        ("P_avg", None, "POP_flux", "mmolP/m^2/day (PITS)"),
    ]
    return _melt(df, specs, "bcodmo_flux")


def _ofp(site_dir, site):
    fs = sorted(glob.glob(os.path.join(site_dir, "raw", "ncei_ofp", "*primary-particle-flux*.tsv")))
    if not fs:
        return pd.DataFrame()
    df = pd.read_csv(fs[-1], sep="\t", na_values=["nd", "", "NA"], low_memory=False)
    t0 = pd.to_datetime(dict(year=df["StartYr"], month=df["StartMo"], day=df["StartDay"]),
                        utc=True, errors="coerce")
    df["time"] = (t0 + pd.to_timedelta(df["Duration"] / 2.0, unit="D")).dt.strftime("%Y-%m-%dT%H:%M:%SZ")
    df["depth_m"] = pd.to_numeric(df["depth"].astype(str).str.replace("m", ""), errors="coerce")
    df["lat"], df["lon"] = 31.83, -64.17  # nominal OFP mooring (~31 50'N 64 10'W)
    u = " (OFP deep trap; sample duration ~2 wk-2 mo; contact PI before use)"
    specs = [
        ("MassFlux", None, "mass_flux", "mg/m^2/day" + u),
        ("CorgFlux", None, "POC_flux", "mgC/m^2/day" + u),
        ("Nflux", None, "PON_flux", "mgN/m^2/day" + u),
        ("CarbFlux", None, "PIC_flux", "mg CaCO3/m^2/day" + u),
        ("PtotalFlux", None, "POP_flux", "mgP/m^2/day (total P)" + u),
        ("OpalFlux", None, "bSi_flux", "mg opal/m^2/day" + u),
        ("LithFlux", None, "lith_flux", "mg/m^2/day (lithogenic)" + u),
    ]
    return _melt(df, specs, "ncei_ofp")


def _btm(site_dir, site):
    fs = sorted(glob.glob(os.path.join(site_dir, "raw", "ncei_btm_mooring", "BTM_*.csv")))
    if not fs:
        return pd.DataFrame()
    parts = []
    for f in fs:
        d = pd.read_csv(f, skiprows=[1], low_memory=False)  # second row = units
        d.columns = [c.strip() for c in d.columns]
        parts.append(d)
    df = pd.concat(parts, ignore_index=True)
    df["time"] = pd.to_datetime(df["Date"].astype(str) + " " + df["Time"].astype(str),
                                format="%m/%d/%Y %H:%M", errors="coerce", utc=True).dt.strftime("%Y-%m-%dT%H:%M:%SZ")
    df["depth_m"] = 0.5  # MAPCO2 equilibrator / SBE at ~0.5-1 m
    df["lat"] = _num(df["Latitude"])
    df["lon"] = _num(df["Longitude"])
    specs = [
        ("fCO2_SW_sat", "xCO2_SW_QF", "fCO2", "uatm (MAPCO2, 100% humidity, SST)"),
        ("fCO2_Air_sat", "xCO2_Air_QF", "pCO2_air", "uatm (fCO2 air at SST, saturated)"),
        ("SST", None, "temp", "degC (mooring SST)"),
        ("SSS", None, "salt", "PSS (mooring SSS)"),
    ]
    # MAPCO2 flag: 2 good only
    global GOOD
    saved = GOOD
    GOOD = {2}
    try:
        out = _melt(df, specs, "ncei_btm_mooring")
    finally:
        GOOD = saved
    return out


def _scripps(site_dir, site):
    """Scripps CO2 program surface DIC/ALK/d13C (Keeling/Dickson lab), BERM.csv or BATS.csv."""
    fs = sorted(glob.glob(os.path.join(site_dir, "raw", "scripps_co2", "*.csv")))
    if not fs:
        return pd.DataFrame()
    parts = []
    for f in fs:
        with open(f, encoding="utf-8-sig", errors="replace") as fh:
            lines = fh.readlines()
        h = next(i for i, l in enumerate(lines) if l.startswith("SAMPDATE"))
        d = pd.read_csv(f, skiprows=h, encoding="utf-8-sig", low_memory=False)
        parts.append(d)
    df = pd.concat(parts, ignore_index=True)
    sd = df["SAMPDATE"].astype(str).str.strip().str.replace(r"\.0$", "", regex=True)
    t = pd.to_datetime(sd, format="%Y%m%d", errors="coerce", utc=True)
    t = t.fillna(pd.to_datetime(sd, format="%m/%d/%y", errors="coerce", utc=True))  # some rows mm/dd/yy
    df["time"] = t.dt.strftime("%Y-%m-%dT%H:%M:%SZ")
    df["depth_m"] = _num(df["Z"])
    df["lat"] = _num(df["LAT"])
    df["lon"] = -_num(df["LONG"]).abs()  # file gives degrees W as positive
    df["_mask"] = ~(_dist_km(df["lat"], df["lon"], site) > MAX_KM[site]) & df["depth_m"].notna()
    specs = [
        ("AVGDIC", None, "DIC", "umol/kg (Scripps, date only - time set 00:00 UTC)"),
        ("AVGALK", None, "ALK", "umol/kg (Scripps, date only - time set 00:00 UTC)"),
        ("ISO", None, "d13C_DIC", "permil VPDB"),
    ]
    return _melt(df, specs, "scripps_co2")


def load(site, site_dir):
    if site == "BATS":
        readers = [_bats_bottle, _bats_pp, _bats_pigments, _bats_flux, _ofp, _btm, _scripps]
    elif site == "HydroS":
        readers = [_hydro_bottle, _scripps]
    else:
        return pd.DataFrame()
    parts = [r(site_dir, site) for r in readers]
    parts = [p for p in parts if len(p)]
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
