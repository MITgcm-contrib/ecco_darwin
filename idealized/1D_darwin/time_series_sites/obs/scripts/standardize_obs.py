#!/usr/bin/env python3
"""Convert raw downloaded observations into one tidy CSV per site.

Usage (always run isolated, paths as arguments):
    python3 -I standardize_obs.py [--obs-root OBS_ROOT] [--site SITE ...] [--summary]

Layout:
    OBS_ROOT/<SITE>/raw/<source>/...   raw downloads (untrusted data, never executed)
    OBS_ROOT/<SITE>/<SITE>_obs.csv     tidy output written here
    OBS_ROOT/scripts/std_<module>.py   per-source readers (this directory; std_tm = trace metals)

Each reader module exposes  load(site, site_dir) -> pandas.DataFrame  with the
columns in COLUMNS (lat/lon optional; NaN means "at the station").
Readers return an empty DataFrame when their raw files are absent.

Tidy columns:
    time     ISO-8601 UTC (YYYY-MM-DDTHH:MM:SSZ)
    depth_m  positive down, metres (pressure in dbar treated ~ m unless noted)
    variable see VARIABLES
    value    as reported (no unit conversion)
    units    original units
    source   short source tag (matches raw/<source>)
    qc       original flag value of the kept sample (as string), or '' if none
    lat, lon sample position (decimal deg, lon -180..180) if known
"""
import argparse
import importlib
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:  # -I drops the script dir from sys.path; add only ours
    sys.path.insert(0, HERE)
# -I also disables user site-packages, where pandas/numpy/netCDF4 live on this Mac.
# Re-add only that (trusted) dir; script dir of raw data and cwd stay excluded.
import site as _site
_usp = _site.getusersitepackages()
if os.path.isdir(_usp) and _usp not in sys.path:
    sys.path.append(_usp)

import numpy as np
import pandas as pd

COLUMNS = ["time", "depth_m", "variable", "value", "units", "source", "qc", "lat", "lon"]

SITES = ["HOT", "BATS", "HydroS", "Papa", "PAP"]

# Reader modules run for each site (each decides whether it applies to that site).
READERS = {
    "HOT": ["std_hot", "std_synth", "std_tm", "std_geotraces"],
    "BATS": ["std_bios", "std_synth", "std_tm", "std_geotraces"],
    "HydroS": ["std_bios", "std_synth", "std_geotraces"],
    "Papa": ["std_papa", "std_synth", "std_geotraces"],
    "PAP": ["std_pap", "std_synth", "std_geotraces", "std_bodc"],
}

# Standard variable names. Core list from the request + extensions.
VARIABLES = {
    "temp": "in-situ temperature",
    "theta": "potential temperature",
    "salt": "salinity (practical, PSS-78 unless noted)",
    "NO3": "nitrate (or NO3+NO2; see README per source)",
    "NO2": "nitrite",
    "NH4": "ammonium",
    "PO4": "phosphate (SRP)",
    "SiO2": "silicate",
    "DIC": "dissolved inorganic carbon",
    "ALK": "total alkalinity",
    "O2": "dissolved oxygen",
    "Chl": "chlorophyll-a (fluorometric extracted, or sensor fluorescence-derived; see source)",
    "Chl_HPLC": "chlorophyll-a (HPLC, total chl a)",
    "Chl_sat": "satellite surface chlorophyll-a",
    "PP": "primary production (14C unless noted)",
    "pCO2": "seawater pCO2 (or fCO2 if variable fCO2)",
    "fCO2": "seawater fCO2",
    "pCO2_air": "atmospheric pCO2/xCO2",
    "pH": "pH (scale in units column)",
    "MLD": "mixed-layer depth",
    "POC": "particulate organic carbon", "PON": "particulate organic nitrogen",
    "POP": "particulate organic phosphorus", "PIC": "particulate inorganic carbon",
    "DOC": "dissolved organic carbon", "DON": "dissolved organic nitrogen",
    "DOP": "dissolved organic phosphorus", "TDN": "total dissolved nitrogen",
    "TDP": "total dissolved phosphorus",
    "bSi": "biogenic silica", "Fe": "dissolved iron",
    "Mn": "dissolved manganese", "Co": "dissolved cobalt (labile where stated)",
    "Ni": "dissolved nickel", "Cu": "dissolved copper", "Zn": "dissolved zinc",
    "Cd": "dissolved cadmium", "Pb": "dissolved lead", "Al": "dissolved aluminium",
    "Ti": "dissolved titanium",
    "POC_flux": "sinking POC flux", "PON_flux": "sinking PON flux",
    "POP_flux": "sinking POP flux", "PIC_flux": "sinking PIC flux",
    "bSi_flux": "sinking biogenic silica flux", "mass_flux": "total mass flux",
    "N2O": "nitrous oxide", "BactProd": "bacterial production",
    "BactAbund": "bacterial abundance", "CO2_flux": "air-sea CO2 flux",
    "bbp": "particulate backscatter", "PAR": "photosynthetically available radiation",
}
# HPLC pigments use the prefix "HPLC_<pigment>" (e.g. HPLC_fucox); allowed generically.


def iso(ts):
    """pandas datetime-like -> ISO UTC strings."""
    t = pd.to_datetime(ts, utc=True, errors="coerce")
    return t.dt.strftime("%Y-%m-%dT%H:%M:%SZ")


def finalize(df):
    for c in COLUMNS:
        if c not in df.columns:
            df[c] = np.nan if c in ("lat", "lon") else ""
    df = df[COLUMNS].copy()
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    df["depth_m"] = pd.to_numeric(df["depth_m"], errors="coerce")
    df = df.dropna(subset=["time", "value"])
    df = df[df["time"].astype(str).str.len() > 0]
    df["qc"] = df["qc"].fillna("").astype(str)
    df = df.drop_duplicates()
    return df.sort_values(["variable", "time", "depth_m"], kind="mergesort")


def summarize(df):
    if df.empty:
        return pd.DataFrame()
    g = df.groupby(["variable", "source"])
    s = g.agg(n=("value", "size"), t0=("time", "min"), t1=("time", "max"),
              zmin=("depth_m", "min"), zmax=("depth_m", "max"))
    s["t0"] = s["t0"].str[:10]
    s["t1"] = s["t1"].str[:10]
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--obs-root", default=os.path.dirname(HERE))
    ap.add_argument("--site", nargs="*", default=SITES)
    ap.add_argument("--summary", action="store_true", help="print per-variable summary")
    ap.add_argument("--summary-csv", default=None, help="write summary table to this CSV")
    a = ap.parse_args()
    allsum = []
    for site in a.site:
        site_dir = os.path.join(a.obs_root, site)
        parts = []
        for modname in READERS[site]:
            try:
                mod = importlib.import_module(modname)
            except ModuleNotFoundError as e:
                print(f"[{site}] reader {modname} not available: {e}")
                continue
            d = mod.load(site, site_dir)
            if d is not None and len(d):
                print(f"[{site}] {modname}: {len(d)} rows")
                parts.append(d)
        if not parts:
            print(f"[{site}] no data")
            continue
        df = finalize(pd.concat(parts, ignore_index=True))
        out = os.path.join(site_dir, f"{site}_obs.csv")
        df.to_csv(out, index=False, float_format="%.6g")
        print(f"[{site}] wrote {out}: {len(df)} rows")
        if a.summary or a.summary_csv:
            s = summarize(df)
            s.insert(0, "site", site)
            allsum.append(s)
            if a.summary:
                with pd.option_context("display.width", 200, "display.max_rows", 500):
                    print(s)
    if a.summary_csv and allsum:
        pd.concat(allsum).to_csv(a.summary_csv)


if __name__ == "__main__":
    main()
