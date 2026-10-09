"""Reader for synthesis products subset around each site.

Raw dirs (per site):  raw/glodap/  raw/socat/  raw/bgcargo/  raw/satchl/
Flags kept:
  GLODAPv2.2023 : per-variable flag f == 2 (acceptable, measured); temperature has
                  no flag (kept if not -9999). qc column = GLODAP flag.
  SOCAT v2026   : WOCE_CO2_water == 2 (query) and dataset QC flag A-D (E dropped).
                  qc column = "<datasetQC>/<WOCE>".
  BGC-Argo      : *_ADJUSTED values only, adjusted QC in {1,2,5,8}; position QC in
                  {1,2,5,8}; profile within 100 km of the site. qc = adjusted QC.
  OC-CCI v6.0   : monthly chlor_a, median of valid 4-km pixels in a ~4x4 pixel box
                  (~17 km) centred on the site; qc = "npix=<n valid>".
"""
import glob
import os

import numpy as np
import pandas as pd

from standardize_obs import iso

SITES = {"HOT": (22.75, -158.0), "BATS": (31.67, -64.17), "HydroS": (32.17, -64.5),
         "Papa": (50.1, -144.9), "PAP": (49.0, -16.5)}
ARGO_GOOD = {"1", "2", "5", "8"}


def _empty():
    return pd.DataFrame()


def _hav_km(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp, dl = p2 - p1, np.radians(lon2 - lon1)
    a = np.sin(dp / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 6371.0 * 2 * np.arcsin(np.sqrt(a))


# ---------------------------------------------------------------- GLODAP
GLODAP_VARS = [  # (column, flag column or None, variable, units)
    ("G2temperature", None, "temp", "degC (in-situ, ITS-90)"),
    ("G2salinity", "G2salinityf", "salt", "PSS-78"),
    ("G2oxygen", "G2oxygenf", "O2", "umol/kg"),
    ("G2nitrate", "G2nitratef", "NO3", "umol/kg"),
    ("G2nitrite", "G2nitritef", "NO2", "umol/kg"),
    ("G2phosphate", "G2phosphatef", "PO4", "umol/kg"),
    ("G2silicate", "G2silicatef", "SiO2", "umol/kg"),
    ("G2tco2", "G2tco2f", "DIC", "umol/kg"),
    ("G2talk", "G2talkf", "ALK", "umol/kg"),
    ("G2phts25p0", "G2phts25p0f", "pH", "pH total scale at 25 degC, 0 dbar"),
    ("G2phtsinsitutp", "G2phtsinsitutpf", "pH", "pH total scale in situ T,P"),
    ("G2doc", "G2docf", "DOC", "umol/kg"),
    ("G2don", "G2donf", "DON", "umol/kg"),
    ("G2tdn", "G2tdnf", "TDN", "umol/kg"),
    ("G2chla", "G2chlaf", "Chl", "ug/kg"),
]


def _glodap(site_dir):
    files = sorted(glob.glob(os.path.join(site_dir, "raw", "glodap", "*.csv")))
    if not files:
        return _empty()
    g = pd.concat([pd.read_csv(f, low_memory=False) for f in files], ignore_index=True)
    g = g.replace(-9999, np.nan)
    hr = g["G2hour"].fillna(0).astype(int)
    mi = g["G2minute"].fillna(0).astype(int)
    t = pd.to_datetime(dict(year=g["G2year"].astype(int), month=g["G2month"].astype(int),
                            day=g["G2day"].astype(int), hour=hr, minute=mi), errors="coerce")
    tiso = iso(t)
    out = []
    for col, fcol, var, units in GLODAP_VARS:
        if col not in g:
            continue
        m = g[col].notna()
        qc = pd.Series("", index=g.index)
        if fcol is not None:
            m &= g[fcol] == 2
            qc = g[fcol].astype("Int64").astype(str)
        if not m.any():
            continue
        out.append(pd.DataFrame({
            "time": tiso[m], "depth_m": g.loc[m, "G2depth"], "variable": var,
            "value": g.loc[m, col], "units": units, "source": "glodap",
            "qc": qc[m], "lat": g.loc[m, "G2latitude"], "lon": g.loc[m, "G2longitude"]}))
    return pd.concat(out, ignore_index=True) if out else _empty()


# ---------------------------------------------------------------- SOCAT
def _socat(site_dir):
    files = sorted(glob.glob(os.path.join(site_dir, "raw", "socat", "socat_*.csv")))
    parts = [pd.read_csv(f, skiprows=[1], low_memory=False) for f in files]
    parts = [p for p in parts if len(p)]
    if not parts:
        return _empty()
    s = pd.concat(parts, ignore_index=True)
    s = s[s["qc_flag"].astype(str).isin(list("ABCD")) &
          (pd.to_numeric(s["WOCE_CO2_water"], errors="coerce") == 2)]
    dep = pd.to_numeric(s["depth"], errors="coerce").fillna(5.0)
    qc = s["qc_flag"].astype(str) + "/" + s["WOCE_CO2_water"].astype(str)
    out = []
    for col, var, units in [("fCO2_recommended", "fCO2", "uatm (in situ SST)"),
                            ("temp", "temp", "degC (SST, intake)"),
                            ("sal", "salt", "PSS-78 (underway/mooring)")]:
        v = pd.to_numeric(s[col], errors="coerce")
        m = v.notna()
        out.append(pd.DataFrame({
            "time": s.loc[m, "time"], "depth_m": dep[m], "variable": var, "value": v[m],
            "units": units, "source": "socat", "qc": qc[m],
            "lat": s.loc[m, "latitude"], "lon": s.loc[m, "longitude"]}))
    return pd.concat(out, ignore_index=True)


# ---------------------------------------------------------------- BGC-Argo
ARGO_VARS = [("temp", "temp", "degC (in-situ)"), ("psal", "salt", "PSS-78"),
             ("doxy", "O2", "umol/kg"), ("nitrate", "NO3", "umol/kg"),
             ("chla", "Chl", "mg/m3 (fluorescence-derived)"),
             ("ph_in_situ_total", "pH", "pH total scale in situ"),
             ("bbp700", "bbp", "m-1 (700 nm)")]


def _argo(site, site_dir):
    files = sorted(glob.glob(os.path.join(site_dir, "raw", "bgcargo", "argo_*.csv")))
    parts = [pd.read_csv(f, skiprows=[1], dtype=str) for f in files]
    parts = [p for p in parts if len(p)]
    if not parts:
        return _empty()
    a = pd.concat(parts, ignore_index=True)
    lat = pd.to_numeric(a["latitude"], errors="coerce")
    lon = pd.to_numeric(a["longitude"], errors="coerce")
    la0, lo0 = SITES[site]
    keep = (_hav_km(la0, lo0, lat, lon) <= 100.0) & \
        a["position_qc"].str.strip().isin(ARGO_GOOD) & \
        a["pres_adjusted_qc"].str.strip().isin(ARGO_GOOD)
    a, lat, lon = a[keep], lat[keep], lon[keep]
    pres = pd.to_numeric(a["pres_adjusted"], errors="coerce")
    good = {}
    for v, var, units in ARGO_VARS:
        val = pd.to_numeric(a[f"{v}_adjusted"], errors="coerce")
        q = a[f"{v}_adjusted_qc"].fillna("").str.strip()
        good[v] = (val, q, val.notna() & q.isin(ARGO_GOOD) & (val < 99999))
    # Keep T/S only on levels that carry at least one good BGC value (the
    # synthetic profiles have ~1-2 dbar CTD resolution, which would swamp the file).
    bgc_level = np.zeros(len(a), dtype=bool)
    for v in ("doxy", "nitrate", "chla", "ph_in_situ_total", "bbp700"):
        bgc_level |= good[v][2].to_numpy()
    out = []
    for v, var, units in ARGO_VARS:
        val, q, m = good[v]
        if v in ("temp", "psal"):
            m = m & bgc_level
        if not m.any():
            continue
        out.append(pd.DataFrame({
            "time": a.loc[m, "time"], "depth_m": pres[m], "variable": var,
            "value": val[m], "units": units, "source": "bgcargo",
            "qc": q[m], "lat": lat[m], "lon": lon[m]}))
    return pd.concat(out, ignore_index=True) if out else _empty()


# ---------------------------------------------------------------- OC-CCI
def _satchl(site, site_dir):
    files = sorted(glob.glob(os.path.join(site_dir, "raw", "satchl", "occci_*.csv")))
    if not files:
        return _empty()
    c = pd.concat([pd.read_csv(f, skiprows=[1]) for f in files], ignore_index=True)
    c["chlor_a"] = pd.to_numeric(c["chlor_a"], errors="coerce")
    g = c.groupby("time")["chlor_a"]
    s = pd.DataFrame({"value": g.median(), "n": g.count()}).reset_index()
    s = s[s["n"] > 0]
    la0, lo0 = SITES[site]
    return pd.DataFrame({
        "time": s["time"], "depth_m": 0.0, "variable": "Chl_sat", "value": s["value"],
        "units": "mg m-3 (OC-CCI v6.0 monthly, OCx blend)", "source": "satchl",
        "qc": "npix=" + s["n"].astype(str), "lat": la0, "lon": lo0})


def load(site, site_dir):
    parts = [_glodap(site_dir), _socat(site_dir), _argo(site, site_dir), _satchl(site, site_dir)]
    parts = [p for p in parts if len(p)]
    return pd.concat(parts, ignore_index=True) if parts else _empty()
