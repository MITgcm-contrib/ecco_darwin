"""Reader for SOTS (Southern Ocean Time Series, IMOS) raw downloads.

load(site, site_dir) -> tidy DataFrame (see standardize_obs.COLUMNS).

Model column: 46.94S 142.28E. All files are IMOS FV01 (quality-controlled) NetCDF from the AODN THREDDS
server, https://thredds.aodn.org.au/thredds/catalog/IMOS/DWM/SOTS/catalog.html (CC BY 4.0).

Sources (all under <site_dir>/raw/):
  imos_mooring/  SOFS, Pulse, FluxPulse and SAZ47 mooring sensors: T/S/P loggers and CTDs, O2 (SBE43, optodes,
                 SBE37-ODO), FLNTUS/ECO chl fluorescence, ECO-PARS PAR, ISUS nitrate, SeaPHOx pH, IMOS mixed-layer
                 depth product. High-frequency records are reduced to DAILY MEANS (time = day 12:00:00Z).
  imos_ras/      RAS (Remote Access Sampler) water samples at ~30 m (Pulse) / ~5 m (SOFS): DIC, ALK, NO3+NO2,
                 PO4, SiO4 (+ T/S of the sample). One row per sample, as reported.
  imos_traps/    SAZ47 McLane PARFLUX sediment traps (~1000, 2000, 3800 m): mass, POC, PIC, PN, BSi fluxes.
                 One row per cup at the collection-interval midpoint.

QC conventions (rows failing these are dropped):
  IMOS FV01 flags (IODE): keep 1 (good) and 2 (probably good); drop 0 (no QC), 3, 4, 6 (not deployed), 7, 8, 9.
  Exceptions: SeaPHOx pHt (2021-2025) carries only flag 0 in every file, so flag 0 + range check is kept (qc '0');
  PRES (used only to set depth_m) is mostly flag 0, so 0/1/2 are accepted for depth, else nominal depth.
  Variables without a QC variable (ancillary T/S copies in FLNTUS files) are not read.
  Sensors with negative nominal depth (on the buoy, in air) are skipped.
  Daily means use only kept samples; `qc` of a daily mean is the worst kept flag of that day.
  Files whose anchor lies > MAX_KM from the model column are skipped (SAZ45/51/54/61 are never read).
  Gross plausibility ranges (RANGES) are applied after the flag filter.
"""
import glob
import math
import os
import re
import site
import sys

sys.path.append(site.getusersitepackages())  # numpy/pandas live in the user site (dropped by -I)

import numpy as np
import pandas as pd

SOTS_LAT, SOTS_LON = -46.94, 142.28
MAX_KM = 150.0  # Pulse-6 (2009-10) sat at 46.32S 140.68E, ~140 km away; all other deployments < 60 km

COLS = ["time", "depth_m", "variable", "value", "units", "source", "qc", "lat", "lon"]
GOOD = (1, 2)

RANGES = {"temp": (-2, 30), "salt": (30, 40), "O2": (50, 450), "Chl": (-0.5, 30), "NO3": (-1, 50),
          "pH": (7.3, 8.5), "PAR": (-5, 3000), "MLD": (0, 1000), "DIC": (1800, 2500),
          "ALK": (2100, 2500), "PO4": (0, 5), "SiO2": (0, 200)}

# mooring variable -> (standard name, units override or None). First match per standard name wins.
MOOR_VARS = [
    ("TEMP", "temp", None), ("PSAL", "salt", "1"),  # some CST files label PSAL 'dbar'; it is PSS-78 salinity
    ("DOX2", "O2", None), ("CPHL", "Chl", None), ("PAR", "PAR", None),
    ("NTRI", "NO3", None), ("NTRI_RAW", "NO3", None), ("MLD", "MLD", None),
    ("pHt", "pH", None),
]
SITES_OK = {"SOFS", "Pulse", "FluxPulse", "SAZ47"}

# Instrument files carry ancillary T/S/P copies from a co-located CTD (often without QC): take only the
# instrument's own product (first match on the file name). Others (CTDs, T loggers, SeaPHOx w/ SBE37-ODO): all.
OWN_PRODUCT = [("Mix-Layer-Depth", {"MLD"}), ("FLNTUS", {"Chl"}), ("PARS", {"PAR"}), ("ISUS", {"NO3"}),
               ("Optode", {"O2"}), ("GTD", set())]

RAS_VARS = {"TCO2": "DIC", "ALKA_CONC": "ALK", "NTRI_CONC": "NO3", "PHOS_CONC": "PO4", "SLCA_CONC": "SiO2",
            "TALK": "ALK", "NTRI": "NO3", "PHOS": "PO4", "SLCA": "SiO2",  # 2019+ KPSRT files
            "TEMP": "temp", "PSAL": "salt"}
TRAP_VARS = {"mass_flux": "mass_flux", "POC_mass_flux": "POC_flux", "PIC_mass_flux": "PIC_flux",
             "PN_mass_flux": "PON_flux", "BSi_mass_flux": "bSi_flux"}  # PC_mass_flux (total C = POC + PIC) not kept


def _empty():
    return pd.DataFrame(columns=COLS)


def _dist_km(lat, lon, lat0=SOTS_LAT, lon0=SOTS_LON):
    la, lo, la0, lo0 = map(math.radians, (lat, lon, lat0, lon0))
    a = math.sin((la - la0) / 2) ** 2 + math.cos(la) * math.cos(la0) * math.sin((lo - lo0) / 2) ** 2
    return 6371.0 * 2 * math.asin(math.sqrt(a))


def _f(x):
    return np.ma.filled(np.ma.masked_invalid(np.ma.asarray(x, dtype=float)), np.nan).ravel()


def _times(d):
    import netCDF4
    tv = d.variables["TIME"]
    t = netCDF4.num2date(_f(tv[:]), tv.units, only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    return pd.DatetimeIndex(pd.to_datetime(list(t), utc=True))


def _scalar(d, names, default=np.nan):
    for n in names:
        if n in d.variables:
            v = _f(d.variables[n][:])
            if v.size and np.isfinite(v[0]):
                return float(v[0])
    return default


def _position(d):
    lat = _scalar(d, ["LATITUDE"], getattr(d, "geospatial_lat_min", np.nan))
    lon = _scalar(d, ["LONGITUDE"], getattr(d, "geospatial_lon_min", np.nan))
    return lat, lon


def _nominal(d, var):
    for n in ("DEPTH_" + var, "NOMINAL_DEPTH", "DEPTH"):
        if n in d.variables and d.variables[n].ndim <= 1:
            v = _f(d.variables[n][:])
            if v.size and np.isfinite(v[0]):
                return float(v[0])
    v = getattr(d, "instrument_nominal_depth", None)
    if v is not None:
        return float(v)
    m = re.search(r"-(\d+)m_END", os.path.basename(d.filepath()))
    return float(m.group(1)) if m else np.nan


def _flags(d, var, n):
    q = var + "_quality_control"
    if q in d.variables:
        return _f(d.variables[q][:])
    return np.full(n, np.nan)


def _file_ok(path):
    m = re.match(r"IMOS_DWM-SOTS_[A-Z]+_\d{8}\S*?_([A-Za-z0-9.]+)_FV01_", os.path.basename(path))
    return bool(m) and m.group(1) in SITES_OK


def _latest(files):
    """Several files differ only by creation date (re-processing): keep the newest."""
    best = {}
    for f in files:
        k = re.sub(r"_C[-_]\d{8}.*$", "", os.path.basename(f))
        if k not in best or os.path.basename(f) > os.path.basename(best[k]):
            best[k] = f
    return sorted(best.values())


# ---------------------------------------------------------------- moored sensors (daily means)
def _read_mooring(path, src="imos_mooring"):
    import netCDF4
    out = []
    with netCDF4.Dataset(path) as d:
        if "TIME" not in d.variables or d.variables["TIME"].shape[0] == 0:
            return out
        lat, lon = _position(d)
        if not (np.isfinite(lat) and _dist_km(lat, lon) <= MAX_KM):
            return out
        t = _times(d)
        n = len(t)
        name = os.path.basename(path)
        allowed = next((a for k, a in OWN_PRODUCT if k in name), None)
        # measured sensor depth from good pressure, else nominal depth
        pres = None
        if "PRES" in d.variables and d.variables["PRES"].shape == (n,):
            p = _f(d.variables["PRES"][:])
            pq = _flags(d, "PRES", n)
            ok = np.isfinite(p) & (np.isin(pq, (0,) + GOOD) | np.isnan(pq)) & (p > 0.3) & (p < 6000)
            pres = np.where(ok, p, np.nan)  # PRES is mostly QC 0 (not checked); used for depth only
        done = set()
        for vname, std, uo in MOOR_VARS:
            if std in done or vname not in d.variables:
                continue
            if (allowed is not None and std not in allowed) or (allowed is None and std == "MLD"):
                continue
            var = d.variables[vname]
            if var.dimensions != ("TIME",):
                continue
            v = _f(var[:])
            q = _flags(d, vname, n)
            keep = (0,) + GOOD if std == "pH" else GOOD  # SeaPHOx pHt carries only QC 0 (not checked)
            good = np.isfinite(v) & np.isin(q, keep)
            lo, hi = RANGES.get(std, (-np.inf, np.inf))
            good &= (v >= lo) & (v <= hi)
            if not good.any():
                continue
            done.add(std)
            units = uo or getattr(var, "units", "")
            if std == "pH":
                units = "pH total scale (SeaPHOx, in situ T)"
            if std == "MLD":
                units = units + " (0.3 degC temperature anomaly from mooring T chain)"
            nom = 0.0 if std == "MLD" else _nominal(d, vname)
            if nom < 0:  # sensor on the buoy above the waterline (e.g. SBE37-3141 2018, PSAL ~0): in air
                continue
            dep = np.full(n, nom) if pres is None else np.where(np.isfinite(pres), pres, nom)
            df = pd.DataFrame({"day": t[good].floor("D"), "value": v[good], "depth": dep[good], "q": q[good]})
            g = df.groupby("day").agg(value=("value", "mean"), depth=("depth", "mean"), q=("q", "max"))
            out.append(pd.DataFrame({
                "time": (g.index + pd.Timedelta(hours=12)).strftime("%Y-%m-%dT%H:%M:%SZ"),
                "depth_m": np.round(g["depth"].values, 1), "variable": std, "value": g["value"].values,
                "units": units, "source": src, "qc": g["q"].astype(int).astype(str).values,
                "lat": round(lat, 4), "lon": round(lon, 4)}))
    return out


# ---------------------------------------------------------------- RAS water samples
def _read_ras(path, src="imos_ras"):
    import netCDF4
    out = []
    with netCDF4.Dataset(path) as d:
        lat, lon = _position(d)
        if not (np.isfinite(lat) and _dist_km(lat, lon) <= MAX_KM):
            return out
        t = _times(d)
        n = len(t)
        times = np.asarray(t.round("s").strftime("%Y-%m-%dT%H:%M:%SZ"))
        nom = _nominal(d, "")
        dep = np.full(n, nom)
        if "PRES" in d.variables and d.variables["PRES"].shape == (n,):
            p = _f(d.variables["PRES"][:])
            ok = np.isfinite(p) & np.isin(_flags(d, "PRES", n), GOOD) & (p > 0.3)
            dep = np.where(ok, p, nom)
        for vname, std in RAS_VARS.items():
            if vname not in d.variables:
                continue
            var = d.variables[vname]
            v = _f(var[:])
            q = _flags(d, vname, n)
            good = np.isfinite(v) & np.isin(q, GOOD)
            lo, hi = RANGES.get(std, (-np.inf, np.inf))
            good &= (v >= lo) & (v <= hi)
            if not good.any():
                continue
            out.append(pd.DataFrame({"time": times[good], "depth_m": np.round(dep[good], 1), "variable": std,
                                     "value": v[good], "units": getattr(var, "units", ""), "source": src,
                                     "qc": q[good].astype(int).astype(str), "lat": round(lat, 4),
                                     "lon": round(lon, 4)}))
    return out


# ---------------------------------------------------------------- sediment traps
def _read_trap(path, src="imos_traps"):
    import netCDF4
    out = []
    with netCDF4.Dataset(path) as d:
        lat, lon = _position(d)
        if not (np.isfinite(lat) and _dist_km(lat, lon) <= MAX_KM):
            return out
        t = _times(d)  # TIME = sample (cup) midpoint
        n = len(t)
        times = np.asarray(t.round("s").strftime("%Y-%m-%dT%H:%M:%SZ"))
        nom = _nominal(d, "")
        dep = np.full(n, nom)
        if "pressure_actual" in d.variables and d.variables["pressure_actual"].shape == (n,):
            p = _f(d.variables["pressure_actual"][:])
            dep = np.where(np.isfinite(p) & (p > 1), p, nom)
        sq = _flags(d, "sample", n)
        for vname, std in TRAP_VARS.items():
            if vname not in d.variables:
                continue
            var = d.variables[vname]
            v = _f(var[:])
            q = _flags(d, vname, n)
            good = np.isfinite(v) & np.isin(q, GOOD) & ~np.isin(sq, (3, 4))
            if not good.any():
                continue
            out.append(pd.DataFrame({"time": times[good], "depth_m": np.round(dep[good], 0), "variable": std,
                                     "value": v[good], "units": getattr(var, "units", ""), "source": src,
                                     "qc": q[good].astype(int).astype(str), "lat": round(lat, 4),
                                     "lon": round(lon, 4)}))
    return out


# ---------------------------------------------------------------- main entry
def load(site, site_dir):
    if site != "SOTS":
        return _empty()
    raw = os.path.join(site_dir, "raw")
    parts = []
    for sub, fn in (("imos_mooring", _read_mooring), ("imos_ras", _read_ras), ("imos_traps", _read_trap)):
        files = [f for f in glob.glob(os.path.join(raw, sub, "IMOS_DWM-SOTS_*.nc")) if _file_ok(f)]
        for f in _latest(files):
            try:
                parts += fn(f)
            except Exception as e:  # corrupt/odd file: report and continue
                print(f"[SOTS] {sub}/{os.path.basename(f)} skipped ({e})")
    parts = [p for p in parts if len(p)]
    if not parts:
        return _empty()
    return pd.concat(parts, ignore_index=True)[COLS]
