"""Reader for PAP (Porcupine Abyssal Plain Sustained Observatory) raw downloads.

load(site, site_dir) -> tidy DataFrame (see standardize_obs.COLUMNS).

Sources (all under <site_dir>/raw/):
  oceansites_gdac_ifremer/  OceanSITES GDAC NetCDF (PAP-1 mooring 0-1000 m, PAP-2, PAP-3 deep T/S)
  oceansites_ndbc/          same file set from the NDBC OceanSITES THREDDS mirror; used only for
                            files missing from the Ifremer copy (Ifremer copy has the 2023 re-issues)
  ncei_ocads_mooring/       NCEI OCADS 0312034: PAP-SO 1 m pCO2 (2013-04..2013-12, 2015-07..2016-03)
  ncei_ocads_jc278/         NCEI OCADS 0312066: RRS James Cook JC278 underway pCO2/fCO2 (2025),
                            kept only within 100 km of PAP-SO
  pangaea_traps/            PANGAEA Lampitt et al. sediment trap fluxes (1989-2005; 1000/3000/~4700 m)
  pangaea_jc087/            PANGAEA JC087 (Jun 2013) bottle nutrients + Winkler O2 at PAP

QC conventions (rows failing these are dropped):
  OceanSITES  keep QC in {0 (no QC performed), 1 (good), 2 (probably good)}; drop 3,4,5,8,9 and fills.
              Plus gross range checks (see RANGES) because many 2014-2019 files carry QC=0 only.
  NCEI mooring pCO2: keep flag 1 (good); 5 bad / 9 missing dropped.
  JC278 underway (SOCAT/WOCE flags): keep 2 (good); drop 3, 4 and missing.
  PANGAEA: no flags; values reported as '<x' (below detection) dropped.
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

PAP_LAT, PAP_LON = 49.0, -16.5

COLS = ["time", "depth_m", "variable", "value", "units", "source", "qc", "lat", "lon"]

# gross plausibility limits applied after flag filtering
RANGES = {"temp": (-2, 35), "salt": (30, 40), "O2": (50, 450), "pCO2": (150, 700),
          "Chl": (0, 30), "NO3": (-0.5, 40), "pH": (7.5, 8.5)}

OS_VARMAP = {  # OceanSITES variable -> standard name
    "TEMP": "temp", "PSAL": "salt", "DOXY": "O2", "DOXM": "O2",
    "PCO2XXXX": "pCO2", "PCO2": "pCO2",
    "CPHLPS01": "Chl", "CPHLPM01": "Chl", "CPHL": "Chl", "CHL": "Chl", "Chl": "Chl",
    "NTRA": "NO3", "NTRZ": "NO3", "NO3": "NO3", "NITRATE": "NO3",
    "PHPH": "pH", "PH": "pH",
}
GOOD_OS = {0, 1, 2}


def _empty():
    return pd.DataFrame(columns=COLS)


def _iso(t):
    t = pd.to_datetime(t, utc=True, errors="coerce")
    if isinstance(t, pd.Series):
        return t.dt.round("s").dt.strftime("%Y-%m-%dT%H:%M:%SZ")
    return np.asarray(pd.DatetimeIndex(t).round("s").strftime("%Y-%m-%dT%H:%M:%SZ"))


def _dist_km(lat, lon, lat0=PAP_LAT, lon0=PAP_LON):
    lat, lon = np.radians(lat), np.radians(lon)
    la0, lo0 = math.radians(lat0), math.radians(lon0)
    a = np.sin((lat - la0) / 2) ** 2 + np.cos(lat) * math.cos(la0) * np.sin((lon - lo0) / 2) ** 2
    return 6371.0 * 2 * np.arcsin(np.sqrt(a))


# ---------------------------------------------------------------- OceanSITES
def _os_select_files(raw):
    """Pick one file per (platform, deployment, product); prefer Ifremer copy, D > R > P mode."""
    files = {}
    for src in ("oceansites_ndbc", "oceansites_gdac_ifremer"):  # later overrides earlier
        for f in glob.glob(os.path.join(raw, src, "OS_PAP*.nc")):
            files[os.path.basename(f)] = (src, f)
    rank = {"D": 0, "R": 1, "P": 2}
    best = {}
    for name, (src, f) in files.items():
        m = re.match(r"OS_(PAP-\d)_(\d{6})_([DRP])_(.*)\.nc$", name)
        if not m:
            continue
        plat, dep, mode, prod = m.groups()
        prod = re.sub(r"^P?_", "", prod)
        key = (plat, dep, prod)
        if key not in best or rank[mode] < rank[best[key][0]]:
            best[key] = (mode, src, f)
    return [(v[1], v[2]) for v in best.values()]


def _as2d(arr, nt):
    a = np.ma.filled(np.ma.masked_invalid(np.ma.asarray(arr, dtype=float)), np.nan)
    return a.reshape(nt, -1)


def _read_oceansites(src, path):
    import netCDF4
    out = []
    with netCDF4.Dataset(path) as d:
        tname = next((k for k in d.variables if k.upper() == "TIME"), None)
        if tname is None or "DEPTH" not in d.variables:
            return out
        tv = d.variables[tname]
        nt = tv.shape[0]
        if nt == 0:
            return out
        tt = netCDF4.num2date(np.ma.filled(np.ma.asarray(tv[:], dtype=float), np.nan), tv.units,
                              only_use_cftime_datetimes=False, only_use_python_datetimes=True)
        times = _iso(pd.DatetimeIndex(list(tt)))
        nom = np.ma.filled(np.ma.asarray(d.variables["DEPTH"][:], dtype=float), np.nan).ravel()
        lat = float(np.ma.filled(d.variables["LATITUDE"][:], np.nan).ravel()[0]) if "LATITUDE" in d.variables else np.nan
        lon = float(np.ma.filled(d.variables["LONGITUDE"][:], np.nan).ravel()[0]) if "LONGITUDE" in d.variables else np.nan
        if lon > 0:  # some early files store 16.4 (positive) for 16.4W
            lon = -lon
        if not _dist_km(np.array([lat]), np.array([lon]))[0] <= 50.0:
            lat = lon = np.nan  # bad position metadata (e.g. 200307 CTD lon -0.69); mooring is at PAP
        # measured sensor depth from pressure where usable, else nominal depth
        depth = np.tile(nom, (nt, 1))
        if "PRES" in d.variables:
            p = _as2d(d.variables["PRES"][:], nt)
            if p.shape == depth.shape:
                ok = np.isfinite(p) & (p > 0.3) & (p < 6000)
                if "PRES_QC" in d.variables:
                    pq = _as2d(d.variables["PRES_QC"][:], nt)
                    ok &= np.isin(pq, list(GOOD_OS))
                depth = np.where(ok, p, depth)
        # each product file carries ancillary CTD/O2 copies; take only its own variables
        prod = os.path.basename(path)
        allowed = ({"pCO2"} if "PCO2" in prod else {"Chl"} if "Chl" in prod else
                   {"NO3"} if "ISUS" in prod else {"O2"} if prod.endswith("_O.nc") else
                   {"temp", "salt", "O2"})
        for vname, var in d.variables.items():
            std = OS_VARMAP.get(vname)
            if std is None or std not in allowed:
                continue
            units = getattr(var, "units", "")
            if std == "O2" and "bar" in units:  # 2014 CTDO file: DOXY is a copy of PRES
                continue
            v = _as2d(var[:], nt)
            if v.shape[1] != len(nom):
                continue
            if vname + "_QC" in d.variables:
                q = _as2d(d.variables[vname + "_QC"][:], nt)
                good = np.isin(q, list(GOOD_OS))
            else:
                q = np.full(v.shape, np.nan)
                good = np.ones(v.shape, bool)
            good &= np.isfinite(v)
            if std == "pCO2":
                med = np.nanmedian(np.where(good, v, np.nan)) if good.any() else np.nan
                if units in ("microPa", "Pa") and med > 150:
                    units = "uatm"  # mislabelled in file; magnitudes are uatm (see SOURCES_PAP.md)
                elif units == "microAtmospheres":
                    units = "uatm"
                lo, hi = (15.0, 71.0) if units == "Pa" else RANGES["pCO2"]
            else:
                lo, hi = RANGES.get(std, (-np.inf, np.inf))
            good &= (v >= lo) & (v <= hi)
            ii, jj = np.nonzero(good)
            if len(ii) == 0:
                continue
            out.append(pd.DataFrame({
                "time": times[ii], "depth_m": np.round(depth[ii, jj], 1), "variable": std,
                "value": v[ii, jj], "units": units, "source": src,
                "qc": [("" if np.isnan(x) else str(int(x))) for x in q[ii, jj]],
                "lat": lat, "lon": lon, "_file": os.path.basename(path), "_nom": nom[jj]}))
    return out


# ---------------------------------------------------------------- NCEI OCADS
def _read_ncei_mooring(raw):
    d = os.path.join(raw, "ncei_ocads_mooring")
    out = []
    f = os.path.join(d, "747F20130425.csv")
    if os.path.exists(f):
        x = pd.read_csv(f, skiprows=4, encoding="latin-1")
        x.columns = [re.sub(r".*\[(.*)\]$", r"\1", c).strip() for c in x.columns]
        t = pd.to_datetime(x["date"].astype(str) + " " + x["time"].astype(str), dayfirst=True,
                           format="mixed", errors="coerce", utc=True)
        v = pd.to_numeric(x["pCO2_(uatm)"], errors="coerce")
        q = pd.to_numeric(x["quality_flag"], errors="coerce")
        ok = (q == 1) & v.between(*RANGES["pCO2"]) & t.notna()
        out.append(pd.DataFrame({"time": _iso(t[ok]), "depth_m": 1.0, "variable": "pCO2",
                                 "value": v[ok], "units": "uatm", "source": "ncei_ocads_mooring",
                                 "qc": "1", "lat": 49.0, "lon": -16.5}))
    f = os.path.join(d, "747F20150702.csv")
    if os.path.exists(f):
        x = pd.read_csv(f, skiprows=4, encoding="latin-1")
        t = pd.to_datetime(x["Date_Time"], dayfirst=True, format="%d/%m/%Y %H:%M", errors="coerce", utc=True)
        v = pd.to_numeric(x["1m_pro_CO2_conc_uatm_SN299745"], errors="coerce")
        ok = v.between(*RANGES["pCO2"]) & t.notna()
        out.append(pd.DataFrame({"time": _iso(t[ok]), "depth_m": 1.0, "variable": "pCO2",
                                 "value": v[ok], "units": "uatm", "source": "ncei_ocads_mooring",
                                 "qc": "", "lat": 49.0, "lon": -16.5}))
    return out


def _read_jc278(raw, max_km=100.0):
    f = os.path.join(raw, "ncei_ocads_jc278", "740H20250530.csv")
    if not os.path.exists(f):
        return []
    with open(f, encoding="latin-1") as fh:  # header row is the first line starting 'date time'
        skip = next(i for i, ln in enumerate(fh) if ln.lower().startswith("date time"))
    x = pd.read_csv(f, skiprows=skip, encoding="latin-1", low_memory=False)
    x.columns = [re.sub(r"\s+", " ", c.strip()) for c in x.columns]
    lat = pd.to_numeric(x["latitude [ deg N ]"], errors="coerce")
    lon = pd.to_numeric(x["longitude [ deg E ]"], errors="coerce")
    near = _dist_km(lat.values, lon.values) <= max_km
    x, lat, lon = x[near], lat[near], lon[near]
    t = _iso(pd.to_datetime(x["date time"], errors="coerce", utc=True))
    dep = pd.to_numeric(x["Depth [m]"], errors="coerce")
    spec = [("SST [ deg C ]", "temp", "degC"), ("Salinity [ PSU ]", "salt", "PSU"),
            ("pCO2_water_SST_wet [ uatm ]", "pCO2", "uatm"),
            ("fCO2_water_SST_wet [ uatm ]", "fCO2", "uatm"),
            ("xCO2_atm_dry_actual [ umol/mol ]", "pCO2_air", "umol/mol (xCO2 dry air)")]
    out = []
    for col, std, units in spec:
        v = pd.to_numeric(x[col], errors="coerce")
        q = pd.to_numeric(x[col + " QC Flag"], errors="coerce")
        ok = (q == 2) & v.notna()
        if std in RANGES:
            ok &= v.between(*RANGES[std])
        out.append(pd.DataFrame({"time": t[ok], "depth_m": dep[ok].fillna(5.0) if std != "pCO2_air" else 0.0,
                                 "variable": std, "value": v[ok], "units": units,
                                 "source": "ncei_ocads_jc278", "qc": "2", "lat": lat[ok], "lon": lon[ok]}))
    return out


# ---------------------------------------------------------------- PANGAEA
def _read_pangaea_tab(f):
    with open(f, encoding="utf-8") as fh:
        txt = fh.read()
    head, body = txt.split("*/\n", 1)
    from io import StringIO
    df = pd.read_csv(StringIO(body), sep="\t", dtype=str)
    m = re.search(r"LATITUDE(?: START)?: ([-\d.]+) \* LONGITUDE(?: START)?: ([-\d.]+)", head)
    lat, lon = (float(m.group(1)), float(m.group(2))) if m else (np.nan, np.nan)
    return head, df, lat, lon


TRAP_MAP = [  # (column regex, std var, units)
    (r"^(Tot mass flux|Flux tot) \[mg/m\*\*2/day\]$", "mass_flux", "mg/m2/d"),
    (r"^TOC flux \[mmol/m\*\*2/day\]$", "POC_flux", "mmol C/m2/d"),
    (r"^POC flux \[mg/m\*\*2/day\]$", "POC_flux", "mg C/m2/d"),
    (r"^IC flux \[mmol/m\*\*2/day\]$", "PIC_flux", "mmol C/m2/d"),
    (r"^PIC flux \[mg/m\*\*2/day\]$", "PIC_flux", "mg C/m2/d"),
    (r"^(bSiO2|PSiO2) flux \[mg/m\*\*2/day\]$", "bSi_flux", "mg SiO2(opal)/m2/d"),
    (r"^N flux \[mg/m\*\*2/day\]$", "PON_flux", "mg N/m2/d (total N)"),
    (r"^PN flux \[mg/m\*\*2/day\]$", "PON_flux", "mg N/m2/d"),
]


def _read_traps(raw):
    out = []
    for f in sorted(glob.glob(os.path.join(raw, "pangaea_traps", "PANGAEA.*.tab"))):
        head, df, lat, lon = _read_pangaea_tab(f)
        if "Lampitt" not in head.split("\n", 2)[1]:
            continue  # skip the Torres-Valdes compilation (duplicates these deployments)
        if not _dist_km(np.array([lat]), np.array([lon]))[0] <= 100.0:
            continue  # 1989-90 deployments (PAP-I/III/V) were at ~47.8N 19.5W, ~260 km away
        t0 = pd.to_datetime(df["Date/Time"], errors="coerce", utc=True)
        t1 = pd.to_datetime(df["Date/time end"], errors="coerce", utc=True)
        tm = t0 + (t1 - t0) / 2  # collection-interval midpoint
        dep = pd.to_numeric(df["Depth water [m]"], errors="coerce")
        tag = "pangaea_traps"
        for col in df.columns:
            for rx, std, units in TRAP_MAP:
                if re.match(rx, col):
                    v = pd.to_numeric(df[col], errors="coerce")
                    ok = v.notna() & tm.notna()
                    out.append(pd.DataFrame({"time": _iso(tm[ok]), "depth_m": dep[ok], "variable": std,
                                             "value": v[ok], "units": units, "source": tag, "qc": "",
                                             "lat": lat, "lon": lon}))
    return out


def _read_jc087(raw):
    f = os.path.join(raw, "pangaea_jc087", "PANGAEA.832864.tab")
    if not os.path.exists(f):
        return []
    head, df, _, _ = _read_pangaea_tab(f)
    t = _iso(pd.to_datetime(df["Date/Time"], errors="coerce", utc=True))
    dep = pd.to_numeric(df["Depth water [m]"], errors="coerce")
    lat = pd.to_numeric(df["Latitude"], errors="coerce")
    lon = pd.to_numeric(df["Longitude"], errors="coerce")
    out = []
    spec = [("Si(OH)4", "SiO2"), ("[PO4]3-", "PO4"), ("[NO3]- + [NO2]-", "NO3"),
            ("[NO2]-", "NO2"), ("[NH4]+", "NH4")]
    for key, std in spec:
        col = next(c for c in df.columns if c.startswith(key + " "))
        v = pd.to_numeric(df[col], errors="coerce")  # '<0.01' -> NaN (below detection, dropped)
        ok = v.notna()
        out.append(pd.DataFrame({"time": t[ok], "depth_m": dep[ok], "variable": std, "value": v[ok],
                                 "units": "umol/L", "source": "pangaea_jc087", "qc": "",
                                 "lat": lat[ok], "lon": lon[ok]}))
    oc = [c for c in df.columns if c.startswith("O2 ")]
    o = pd.concat([pd.to_numeric(df[c], errors="coerce") for c in oc], axis=1).mean(axis=1)
    ok = o.notna()
    out.append(pd.DataFrame({"time": t[ok], "depth_m": dep[ok], "variable": "O2", "value": o[ok],
                             "units": "umol/L", "source": "pangaea_jc087", "qc": "",
                             "lat": lat[ok], "lon": lon[ok]}))
    return out


# ---------------------------------------------------------------- main entry
def load(site, site_dir):
    if site != "PAP":
        return _empty()
    raw = os.path.join(site_dir, "raw")
    parts = []
    for src, f in sorted(_os_select_files(raw), key=lambda x: os.path.basename(x[1])):
        try:
            parts += _read_oceansites(src, f)
        except Exception as e:  # truncated/corrupt copy: fall back to the other mirror
            alt = os.path.join(raw, "oceansites_ndbc", os.path.basename(f))
            print(f"[PAP] {src}/{os.path.basename(f)} unreadable ({e}); trying NDBC copy")
            if src != "oceansites_ndbc" and os.path.exists(alt):
                try:
                    parts += _read_oceansites("oceansites_ndbc", alt)
                except Exception as e2:
                    print(f"[PAP] NDBC copy also unreadable ({e2}); skipped")
    os_df = pd.concat(parts, ignore_index=True) if parts else _empty()
    ncei = _read_ncei_mooring(raw)
    if ncei and len(os_df):
        # NCEI OCADS 1 m pCO2 is the archived carbon-QC'd version of the OceanSITES 1 m pCO2
        # for the 2013 and 2015 deployments: drop the OceanSITES duplicate in those windows.
        nc = pd.concat(ncei)
        tn = pd.to_datetime(nc["time"])
        for t0, t1 in [(tn[tn.dt.year == 2013].min(), tn[tn.dt.year == 2013].max()),
                       (tn[tn.dt.year >= 2015].min(), tn[tn.dt.year >= 2015].max())]:
            to = pd.to_datetime(os_df["time"])
            dup = (os_df["variable"] == "pCO2") & (os_df["_nom"] <= 5) & (to >= t0) & (to <= t1)
            os_df = os_df[~dup]
    parts = [os_df.drop(columns=["_file", "_nom"], errors="ignore")] + ncei
    parts += _read_jc278(raw) + _read_traps(raw) + _read_jc087(raw)
    parts = [p for p in parts if len(p)]
    if not parts:
        return _empty()
    return pd.concat(parts, ignore_index=True)[COLS]
