"""Reader for DYFAMED (Ligurian Sea, MOOSE/EMSO-Ligure) raw downloads.

load(site, site_dir) -> tidy DataFrame (see standardize_obs.COLUMNS).
Model column: 43.4167N 7.8667E. Details, URLs, licenses: DYFAMED/SOURCES_DYFAMED.md.

Sources (all under <site_dir>/raw/):
  seanoe/  SEANOE "Dyfamed observatory data" (doi:10.17882/43749), open files only:
    89000.mat         monthly bottle file 1994-2019 (MATLAB table object; read via
                      csv_matlab/89000_DYFAMED_BTLE.csv made by scripts/dyfamed_mat2csv.m)
                      -> tag seanoe_bottle
    105944.txt        monthly AT/CT (+ derived carbonate system) 1998-2019 -> seanoe_carb
                      (ALK/DIC only, and only > 4 d from any good bottle-file ALK/DIC)
    58203.mat.mat     CTD T/S/O2 profiles 1994-2014, O2 and S adjusted to Winkler/post-cal
                      (Coppola et al. 2018) -> seanoe_tso2
    134464.tgz -> x134464/  1988-2025 NetCDF set:
      dyfamed-vessel/*.ctd      MEDATLAS CTD casts 1990-2014 -> seanoe_vessel_ctd
                                (T/S/O2 only for casts > 4 d from a seanoe_tso2 cast; FLU2 fluorescence
                                 chl from all casts)
      dyfamed-mooring/*TSOF.nc  OceanSITES CTD casts 1998-2007 -> seanoe_tsof (only casts > 4 d
                                from seanoe_tso2 / seanoe_vessel_ctd casts)
      dyfamed-mooring/*TSCTD.nc deep mooring CTD/O2 2009-2026 -> seanoe_mooring (daily means)
      dyfamed-mooring/*TScurrents.nc  current-meter T 2009-2016 -> seanoe_mooring_rcm (daily)
      dyfamed-mooring/*FCO2TW.nc     CARIOCA buoy fCO2/SST 1995-1997 -> seanoe_carioca (daily)
      dyfamed-buoy/*_R.nc            Meteo-France buoy SST 2007-2016 -> seanoe_buoy (daily)
      dyfamed-mooring/DYFAMED_trap_fluxes_1988_2017.xls  (via csv_matlab) -> seanoe_traps
  obsvlfr_sodyf/  obs-vlfr SO-DYFAMED (JGOFS-France/PROOF) CD-ROM pages (csv via MATLAB):
    sodyf_pigments_marty.xls  HPLC pigments 1990-2001 -> Chl_HPLC + HPLC_<pig>
    sodyf_primprodC_93_99_chiaverini.xls + pp-2000.htm, pp-2001.htm -> PP 1993-2001
    sodyf_POC_PON_97_01_marty.xls -> POC, PON 1997-2001
    sodyf_nut_nutriments.xls  nutrients 1991-2002 (only > 4 d from seanoe_bottle nutrients)

QC conventions (rows failing these are dropped):
  seanoe_bottle  own scheme (QCcodes: 2 probably good, 3 probably bad, 4 bad, 9 missing; no 1);
                 keep 2 only (nitrite flag 0, undefined, dropped). CTD temp has no flag.
  OceanSITES NetCDF (mooring, rcm, carioca, buoy, tsof): keep QC 1, 2; drop 0, 3, 4, 7, 8, 9
                 and fills, plus gross range checks (RANGES) and, for moored sensors, samples
                 whose pressure is < 50% of the nominal depth (deployment/recovery at surface).
                 Daily means: qc = worst kept flag of the day.
  MEDATLAS .ctd  per-value QC character: keep 1, 2; drop default (missing) values.
  seanoe_tso2, seanoe_carb, obs-vlfr files: no flags (qc ''); -999 / blanks dropped.
  Sediment-trap NetCDF (OS_DYFAMED_201[0-3]_D_SedimentTrap.nc) NOT used: their POCF equals
  the xls nitrogen flux (mislabeled); the xls covers the same deployments.
"""
import glob
import math
import os
import re
import site
import sys
import warnings

sys.path.append(site.getusersitepackages())  # numpy/pandas live in the user site (dropped by -I)

import numpy as np
import pandas as pd

DYF_LAT, DYF_LON = 43.4167, 7.8667
MAX_KM = 60.0  # casts kept within this distance (some headers carry 7deg15' instead of 7deg51')

COLS = ["time", "depth_m", "variable", "value", "units", "source", "qc", "lat", "lon"]

RANGES = {"temp": (10, 30), "salt": (36, 40), "O2": (100, 320), "O2_CTD": (100, 320),
          "fCO2": (150, 700), "Chl": (-0.1, 20)}


def _empty():
    return pd.DataFrame(columns=COLS)


def _iso(t):
    t = pd.to_datetime(t, utc=True, errors="coerce")
    if isinstance(t, pd.Series):
        return t.dt.round("s").dt.strftime("%Y-%m-%dT%H:%M:%SZ")
    return np.asarray(pd.DatetimeIndex(t).round("s").strftime("%Y-%m-%dT%H:%M:%SZ"))


def _dist_km(lat, lon, lat0=DYF_LAT, lon0=DYF_LON):
    lat, lon = np.radians(np.asarray(lat, float)), np.radians(np.asarray(lon, float))
    la0, lo0 = math.radians(lat0), math.radians(lon0)
    a = np.sin((lat - la0) / 2) ** 2 + np.cos(lat) * math.cos(la0) * np.sin((lon - lo0) / 2) ** 2
    return 6371.0 * 2 * np.arcsin(np.sqrt(a))


def _frame(t, z, var, v, units, src, qc="", lat=DYF_LAT, lon=DYF_LON):
    return pd.DataFrame({"time": t, "depth_m": z, "variable": var, "value": v, "units": units,
                         "source": src, "qc": qc, "lat": lat, "lon": lon})


def _rng(var, v):
    lo, hi = RANGES.get(var, (-np.inf, np.inf))
    return (v >= lo) & (v <= hi)


# ---------------------------------------------------------------- SEANOE bottle (89000)
BTL_SPEC = [  # column, flag column, variable, units
    ("CTDtemp", None, "temp", "degC (in-situ, CTD at bottle)"),
    ("CTDsal", "CTDsal_q", "salt", "PSS-78 (CTD at bottle)"),
    ("oxygen", "oxygen_q", "O2", "umol/kg (Winkler)"),
    ("nitrat", "nitrat_q", "NO3", "umol/kg"),
    ("nitrit", "nitrit_q", "NO2", "umol/kg"),
    ("phosphat", "phosphat_q", "PO4", "umol/kg"),
    ("silcat", "silcat_q", "SiO2", "umol/kg"),
    ("alkali", "alkali_q", "ALK", "umol/kg"),
    ("tcarbn", "tcarbn_q", "DIC", "umol/kg"),
    ("phts25p0", "phts25p0_q", "pH", "total scale, 25 degC, 0 dbar (measured)"),
]


def _read_bottle(raw):
    f = os.path.join(raw, "seanoe", "csv_matlab", "89000_DYFAMED_BTLE.csv")
    if not os.path.exists(f):
        return []
    x = pd.read_csv(f)
    t = _iso(pd.to_datetime(x["date_s"].astype(int).astype(str), format="%Y%m%d", utc=True))
    z = pd.to_numeric(x["CTDpres"], errors="coerce")
    lat = pd.to_numeric(x["latitude"], errors="coerce")
    lon = pd.to_numeric(x["longitude"], errors="coerce")
    out = []
    for col, qcol, var, units in BTL_SPEC:
        v = pd.to_numeric(x[col], errors="coerce")
        ok = v.notna() & z.notna()
        if qcol:
            q = pd.to_numeric(x[qcol], errors="coerce")
            ok &= q == 2
            qs = "2"
        else:
            qs = ""
        out.append(_frame(t[ok], z[ok], var, v[ok], units, "seanoe_bottle", qs, lat[ok], lon[ok]))
    return out


# ---------------------------------------------------------------- SEANOE AT/CT (105944)
def _read_carb(raw, btl):
    f = os.path.join(raw, "seanoe", "105944.txt")
    if not os.path.exists(f):
        return []
    from io import StringIO
    with open(f, encoding="latin-1") as fh:
        txt = fh.read().replace("\r\n", "\n").replace("\r", "\n")
    txt = "\n".join(ln for ln in txt.split("\n") if ln.strip())
    x = pd.read_csv(StringIO(txt), sep="\t", na_values=["NA"], dtype=str)
    x = x.loc[:, ~x.columns.duplicated()]
    t = pd.to_datetime(x["Date"], format="%d/%m/%Y", errors="coerce", utc=True)
    z = pd.to_numeric(x["PRESS"], errors="coerce")
    out = []
    for col, var in (("ALK", "ALK"), ("TIC", "DIC")):
        v = pd.to_numeric(x[col], errors="coerce")
        ok = v.notna() & t.notna() & z.notna()
        have = set()
        if btl is not None and len(btl):
            have = set(btl.loc[btl["variable"] == var, "time"].str[:10])
        ok &= ~t.dt.strftime("%Y-%m-%d").map(lambda dd: _near(dd, have) if isinstance(dd, str) else True).astype(bool)
        out.append(_frame(_iso(t[ok]), z[ok], var, v[ok], "umol/kg (units not stated in file)",
                          "seanoe_carb"))
    return out


# ---------------------------------------------------------------- SEANOE TSO2 (58203, v7.3)
def _read_tso2(raw):
    f = os.path.join(raw, "seanoe", "58203.mat.mat")
    if not os.path.exists(f):
        return []
    import h5py
    out = []
    spec = [(2, "temp", "degC (CTD; in-situ assumed, not stated in file)"), (3, "salt", "PSS-78 (CTD, adjusted)"),
            (4, "O2_CTD", "umol/kg (CTD sensor, adjusted to Winkler)")]
    with h5py.File(f, "r") as h:
        data = h["DYFAMED"]["Data"]
        for j in range(data.shape[1]):
            g = h[data[0, j]]
            if not hasattr(g, "keys") or "TSO2" not in g:
                continue
            prof = g["TSO2"]
            for i in range(prof.shape[0]):
                a = np.asarray(h[prof[i, 0]], dtype=float)
                if a.ndim != 2 or a.shape[0] != 5 or a.shape[1] == 0:
                    continue
                dn = a[0]
                t = pd.to_datetime(dn - 719529.0, unit="D", utc=True)
                ts = _iso(t)
                z = a[1]
                for row, var, units in spec:
                    v = a[row]
                    ok = np.isfinite(v) & np.isfinite(z) & np.isfinite(dn) & _rng(var, v)
                    if ok.any():
                        out.append(_frame(ts[ok], z[ok], var, v[ok], units, "seanoe_tso2"))
    return out


# ---------------------------------------------------------------- MEDATLAS vessel CTD
def _ddm(s):
    m = re.match(r"([NSEW])(\d+)\s+([\d.]+)", s.strip())
    if not m:
        return np.nan
    v = float(m.group(2)) + float(m.group(3)) / 60.0
    return -v if m.group(1) in "SW" else v


def _parse_medatlas(path):
    """Yield (time, lat, lon, DataFrame of values, DataFrame of qc chars, defaults)."""
    with open(path, encoding="latin-1") as fh:
        lines = fh.read().replace("\r\n", "\n").replace("\r", "\n").split("\n")
    i, n = 0, len(lines)
    while i < n:
        ln = lines[i]
        if ln.startswith("*DATE="):
            m = re.match(r"\*DATE=(\d{8}) TIME=(\d{4}) LAT=(\S+\s+\S+) LON=(\S+\s+\S+)", ln)
            if not m:
                i += 1
                continue
            t = pd.to_datetime(m.group(1) + m.group(2), format="%d%m%Y%H%M", utc=True, errors="coerce")
            lat, lon = _ddm(m.group(3)), _ddm(m.group(4))
            defs = {}
            i += 1
            while i < n and not re.match(r"\*[A-Z0-9]{4}( {2,}[A-Z0-9]{4})+\s*$", lines[i]):
                pm = re.match(r"\*([A-Z0-9]{4}) .*def\.=\s*(\S+)", lines[i])
                if pm:
                    defs[pm.group(1)] = float(pm.group(2))
                if lines[i].startswith("*DATE=") or lines[i].startswith("*FI"):
                    break
                i += 1
            if i >= n or not re.match(r"\*[A-Z0-9]{4}( {2,}[A-Z0-9]{4})+\s*$", lines[i]):
                continue
            names = lines[i][1:].split()
            i += 1
            rows = []
            while i < n and lines[i] and not lines[i].startswith("*"):
                p = lines[i].split()
                if len(p) == len(names) + 1 and len(p[-1]) == len(names):
                    rows.append(p)
                i += 1
            if rows:
                vals = pd.DataFrame([r[:-1] for r in rows], columns=names).apply(pd.to_numeric, errors="coerce")
                qcs = pd.DataFrame([list(r[-1]) for r in rows], columns=names)
                yield t, lat, lon, vals, qcs, defs
            continue
        i += 1


DEDUP_DAYS = 4  # a cast within this many days of a higher-priority cast is the same cruise visit


def _near(day, dates):
    if not dates:
        return False
    d = pd.Timestamp(day)
    return any(abs((pd.Timestamp(x) - d).days) <= DEDUP_DAYS for x in dates)


def _read_vessel(raw, skip_dates):
    files = sorted(glob.glob(os.path.join(raw, "seanoe", "x134464", "dyfamed-vessel", "*.ctd")))
    out, dates = [], set()
    for f in files:
        for t, lat, lon, vals, qcs, defs in _parse_medatlas(f):
            if pd.isna(t) or "PRES" not in vals or not (_dist_km(lat, lon) <= MAX_KM):
                continue
            day = t.strftime("%Y-%m-%d")
            z = vals["PRES"]
            zok = z.notna() & qcs["PRES"].isin(["1", "2"]) & (z != defs.get("PRES", -999.9))
            ts = t.strftime("%Y-%m-%dT%H:%M:%SZ")
            spec = []
            if not _near(day, skip_dates):
                dates.add(day)
                spec += [("TEMP", "temp", "degC (in-situ, CTD)"), ("PSAL", "salt", "PSS-78 (CTD)")]
                ox = next((c for c in ("DOX2", "DOXY", "DOX1") if c in vals), None)
                if ox:
                    spec.append((ox, "O2_CTD", {"DOX2": "umol/kg (CTD sensor)",
                                                "DOXY": "mmol/m3 (CTD sensor)",
                                                "DOX1": "mL/L (CTD sensor)"}[ox]))
            if "FLU2" in vals:
                spec.append(("FLU2", "Chl", "mg/m3 (CTD fluorescence, FLU2)"))
            for col, var, units in spec:
                if col not in vals:
                    continue
                v = vals[col]
                ok = zok & v.notna() & qcs[col].isin(["1", "2"]) & (v != defs.get(col, 99.999))
                ok &= _rng(var, v)
                if ok.any():
                    out.append(_frame(ts, z[ok].values, var, v[ok].values, units, "seanoe_vessel_ctd",
                                      qcs[col][ok].values, round(lat, 4), round(lon, 4)))
    return out, dates


# ---------------------------------------------------------------- OceanSITES helpers
def _nc_time(d):
    import netCDF4
    tv = d.variables["TIME"]
    tt = netCDF4.num2date(np.ma.filled(np.ma.asarray(tv[:], dtype=float), np.nan), tv.units,
                          only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    return pd.DatetimeIndex(pd.to_datetime(list(tt), utc=True))


def _a2(d, name, nt):
    a = np.ma.filled(np.ma.masked_invalid(np.ma.asarray(d.variables[name][:], dtype=float)), np.nan)
    return a.reshape(nt, -1)


def _good(d, name, nt):
    if name + "_QC" not in d.variables:
        return np.zeros((nt, 1), bool) | True, np.full((nt, 1), np.nan)
    q = _a2(d, name + "_QC", nt)
    return np.isin(q, [1, 2]), q


def _daily(df, keys):
    """Daily mean of value/depth per keys; qc = worst kept flag."""
    if df.empty:
        return df
    df = df.copy()
    df["day"] = df["t"].dt.floor("D")
    g = df.groupby(keys + ["day"], dropna=False)
    r = g.agg(value=("value", "mean"), depth_m=("depth_m", "mean"), qc=("qcn", "max"),
              lat=("lat", "first"), lon=("lon", "first")).reset_index()
    r["time"] = r["day"].dt.strftime("%Y-%m-%dT12:00:00Z")
    r["qc"] = r["qc"].map(lambda q: "" if pd.isna(q) else str(int(q)))
    r["depth_m"] = r["depth_m"].round(1)
    return r


def _read_moored(path, spec, src, use_pres=True):
    import netCDF4
    rows = []
    with netCDF4.Dataset(path) as d:
        nt = d.dimensions["TIME"].size
        if nt == 0:
            return rows
        t = _nc_time(d)
        if "DEPTH" in d.variables and d.variables["DEPTH"].dimensions == ("DEPTH",):
            nom = np.ma.filled(np.ma.asarray(d.variables["DEPTH"][:], dtype=float), np.nan)
        else:
            nom = np.nanmean(_a2(d, "DEPH", nt), axis=0)
        lat = float(np.nanmean(np.ma.filled(np.ma.asarray(d.variables["LATITUDE"][:], dtype=float), np.nan)))
        lon = float(np.nanmean(np.ma.filled(np.ma.asarray(d.variables["LONGITUDE"][:], dtype=float), np.nan)))
        depth = np.tile(nom, (nt, 1))
        pbad = np.zeros_like(depth, bool)
        if use_pres and "PRES" in d.variables:
            p = _a2(d, "PRES", nt)
            if p.shape == depth.shape:
                pg, _ = _good(d, "PRES", nt)
                pok = np.isfinite(p) & pg
                depth = np.where(pok, p, depth)
                pbad = pok & (p < 0.5 * np.tile(nom, (nt, 1)))  # sensor at surface
        for name, var, units in spec:
            if name not in d.variables:
                continue
            v = _a2(d, name, nt)
            if v.shape != depth.shape:
                continue
            g, q = _good(d, name, nt)
            ok = g & np.isfinite(v) & _rng(var, v) & ~pbad & np.isfinite(depth)
            ii, jj = np.nonzero(ok)
            if len(ii) == 0:
                continue
            rows.append(pd.DataFrame({"t": t[ii], "depth_m": depth[ii, jj], "nom": nom[jj],
                                      "variable": var, "value": v[ii, jj], "units": units,
                                      "qcn": q[ii, jj] if q.shape == v.shape else np.nan,
                                      "lat": round(lat, 4), "lon": round(lon, 4)}))
    return rows


def _moored_daily(raw, pattern, spec, src, use_pres=True):
    parts = []
    for f in sorted(glob.glob(os.path.join(raw, "seanoe", "x134464", pattern))):
        parts += _read_moored(f, spec, src, use_pres)
    if not parts:
        return []
    df = pd.concat(parts, ignore_index=True)
    r = _daily(df, ["variable", "units", "nom"])
    r["source"] = src
    return [r[COLS]]


# ---------------------------------------------------------------- OceanSITES CTD casts (TSOF)
def _read_tsof(raw, skip_dates):
    import netCDF4
    out = []
    for f in sorted(glob.glob(os.path.join(raw, "seanoe", "x134464", "dyfamed-mooring", "*_TSOF.nc"))):
        with netCDF4.Dataset(f) as d:
            nt = d.dimensions["TIME"].size
            t = _nc_time(d)
            lat = np.ma.filled(np.ma.asarray(d.variables["LATITUDE"][:], dtype=float), np.nan)
            lon = np.ma.filled(np.ma.asarray(d.variables["LONGITUDE"][:], dtype=float), np.nan)
            p = _a2(d, "PRES", nt)
            pg, _ = _good(d, "PRES", nt)
            for k in range(nt):
                day = t[k].strftime("%Y-%m-%d")
                if _near(day, skip_dates) or not (_dist_km(lat[k], lon[k]) <= MAX_KM):
                    continue
                for name, var, units in (("TEMP", "temp", "degC (in-situ, CTD)"), ("PSAL", "salt", "PSS-78 (CTD)")):
                    v = _a2(d, name, nt)[k]
                    g, q = _good(d, name, nt)
                    ok = g[k] & pg[k] & np.isfinite(v) & np.isfinite(p[k]) & _rng(var, v)
                    if ok.any():
                        out.append(_frame(t[k].strftime("%Y-%m-%dT%H:%M:%SZ"), p[k][ok], var, v[ok], units,
                                          "seanoe_tsof", [str(int(x)) for x in q[k][ok]],
                                          round(float(lat[k]), 4), round(float(lon[k]), 4)))
    return out


# ---------------------------------------------------------------- sediment traps (xls)
def _read_traps(raw):
    f = os.path.join(raw, "seanoe", "csv_matlab", "DYFAMED_trap_fluxes_1988_2017__01_Feuil1.csv")
    if not os.path.exists(f):
        return []
    x = pd.read_csv(f, header=None, skiprows=7, dtype=str)
    x = x[x[0].fillna("").str.startswith("DYF")]
    t0 = pd.to_datetime(x[3], errors="coerce", utc=True)
    t1 = pd.to_datetime(x[4], errors="coerce", utc=True)
    tm = _iso(t0 + (t1 - t0) / 2)
    out = []
    for depth, cols in ((200.0, (7, 8, 10)), (1000.0, (12, 13, 15))):
        for c, var, units in zip(cols, ("mass_flux", "PON_flux", "POC_flux"),
                                 ("mg/m2/d", "mg N/m2/d (total N)", "mg C/m2/d (organic C)")):
            v = pd.to_numeric(x[c], errors="coerce")
            ok = v.notna() & t0.notna() & t1.notna()
            out.append(_frame(tm[ok.values], depth, var, v[ok].values, units, "seanoe_traps"))
    return out


# ---------------------------------------------------------------- obs-vlfr SO-DYFAMED
PIG = {"chl c3": "chl_c3", "chl c12": "chl_c1c2", "peri": "perid", "19'BF": "but_fucox",
       "fuco": "fucox", "9'-cis-neo": "neo", "19'HF": "hex_fucox", "prasi": "prasino",
       "diadi": "diadino", "allo": "allo", "diato": "diato", "ze/lu": "zea_lut", "lut": "lut",
       "dvchlb": "dv_chl_b", "chl b": "chl_b", "achl a": "chlide_a", "dv-chl a": "dv_chl_a",
       "chl a": "mv_chl_a", "a+b-car": "ab_carot", "Tchlb": "tot_chl_b"}


def _read_pigments(raw):
    f = os.path.join(raw, "obsvlfr_sodyf", "csv_matlab", "sodyf_pigments_marty__01_Feuil1.csv")
    if not os.path.exists(f):
        return []
    x = pd.read_csv(f, dtype=str)
    t = _iso(pd.to_datetime(x["date"], errors="coerce", utc=True))
    z = pd.to_numeric(x["Profondeur"], errors="coerce")
    out = []
    for col, var in [("Tchla", "Chl_HPLC")] + [(k, "HPLC_" + v) for k, v in PIG.items()]:
        v = pd.to_numeric(x[col], errors="coerce")
        ok = v.notna() & (v > -900) & z.notna() & pd.notna(t)
        units = "ng/L (HPLC total chl a = MV + DV chl a)" if var == "Chl_HPLC" else "ng/L"
        out.append(_frame(t[ok.values], z[ok], var, v[ok], units, "obsvlfr_sodyf"))
    return out


PP_UNITS = ("mgC/m3 (14C in-situ 4-h midday incubation, light bottles, labelled "
            "'C assimile corrige Dtot': probably scaled to daily total, i.e. per day; see SOURCES)")


def _read_pp(raw):
    out = []
    f = os.path.join(raw, "obsvlfr_sodyf", "csv_matlab", "sodyf_primprodC_93_99_chiaverini__01_Feuil1.csv")
    if os.path.exists(f):
        x = pd.read_csv(f, header=None, skiprows=3, dtype=str)
        d = pd.to_datetime(x[0], errors="coerce", utc=True).ffill()
        z = pd.to_numeric(x[4], errors="coerce")
        v = pd.to_numeric(x[5], errors="coerce")
        ok = d.notna() & z.notna() & v.notna()
        out.append(_frame(_iso(d[ok]), z[ok], "PP", v[ok], PP_UNITS, "obsvlfr_sodyf"))
    import html as _html
    for y in (2000, 2001):
        f = os.path.join(raw, "obsvlfr_sodyf", f"pp-{y}.htm")
        if not os.path.exists(f):
            continue
        with open(f, encoding="latin-1") as fh:
            t = fh.read()
        recs = []
        for r in re.findall(r"<tr[^>]*>(.*?)</tr>", t, re.S | re.I):
            c = [re.sub(r"\s+", " ", _html.unescape(re.sub(r"<[^>]+>", " ", s))).strip()
                 for s in re.findall(r"<td[^>]*>(.*?)</td>", r, re.S | re.I)]
            if len(c) >= 3 and re.match(r"\d\d/\d\d/\d\d$", c[0]):
                recs.append((c[0], c[1], c[2]))
        if recs:
            x = pd.DataFrame(recs)
            d = pd.to_datetime(x[0], format="%d/%m/%y", errors="coerce", utc=True)
            z = pd.to_numeric(x[1], errors="coerce")
            v = pd.to_numeric(x[2].str.replace(",", "."), errors="coerce")
            ok = d.notna() & z.notna() & v.notna()
            out.append(_frame(_iso(d[ok]), z[ok], "PP", v[ok], PP_UNITS, "obsvlfr_sodyf"))
    return out


def _read_pocpon(raw):
    out = []
    for f in sorted(glob.glob(os.path.join(raw, "obsvlfr_sodyf", "csv_matlab", "sodyf_POC_PON_97_01_marty__*.csv"))):
        x = pd.read_csv(f, dtype=str)
        d = pd.to_datetime(x.iloc[:, 0], errors="coerce", utc=True)
        z = pd.to_numeric(x.iloc[:, 1], errors="coerce").abs()  # CN2000 sheet stores depth as negative
        for k, var in ((2, "POC"), (3, "PON")):
            v = pd.to_numeric(x.iloc[:, k], errors="coerce")
            ok = d.notna() & z.notna() & v.notna()
            out.append(_frame(_iso(d[ok]), z[ok], var, v[ok], "ug/L", "obsvlfr_sodyf"))
    return out


def _read_vlfr_nut(raw, btl):
    f = os.path.join(raw, "obsvlfr_sodyf", "csv_matlab", "sodyf_nut_nutriments__01_Feuil1.csv")
    if not os.path.exists(f):
        return []
    x = pd.read_csv(f, header=None, skiprows=2, dtype=str)
    d = pd.to_datetime(x[0], errors="coerce", utc=True)
    z = pd.to_numeric(x[1], errors="coerce")
    have = set()
    if btl is not None and len(btl):
        have = set(btl.loc[btl["variable"].isin(["NO3", "PO4", "SiO2", "NO2"]), "time"].str[:10])
    keep = ~d.dt.strftime("%Y-%m-%d").map(lambda dd: _near(dd, have) if isinstance(dd, str) else True).astype(bool)
    out = []
    for k, var in ((2, "NO2"), (3, "NO3"), (4, "PO4"), (5, "SiO2")):
        v = pd.to_numeric(x[k], errors="coerce")
        ok = d.notna() & z.notna() & v.notna() & (v >= 0) & keep
        out.append(_frame(_iso(d[ok]), z[ok], var, v[ok], "umol/L", "obsvlfr_sodyf"))
    return out


# ---------------------------------------------------------------- main entry
def load(site, site_dir):
    if site != "DYFAMED":
        return _empty()
    warnings.filterwarnings("ignore", category=UserWarning)
    raw = os.path.join(site_dir, "raw")
    parts = []
    btl_parts = _read_bottle(raw)
    btl = pd.concat(btl_parts, ignore_index=True) if btl_parts else None
    parts += btl_parts + _read_carb(raw, btl)
    tso2 = _read_tso2(raw)
    parts += tso2
    tso2_dates = set(pd.concat(tso2)["time"].str[:10]) if tso2 else set()
    ves, ves_dates = _read_vessel(raw, tso2_dates)
    parts += ves
    parts += _read_tsof(raw, tso2_dates | ves_dates)
    parts += _moored_daily(raw, "dyfamed-mooring/OS_DYFAMED_*_D_TSCTD.nc",
                           [("TEMP", "temp", "degC (in-situ, moored CTD)"), ("PSAL", "salt", "PSS-78 (moored CTD)"),
                            ("DOX2", "O2", "umol/kg (moored sensor)")], "seanoe_mooring")
    parts += _moored_daily(raw, "dyfamed-mooring/OS_DYFAMED_*_D_TScurrents.nc",
                           [("TEMP", "temp", "degC (in-situ, current-meter thermistor)")], "seanoe_mooring_rcm")
    parts += _moored_daily(raw, "dyfamed-mooring/OS_DYFAMED_*_FCO2TW.nc",
                           [("FCO2", "fCO2", "uatm"), ("TEMP", "temp", "degC (buoy SST)")],
                           "seanoe_carioca", use_pres=False)
    parts += _moored_daily(raw, "dyfamed-buoy/OS_DYFAMED_*_R.nc",
                           [("TEMP", "temp", "degC (buoy SST, GTS real time)")], "seanoe_buoy", use_pres=False)
    parts += _read_traps(raw)
    parts += _read_pigments(raw) + _read_pp(raw) + _read_pocpon(raw) + _read_vlfr_nut(raw, btl)
    parts = [p for p in parts if len(p)]
    if not parts:
        return _empty()
    return pd.concat(parts, ignore_index=True)[COLS]
