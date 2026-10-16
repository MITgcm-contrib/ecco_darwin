"""Reader for JAMSTEC time-series station K2 (47.0N 160.0E, western subarctic Pacific) -> tidy rows.

Sources (raw dirs under <site_dir>/raw/; see K2/SOURCES_K2.md):
  ncei_ocads_k2/        NCEI OCADS 0100115, K2 bottle data 1999-2008 (Wakita et al. 2010), WHP-style CSV
  jamstec_k2s1_bottle/  JAMSTEC K2S1 Mirai routine bottle data 2010-2012 (ver2 xlsx)
  jamstec_k2s1_ctd/     JAMSTEC K2S1 CTD: WHP .wct casts (MR10-01, MR11-02/03/05, MR12-02), ODV .o4x (MR10-06)
  jamstec_k2s1_pp/      JAMSTEC K2S1 13C primary production (open files MR10-01, MR10-06 only)
  jamstec_k2s1_trap/    JAMSTEC K2S1 moored sediment traps 2010-2014 (MSTdata_K2.xlsx; reported fluxes)
  oceansites_k2_trap/   JAMSTEC OceanSITES K2 sediment-trap NetCDF 2005-2025

Only K2 casts are used (station K02 / K2); KNOT, S1 and transit stations are dropped.

Flag conventions (rows failing are dropped):
  NCEI bottle: WOCE *_FLAG_W, keep 2 (good); 3/4/9 and -999 dropped. CTDTMP has no flag.
  JAMSTEC bottle xlsx: own QF scheme 0 accepted, 1 unknown, 4 questionable, 8 bad/missing.
      Keep 0 and 1 (1 = not assessed, as HOT/BIOS 'unchecked'; all chlorophyll is QF 1); drop 4, 8, -9.
  JAMSTEC CTD .wct: WHP QUALT1 digit per parameter, keep 2 (good). .o4x: ODV flags, keep 0 (good).
  JAMSTEC PP xlsx: WOCE-style 2 acceptable kept; 3/4/9 dropped.
  JAMSTEC trap xlsx: no flags; 999 = not available (dropped); fluxes of exactly 0 whose concentration
      is 999 (not measured) are dropped.
  OceanSITES traps: <VAR>_QC, keep 1 (good) and 2 (probably good); 4 dropped. Fluxes reported as exactly
      0 while the underlying concentration is flagged bad (sample too small to analyze) are dropped.
"""
import glob
import os
import re

import numpy as np
import pandas as pd

try:
    from standardize_obs import iso
except Exception:  # pragma: no cover
    def iso(ts):
        return pd.to_datetime(ts, utc=True, errors="coerce").dt.strftime("%Y-%m-%dT%H:%M:%SZ")

LAT0, LON0 = 47.0, 160.0
CTD_FULLRES_MAX_DBAR = 1000.0  # CTD: keep every level <= 1000 dbar, then every 10 dbar below
CTD_STRIDE_DBAR = 10.0


def _rows(time, depth, var, value, units, source, qc="", lat=np.nan, lon=np.nan):
    d = pd.DataFrame({"time": time, "depth_m": depth, "variable": var, "value": value,
                      "units": units, "source": source, "qc": qc, "lat": lat, "lon": lon})
    d["value"] = pd.to_numeric(d["value"], errors="coerce")
    return d.dropna(subset=["value"])


def p2z(p, lat=LAT0):
    """UNESCO 1983 (Fofonoff & Millard) pressure (dbar) -> depth (m)."""
    p = np.asarray(p, dtype=float)
    x = np.sin(np.deg2rad(lat)) ** 2
    g = 9.780318 * (1.0 + (5.2788e-3 + 2.36e-5 * x) * x) + 1.092e-6 * p
    return ((((-1.82e-15 * p + 2.279e-10) * p - 2.2512e-5) * p + 9.72659) * p) / g


# ---------------------------------------------------------------- NCEI OCADS bottle 1999-2008
def _ncei_bottle(raw):
    f = os.path.join(raw, "ncei_ocads_k2", "K2_bottle_data.csv")
    if not os.path.exists(f):
        return []
    d = pd.read_csv(f, skiprows=[1], low_memory=False)
    d = d[d["SECT_ID"].astype(str).str.strip() == "K2"]
    tt = pd.to_datetime(d["DATE"].astype(str) + " " + d["TIME"].astype(str),
                        format="%m/%d/%Y %H:%M", errors="coerce")
    t = iso(tt)
    p = pd.to_numeric(d["CTDPRS"], errors="coerce")
    z = pd.Series(p2z(p.values, LAT0), index=d.index)
    lat = pd.to_numeric(d["LATITUDE"], errors="coerce")
    lon = pd.to_numeric(d["LONGITUDE"], errors="coerce")
    src = "ncei_ocads_k2"
    out = []
    v = pd.to_numeric(d["CTDTMP"], errors="coerce").where(lambda x: x > -99)
    out.append(_rows(t, z, "temp", v, "degC (in-situ, CTD; ITS-90 not stated)", src, lat=lat, lon=lon))
    # salinity: bottle (SALNTY) where flag 2, else CTDSAL flag 2
    sb = pd.to_numeric(d["SALNTY"], errors="coerce")
    fb = pd.to_numeric(d["SALNTY_FLAG_W"], errors="coerce")
    sc = pd.to_numeric(d["CTDSAL"], errors="coerce")
    fc = pd.to_numeric(d["CTDSAL_FLAG_W"], errors="coerce")
    mb = (fb == 2) & (sb > 0)
    mc = (~mb) & (fc == 2) & (sc > 0)
    out.append(_rows(t[mb], z[mb], "salt", sb[mb], "PSS-78 (bottle)", src, qc="2", lat=lat[mb], lon=lon[mb]))
    out.append(_rows(t[mc], z[mc], "salt", sc[mc], "PSS-78 (CTD)", src, qc="2", lat=lat[mc], lon=lon[mc]))
    for c, var in [("OXYGEN", "O2"), ("NITRAT", "NO3"), ("NITRIT", "NO2"), ("SILCAT", "SiO2"),
                   ("PHSPHT", "PO4"), ("AMMONIUM", "NH4"), ("ALKALI", "ALK"), ("TCARBN", "DIC")]:
        v = pd.to_numeric(d[c], errors="coerce")
        fl = pd.to_numeric(d[c + "_FLAG_W"], errors="coerce")
        m = (fl == 2) & (v > -99)
        out.append(_rows(t[m], z[m], var, v[m], "umol/kg", src, qc="2", lat=lat[m], lon=lon[m]))
    return out


# ---------------------------------------------------------------- JAMSTEC bottle 2010-2012
BOT_GOOD = {0, 1}
BOT_VARS = [  # (column, variable, units)
    ("DO [~$m~#mol/kg]", "O2", "umol/kg"),
    ("Silicate [~$m~#mol/kg]", "SiO2", "umol/kg"),
    ("Nitrate [~$m~#mol/kg]", "NO3", "umol/kg"),
    ("Nitrite [~$m~#mol/kg]", "NO2", "umol/kg"),
    ("Phosphate [~$m~#mol/kg]", "PO4", "umol/kg"),
    ("Ammonia [~$m~#mol/kg]", "NH4", "umol/kg"),
    ("TA [~$m~#mol/kg]", "ALK", "umol/kg"),
    ("pHt (25C)", "pH", "total scale at 25 C"),
    ("DIC [~$m~#mol/kg]", "DIC", "umol/kg"),
    ("DOC [~$m~#mol/kg]", "DOC", "umol/kg"),
    ("CHLWEL [~$m~#g/l]", "Chl", "ug/L (fluorometric, non-acidification/Welschmeyer)"),
    ("CHLHOL [~$m~#g/l]", "Chl_acid", "ug/L (fluorometric, acidification/Holm-Hansen)"),
    ("CHLHPLC [~$m~#g/l]", "Chl_HPLC", "ug/L"),
    ("SIZECHL [~$m~#g/l]", "Chl_gt10um", "ug/L (size-fractionated, >10 um, non-acidification)"),
    ("1SIZECHL [~$m~#g/l]", "Chl_3to10um", "ug/L (size-fractionated, 3-10 um, non-acidification)"),
    ("2SIZECHL [~$m~#g/l]", "Chl_1to3um", "ug/L (size-fractionated, 1-3 um, non-acidification)"),
    ("3SIZECHL [~$m~#g/l]", "Chl_lt1um", "ug/L (size-fractionated, <1 um, non-acidification)"),
]


def _bottle_table(raw):
    f = os.path.join(raw, "jamstec_k2s1_bottle", "Mirai_routinedata_2010-2012_ver2.xlsx")
    if not os.path.exists(f):
        return None
    d = pd.read_excel(f, sheet_name="Hydrographic data")
    d = d[d["Station"].astype(str).str.match(r"K02_\d+$")].copy()
    tm = d["hh:mm"].astype(str).str.slice(0, 8)
    d["_t"] = pd.to_datetime(d["mon/day/yr"].astype(str) + " " + tm, format="%m/%d/%Y %H:%M:%S",
                             errors="coerce")
    d["_cast"] = d["Station"].astype(str).str.extract(r"K02_(\d+)")[0].astype(int)
    return d


def _jam_bottle(raw, d):
    if d is None:
        return []
    cols = list(d.columns)
    qf = {c: cols[i + 1] for i, c in enumerate(cols[:-1]) if str(cols[i + 1]).startswith("QF")}
    t = iso(d["_t"])
    z = pd.to_numeric(d["Depth [m]"], errors="coerce")
    lat = pd.to_numeric(d["Latitude"], errors="coerce")
    lon = pd.to_numeric(d["Longitude"], errors="coerce")
    src = "jamstec_k2s1_bottle"
    out = []

    def good(c):
        v = pd.to_numeric(d[c], errors="coerce")
        q = pd.to_numeric(d[qf[c]], errors="coerce")
        return v, q, q.isin(BOT_GOOD) & v.notna() & (v != -9)

    v, q, m = good("Temp. [C]")
    out.append(_rows(t[m], z[m], "temp", v[m], "degC (in-situ, CTD)", src, qc=q[m].astype(int).astype(str),
                     lat=lat[m], lon=lon[m]))
    vb, qb, mb = good("Salinity-Bottle")
    vc, qc_, mc = good("Salinity-CTD")
    mc = mc & ~mb
    out.append(_rows(t[mb], z[mb], "salt", vb[mb], "PSS-78 (bottle)", src, qc=qb[mb].astype(int).astype(str),
                     lat=lat[mb], lon=lon[mb]))
    out.append(_rows(t[mc], z[mc], "salt", vc[mc], "PSS-78 (CTD)", src, qc=qc_[mc].astype(int).astype(str),
                     lat=lat[mc], lon=lon[mc]))
    for c, var, u in BOT_VARS:
        v, q, m = good(c)
        out.append(_rows(t[m], z[m], var, v[m], u, src, qc=q[m].astype(int).astype(str),
                         lat=lat[m], lon=lon[m]))
    return out


# ---------------------------------------------------------------- JAMSTEC CTD
def _ctd_keep(p):
    p = np.asarray(p, dtype=float)
    return (p <= CTD_FULLRES_MAX_DBAR) | (np.isclose(np.mod(p, CTD_STRIDE_DBAR), 0.0))


def _read_wct(f):
    with open(f, encoding="latin-1") as fh:
        lines = fh.read().splitlines()
    h0 = lines[0]
    date = re.search(r"DATE\s+(\d{6})", h0).group(1)
    la = re.search(r"LATITUDE\s+([-\d.]+)", h0)
    lo = re.search(r"LONGITUDE\s+([-\d.]+)", h0)
    cast = int(re.search(r"CASTNO\s+(\d+)", lines[1]).group(1))
    names = lines[3].split()
    params = [n for n in names if n not in ("NUMBER", "QUALT1")]
    rows = []
    for s in lines[6:]:
        tok = s.split()
        if len(tok) != len(params) + 2:
            continue
        rows.append(tok)
    a = pd.DataFrame(rows, columns=params + ["NUMBER", "QUALT1"])
    for c in params:
        a[c] = pd.to_numeric(a[c], errors="coerce")
    for i, c in enumerate(params):
        a[c + "_Q"] = a["QUALT1"].str[i]
    meta = dict(date=pd.to_datetime(date, format="%m%d%y"), cast=cast,
                lat=float(la.group(1)) if la else np.nan, lon=float(lo.group(1)) if lo else np.nan)
    return a, meta


def _read_o4x(f):
    """ODV4 listing: '#' station header then value/flag pairs per row."""
    with open(f, encoding="latin-1") as fh:
        lines = fh.read().splitlines()
    names = []
    for s in lines:
        m = re.match(r"^(\S.*?\[[^\]]*\])\s+[\d.]+\s*$", s)
        if m:
            names.append(m.group(1).strip())
        if s.startswith("#"):
            break
    casts, cur, rows = [], None, []
    for s in lines:
        if s.startswith("#"):
            if cur is not None:
                casts.append((cur, rows))
            tok = s[1:].split()
            cur = dict(station=tok[1], date=pd.to_datetime(tok[3], format="%m/%d/%Y"),
                       lon=float(tok[4]), lat=float(tok[5]))
            rows = []
            continue
        if cur is None:
            continue
        tok = s.split()
        if len(tok) == 2 * len(names):
            rows.append(tok)
    if cur is not None:
        casts.append((cur, rows))
    return names, casts


def _ctd(raw, bot):
    rdir = os.path.join(raw, "jamstec_k2s1_ctd")
    if not os.path.isdir(rdir):
        return []
    src = "jamstec_k2s1_ctd"
    # cast start time/position from the bottle file (same cruise + K02 cast number), else header date 00:00
    lut = {}
    if bot is not None:
        g = bot.groupby(["Cruise", "_cast"]).first()
        lut = {k: (r["_t"], r["Latitude"], r["Longitude"]) for k, r in g.iterrows()}
    out = []

    def emit(cruise, cast, date, hlat, hlon, p, z, cols, units):
        t0, la, lo = lut.get((cruise, cast), (None, hlat, hlon))
        if t0 is None or pd.isna(t0) or t0.normalize() != date.normalize():
            t0 = date  # date only -> 00:00 UTC
            tq = "time=date_only"
        else:
            tq = "time=bottle_cast"
        keep = _ctd_keep(p)
        for var, (v, ok, flag) in cols.items():
            m = keep & ok & np.isfinite(v) & (v > -99)
            n = int(m.sum())
            if n == 0:
                continue
            out.append(_rows([iso(pd.Series([t0])).iloc[0]] * n, z[m], var, v[m], units[var], src,
                             qc=f"{flag}|{tq}", lat=la, lon=lo))

    for f in sorted(glob.glob(os.path.join(rdir, "MR*", "*.wct"))):
        cruise = os.path.basename(os.path.dirname(f))
        try:
            a, meta = _read_wct(f)
        except Exception as e:
            print(f"  std_k2: cannot read {f}: {e!r}")
            continue
        p = a["CTDPRS"].values
        z = a["CTDDEP"].values if "CTDDEP" in a else p2z(p)
        cols = {}
        for c, var in [("CTDTMP", "temp"), ("CTDSAL", "salt"), ("CTDOXY", "O2_CTD")]:
            if c in a:
                cols[var] = (a[c].values, (a[c + "_Q"] == "2").values, "2")
        units = {"temp": "degC (in-situ, CTD, ITS-90)", "salt": "PSS-78 (CTD)", "O2_CTD": "umol/kg (CTD sensor)"}
        emit(cruise, meta["cast"], meta["date"], meta["lat"], meta["lon"], p, z, cols, units)
    for f in sorted(glob.glob(os.path.join(rdir, "MR*", "*.o4x"))):
        cruise = os.path.basename(os.path.dirname(f))
        names, casts = _read_o4x(f)
        idx = {n.split(" [")[0]: i for i, n in enumerate(names)}
        for meta, rows in casts:
            if not meta["station"].startswith("K02") or not rows:
                continue
            arr = np.array(rows, dtype=float)
            p = arr[:, 2 * idx["Pressure"]]
            z = p2z(p, meta["lat"])
            cols = {}
            for c, var in [("Temperature", "temp"), ("Salinity", "salt"), ("Oxygen", "O2_CTD")]:
                cols[var] = (arr[:, 2 * idx[c]], arr[:, 2 * idx[c] + 1] == 0, "0")
            units = {"temp": "degC (in-situ, CTD)", "salt": "psu (CTD)", "O2_CTD": "umol/kg (CTD sensor)"}
            cast = int(meta["station"].split("_")[1])
            emit(cruise, cast, meta["date"], meta["lat"], meta["lon"], p, z, cols, units)
    return out


# ---------------------------------------------------------------- JAMSTEC primary production
def _pp(raw):
    out = []
    src = "jamstec_k2s1_pp"
    for f in sorted(glob.glob(os.path.join(raw, "jamstec_k2s1_pp", "*_PP.xlsx"))):
        d = pd.read_excel(f, header=1)
        d.columns = [re.sub(r"\s+", " ", str(c)).strip() for c in d.columns]
        d = d[d["Stn."].astype(str).str.strip() == "K2"]
        if d.empty:
            continue
        t = iso(pd.to_datetime(d["Date"].astype(str), format="%m/%d/%Y", errors="coerce"))
        z = pd.to_numeric(d["Inc. depth"], errors="coerce")
        lat = pd.to_numeric(d["Latitude"].astype(str).str.rstrip("N"), errors="coerce")
        lon = pd.to_numeric(d["Longitude"].astype(str).str.rstrip("E"), errors="coerce")
        inc = d["Inc. type"].astype(str).str.strip()
        for c, fc, var, rep in [("PP-1", "PP-1 _FLAG", "PP", "1"), ("PP-2", "PP-2 _FLAG", "PP", "2"),
                                ("PP (dark)", "PP (dark) _FLAG", "PP_dark", "1")]:
            v = pd.to_numeric(d[c], errors="coerce")
            fl = pd.to_numeric(d[fc], errors="coerce")
            m = (fl == 2) & v.notna()
            q = "2|rep=" + rep + "|inc=" + inc[m] + "|light_pct=" + d.loc[m, "Light Intensity"].astype(str)
            out.append(_rows(t[m], z[m], var, v[m], "mgC/m3/d (13C, dawn-to-dawn 24 h)", src, qc=q,
                             lat=lat[m], lon=lon[m]))
    return out


# ---------------------------------------------------------------- sediment traps
XLSX_TRAP_DEPLOYMENTS = {"20100222", "20101101", "20110705", "20120616", "20130725"}  # covered by MSTdata_K2.xlsx


def _trap_xlsx(raw):
    f = os.path.join(raw, "jamstec_k2s1_trap", "MSTdata_K2.xlsx")
    if not os.path.exists(f):
        return []
    src = "jamstec_k2s1_trap"
    out = []
    xl = pd.ExcelFile(f)
    for s in xl.sheet_names:
        d = pd.read_excel(xl, sheet_name=s, header=None)
        # data rows: col2 = cup open date (datetime), col3 = days, col5 = trap depth
        op = pd.to_datetime(d[2], errors="coerce", format="mixed")
        dur = pd.to_numeric(d[3], errors="coerce")
        dep = pd.to_numeric(d[5], errors="coerce")
        m = op.notna() & dur.notna() & dep.notna()
        d, op, dur, dep = d[m], op[m], dur[m], dep[m]
        mid = op + pd.to_timedelta(dur / 2.0, unit="D")
        q = "open=" + op.dt.strftime("%Y-%m-%d") + "|dur_d=" + dur.astype(int).astype(str)
        # (flux col, variable, units, concentration col: flux == 0 with conc 999 = not measured)
        for col, var, u, cc in [(7, "mass_flux", "mg/m2/d", None), (21, "POC_flux", "mg C/m2/d", 8),
                                (22, "bSi_flux", "mg opal/m2/d (Si-based, (Si-3.42Al)x67.2/28)", 12),
                                (23, "PIC_flux", "mg CaCO3/m2/d (Ca-based)", 13),
                                (24, "lith_flux", "mg/m2/d (Al-based, Al/8%)", 11)]:
            v = pd.to_numeric(d[col], errors="coerce")
            ok = v.notna() & (v != 999)
            if cc is not None:
                ok &= ~((v == 0) & (pd.to_numeric(d[cc], errors="coerce") == 999))
            out.append(_rows(iso(mid[ok]), dep[ok], var, v[ok], u, src, qc=q[ok], lat=LAT0, lon=LON0))
    return out


OS_VARS = [  # (nc var, tidy variable, units, concentration var whose QC must be good when flux == 0)
    ("TMF", "mass_flux", "mg/m2/d", None),
    ("Org-Cflux", "POC_flux", "mg C/m2/d", "Org-C"),
    ("Opalflux", "bSi_flux", "mg opal/m2/d (Si-based)", "Si"),
    ("CabasedCaCO3flux", "PIC_flux", "mg CaCO3/m2/d (Ca-based)", "Ca"),
    ("Lithogenicflux", "lith_flux", "mg/m2/d (Al-based)", "Al"),
]


def _trap_oceansites(raw):
    files = sorted(glob.glob(os.path.join(raw, "oceansites_k2_trap", "OS_JAMSTECK2_*_D_ST*.nc")))
    if not files:
        return []
    try:
        import netCDF4
    except Exception as e:
        print(f"  std_k2: netCDF4 unavailable: {e}")
        return []
    src = "oceansites_k2_trap"
    out = []
    for f in files:
        dep_id = os.path.basename(f).split("_")[2]
        if dep_id in XLSX_TRAP_DEPLOYMENTS:
            continue  # same deployments with reported fluxes in jamstec_k2s1_trap
        with netCDF4.Dataset(f) as ds:
            tv = ds.variables["TIME"]
            op = pd.to_datetime(netCDF4.num2date(tv[:], tv.units, only_use_cftime_datetimes=False,
                                                 only_use_python_datetimes=True))
            dur = np.ma.filled(ds.variables["SamplingInterval"][:, 0].astype(float), np.nan)
            mid = op + pd.to_timedelta(dur / 2.0, unit="D")
            dep = float(ds.variables["DEPTH"][0])
            q0 = ("open=" + pd.Series(op).dt.strftime("%Y-%m-%d") + "|dur_d="
                  + pd.Series(dur).map(lambda x: f"{x:g}"))

            def get(n):
                v = np.ma.filled(ds.variables[n][:, 0].astype(float), np.nan)
                qn = n + "_QC"
                qq = (np.ma.filled(ds.variables[qn][:, 0], -1).astype(int) if qn in ds.variables
                      else np.zeros(len(v), int))
                return v, qq

            for n, var, u, conc in OS_VARS:
                if n not in ds.variables:
                    continue
                v, qq = get(n)
                ok = np.isfinite(v) & np.isin(qq, [1, 2]) & (v > -99) & np.isfinite(dur)
                if conc is not None and conc in ds.variables:
                    _, cq = get(conc)
                    ok &= ~((v == 0) & ~np.isin(cq, [1, 2]))
                if ok.any():
                    out.append(_rows(iso(pd.Series(mid[ok])), dep, var, v[ok], u, src,
                                     qc=(pd.Series(qq[ok]).astype(str) + "|" + q0[ok].values).values,
                                     lat=LAT0, lon=LON0))
    return out


def load(site, site_dir):
    if site != "K2":
        return pd.DataFrame()
    raw = os.path.join(site_dir, "raw")
    parts = []
    try:
        bot = _bottle_table(raw)
    except Exception as e:
        print(f"  std_k2: bottle xlsx failed: {e!r}")
        bot = None
    for fn, args in ((_ncei_bottle, (raw,)), (_jam_bottle, (raw, bot)), (_ctd, (raw, bot)),
                     (_pp, (raw,)), (_trap_xlsx, (raw,)), (_trap_oceansites, (raw,))):
        try:
            parts += [p for p in fn(*args) if len(p)]
        except Exception as e:
            print(f"  std_k2: {fn.__name__} failed: {e!r}")
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
