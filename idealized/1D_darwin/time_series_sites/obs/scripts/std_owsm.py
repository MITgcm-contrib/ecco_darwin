"""Reader for Ocean Weather Station M ("Mike", OWSM; 66.0N 2.0E, Norwegian Sea) -> tidy rows.

Sources (raw dirs under <site_dir>/raw/; see OWSM/SOURCES_OWSM.md):
  oceansites_stationm/    OceanSITES STATION-M (UiB/Bjerknes, NCEI 0130051): weather-ship hydrography 1948-2009
                          (CTD.nc -> T/S; OCTD.nc -> Winkler O2 in mL/L). Only >= 1992 is kept.
  ncei_ocads_owsm/        NCEI OCADS 0112884 (Skjelvan): bottle DIC/ALK + CTD T/S 2001-2021, nutrients/O2 2001-2007
  nmdc_owsm_sdg/          NMDC (NORCE, Skjelvan) SDG 14.3.1 xlsx: bottle DIC/ALK/nutrients + CTD T/S 2022-2025
  nmdc_imr_nutrients/     NMDC IMR Plankton Chemistry Lab (Gundersen et al.): nutrients, Chl a, Phaeo, Winkler O2,
                          Nordic Seas 1990-2024 (subset to <= MAX_KM of the station)
  nmdc_noremso_mooring/   NMDC NorEMSO mooring (NORCE): hourly T/S at 500-2000 m 2020-2023, HydroC fCO2 at 500 m
                          2021-2023
  bcodmo_spots/           SPOTS OWSM subset (same data as ncei_ocads_owsm; downloaded for reference, NOT read)

Only samples within MAX_KM (30 km) of 66.0N 2.0E are kept (lat/lon of each sample written to the tidy file).

Flag conventions (rows failing are dropped):
  OceanSITES STATION-M: TEMP_QC/PSAL_QC keep 1 (good; values present are 1, 4, 9). POSITION_QC keep 1.
      DOXYWITX_QC is 0 (no QC performed) for every O2 value; kept with a 3-10 mL/L range check (as PAP OceanSITES 0).
      DOXY (umol/kg) is NOT used: its values are ~1e-12 to 0.07 (broken conversion in the file).
  NCEI OCADS / SPOTS-style WOCE *_FLAG_W: DIC, ALK, CTDSAL keep 2. CTDTMP has no flag (-999 dropped).
      OXYGEN/NITRAT/PHSPHT/SILCAT carry WOCE 0 in the NCEI file (PI note: "values interpolated to match carbon
      sampling times"; SPOTS re-flags the same values as 2). Kept with qc "0"; flags 3/4/9 dropped.
  NMDC SDG xlsx: *_flag keep 2 (values present: 2 and 9 with -999).
  IMR nutrients: SeaDataNet flags; keep 1 (good), 2 (unexpected but not discarded) and 5 (corrected/relocated),
      the three flags the provider publishes. 0 with a value does not occur in the subset.
  NorEMSO mooring: EuroSITES/OceanSITES flags, keep 1 and 2; drop 0 (raw, uncalibrated), 3, 4.

De-duplication: SPOTS is not read (identical to NCEI 0112884). SDG xlsx nutrients (2022-2024) are dropped on dates
where IMR nutrient samples exist within +-1 day (same casts, IMR is the analyzing lab), kept otherwise.
"""
import glob
import os

import numpy as np
import pandas as pd

try:
    from standardize_obs import iso
except Exception:  # pragma: no cover
    def iso(ts):
        return pd.to_datetime(ts, utc=True, errors="coerce").dt.strftime("%Y-%m-%dT%H:%M:%SZ")

LAT0, LON0 = 66.0, 2.0
MAX_KM = 30.0
T0 = pd.Timestamp("1992-01-01")


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


def dist_km(lat, lon):
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)
    a = (np.sin(np.deg2rad(lat - LAT0) / 2) ** 2 + np.cos(np.deg2rad(LAT0)) * np.cos(np.deg2rad(lat))
         * np.sin(np.deg2rad(lon - LON0) / 2) ** 2)
    return 2 * 6371.0 * np.arcsin(np.sqrt(a))


def _num(s):
    return pd.to_numeric(s, errors="coerce")


# ---------------------------------------------------------------- OceanSITES STATION-M 1948-2009
def _oceansites(raw):
    rdir = os.path.join(raw, "oceansites_stationm")
    f_ctd = os.path.join(rdir, "OS_STATION-M-1_194810_D_CTD.nc")
    f_oct = os.path.join(rdir, "OS_STATION-M-1_195304_D_OCTD.nc")
    if not os.path.exists(f_ctd):
        return []
    import netCDF4
    src = "oceansites_stationm"
    out = []

    def table(f, names):
        with netCDF4.Dataset(f) as ds:
            t = pd.Timestamp("1950-01-01") + pd.to_timedelta(np.ma.filled(ds["TIME"][:].astype(float), np.nan),
                                                              unit="D")
            d = pd.DataFrame({"t": t,
                              "lat": np.ma.filled(ds["LATITUDE"][:].astype(float), np.nan),
                              "lon": np.ma.filled(ds["LONGITUDE"][:].astype(float), np.nan),
                              "pq": np.ma.filled(ds["POSITION_QC"][:], -1).astype(int),
                              "z": np.ma.filled(ds["DEPH"][:].astype(float), np.nan).ravel(),
                              "zq": np.ma.filled(ds["DEPH_QC"][:], -1).astype(int).ravel()})
            for n in names:
                d[n] = np.ma.filled(ds[n][:].astype(float), np.nan).ravel()
                d[n + "_QC"] = np.ma.filled(ds[n + "_QC"][:], -1).astype(int).ravel()
        d = d[(d["t"] >= T0) & (d["pq"] == 1) & (d["zq"] == 1)]
        return d[dist_km(d["lat"], d["lon"]) <= MAX_KM]

    d = table(f_ctd, ["TEMP", "PSAL"])
    t = iso(d["t"])
    for n, var, u in [("TEMP", "temp", "degC (in-situ; Nansen/CTD, scale not stated)"),
                      ("PSAL", "salt", "psu (bottle salinometer per qc_manual)")]:
        m = d[n + "_QC"] == 1
        out.append(_rows(t[m], d["z"][m], var, d[n][m], u, src, qc="1", lat=d["lat"][m], lon=d["lon"][m]))
    if os.path.exists(f_oct):
        d = table(f_oct, ["DOXYWITX"])
        m = (d["DOXYWITX_QC"] == 0) & d["DOXYWITX"].between(3.0, 10.0)
        out.append(_rows(iso(d["t"][m]), d["z"][m], "O2", d["DOXYWITX"][m], "mL/L (Winkler)", src, qc="0",
                         lat=d["lat"][m], lon=d["lon"][m]))
    return out


# ---------------------------------------------------------------- NCEI OCADS 0112884, 2001-2021
def _ncei(raw):
    f = os.path.join(raw, "ncei_ocads_owsm", "58P320011031.csv")
    if not os.path.exists(f):
        return []
    with open(f, encoding="latin-1") as fh:
        lines = fh.read().splitlines()
    i0 = next(i for i, s in enumerate(lines) if s.startswith("EXPOCODE,"))
    body = [s for s in lines[i0 + 2:] if s.strip() and not s.startswith("END_DATA")]
    cols = lines[i0].split(",")
    d = pd.DataFrame([s.split(",") for s in body], columns=cols)
    for c in cols[2:]:
        d[c] = _num(d[c])
    d = d[d["SECT_ID"].str.strip() == "OWSM"]
    lat, lon = d["LATITUDE"], d["LONGITUDE"]
    d = d[dist_km(lat, lon) <= MAX_KM]
    lat, lon = d["LATITUDE"], d["LONGITUDE"]
    # TIME is 0 for every row (date only) -> 00:00 UTC
    t = iso(pd.to_datetime(d["DATE"].astype(int).astype(str), format="%Y%m%d", errors="coerce"))
    z = pd.Series(p2z(d["CTDPRS"].values, LAT0), index=d.index)
    src = "ncei_ocads_owsm"
    out = []
    v = d["CTDTMP"].where(d["CTDTMP"] > -99)
    out.append(_rows(t, z, "temp", v, "degC (in-situ, CTD, ITS-90)", src, lat=lat, lon=lon))
    m = (d["CTDSAL_FLAG_W"] == 2) & (d["CTDSAL"] > 0)
    out.append(_rows(t[m], z[m], "salt", d["CTDSAL"][m], "PSS-78 (CTD, adjusted to bottle salinity)", src,
                     qc="2", lat=lat[m], lon=lon[m]))
    for c, var in [("TCARBN", "DIC"), ("ALKALI", "ALK")]:
        m = (d[c + "_FLAG_W"] == 2) & (d[c] > -99)
        out.append(_rows(t[m], z[m], var, d[c][m], "umol/kg", src, qc="2", lat=lat[m], lon=lon[m]))
    for c, var in [("OXYGEN", "O2"), ("NITRAT", "NO3"), ("PHSPHT", "PO4"), ("SILCAT", "SiO2")]:
        fl = d[c + "_FLAG_W"]
        m = fl.isin([0, 2]) & (d[c] > -99)
        out.append(_rows(t[m], z[m], var, d[c][m], "umol/kg (interpolated to carbon sampling times by PI)", src,
                         qc=fl[m].astype(int).astype(str), lat=lat[m], lon=lon[m]))
    return out


# ---------------------------------------------------------------- NMDC SDG 14.3.1 xlsx, 2022-2025
SDG_VARS = [("T_degC", "T_flag", "temp", "degC (in-situ, CTD)"),
            ("S_CTD", "S_flag", "salt", "PSS-78 (CTD, bottle-corrected)"),
            ("DIC_umol/kg", "DIC_flag", "DIC", "umol/kg"),
            ("TA_umol/kg", "TA_flag", "ALK", "umol/kg"),
            ("NO3_umol/kg", "NO3_flag", "NO3", "umol/kg"),
            ("PO4_umol/kg", "PO4_flag", "PO4", "umol/kg"),
            ("Si_umol/kg", "Si_flag", "SiO2", "umol/kg")]


def _sdg(raw, imr_nut_dates):
    files = sorted(glob.glob(os.path.join(raw, "nmdc_owsm_sdg", "OWSM_*_SDG14-3-1*.xlsx")))
    src = "nmdc_owsm_sdg"
    out = []
    for f in files:
        d = pd.read_excel(f, sheet_name=0)
        d.columns = [str(c).strip() for c in d.columns]  # header has "  DIC_umol/kg"
        d = d[_num(d["Lat_degreeN"]).notna()]
        lat, lon = _num(d["Lat_degreeN"]), _num(d["Lon_degreeE"])
        d = d[dist_km(lat, lon) <= MAX_KM]
        lat, lon = _num(d["Lat_degreeN"]), _num(d["Lon_degreeE"])
        tt = pd.to_datetime(d["YYY-MM-DD"].astype(str).str[:10] + " " + d["hh:mm"].astype(str).str[:8],
                            errors="coerce")
        t = iso(tt)
        z = _num(d["Depth_m"])
        day = tt.dt.normalize()
        near_imr = day.apply(lambda x: any(abs((x - y).days) <= 1 for y in imr_nut_dates))
        for c, fc, var, u in SDG_VARS:
            v, fl = _num(d[c]), _num(d[fc])
            m = (fl == 2) & (v > -99)
            if var in ("NO3", "PO4", "SiO2"):
                m &= ~near_imr
            out.append(_rows(t[m], z[m], var, v[m], u, src, qc="2", lat=lat[m], lon=lon[m]))
    return out


# ---------------------------------------------------------------- IMR nutrients / Chl, 1990-2024
IMR_VARS = [("Nitrate", "NO3", "umol/L (nitrate, corrected for nitrite)"), ("Nitrite", "NO2", "umol/L"),
            ("Phosphate", "PO4", "umol/L"), ("Silicate", "SiO2", "umol/L"),
            ("ChlA", "Chl", "mg/m3 (extracted, 90% acetone, fluorometric, acidified)"),
            ("Phaeo", "Phaeo", "mg/m3 (fluorometric)"), ("POC", "POC", "umol/L"), ("PN", "PON", "umol/L"),
            ("DOW", "O2", "mL/L (Winkler)")]
IMR_GOOD = {1, 2, 5}


def _imr_table(raw):
    files = sorted(glob.glob(os.path.join(raw, "nmdc_imr_nutrients", "NMDC_BGC_Norwegian-and-Greenland-Seas_*.csv")))
    parts = []
    for f in files:
        d = pd.read_csv(f, skiprows=[1], low_memory=False)
        lat, lon = _num(d["latitude"]), _num(d["longitude"])
        d = d[dist_km(lat, lon) <= MAX_KM].copy()
        # 2023-2024 file: 'DEPH' = sample depth, 'depth' = bottom depth; older files: 'depth' = sample depth
        d["_z"] = _num(d["DEPH"]) if "DEPH" in d.columns else _num(d["depth"])
        d["_f"] = os.path.basename(f)
        parts.append(d)
    return pd.concat(parts, ignore_index=True) if parts else None


def _imr(d):
    if d is None or d.empty:
        return []
    src = "nmdc_imr_nutrients"
    t = iso(pd.to_datetime(d["time"], utc=True, errors="coerce"))
    lat, lon = _num(d["latitude"]), _num(d["longitude"])
    out = []
    for c, var, u in IMR_VARS:
        if c not in d.columns:
            continue
        v, fl = _num(d[c]), _num(d[c + "_SEADATANET_QC"])
        m = v.notna() & fl.isin(IMR_GOOD)
        out.append(_rows(t[m], d["_z"][m], var, v[m], u, src, qc=fl[m].astype(int).astype(str),
                         lat=lat[m], lon=lon[m]))
    return out


# ---------------------------------------------------------------- NorEMSO mooring 2020-2023
def _mooring(raw):
    rdir = os.path.join(raw, "nmdc_noremso_mooring")
    if not os.path.isdir(rdir):
        return []
    import netCDF4
    src = "nmdc_noremso_mooring"
    out = []

    def grab(ds, n):
        v = np.ma.filled(ds[n][:].astype(float), np.nan)
        q = np.ma.filled(ds[n + "_QC"][:], -1).astype(int)
        if ds[n].dimensions[0] == "DEPTH":  # (DEPTH, TIME) -> (TIME, DEPTH)
            v, q = v.T, q.T
        return v, q

    for f in sorted(glob.glob(os.path.join(rdir, "*_hydrography.nc")) + glob.glob(os.path.join(rdir, "*_carbon.nc"))):
        name = os.path.basename(f)
        with netCDF4.Dataset(f) as ds:
            la, lo = float(ds["LATITUDE"][0]), float(ds["LONGITUDE"][0])
            if dist_km(la, lo) > MAX_KM:
                continue
            tim = pd.Timestamp("1950-01-01") + pd.to_timedelta(np.asarray(ds["TIME"][:], dtype=float), unit="D")
            tim = tim.round("s")
            nom = np.asarray(ds["DEPTH"][:], dtype=float)
            p, _ = grab(ds, "PRES") if "PRES" in ds.variables else (np.full((len(tim), len(nom)), np.nan), None)
            z = np.where(np.isfinite(p), p2z(p, la), nom[None, :])
            if name.endswith("_hydrography.nc"):
                todo = [("TEMP", "temp", "degC (in-situ, SBE37, corrected to ship CTD)"),
                        ("PSAL", "salt", "psu (SBE37, corrected to ship CTD)")]
            else:
                todo = [("FCO2", "fCO2", "uatm (HydroC at in-situ T, adjusted to DIC/TA-derived fCO2)")]
            for n, var, u in todo:
                if n not in ds.variables:
                    continue  # StationM_2021_carbon.nc has only uncalibrated pCO2 (flag 0)
                v, q = grab(ds, n)
                ts = pd.Series(np.repeat(np.asarray(tim), v.shape[1]))
                vv, qq, zz = v.ravel(), q.ravel(), z.ravel()
                m = np.isin(qq, [1, 2]) & np.isfinite(vv)
                if not m.any():
                    continue
                out.append(_rows(iso(ts[m]).values, zz[m], var, vv[m], u, src,
                                 qc=pd.Series(qq[m]).astype(str).values + "|" + name.replace(".nc", ""),
                                 lat=la, lon=lo))
    return out


def load(site, site_dir):
    if site != "OWSM":
        return pd.DataFrame()
    raw = os.path.join(site_dir, "raw")
    parts = []
    try:
        imr = _imr_table(raw)
    except Exception as e:
        print(f"  std_owsm: IMR table failed: {e!r}")
        imr = None
    imr_dates = []
    if imr is not None and len(imr):
        has_nut = imr[["Nitrate", "Phosphate", "Silicate"]].apply(_num).notna().any(axis=1)
        imr_dates = sorted(set(pd.to_datetime(imr.loc[has_nut, "time"], utc=True).dt.tz_localize(None)
                               .dt.normalize()))
    for fn, args in ((_oceansites, (raw,)), (_ncei, (raw,)), (_sdg, (raw, imr_dates)), (_imr, (imr,)),
                     (_mooring, (raw,))):
        try:
            parts += [p for p in fn(*args) if len(p)]
        except Exception as e:
            print(f"  std_owsm: {fn.__name__} failed: {e!r}")
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
