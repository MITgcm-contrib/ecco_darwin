"""Reader for the Munida Time-Series Transect, offshore subantarctic-water (SAW) station -> tidy rows.

Site (model column): MUNIDA, 45.83S 171.53E (offshore end of the Munida transect off Otago, New Zealand,
~63 km east of the Otago Peninsula; bottom ~1180 m, ETOPO1). See MUNIDA/SOURCES_MUNIDA.md.

Sources (raw dirs under <site_dir>/raw/):
  bcodmo_spots/        SPOTS v2 (BCO-DMO 896862 v2, Lange et al. 2024 ESSD), TimeSeriesSite = "Munida":
                       the SAW deep cast (7-1007 dbar), 2001-2019. WOCE-style CSV.
  ncei_ocads_0100218/  NCEI OCADS 0100218, R/V Munida underway SST/SSS/pCO2, Jan 2004 - Jul 2006 (whole transect)
  ncei_ocads_0112883/  NCEI OCADS 0112883, R/V Polaris II underway SST/SSS/pCO2, Aug - Oct 2006 (whole transect)

Station cut: only samples within MAX_KM (15 km) great-circle distance of the site position are kept
(the underway files cover the whole 65 km transect from the coast; 15 km keeps the SAW end only,
>= ~48 km offshore, east of the Subtropical Front). All SPOTS Munida rows are already within 12.6 km.

Flag conventions (rows failing are dropped):
  SPOTS: WOCE *_FLAG_W, keep 2 (good) and 6 (replicate mean); 9 = not measured. CTDTMP has no flag
      (kept if > -99). SALNTY/OXYGEN/NITRIT/NH4/pH/pCO2/particulates/DOC are all flag 9 for Munida.
  NCEI underway: no quality flags; -999 = missing (dropped); pCO2 or SSS <= 0 dropped (one pCO2 = 0 in M050121).
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

LAT0, LON0 = -45.83, 171.53
MAX_KM = 15.0
UNDERWAY_DEPTH_M = 2.0  # intake depth not stated in the metadata; nominal (Munida surface bottles at ~2 m)
SPOTS_KEEP = {2, 6}


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
    p1, p2 = np.deg2rad(LAT0), np.deg2rad(lat)
    a = (np.sin((p2 - p1) / 2) ** 2
         + np.cos(p1) * np.cos(p2) * np.sin(np.deg2rad(lon - LON0) / 2) ** 2)
    return 2 * 6371.0 * np.arcsin(np.sqrt(a))


# ---------------------------------------------------------------- SPOTS (BCO-DMO 896862 v2)
SPOTS_MAP = [  # (column, variable, units, flag column or None)
    ("NITRAT", "NO3", "umol/kg (SPOTS NITRAT; NO3 or NO3+NO2 not stated)", "NITRAT_FLAG_W"),
    ("PHSPHT", "PO4", "umol/kg", "PHSPHT_FLAG_W"),
    ("SILCAT", "SiO2", "umol/kg", "SILCAT_FLAG_W"),
    ("TCARBN", "DIC", "umol/kg", "TCARBN_FLAG_W"),
    ("ALKALI", "ALK", "umol/kg", "ALKALI_FLAG_W"),
    # all flag 9 for Munida in SPOTS v2 (kept here so any future good values come through)
    ("OXYGEN", "O2", "umol/kg", "OXYGEN_FLAG_W"),
    ("NITRIT", "NO2", "umol/kg", "NITRIT_FLAG_W"),
    ("NH4", "NH4", "umol/kg", "NH4_FLAG_W"),
    ("PH_TOT", "pH", "total scale at PH_TMP", "PH_TOT_FLAG_W"),
    ("PCO2", "pCO2", "uatm at PCO2_TMP", "PCO2_FLAG_W"),
    ("DOC", "DOC", "umol/kg", "DOC_FLAG_W"),
]


def _spots(raw):
    files = sorted(glob.glob(os.path.join(raw, "bcodmo_spots", "spots_*Munida*.csv")))
    if not files:
        return []
    d = pd.concat([pd.read_csv(f, skiprows=[1], low_memory=False) for f in files], ignore_index=True)
    d = d[d["TimeSeriesSite"].astype(str).str.strip() == "Munida"]
    lat = pd.to_numeric(d["LATITUDE"], errors="coerce")
    lon = pd.to_numeric(d["longitude"], errors="coerce")
    d = d[dist_km(lat, lon) <= MAX_KM]
    lat, lon = lat[d.index], lon[d.index]
    # TIME is -999 for all Munida rows -> date only, 00:00 UTC
    tt = pd.to_datetime(d["DATE"].astype(int).astype(str), format="%Y%m%d", errors="coerce")
    t = iso(tt)
    p = pd.to_numeric(d["CTDPRS"], errors="coerce")
    z = pd.Series(p2z(p.values, LAT0), index=d.index)
    src = "bcodmo_spots"
    out = []
    v = pd.to_numeric(d["CTDTMP"], errors="coerce")
    m = v > -99
    out.append(_rows(t[m], z[m], "temp", v[m], "degC (in-situ, CTD, ITS-90)", src, lat=lat[m], lon=lon[m]))
    # salinity: bottle (SALNTY) where flag 2/6, else CTDSAL flag 2/6 (all Munida SALNTY are flag 9)
    sb, fb = pd.to_numeric(d["SALNTY"], errors="coerce"), pd.to_numeric(d["SALNTY_FLAG_W"], errors="coerce")
    sc, fc = pd.to_numeric(d["CTDSAL"], errors="coerce"), pd.to_numeric(d["CTDSAL_FLAG_W"], errors="coerce")
    mb = fb.isin(SPOTS_KEEP) & (sb > 0)
    mc = ~mb & fc.isin(SPOTS_KEEP) & (sc > 0)
    out.append(_rows(t[mb], z[mb], "salt", sb[mb], "PSS-78 (bottle)", src, qc=fb[mb].astype(int).astype(str),
                     lat=lat[mb], lon=lon[mb]))
    out.append(_rows(t[mc], z[mc], "salt", sc[mc], "PSS-78 (CTD)", src, qc=fc[mc].astype(int).astype(str),
                     lat=lat[mc], lon=lon[mc]))
    for col, var, units, fcol in SPOTS_MAP:
        if col not in d:
            continue
        v = pd.to_numeric(d[col], errors="coerce")
        f = pd.to_numeric(d[fcol], errors="coerce")
        m = f.isin(SPOTS_KEEP) & (v > -99)
        out.append(_rows(t[m], z[m], var, v[m], units, src, qc=f[m].astype(int).astype(str),
                         lat=lat[m], lon=lon[m]))
    return out


# ---------------------------------------------------------------- NCEI OCADS underway pCO2
def _underway(raw):
    out = []
    for acc in ("0100218", "0112883"):  # R/V Munida 2004-2006, R/V Polaris II 2006
        files = sorted(glob.glob(os.path.join(raw, f"ncei_ocads_{acc}", "*.csv")))
        if not files:
            continue
        src = f"ncei_ocads_{acc}"
        for f in files:
            d = pd.read_csv(f, skiprows=[1], low_memory=False)
            la = pd.to_numeric(d["LATITUDE"], errors="coerce")
            lo = pd.to_numeric(d["LONGITUDE"], errors="coerce")
            d = d[dist_km(la, lo) <= MAX_KM]
            if d.empty:
                continue
            la, lo = la[d.index], lo[d.index]
            # times taken as UTC (not stated; consistent with cruise IDs = next-day NZ local date)
            tt = pd.to_datetime(d["DATE"].astype(str) + " " + d["TIME"].astype(str),
                                format="%m/%d/%Y %H:%M:%S", errors="coerce")
            t = iso(tt)
            q = "cruise=" + d["CRUISE"].astype(str)
            for col, var, units in (("SST", "temp", "degC (underway SBE21, intake depth not stated)"),
                                    ("SSS", "salt", "PSS-78 (underway SBE21)"),
                                    ("PCO2_WAT", "pCO2", "uatm (underway equilibrator, seawater pCO2)")):
                v = pd.to_numeric(d[col], errors="coerce")
                m = (v > 0) if var in ("pCO2", "salt") else (v > -99)  # -999 missing; one pCO2 = 0 (M050121) dropped
                out.append(_rows(t[m], UNDERWAY_DEPTH_M, var, v[m], units, src, qc=q[m],
                                 lat=la[m], lon=lo[m]))
    return out


def load(site, site_dir):
    if site != "MUNIDA":
        return pd.DataFrame()
    raw = os.path.join(site_dir, "raw")
    parts = []
    for fn in (_spots, _underway):
        try:
            parts += [p for p in fn(raw) if len(p)]
        except Exception as e:
            print(f"  std_munida: {fn.__name__} failed: {e!r}")
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
