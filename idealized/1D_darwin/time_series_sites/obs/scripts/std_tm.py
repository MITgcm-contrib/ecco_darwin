"""Reader for dissolved trace-metal data (added 2026-10-08) -> tidy rows.

Sources (raw dirs per site):
  HOT/raw/bcodmo_hot_metals/            BCO-DMO 962986 v2: HOT dissolved (+ total dissolvable) metals at
                                        ALOHA, 21 cruises from HOT-325 (Dec 2020-Nov 2023), 15-900 m. Hawco & Bates.
  BATS/raw/bcodmo_bait_fe/              BCO-DMO 936824 v1: BAIT 2019 dissolved Fe (+ d56Fe), Conway et al.
  BATS/raw/bcodmo_bait_tm_bottle/       BCO-DMO 937302 v1: BAIT 2019 trace-metal rosette bottle data, Sedwick et al.

The GEOTRACES IDP2025 itself (GA03, GP15, Canadian GEOTRACES P26, GEOVIDE) could NOT be fetched
(BODC host unreachable from this machine; webODV is interactive); see README_obs.md. If the IDP
discrete-sample file is later put in <SITE>/raw/geotraces_idp2025/, add a reader here.

Variables: Fe (dissolved Fe), Mn, Co (labile dissolved Co where stated), Ni, Cu, Zn, Cd, Pb, Al, Ti (dissolved).
Total-dissolvable (td*), soluble (s*), isotope ratios and the BAIT macronutrients are not used.
(BCO-DMO 792817, Hayes Lagrangian TMs near ALOHA, was checked: no Fe within 1 deg of the site, not used.)
Values and units as reported (nmol/L, nmol/kg or pmol/kg; see units column).

Flags:
  962986 (HOT): GEOTRACES/SeaDataNet 1 good, 2 probably good, 3 probably bad, 4 bad, 6 below DL, 9 missing.
                Kept 1, 2 (PI recommendation).
  936824 (BAIT Fe, Conway): SeaDataNet; "accurate" data flagged 2; 9 missing. Kept 1, 2.
  937302 (BAIT bottle, Sedwick): 1 good, 2 likely contaminated/questionable, 3 questionable (Sc internal std),
                4 not determined, 5 below detection. Kept 1 only.
"""
import glob
import os

import numpy as np
import pandas as pd

COLS = ["time", "depth_m", "variable", "value", "units", "source", "qc", "lat", "lon"]
SITE_LL = {"HOT": (22.75, -158.0), "BATS": (31.67, -64.17)}


def _iso(t):
    t = pd.to_datetime(t, utc=True, errors="coerce")
    return t.dt.strftime("%Y-%m-%dT%H:%M:%SZ")


def _dist_km(lat, lon, lat0, lon0):
    r = 6371.0
    p1, p2 = np.deg2rad(lat), np.deg2rad(lat0)
    dphi = p2 - p1
    dl = np.deg2rad(lon0 - lon)
    a = np.sin(dphi / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * r * np.arcsin(np.sqrt(a))


def _read(raw, sub, pat="*.csv"):
    files = sorted(glob.glob(os.path.join(raw, sub, pat)))
    if not files:
        return None
    return pd.concat([pd.read_csv(f, skiprows=[1], dtype=str) for f in files], ignore_index=True)


def _melt(d, spec, source, t, z, lat, lon, keep_flags):
    """spec: list of (value col, flag col or None, variable, units)."""
    out = []
    for col, fcol, var, units in spec:
        if col not in d.columns:
            continue
        v = pd.to_numeric(d[col], errors="coerce")
        ok = v.notna() & t.notna()
        if fcol is not None:
            f = pd.to_numeric(d[fcol], errors="coerce")
            ok &= f.isin(keep_flags)
            q = f.astype("Int64").astype(str)
        else:
            q = pd.Series("", index=d.index)
        out.append(pd.DataFrame({"time": t[ok], "depth_m": z[ok], "variable": var, "value": v[ok],
                                 "units": units, "source": source, "qc": q[ok],
                                 "lat": lat[ok], "lon": lon[ok]}))
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame(columns=COLS)


# ------------------------------------------------------------------ HOT
HOT_METALS = [("dFe", "dFe_qc", "Fe", "nmol/L"), ("dMn", "dMn_qc", "Mn", "nmol/L"),
              ("dCo", "dCo_qc", "Co", "nmol/L (labile dissolved Co, not UV-oxidised)"),
              ("dNi", "dNi_qc", "Ni", "nmol/L"), ("dCu", "dCu_qc", "Cu", "nmol/L"),
              ("dZn", "dZn_qc", "Zn", "nmol/L"), ("dCd", "dCd_qc", "Cd", "nmol/L"),
              ("dPb", "dPb_qc", "Pb", "nmol/L"), ("dTi", "dTi_qc", "Ti", "nmol/L")]


def read_hot_metals(raw):
    d = _read(raw, "bcodmo_hot_metals")
    if d is None:
        return pd.DataFrame(columns=COLS)
    t = _iso(d["time"])  # UTC (BCO-DMO converted from ISO_DateTime_Local, HST)
    z = pd.to_numeric(d["depth"], errors="coerce")
    lat = pd.to_numeric(d["latitude"], errors="coerce")
    lon = pd.to_numeric(d["longitude"], errors="coerce")
    return _melt(d, HOT_METALS, "bcodmo_hot_metals", t, z, lat, lon, {1, 2})


# ------------------------------------------------------------------ BATS (BAIT 2019)
def read_bait_fe(raw, site):
    d = _read(raw, "bcodmo_bait_fe")
    if d is None:
        return pd.DataFrame(columns=COLS)
    lat = pd.to_numeric(d["latitude"], errors="coerce")
    lon = pd.to_numeric(d["longitude"], errors="coerce")
    keep = _dist_km(lat, lon, *SITE_LL[site]) <= 100.0  # BATS spatial convention (std_bios.py)
    d, lat, lon = d[keep], lat[keep], lon[keep]
    t = _iso(d["Start_Date_UTC"])  # date only -> 00:00 UTC
    z = pd.to_numeric(d["depth"], errors="coerce")
    spec = [("Fe_D_CONC_BOTTLE", "Flag_Fe_D_CONC_BOTTLE", "Fe", "nmol/kg (GO-Flo bottle)"),
            ("Fe_D_CONC_BOAT_PUMP", "Flag_Fe_D_CONC_BOAT_PUMP", "Fe", "nmol/kg (surface towed/boat pump)")]
    return _melt(d, spec, "bcodmo_bait_fe", t, z, lat, lon, {1, 2})


def read_bait_bottle(raw, site):
    d = _read(raw, "bcodmo_bait_tm_bottle")
    if d is None:
        return pd.DataFrame(columns=COLS)
    lat = pd.to_numeric(d["latitude"], errors="coerce")
    lon = pd.to_numeric(d["longitude"], errors="coerce")
    keep = _dist_km(lat, lon, *SITE_LL[site]) <= 100.0
    d, lat, lon = d[keep], lat[keep], lon[keep]
    t = _iso(d["time"])  # UTC hydrocast time
    z = pd.to_numeric(d["depth"], errors="coerce")
    spec = [("DFe", "DFe_Flag", "Fe", "nmol/L"), ("DMn", "DMn_Flag", "Mn", "nmol/L"),
            ("DAl", "DAl_Flag", "Al", "nmol/L")]
    return _melt(d, spec, "bcodmo_bait_tm_bottle", t, z, lat, lon, {1})


def load(site, site_dir):
    raw = os.path.join(site_dir, "raw")
    if site == "HOT":
        fns = [read_hot_metals]
    elif site == "BATS":
        fns = [lambda r: read_bait_fe(r, site), lambda r: read_bait_bottle(r, site)]
    else:
        return pd.DataFrame(columns=COLS)
    parts = []
    for fn in fns:
        try:
            d = fn(raw)
        except Exception as e:
            print(f"  [std_tm] reader failed: {e!r}")
            continue
        print(f"  [std_tm] {len(d)} rows")
        parts.append(d)
    parts = [p for p in parts if len(p)]
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=COLS)
