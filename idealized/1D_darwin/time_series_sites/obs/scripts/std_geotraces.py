"""Reader for GEOTRACES IDP2025 dissolved Fe near each site (added 2026-10-08) -> tidy rows.

Source: GEOTRACES IDP2025 seawater discrete-sample dataset, extracted with the GEOTRACES webODV
data extractor (https://geotraces.webodv.awi.de/, IDP2025 > seawater; variables DEPTH,
CTDPRS_UP_T_VALUE, Fe_D_CONC; all 1926 stations with dissolved Fe), ASCII spreadsheet.
Raw extract kept once in obs/_geotraces_idp2025_raw/; `split` writes the per-site subsets
(stations within RADIUS_KM of the site) to <SITE>/raw/geotraces_idp2025/.

Citation: GEOTRACES Intermediate Data Product Group (2025). The GEOTRACES Intermediate Data Product
2025 (IDP2025). NERC EDS British Oceanographic Data Centre NOC. doi:10.5285/42c92148-8d03-8be6-e063-7086abc09f0c
(CC-BY 4.0; see the GEOTRACES Fair Data Use Statement).

Variable: Fe (dissolved Fe, Fe_D_CONC), nmol/kg as reported. Flags: SeaDataNet QV, kept 1 (good)
and 2 (probably good).

  python3 -I std_geotraces.py split <raw_extract.txt> <obs_root>     # write per-site subsets
"""
import glob
import os
import sys

import numpy as np
import pandas as pd

COLS = ["time", "depth_m", "variable", "value", "units", "source", "qc", "lat", "lon"]
SITE_LL = {"HOT": (22.75, -158.0), "BATS": (31.67, -64.17), "HydroS": (32.17, -64.5),
           "Papa": (50.1, -144.9), "PAP": (49.0, -16.5)}
RADIUS_KM = 120.
SUB = "geotraces_idp2025"
SKIP_CRUISES = {"GApr13"}


def _dist_km(lat, lon, lat0, lon0):
    r = 6371.
    p = np.radians
    a = np.sin(p(lat - lat0) / 2) ** 2 + np.cos(p(lat)) * np.cos(p(lat0)) * np.sin(p(lon - lon0) / 2) ** 2
    return 2 * r * np.arcsin(np.sqrt(a))


def _read_odv(path):
    with open(path, errors="replace") as fh:                 # skip the ODV '//' comment header
        nskip = next(i for i, line in enumerate(fh) if not line.startswith("//"))
    d = pd.read_csv(path, sep="\t", dtype=str, low_memory=False, skiprows=nskip)
    meta = ["Cruise", "Station", "yyyy-mm-ddThh:mm:ss.sss", "Longitude [degrees_east]",
            "Latitude [degrees_north]"]
    d[meta] = d[meta].replace("", np.nan).ffill()          # ODV may leave metadata blank on later rows
    return d


def split(raw, obs_root):
    d = _read_odv(raw)
    lat = d["Latitude [degrees_north]"].astype(float)
    lon = d["Longitude [degrees_east]"].astype(float)
    lon = np.where(lon > 180, lon - 360, lon)
    for site, (la, lo) in SITE_LL.items():
        k = _dist_km(lat, lon, la, lo) <= RADIUS_KM
        out = os.path.join(obs_root, site, "raw", SUB)
        os.makedirs(out, exist_ok=True)
        d[k].to_csv(os.path.join(out, "idp2025_fe_within_%dkm.tsv" % RADIUS_KM), sep="\t", index=False)
        st = d[k].groupby(["Cruise", "Station"]).size()
        print("%-6s %4d rows, %d stations, cruises %s" % (site, k.sum(), len(st),
                                                         sorted(set(d[k]["Cruise"]))))


def load(site, site_dir):
    files = glob.glob(os.path.join(site_dir, "raw", SUB, "*.tsv"))
    if not files:
        return pd.DataFrame(columns=COLS)
    d = pd.concat([pd.read_csv(f, sep="\t", dtype=str) for f in files], ignore_index=True)
    fe_col = [c for c in d.columns if c.startswith("Fe_D_CONC")][0]
    qv_col = d.columns[list(d.columns).index(fe_col) + 2]           # Fe, STANDARD_DEV, QV:SEADATANET
    v = pd.to_numeric(d[fe_col], errors="coerce")
    q = pd.to_numeric(d[qv_col], errors="coerce")
    # GApr13 (BAIT 2019) is already ingested from BCO-DMO 936824/937302 by std_tm: skip to avoid duplicates
    keep = v.notna() & q.isin([1, 2]) & ~d["Cruise"].isin(SKIP_CRUISES)
    t = pd.to_datetime(d["yyyy-mm-ddThh:mm:ss.sss"], utc=True, errors="coerce")
    out = pd.DataFrame({
        "time": t.dt.strftime("%Y-%m-%dT%H:%M:%SZ"), "depth_m": pd.to_numeric(d["DEPTH [m]"], errors="coerce"),
        "variable": "Fe", "value": v, "units": "nmol/kg (GEOTRACES IDP2025 Fe_D_CONC)",
        "source": "geotraces_idp2025_" + d["Cruise"].astype(str), "qc": q.astype("Int64").astype(str),
        "lat": pd.to_numeric(d["Latitude [degrees_north]"], errors="coerce"),
        "lon": pd.to_numeric(d["Longitude [degrees_east]"], errors="coerce")})[keep]
    print(f"  [std_geotraces] {len(out)} rows")
    return out[COLS]


if __name__ == "__main__":
    if sys.argv[1] == "split":
        split(sys.argv[2], sys.argv[3])
