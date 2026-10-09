"""Reader for BODC NODB discrete water-sample data at PAP (added 2026-10-08) -> tidy rows.

Source: BODC National Oceanographic Database (NODB) request RN-2506 (2026-10-08): the 209
unrestricted series in the BODC search with Site = Porcupine Abyssal Plain (PAP), data types
Water column chemistry / Water sample data / CTD / fluorescence / PAR, 1992 onward
(CLASS, Porcupine Abyssal Plain Observatory and RAGNARoCC collections; cruises 2007-2022).
ODV text exports in PAP/raw/bodc_nodb/unz/*.txt; per-series HTML documentation alongside.
Data supplied by the British Oceanographic Data Centre (BODC), National Oceanography Centre,
under the BODC data licence accepted at download.

Only discrete bottle measurements are ingested here (nutrients, Winkler O2, DIC, ALK, extracted
chl). Moored/CTD sensor series (pCO2, optode O2, ISUS NO3, fluorescence, T/S) duplicate the
OceanSITES PAP files already read by std_pap, so they are skipped.

Flags: SeaDataNet QV; kept 1 (good), 2 (probably good) and 0 (no QC) and blank.
"""
import glob
import os

import numpy as np
import pandas as pd

COLS = ["time", "depth_m", "variable", "value", "units", "source", "qc", "lat", "lon"]
SUB = "bodc_nodb"
# ODV column -> (std variable, units)
VARMAP = {
    "NO3+NO2_Unfilt_ColAA [umol/l]": ("NO3", "umol/L (NO3+NO2)"),
    "NO2_Unfilt_ColAA [umol/l]": ("NO2", "umol/L"),
    "PO4_Unfilt_ColAA [umol/l]": ("PO4", "umol/L"),
    "SiOx_Unfilt_ColAA [umol/l]": ("SiO2", "umol/L"),
    "WC_dissO2_Winkler [umol/l]": ("O2", "umol/L (Winkler)"),
    "TCO2/kg [umol/kg]": ("DIC", "umol/kg"),
    "TotAlk/kg [umol/kg]": ("ALK", "umol/kg"),
    "chl-a_water>GF/F_FA_fluor [mg/m^3]": ("Chl", "mg/m3 (extracted, GF/F)"),
}
KEEP_QV = {"", "0", "1", "2"}
# JC087 (2013) nutrients are already ingested from PANGAEA by std_pap: skip to avoid duplicates
SKIP = {("JC087", "NO3"), ("JC087", "NO2"), ("JC087", "PO4"), ("JC087", "SiO2")}


def _read(path):
    with open(path, errors="replace") as fh:
        lines = [l.rstrip("\n") for l in fh if not l.startswith("//")]
    if len(lines) < 2:
        return None
    hdr = lines[0].split("\t")
    rows = [l.split("\t") for l in lines[1:]]
    rows = [r + [""] * (len(hdr) - len(r)) for r in rows]
    d = pd.DataFrame(rows, columns=range(len(hdr)), dtype=str)
    meta = {"Cruise": 0, "time": 3, "lon": 4, "lat": 5}
    for k, i in meta.items():                       # ODV leaves metadata blank after a station's first row
        d[i] = d[i].replace("", np.nan).ffill()
    return hdr, d, meta


def load(site, site_dir):
    if site != "PAP":
        return pd.DataFrame(columns=COLS)
    files = sorted(glob.glob(os.path.join(site_dir, "raw", SUB, "unz", "*.txt")))
    out = []
    for f in files:
        r = _read(f)
        if r is None:
            continue
        hdr, d, meta = r
        if "DepBelowSurf [m]" not in hdr:
            continue
        iz = hdr.index("DepBelowSurf [m]")
        for col, (std, units) in VARMAP.items():
            if col not in hdr:
                continue
            i = hdr.index(col)
            v = pd.to_numeric(d[i].str.strip(), errors="coerce")
            q = d[i + 1].str.strip() if i + 1 < len(hdr) and hdr[i + 1].startswith("QV") else pd.Series("", index=d.index)
            cru = d[meta["Cruise"]].astype(str).str.strip()
            ok = v.notna() & q.isin(KEEP_QV) & ~cru.map(lambda c: (c, std) in SKIP)
            if not ok.any():
                continue
            out.append(pd.DataFrame({
                "time": pd.to_datetime(d[meta["time"]], utc=True, errors="coerce").dt.strftime("%Y-%m-%dT%H:%M:%SZ"),
                "depth_m": pd.to_numeric(d[iz].str.strip(), errors="coerce"),
                "variable": std, "value": v, "units": units,
                "source": "bodc_nodb_" + d[meta["Cruise"]].astype(str).str.strip(),
                "qc": q, "lat": pd.to_numeric(d[meta["lat"]], errors="coerce"),
                "lon": pd.to_numeric(d[meta["lon"]], errors="coerce")})[ok])
    if not out:
        return pd.DataFrame(columns=COLS)
    df = pd.concat(out, ignore_index=True)
    print(f"  [std_bodc] {len(df)} rows from {len(files)} files")
    return df[COLS]
