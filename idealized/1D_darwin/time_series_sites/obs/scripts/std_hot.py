"""Reader for HOT / Station ALOHA raw downloads -> tidy rows.

Sources (under HOT/raw/):
  bcodmo_bottle/      BCO-DMO dataset 3773 v1 (HOT Niskin bottle file, Station 2 = ALOHA), ERDDAP CSV
  bcodmo_pp/          BCO-DMO dataset 737163 v1 (HOT 14C primary production), ERDDAP CSV
  bcodmo_flux/        BCO-DMO dataset 737393 v1 (HOT sediment-trap particle flux), ERDDAP CSV
  hot_co2/            HOT_surface_CO2.txt (Dore et al. HOT surface carbonate product)
  ncei_whots_pco2/    NCEI OCADS 0100080, PMEL MAPCO2 on WHOTS mooring
  scripps_co2/        Scripps CO2 program HAWI.csv (Keeling surface DIC/ALK)
  bcodmo_pp_v5/       BCO-DMO 737163 v5 PP (HOT-1..355); only cruises absent from bcodmo_pp used
  bcodmo_flux_v4/     BCO-DMO 737393 v4 traps (HOT-2..355); only cruises absent from bcodmo_flux used
  bcodmo_spots/       BCO-DMO 896862 v2 SPOTS, ALOHA subset; only cruises absent from bcodmo_bottle used

Flag conventions (documented in HOT/SOURCES_HOT.md):
  HOT bottle/PP: 1=not QC'd, 2=good, 3=suspect, 4=bad, 5=missing, 9=not measured.
                 Kept: 1 and 2. Dropped: 3,4,5,9.
  PMEL MAPCO2 (WHOTS): WOCE-style 2=good, 3=questionable, 4=bad. Kept: 2 only.
"""
import glob
import os

import numpy as np
import pandas as pd

LAT0, LON0 = 22.75, -158.0

COLS = ["time", "depth_m", "variable", "value", "units", "source", "qc", "lat", "lon"]


def _iso(t):
    t = pd.to_datetime(t, utc=True, errors="coerce")
    return t.dt.strftime("%Y-%m-%dT%H:%M:%SZ")


def p2z(p, lat):
    """UNESCO 1983 (Fofonoff & Millard) pressure (dbar) -> depth (m)."""
    p = np.asarray(p, dtype=float)
    x = np.sin(np.deg2rad(lat)) ** 2
    g = 9.780318 * (1.0 + (5.2788e-3 + 2.36e-5 * x) * x) + 1.092e-6 * p
    return ((((-1.82e-15 * p + 2.279e-10) * p - 2.2512e-5) * p + 9.72659) * p) / g


# ------------------------------------------------------------------ bottle
# QUALTn digit order (BCO-DMO 3773 dataset description)
QGROUPS = {
    "QUALT1": ["CTDSAL", "CTDOXY", "SALNITY", "OXYGEN", "DIC", "PH", "ALKALIN"],
    "QUALT2": ["pCO2", "PHSPHT", "NO2_NO3", "SILCAT", "DOP", "DON", "DOC"],
    "QUALT3": ["TDP", "TDN", "PC", "PN", "PP", "LLN", "LLP"],
    "QUALT4": ["LLSi", "CHL_A", "PHEO", "CHL_C3", "CHLC1_2", "CHL_Plus", "PERID"],
    "QUALT5": ["BUT_19", "FUCO", "HEX_19", "PRASINO", "DIADINO", "ZEAXAN", "CHL_B"],
    "QUALT6": ["HPLCchl", "CHL_C4", "A_CAR", "B_CAR", "CAROTEN", "CHLDA_A"],
    "QUALT7": ["VIOL", "LUTEIN", "MV_CHLA", "DV_CHLA", "H_BACT", "P_BACT"],
    "QUALT8": ["S_BACT", "E_BACT", "ATP", "GTP", "H2O2", "N2O"],
    "QUALT9": ["PSi", "PIC", "PE_pt4u", "PE_5u", "PE_10u", "P15N"],
    "QUAL10": ["P13C", "TD700A", "TD700B", "TD700C", "NO2", "SPEC_SI"],
}
FLAGPOS = {v: (g, i, len(vs)) for g, vs in QGROUPS.items() for i, v in enumerate(vs)}

# column -> (tidy variable, units)
BOTTLE_MAP = {
    "CTDTMP": ("temp", "degC (ITS-90, CTD in-situ)"),
    "THETA": ("theta", "degC (ITS-90, potential, CTD)"),
    "CTDSAL": ("salt", "PSS-78 (CTD)"),
    "OXYGEN": ("O2", "umol/kg"),
    "CTDOXY": ("O2_CTD", "umol/kg"),
    "NO2_NO3": ("NO3", "umol/kg (NO3+NO2)"),
    "NO2": ("NO2", "nmol/kg"),
    "LLN": ("NO3_LL", "nmol/kg (low-level NO3+NO2)"),
    "PHSPHT": ("PO4", "umol/kg"),
    "LLP": ("PO4_LL", "nmol/kg (low-level SRP)"),
    "SILCAT": ("SiO2", "umol/kg"),
    "DIC": ("DIC", "umol/kg"),
    "ALKALIN": ("ALK", "ueq/kg"),
    "PH": ("pH", "pH total scale @25C (TOT25; early data may be NBS25)"),
    "pCO2": ("pCO2", "uatm"),
    "DOC": ("DOC", "umol/kg"),
    "DON": ("DON", "umol/kg"),
    "DOP": ("DOP", "umol/kg"),
    "TDN": ("TDN", "umol/kg"),
    "TDP": ("TDP", "umol/kg"),
    "PC": ("POC", "umol/kg (HOT PC, total particulate C)"),
    "PN": ("PON", "umol/kg"),
    "PP": ("POP", "nmol/kg"),
    "PIC": ("PIC", "umol/kg"),
    "PSi": ("bSi", "nmol/kg"),
    "N2O": ("N2O", "nmol/kg"),
    "CHL_A": ("Chl", "ug/L (fluorometric)"),
    "HPLCchl": ("Chl_HPLC", "ng/L"),
    "H_BACT": ("BactAbund", "1e5 cells/mL (heterotrophic bacteria)"),
    "P_BACT": ("ProAbund", "1e5 cells/mL (Prochlorococcus)"),
    "S_BACT": ("SynAbund", "1e5 cells/mL (Synechococcus)"),
    "E_BACT": ("PicoEukAbund", "1e5 cells/mL (picoeukaryotes)"),
}
for _p in ["MV_CHLA", "DV_CHLA", "CHL_B", "CHL_C3", "CHLC1_2", "CHL_Plus", "CHL_C4",
           "CHLDA_A", "PERID", "BUT_19", "FUCO", "HEX_19", "PRASINO", "DIADINO",
           "ZEAXAN", "VIOL", "LUTEIN", "A_CAR", "B_CAR", "CAROTEN"]:
    BOTTLE_MAP[_p] = ("HPLC_" + _p, "ng/L")
KEEP_FLAGS = {"1", "2"}


def _flag_digit(series, pos, n):
    s = series.astype("Int64").astype(str)
    s = s.where(s != "<NA>", "")
    return s.map(lambda x: x[pos] if len(x) == n else "")


def read_bottle(raw):
    files = sorted(glob.glob(os.path.join(raw, "bcodmo_bottle", "*.csv")))
    if not files:
        return pd.DataFrame(columns=COLS)
    parts = []
    for f in files:
        try:
            d = pd.read_csv(f, skiprows=[1], low_memory=False, on_bad_lines="skip")
        except Exception as e:  # e.g. an HTML error page
            print(f"  skip {f}: {e}")
            continue
        if "QUAL10" not in d.columns:
            continue
        d = d[pd.to_numeric(d["QUAL10"], errors="coerce").notna()]  # drop truncated rows
        parts.append(d)
    if not parts:
        return pd.DataFrame(columns=COLS)
    d = pd.concat(parts, ignore_index=True)
    d = d[pd.to_numeric(d["STNNBR"], errors="coerce") == 2]
    d = d.drop_duplicates(subset=["EXPOCODE", "CASTNO", "ROSETTE", "CTDPRS", "time"])
    t = _iso(d["time"])
    lat = pd.to_numeric(d["latitude"], errors="coerce")
    lon = pd.to_numeric(d["longitude"], errors="coerce")
    z = p2z(pd.to_numeric(d["CTDPRS"], errors="coerce"), LAT0)
    for g in QGROUPS:
        d[g] = pd.to_numeric(d[g], errors="coerce")
    out = []
    for col, (var, units) in BOTTLE_MAP.items():
        if col not in d.columns:
            continue
        v = pd.to_numeric(d[col], errors="coerce")
        if col in FLAGPOS:
            g, i, n = FLAGPOS[col]
            q = _flag_digit(d[g], i, n)
            ok = v.notna() & q.isin(KEEP_FLAGS)
        else:
            q = pd.Series("", index=d.index)
            ok = v.notna()
        out.append(pd.DataFrame({"time": t[ok], "depth_m": z[ok.values], "variable": var,
                                 "value": v[ok], "units": units, "source": "bcodmo_bottle",
                                 "qc": q[ok], "lat": lat[ok], "lon": lon[ok]}))
    return pd.concat(out, ignore_index=True)


# ------------------------------------------------------------------ primary production
def read_pp(raw):
    files = sorted(glob.glob(os.path.join(raw, "bcodmo_pp", "*.csv")))
    if not files:
        return pd.DataFrame(columns=COLS)
    d = pd.read_csv(files[0], skiprows=[1], dtype={"Flag": str, "Date": str, "Start_time": str})
    # Times are local HST (UTC-10) despite the 'Z' BCO-DMO appended to `time`;
    # use incubation start (dawn) and convert to UTC. Missing start time -> 06:00 HST.
    date = d["Date"].str.replace(r"\.0$", "", regex=True).str.zfill(6)
    st = d["Start_time"].str.replace(r"\.0$", "", regex=True)
    st = st.where(st.notna() & (st != "nan"), "600").str.zfill(4)
    local = pd.to_datetime(date + st, format="%y%m%d%H%M", errors="coerce")
    t = _iso((local + pd.Timedelta(hours=10)).dt.tz_localize("UTC"))
    flag = d["Flag"].fillna("").str.strip()
    ok_bottle = flag.str[0].isin(KEEP_FLAGS)

    def fl(i):
        return flag.str[i].fillna("")

    out = []
    light = []
    for k, i in zip([1, 2, 3], [3, 4, 5]):
        v = pd.to_numeric(d[f"Light_rep{k}"], errors="coerce")
        light.append(v.where(fl(i).isin(KEEP_FLAGS) & ok_bottle))
    L = pd.concat(light, axis=1)
    pp = L.mean(axis=1)
    nrep = L.notna().sum(axis=1)
    ok = pp.notna()
    qc = flag.str[3:6] + "|inc=" + d["Incubation_type"].fillna("") + "|n=" + nrep.astype(str)
    out.append(pd.DataFrame({"time": t[ok], "depth_m": d["depth"][ok], "variable": "PP",
                             "value": pp[ok],
                             "units": "mg C/m3 per dawn-dusk incubation (~12-15 h; light mean, dark NOT subtracted)",
                             "source": "bcodmo_pp", "qc": qc[ok], "lat": LAT0, "lon": LON0}))
    dark = []
    for k, i in zip([1, 2, 3], [6, 7, 8]):
        v = pd.to_numeric(d[f"Dark_rep{k}"], errors="coerce")
        dark.append(v.where(fl(i).isin(KEEP_FLAGS) & ok_bottle))
    D = pd.concat(dark, axis=1).mean(axis=1)
    ok = D.notna()
    out.append(pd.DataFrame({"time": t[ok], "depth_m": d["depth"][ok], "variable": "PP_dark",
                             "value": D[ok], "units": "mg C/m3 per incubation (dark bottle)",
                             "source": "bcodmo_pp", "qc": flag.str[6:9][ok], "lat": LAT0, "lon": LON0}))
    return pd.concat(out, ignore_index=True)


# ------------------------------------------------------------------ sediment traps
FLUX_MAP = {"Carbon": ("POC_flux", "mg C/m2/d (HOT trap particulate C)"),
            "Nitrogen": ("PON_flux", "mg N/m2/d"),
            "Phosphorus": ("POP_flux", "mg P/m2/d"),
            "Silica": ("bSi_flux", "mg Si/m2/d"),
            "Mass": ("mass_flux", "mg/m2/d"),
            "PIC": ("PIC_flux", "mg C/m2/d")}


def _cruise_dates(raw):
    """HOT cruise number -> mid-cruise date, from HOT_surface_CO2.txt."""
    f = os.path.join(raw, "hot_co2", "HOT_surface_CO2.txt")
    if not os.path.exists(f):
        return {}
    d = _read_hot_co2_table(f)
    return dict(zip(d["cruise"].astype(int), pd.to_datetime(d["date"], format="%d-%b-%y")))


def read_flux(raw):
    files = sorted(glob.glob(os.path.join(raw, "bcodmo_flux", "*.csv")))
    if not files:
        return pd.DataFrame(columns=COLS)
    d = pd.read_csv(files[0], skiprows=[1])
    cd = _cruise_dates(raw)
    t = pd.to_datetime(d["Cruise"].map(cd))
    t = _iso(t.dt.tz_localize("UTC"))
    out = []
    for col, (var, units) in FLUX_MAP.items():
        v = pd.to_numeric(d[col], errors="coerce")
        ok = v.notna() & t.notna()
        out.append(pd.DataFrame({"time": t[ok], "depth_m": d["depth"][ok], "variable": var,
                                 "value": v[ok], "units": units, "source": "bcodmo_flux",
                                 "qc": "trt=" + d["Treatment"][ok].astype(str) + "|n=" + d[col + "_n"][ok].astype(str),
                                 "lat": LAT0, "lon": LON0}))
    return pd.concat(out, ignore_index=True)


# ------------------------------------------------------------------ HOT surface CO2 (Dore)
def _read_hot_co2_table(f):
    with open(f, encoding="latin-1") as fh:
        lines = fh.read().splitlines()
    i = next(k for k, l in enumerate(lines) if l.startswith("cruise\t"))
    rows = [l.split("\t") for l in lines[i:] if l.strip()]
    hdr = rows[0]
    rows = [r + [""] * (len(hdr) - len(r)) for r in rows[1:]]
    return pd.DataFrame([r[:len(hdr)] for r in rows], columns=hdr)


HOTCO2_MAP = {"temp": ("temp", "degC (in-situ, mean 0-30 dbar)"),
              "sal": ("salt", "PSS-78 (mean 0-30 dbar)"),
              "phos": ("PO4", "umol/kg"),
              "sil": ("SiO2", "umol/kg"),
              "DIC": ("DIC", "umol/kg"),
              "TA": ("ALK", "ueq/kg"),
              "pHmeas_insitu": ("pH", "pH total scale, measured (25C) adjusted to in-situ T"),
              "pCO2calc_insitu": ("pCO2", "uatm (calculated from DIC+TA, CO2SYS)")}


def read_hot_co2(raw):
    f = os.path.join(raw, "hot_co2", "HOT_surface_CO2.txt")
    if not os.path.exists(f):
        return pd.DataFrame(columns=COLS)
    d = _read_hot_co2_table(f)
    t = _iso(pd.to_datetime(d["date"], format="%d-%b-%y", errors="coerce"))
    out = []
    for col, (var, units) in HOTCO2_MAP.items():
        v = pd.to_numeric(d[col], errors="coerce")
        ok = v.notna() & (v > -998)
        out.append(pd.DataFrame({"time": t[ok], "depth_m": 15.0, "variable": var, "value": v[ok],
                                 "units": units, "source": "hot_co2",
                                 "qc": "notes=" + d["notes"][ok].astype(str),
                                 "lat": LAT0, "lon": LON0}))
    return pd.concat(out, ignore_index=True)


# ------------------------------------------------------------------ WHOTS MAPCO2 mooring
def _norm(s):
    return " ".join(s.lower().replace('"', "").split())


def read_whots(raw):
    files = sorted(f for f in glob.glob(os.path.join(raw, "ncei_whots_pco2", "WHOTS_*.csv"))
                   if "qflog" not in f.lower())
    out = []
    for f in files:
        with open(f, encoding="latin-1", newline=None) as fh:
            lines = fh.read().splitlines()
        hi = next((k for k, l in enumerate(lines) if "Date" in l and "Time" in l and "Lat" in l), None)
        if hi is None:
            continue
        import csv
        rows = list(csv.reader(lines[hi:]))
        hdr = [_norm(h) for h in rows[0]]
        body = [r for r in rows[1:] if len(r) >= 6 and r[3][:1].isdigit()]
        if not body:
            continue
        n = len(hdr)
        body = [(r + [""] * n)[:n] for r in body]
        d = pd.DataFrame(body, columns=range(n))
        dt = d[3].str.strip() + " " + d[4].str.strip()
        t = pd.to_datetime(dt, format="%m/%d/%Y %H:%M", errors="coerce")
        t2 = pd.to_datetime(dt, format="%m/%d/%y %H:%M", errors="coerce")
        t = t.fillna(t2)
        t = _iso(t.dt.tz_localize("UTC"))
        lat = pd.to_numeric(d[1], errors="coerce")
        lon = pd.to_numeric(d[2], errors="coerce")

        def find(*keys, start=0):
            for k in keys:
                for j in range(start, n):
                    if hdr[j].startswith(k):
                        return j
            return None

        jsw = find("xco2 sw (wet)", "xco2_sw_wet")
        qsw = jsw + 1 if jsw is not None else None  # QF follows xCO2 SW (wet)
        jair = find("xco2 air (wet)", "xco2_air_wet")
        qair = jair + 1 if jair is not None else None
        spec = [  # (col idx, flag idx, var, units, depth)
            (find("pco2 sw (sat)"), qsw, "pCO2", "uatm (MAPCO2, sat. at SST)", 0.5),
            (find("fco2 sw (sat)", "fco2_sw_sat"), qsw, "fCO2", "uatm (MAPCO2)", 0.5),
            (find("pco2 air (sat)"), qair, "pCO2_air", "uatm (MAPCO2 air, sat.)", 0.0),
            (find("sst"), None, "temp", "degC (mooring SST)", 0.5),
            (find("salinity", "sss"), None, "salt", "PSS-78 (mooring SSS)", 0.5),
            (find("ph (total scale)"), find("ph qf"), "pH", "pH total scale (SAMI-pH)", 0.5),
            (find("chl"), find("chl qf"), "Chl", "ug/L (fluorometer)", 0.5),
            (find("doxy"), find("doxy qf"), "O2", "umol/kg (optode)", 0.5),
        ]
        for j, q, var, units, z in spec:
            if j is None or (var == "Chl" and "qf" in hdr[j]):
                continue
            v = pd.to_numeric(d[j], errors="coerce")
            ok = v.notna() & (v > -900)
            if q is not None:
                qq = pd.to_numeric(d[q], errors="coerce")
                qs = qq.astype("Int64").astype(str)
                ok &= (qq == 2)
            else:
                qs = pd.Series("", index=d.index)
            out.append(pd.DataFrame({"time": t[ok], "depth_m": z, "variable": var, "value": v[ok],
                                     "units": units, "source": "ncei_whots_pco2", "qc": qs[ok],
                                     "lat": lat[ok], "lon": lon[ok]}))
    if not out:
        return pd.DataFrame(columns=COLS)
    return pd.concat(out, ignore_index=True)


# ------------------------------------------------------------------ Scripps (Keeling) HAWI
def read_scripps(raw):
    f = os.path.join(raw, "scripps_co2", "HAWI.csv")
    if not os.path.exists(f):
        return pd.DataFrame(columns=COLS)
    with open(f, encoding="utf-8-sig", newline=None) as fh:
        lines = fh.read().splitlines()
    i = next(k for k, l in enumerate(lines) if l.startswith("DATE,"))
    import io
    d = pd.read_csv(io.StringIO("\n".join(lines[i:])), dtype={"DATE": str})
    t = _iso(pd.to_datetime(d["DATE"].str.strip(), format="%Y%m%d", errors="coerce"))
    lon = -pd.to_numeric(d["LONG"], errors="coerce")  # file gives degrees West as positive
    out = []
    for col, var, units in [("AVGDIC", "DIC", "umol/kg"), ("AVGALK", "ALK", "umol/kg"),
                            ("TEMP", "temp", "degC (CTD)"), ("SAL", "salt", "PSS-78")]:
        v = pd.to_numeric(d[col], errors="coerce")
        ok = v.notna()
        out.append(pd.DataFrame({"time": t[ok], "depth_m": pd.to_numeric(d["Z"], errors="coerce")[ok],
                                 "variable": var, "value": v[ok], "units": units,
                                 "source": "scripps_co2", "qc": "", "lat": d["LAT"][ok], "lon": lon[ok]}))
    return pd.concat(out, ignore_index=True)


# ------------------------------------------------------------------ gap fill (added 2026-10-08)
# Newer BCO-DMO releases extend PP/flux to HOT-355 (Dec 2024); SPOTS extends bottle chemistry to
# HOT-317 (Dec 2019). De-duplication against the existing v1 sources is BY CRUISE: a new source only
# contributes cruises that are absent from the corresponding v1 file, so v1 rows are never replaced.
def _v1_cruises(raw, sub, col):
    files = sorted(glob.glob(os.path.join(raw, sub, "*.csv")))
    out = set()
    for f in files:
        try:
            d = pd.read_csv(f, skiprows=[1], usecols=[col], low_memory=False, on_bad_lines="skip")
        except Exception:
            continue
        out |= set(pd.to_numeric(d[col], errors="coerce").dropna().astype(int))
    return out


def read_pp_v5(raw):
    """BCO-DMO 737163 v5 (HOT-1..355). Only cruises missing from bcodmo_pp (v1) are used."""
    files = sorted(glob.glob(os.path.join(raw, "bcodmo_pp_v5", "737163_v5*.csv")))
    if not files:
        return pd.DataFrame(columns=COLS)
    d = pd.read_csv(files[0], dtype=str)
    old = _v1_cruises(raw, "bcodmo_pp", "Cruise")
    cr = pd.to_numeric(d["Cruise_num"], errors="coerce")
    d = d[cr.notna() & ~cr.isin(old)].copy()
    # v5 has one flag per light set (Flag_Light) and per dark set (Flag_Dark), plus Flag_Bottle;
    # same HOT scheme (1 not QC'd, 2 good, 3 suspect, 4 bad, 5 missing, 9 not measured). Keep 1/2.
    okb = d["Flag_Bottle"].str.strip().isin(KEEP_FLAGS)
    okl = d["Flag_Light"].str.strip().isin(KEEP_FLAGS)
    okd = d["Flag_Dark"].str.strip().isin(KEEP_FLAGS)
    t = _iso(d["Start_ISO_DateTime_UTC"])  # true UTC incubation start in v5
    lat = pd.to_numeric(d["Latitude"], errors="coerce")
    lon = pd.to_numeric(d["Longitude"], errors="coerce")
    z = pd.to_numeric(d["Depth"], errors="coerce")
    L = pd.concat([pd.to_numeric(d[f"Light_rep{k}"], errors="coerce") for k in (1, 2, 3)], axis=1)
    pp = L.mean(axis=1).where(okb & okl)
    nrep = L.notna().sum(axis=1)
    ok = pp.notna()
    qc = (d["Flag_Light"].fillna("") + "|inc=" + d["Incubation_type"].fillna("") + "|n=" + nrep.astype(str))
    out = [pd.DataFrame({"time": t[ok], "depth_m": z[ok], "variable": "PP", "value": pp[ok],
                         "units": "mg C/m3 per dawn-dusk incubation (~12-15 h; light mean, dark NOT subtracted)",
                         "source": "bcodmo_pp_v5", "qc": qc[ok], "lat": lat[ok], "lon": lon[ok]})]
    D = pd.concat([pd.to_numeric(d[f"Dark_rep{k}"], errors="coerce") for k in (1, 2, 3)], axis=1).mean(axis=1)
    D = D.where(okb & okd)
    ok = D.notna()
    out.append(pd.DataFrame({"time": t[ok], "depth_m": z[ok], "variable": "PP_dark", "value": D[ok],
                             "units": "mg C/m3 per incubation (dark bottle)", "source": "bcodmo_pp_v5",
                             "qc": d["Flag_Dark"][ok], "lat": lat[ok], "lon": lon[ok]}))
    return pd.concat(out, ignore_index=True)


FLUX_V4_MAP = {"Carbon_flux": ("POC_flux", "mg C/m2/d (HOT trap particulate C)", "Carbon_numreps"),
               "Nitrogen_flux": ("PON_flux", "mg N/m2/d", "Nitrogen_numreps"),
               "Phosphorus_flux": ("POP_flux", "mg P/m2/d", "Phosphorus_numreps"),
               "Silica_flux": ("bSi_flux", "mg Si/m2/d", "Silica_numreps"),
               "Mass_flux": ("mass_flux", "mg/m2/d", "Mass_numreps"),
               "PIC_flux": ("PIC_flux", "mg C/m2/d", "PIC_numreps")}


def read_flux_v4(raw):
    """BCO-DMO 737393 v4 (HOT-2..355). Only cruises missing from bcodmo_flux (v1) are used.
    v4 carries real deployment start/end (UTC); time = deployment midpoint."""
    files = sorted(glob.glob(os.path.join(raw, "bcodmo_flux_v4", "737393_v4*.csv")))
    if not files:
        return pd.DataFrame(columns=COLS)
    d = pd.read_csv(files[0], dtype=str)
    old = _v1_cruises(raw, "bcodmo_flux", "Cruise")
    cr = pd.to_numeric(d["Cruise_num"], errors="coerce")
    d = d[cr.notna() & ~cr.isin(old)].copy()
    t0 = pd.to_datetime(d["Start_ISO_DateTime_UTC"], utc=True, errors="coerce")
    t1 = pd.to_datetime(d["End_ISO_DateTime_UTC"], utc=True, errors="coerce")
    t = _iso(t0 + (t1 - t0) / 2)
    z = pd.to_numeric(d["Depth"], errors="coerce")
    lat = pd.to_numeric(d["Latitude"], errors="coerce")
    lon = pd.to_numeric(d["Longitude"], errors="coerce")
    out = []
    for col, (var, units, ncol) in FLUX_V4_MAP.items():
        v = pd.to_numeric(d[col], errors="coerce")
        ok = v.notna() & t.notna()
        out.append(pd.DataFrame({"time": t[ok], "depth_m": z[ok], "variable": var, "value": v[ok],
                                 "units": units, "source": "bcodmo_flux_v4",
                                 "qc": "trt=" + d["Treatment"][ok].astype(str) + "|n=" + d[ncol][ok].fillna("").astype(str),
                                 "lat": lat[ok], "lon": lon[ok]}))
    return pd.concat(out, ignore_index=True)


# SPOTS (BCO-DMO 896862 v2) column -> (tidy variable, units, WOCE flag column or None)
SPOTS_MAP = {
    "CTDTMP": ("temp", "degC (ITS-90, CTD in-situ)", None),
    "CTDSAL": ("salt", "PSS-78 (CTD)", "CTDSAL_FLAG_W"),
    "OXYGEN": ("O2", "umol/kg", "OXYGEN_FLAG_W"),
    "CTDOXY": ("O2_CTD", "umol/kg", "CTDOXY_FLAG_W"),
    "NITRAT": ("NO3", "umol/kg (HOT NO3+NO2; SPOTS NITRAT)", "NITRAT_FLAG_W"),
    "NITRIT": ("NO2", "umol/kg", "NITRIT_FLAG_W"),
    "PHSPHT": ("PO4", "umol/kg", "PHSPHT_FLAG_W"),
    "SILCAT": ("SiO2", "umol/kg", "SILCAT_FLAG_W"),
    "NH4": ("NH4", "umol/kg", "NH4_FLAG_W"),
    "TCARBN": ("DIC", "umol/kg", "TCARBN_FLAG_W"),
    "ALKALI": ("ALK", "umol/kg", "ALKALI_FLAG_W"),
    "PH_TOT": ("pH", "pH total scale @25C (SPOTS PH_TOT, PH_TMP=25)", "PH_TOT_FLAG_W"),
    "DOC": ("DOC", "umol/kg", "DOC_FLAG_W"),
    "TPC": ("POC", "umol/kg (HOT PC, total particulate C; SPOTS TPC)", "TPC_FLAG_W"),
    "TPN": ("PON", "umol/kg (SPOTS TPN)", "TPN_FLAG_W"),
    "TPP": ("POP", "umol/kg (SPOTS TPP)", "TPP_FLAG_W"),
}
SPOTS_KEEP = {2, 6}  # WOCE: 2 acceptable, 6 mean of replicates


def read_spots(raw):
    """SPOTS v2 ALOHA subset. STNNBR in SPOTS = HOT cruise number. Only cruises missing from
    bcodmo_bottle (3773 v1) are used (in practice HOT-289..317, Dec 2016-Dec 2019)."""
    files = sorted(glob.glob(os.path.join(raw, "bcodmo_spots", "spots_*.csv")))
    if not files:
        return pd.DataFrame(columns=COLS)
    d = pd.concat([pd.read_csv(f, skiprows=[1], dtype=str) for f in files], ignore_index=True)
    d = d[d["TimeSeriesSite"].str.strip() == "ALOHA"]
    old = _v1_cruises(raw, "bcodmo_bottle", "cruise_name")
    cr = pd.to_numeric(d["STNNBR"], errors="coerce")
    d = d[cr.notna() & ~cr.isin(old)].copy()
    d = d.drop_duplicates(subset=["CRUISE", "STNNBR", "CASTNO", "BTLNBR", "CTDPRS", "DATE", "TIME"])
    tt = pd.to_datetime(d["DATE"].str.strip() + d["TIME"].str.strip().str.zfill(4),
                        format="%Y%m%d%H%M", errors="coerce")
    # SPOTS bug for January HOT cruises: the HOT mmddyy date lost its leading zero, e.g. 012317
    # (23 Jan 2017) -> "12317" -> 2017-12-03. Repair: month 1, day = <last digit of wrong month><wrong day>,
    # accepted only if within 10 days of the HOT cruise date (HOT_surface_CO2.txt); otherwise dropped.
    cd = _cruise_dates(raw)
    if cd:
        ref = pd.to_datetime(pd.to_numeric(d["STNNBR"], errors="coerce").map(cd))
        bad = ref.notna() & ((tt - ref).abs() > pd.Timedelta(days=20))
        if bad.any():
            wm, wd = tt[bad].dt.month, tt[bad].dt.day
            day = ((wm % 10).astype(str) + wd.astype(str)).astype(int)
            day = day.where((wm >= 10) & (day >= 1) & (day <= 31))
            fixed = pd.to_datetime(dict(year=tt[bad].dt.year, month=1, day=day,
                                        hour=tt[bad].dt.hour, minute=tt[bad].dt.minute), errors="coerce")
            fixed = fixed.where((fixed - ref[bad]).abs() <= pd.Timedelta(days=10))
            tt.loc[bad] = fixed
            print(f"  [std_hot] SPOTS: repaired {int(fixed.notna().sum())} January dates, "
                  f"dropped {int(fixed.isna().sum())}")
    t = _iso(tt)
    lat = pd.to_numeric(d["LATITUDE"], errors="coerce")
    lon = pd.to_numeric(d["longitude"], errors="coerce")
    p = pd.to_numeric(d["CTDPRS"], errors="coerce")
    z = pd.Series(p2z(p, LAT0), index=d.index)
    out = []
    for col, (var, units, fcol) in SPOTS_MAP.items():
        v = pd.to_numeric(d[col], errors="coerce")
        ok = v.notna() & (v > -998) & p.notna() & (p > -998)
        if fcol is not None:
            f = pd.to_numeric(d[fcol], errors="coerce")
            ok &= f.isin(SPOTS_KEEP)
            q = f.astype("Int64").astype(str)
        else:
            q = pd.Series("", index=d.index)
        out.append(pd.DataFrame({"time": t[ok], "depth_m": z[ok], "variable": var, "value": v[ok],
                                 "units": units, "source": "bcodmo_spots", "qc": q[ok],
                                 "lat": lat[ok], "lon": lon[ok]}))
    return pd.concat(out, ignore_index=True)


def load(site, site_dir):
    if site != "HOT":
        return pd.DataFrame(columns=COLS)
    raw = os.path.join(site_dir, "raw")
    parts = []
    for fn in (read_bottle, read_pp, read_flux, read_hot_co2, read_whots, read_scripps,
               read_pp_v5, read_flux_v4, read_spots):
        try:
            d = fn(raw)
        except Exception as e:
            print(f"  [std_hot] {fn.__name__} failed: {e!r}")
            continue
        print(f"  [std_hot] {fn.__name__}: {len(d)} rows")
        parts.append(d)
    parts = [p for p in parts if len(p)]
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=COLS)
