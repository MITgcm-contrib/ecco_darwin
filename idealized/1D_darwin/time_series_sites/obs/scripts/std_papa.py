"""Reader for Ocean Station Papa / Line P P26 raw downloads -> tidy rows.

Sources (raw dirs under <site_dir>/raw/):
  ios_bot_erddap/          DFO-IOS bottle profiles (CIOOS Pacific ERDDAP IOS_BOT_Profiles),
                           box 49.9-50.3N, 145.2-144.6W (P26 = 50.0N 145.0W)
  ncei_linep_carbon/       Line P DIC/TA (NCEI OCADS 0234342 synthesis 1990-2019, 0300980 2021,
                           0310718 2024); only DIC and ALK taken (nutrients/O2 come from IOS)
  pmel_co2_mooring/        PMEL Papa CO2/OA mooring, 3-hourly, 0.5 m (PMEL ERDDAP)
  oceansites_pmel_daily/   PMEL OCS Papa daily gridded T/S (OceanSITES OS_PAPA_*_M_TSVM_*_dy.nc)
  ooi_papa_flanking/       OOI Global Station Papa flanking moorings A/B riser O2, pH, chl
                           (daily means of QARTOD qc_agg<=2 samples, ~30-40 m)
  pangaea_osp_trap/        Wong et al. OSP sediment traps (PANGAEA.92552), annual-mean fluxes

Flag conventions (rows failing are dropped):
  IOS_BOT ERDDAP: no flag channels exposed; values used as archived by IOS (qc='').
  Line P carbon: WOCE bottle flags; keep 2 (good) and 6 (replicate mean); drop 3,4,9 and -999.
  PMEL CO2 mooring: ERDDAP file already final QC'd (only flag-2 data published); NaN dropped.
  OceanSITES PMEL daily: OceanSITES <VAR>_QC; keep 1 (good) and 2 (probably good); drop 0 (no QC),
      3/4 (bad), 7-9; also drop records with fill TIME (<2007).
  OOI: server-side filter qc_agg in {1 pass, 2 not evaluated}; suspect(3)/fail(4) excluded.
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

BOX = dict(lat0=49.9, lat1=50.3, lon0=-145.2, lon1=-144.6)


def _rows(time, depth, var, value, units, source, qc="", lat=np.nan, lon=np.nan):
    d = pd.DataFrame({"time": time, "depth_m": depth, "variable": var, "value": value,
                      "units": units, "source": source, "qc": qc, "lat": lat, "lon": lon})
    d["value"] = pd.to_numeric(d["value"], errors="coerce")
    return d.dropna(subset=["value"])


# ---------------------------------------------------------------- IOS bottle
def _ios(raw):
    f = os.path.join(raw, "ios_bot_erddap", "IOS_BOT_Profiles_P26.csv")
    if not os.path.exists(f):
        return []
    hdr = pd.read_csv(f, nrows=1)  # row 0 after header = units
    units = hdr.iloc[0].to_dict()
    d = pd.read_csv(f, skiprows=[1], low_memory=False)
    d = d[(d.latitude.between(BOX["lat0"], BOX["lat1"])) &
          (d.longitude.between(BOX["lon0"], BOX["lon1"]))]
    t = iso(d["time"])
    z = d["depth"].where(d["depth"].notna(), d.get("PRESPR01"))
    src = "ios_bot_erddap"
    out = []

    def first(cols):
        v = pd.Series(np.nan, index=d.index)
        u = pd.Series("", index=d.index, dtype=object)
        for c in cols:
            if c in d.columns:
                m = v.isna() & d[c].notna()
                v[m] = d.loc[m, c]
                u[m] = str(units.get(c, ""))
        return v, u

    tcols = ["TEMPS901", "TEMPS902", "TEMPST01", "TEMPS601", "TEMPS602", "TEMPRTN1",
             "sea_water_temperature"]
    v, u = first(tcols)
    out.append(_rows(t, z, "temp", v, u.replace({"degC": "degC (in-situ)"}), src,
                     lat=d.latitude, lon=d.longitude))
    scols = ["PSALBST01", "PSALST01", "PSALST02", "SSALBST01", "SSALST01",
             "sea_water_practical_salinity"]
    v, u = first(scols)
    out.append(_rows(t, z, "salt", v, u, src, lat=d.latitude, lon=d.longitude))
    v, u = first(["DOXMZZ01"])
    m = v.notna()
    out.append(_rows(t[m], z[m], "O2", v[m], u[m], src, lat=d.latitude[m], lon=d.longitude[m]))
    # oxygen only reported in mL/L for some samples -> keep as mL/L
    if "DOXYZZ01" in d.columns:
        m2 = (~m) & d["DOXYZZ01"].notna()
        out.append(_rows(t[m2], z[m2], "O2", d.loc[m2, "DOXYZZ01"], units.get("DOXYZZ01", "mL/L"),
                         src, lat=d.latitude[m2], lon=d.longitude[m2]))
    for c, var in [("NTRZAAZ1", "NO3"), ("SLCAAAZ1", "SiO2"), ("PHOSAAZ1", "PO4"),
                   ("CPHLFLP1", "Chl")]:
        if c in d.columns:
            out.append(_rows(t, z, var, d[c], units.get(c, ""), src,
                             lat=d.latitude, lon=d.longitude))
    return out


# ---------------------------------------------------------------- Line P carbon
CARB_GOOD = {2, 6}


def _carbon_frame(d, latc, lonc, timec, zc, src):
    d = d.assign(_t=timec.values)
    d = d[(pd.to_numeric(d[latc], errors="coerce").between(BOX["lat0"], BOX["lat1"])) &
          (pd.to_numeric(d[lonc], errors="coerce").between(BOX["lon0"], BOX["lon1"]))]
    z = pd.to_numeric(d[zc], errors="coerce").where(lambda x: x >= 0)
    if "CTDPRS_DBAR" in d.columns:  # fall back to pressure (dbar ~ m) when depth is -999
        z = z.fillna(pd.to_numeric(d["CTDPRS_DBAR"], errors="coerce").where(lambda x: x >= 0))
    d = d.assign(_z=z)
    out = []
    for c, fc, var in [("DIC_UMOL_KG", "DIC_FLAG_W", "DIC"), ("TA_UMOL_KG", "TA_FLAG_W", "ALK")]:
        v = pd.to_numeric(d[c], errors="coerce")
        fl = pd.to_numeric(d[fc], errors="coerce")
        m = fl.isin(CARB_GOOD) & (v > 0)
        out.append(_rows(d.loc[m, "_t"], d.loc[m, "_z"], var,
                         v[m], "umol/kg", src, qc=fl[m].astype(int).astype(str),
                         lat=pd.to_numeric(d.loc[m, latc]), lon=pd.to_numeric(d.loc[m, lonc])))
    return out


def _linep_carbon(raw):
    out = []
    rdir = os.path.join(raw, "ncei_linep_carbon")
    src = "ncei_linep_carbon"
    f = os.path.join(rdir, "LineP_for_Data_Synthesis_1990-2019_v1.csv")
    if os.path.exists(f):
        d = pd.read_csv(f, low_memory=False)
        d = d.loc[:, [c for c in d.columns if not c.startswith("Unnamed")]]
        tt = pd.to_datetime(dict(year=d.YEAR_UTC, month=d.MONTH_UTC, day=d.DAY_UTC), errors="coerce")
        hm = d.TIME_UTC.astype(str).str.extract(r"(\d+):(\d+)").astype(float).fillna(0)
        tt = tt + pd.to_timedelta(hm[0], unit="h") + pd.to_timedelta(hm[1], unit="m")
        out += _carbon_frame(d, "LATITUDE_DEC", "LONGITUDE_DEC", iso(tt), "DEPTH_METER", src)
    for f, latc, lonc in [(os.path.join(rdir, "0300980_2021-008_data.xlsx"), "LATITUDE_DEC", "LONGITUDE_DEC"),
                          (os.path.join(rdir, "0310718_18DD20240124_Data.xlsx"), "LATITUDE", "LONGITUDE")]:
        if not os.path.exists(f):
            continue
        try:
            d = pd.read_excel(f)
        except Exception as e:  # openpyxl missing etc.
            print(f"  std_papa: cannot read {os.path.basename(f)}: {e}")
            continue
        tt = pd.to_datetime(dict(year=d.YEAR_UTC, month=d.MONTH_UTC, day=d.DAY_UTC), errors="coerce")
        hm = d.TIME_UTC.astype(str).str.extract(r"(\d+):(\d+)").astype(float).fillna(0)
        tt = tt + pd.to_timedelta(hm[0], unit="h") + pd.to_timedelta(hm[1], unit="m")
        out += _carbon_frame(d, latc, lonc, iso(tt), "DEPTH_M", src)
    return out


# ---------------------------------------------------------------- PMEL CO2 mooring
def _pmel_co2(raw):
    f = os.path.join(raw, "pmel_co2_mooring", "papa_co2_mooring.csv")
    if not os.path.exists(f):
        return []
    d = pd.read_csv(f, skiprows=[1])
    t = iso(d["time"])
    src = "pmel_co2_mooring"
    out = []
    for c, var, u in [("SST", "temp", "degC (in-situ)"), ("SSS", "salt", "PSU"),
                      ("pCO2_sw", "pCO2", "uatm"), ("pCO2_air", "pCO2_air", "uatm"),
                      ("pH_sw", "pH", "total scale"), ("DOXY", "O2", "umol/kg"),
                      ("CHL", "Chl", "ug/L (fluorescence, nighttime, cal. bias 2 applied)")]:
        if c in d.columns:
            z = 0.0 if var == "pCO2_air" else 0.5  # air: marine boundary layer, depth 0 by convention
            out.append(_rows(t, z, var, d[c], u, src, lat=d.latitude, lon=d.longitude))
    return out


# ---------------------------------------------------------------- OceanSITES daily T/S
def _oceansites_daily(raw):
    out = []
    files = sorted(glob.glob(os.path.join(raw, "oceansites_pmel_daily", "*.nc")))
    if not files:
        return out
    try:
        import netCDF4
    except Exception as e:
        print(f"  std_papa: netCDF4 unavailable: {e}")
        return out
    src = "oceansites_pmel_daily"
    for f in files:
        try:
            ds = netCDF4.Dataset(f)
        except Exception as e:
            print(f"  std_papa: cannot open {os.path.basename(f)}: {e}")
            continue
        with ds:
            tv = ds.variables["TIME"]
            times = pd.to_datetime(netCDF4.num2date(tv[:], tv.units, only_use_cftime_datetimes=False,
                                                    only_use_python_datetimes=True))
            for vname, var, u in [("TEMP", "temp", "degC (in-situ)"), ("PSAL", "salt", "PSU")]:
                if vname not in ds.variables:
                    continue
                v = ds.variables[vname]
                dname = [dn for dn in v.dimensions if dn.startswith("DEP")][0]
                dep = np.asarray(ds.variables[dname][:], dtype=float)
                a = np.ma.filled(np.squeeze(v[:]).astype(float), np.nan)
                qn = f"QUALITY_{vname}" if f"QUALITY_{vname}" in ds.variables else f"{vname}_QC"
                q = None
                if qn in ds.variables:
                    q = np.ma.filled(np.squeeze(ds.variables[qn][:]), 0).astype(int)
                tt = np.repeat(np.asarray(times), a.shape[1])
                zz = np.tile(dep, a.shape[0])
                vv = a.ravel()
                qq = q.ravel() if q is not None else np.zeros_like(vv, dtype=int)
                keep = (np.isfinite(vv) & (vv < 1e30) & np.isin(qq, [1, 2])
                        & (pd.DatetimeIndex(tt) >= pd.Timestamp("2007-01-01")))
                out.append(_rows(iso(pd.Series(tt[keep])), zz[keep], var, vv[keep], u, src,
                                 qc=pd.Series(qq[keep]).astype(str).values))
    return out


# ---------------------------------------------------------------- OOI flanking moorings
def _ooi(raw):
    out = []
    rdir = os.path.join(raw, "ooi_papa_flanking")
    if not os.path.isdir(rdir):
        return out
    src = "ooi_papa_flanking"
    for m in ["flma", "flmb"]:
        pf = os.path.join(rdir, f"ooi_gp03{m}_dostad_pressure_daily.csv")
        pres = None
        if os.path.exists(pf):
            p = pd.read_csv(pf, skiprows=[1])
            pres = p.set_index(p["time"].str[:10])["sea_water_pressure"]
        for fn, col, var, u in [
                (f"ooi_gp03{m}_dostad_O2_daily.csv", "moles_of_oxygen_per_unit_mass_in_sea_water", "O2", "umol/kg"),
                (f"ooi_gp03{m}_phsen_pH_daily.csv", "sea_water_ph_reported_on_total_scale", "pH", "total scale"),
                (f"ooi_gp03{m}_flort_chl_daily.csv", "mass_concentration_of_chlorophyll_a_in_sea_water", "Chl",
                 "ug/L (fluorescence, factory cal)")]:
            f = os.path.join(rdir, fn)
            if not os.path.exists(f):
                continue
            d = pd.read_csv(f, skiprows=[1])
            day = d["time"].str[:10]
            z = day.map(pres) if pres is not None else pd.Series(np.nan, index=d.index)
            z = z.fillna(30.0)  # nominal riser instrument depth when no daily pressure
            out.append(_rows(iso(d["time"]), z.values, var, d[col], u, f"{src}_{m}",
                             qc="daily_mean_qc_agg<=2", lat=d.latitude, lon=d.longitude))
    return out


# ---------------------------------------------------------------- PANGAEA OSP traps
def _traps(raw):
    f = os.path.join(raw, "pangaea_osp_trap", "PANGAEA_92552.tab")
    if not os.path.exists(f):
        return []
    with open(f, encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()
    i = [k for k, s in enumerate(lines) if s.strip() == "*/"][0]
    from io import StringIO
    d = pd.read_csv(StringIO("\n".join(lines[i + 1:])), sep="\t")
    d = d[d["Duration [days]"] < 400]  # drop multi-year (1834-4253 d) composite means
    t = iso(d["Date/Time"])
    z = d["Depth water [m]"]
    src = "pangaea_osp_trap"
    q = "annual-mean flux; period start; duration_d=" + d["Duration [days]"].astype(str)
    out = []
    for c, var, u in [("Ann mass flux [g/m**2/a]", "mass_flux", "g/m2/yr"),
                      ("POC flux [mol/m**2/a]", "POC_flux", "mol/m2/yr"),
                      ("bSiO2 flux [mol/m**2/a]", "bSi_flux", "mol/m2/yr"),
                      ("PIC flux [mol/m**2/a]", "PIC_flux", "mol/m2/yr"),
                      ("PON flux [g/m**2/a]", "PON_flux", "g/m2/yr")]:
        out.append(_rows(t, z, var, d[c], u, src, qc=q, lat=50.0, lon=-145.0))
    return out


def load(site, site_dir):
    if site != "Papa":
        return pd.DataFrame()
    raw = os.path.join(site_dir, "raw")
    parts = []
    for fn in (_ios, _linep_carbon, _pmel_co2, _oceansites_daily, _ooi, _traps):
        try:
            parts += [p for p in fn(raw) if len(p)]
        except Exception as e:
            print(f"  std_papa: {fn.__name__} failed: {e!r}")
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
