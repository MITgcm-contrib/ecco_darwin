"""Global-mean forcing time series/climatology for apCO2 and iron dust
deposition, both read directly from the real delivered files (no
regridding -- both are read via pkg/darwin's own runtime horizontal
interpolation, so these plots are of the RAW native-grid files exactly as
data.darwin references them, not anything pre-processed by this repo).

apCO2 (real multi-decade daily time series): steph_input/v06/apCO2/apCO2_<year>,
one file per year, (ndays, 256, 2) big-endian float32, matching
data.darwin's pCO2_lon0/lat0/nlon/nlat/lat_inc exactly (confirmed:
ndays*256*2*4 bytes matches every file's size for both leap and non-leap
years; the 2 "longitude" points are bit-identical at every sample checked,
consistent with pCO2_nlon=2 being a no-real-zonal-structure trick for the
generic reader, not a real second field). Global mean = cos(latitude)-
area-weighted mean across the 256 latitude bands (NOT a plain arithmetic
mean, which would over-weight the poles).

Iron dust (a repeating 12-month climatology, NOT a real time series --
matches ironPeriod=-12 in data.darwin): steph_input/v06/forcing/
CESM-MIMI_..._144x96x12.bin, delivered. Confirmed via byte
count (144*96*12*4 = 663552, matching exactly) that the real grid is
144 lon x 96 lat -- data.darwin's iron_nlat was 95 (copied verbatim from
ecco_darwin's own 1-degree config) and has been corrected to 96 as part of
this check (see data.darwin's own comment for the full derivation). Global
mean here = cos(latitude)-area-weighted mean across the 96 latitude bands,
using iron_lat0=-90/iron_lat_inc=1.894736842105263 (uniform spacing,
unlike apCO2's irregular bands).
"""
import glob
import os
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

APCO2_DIR = "steph_input/v06/apCO2"
IRON_PATH = ("steph_input/v06/forcing/"
             "CESM-MIMI_1980-2015_CAM4-6MEAN_MonthlyDep_Hamiltonetal2020_clim.R4_144x96x12.bin")
OUT_PDF = "figures/forcing_timeseries.pdf"

# Verbatim from input_darwin_offline/data.darwin's pCO2_lat0/pCO2_lat_inc.
PCO2_LAT0 = -89.4628220
PCO2_LAT_INC = np.array(
    [0.6958694, 0.6999817, 0.7009048, 0.7012634, 0.7014313]
    + [0.7017418] * 245
    + [0.7014313, 0.7012634, 0.7009048, 0.6999817, 0.6958694]
)
assert PCO2_LAT_INC.size == 255, PCO2_LAT_INC.size  # 256 points -> 255 steps

# Verbatim from input_darwin_offline/data.darwin's iron_lat0/iron_lat_inc
# (post-fix: 96 points, 95 uniform increments spanning -90 to +90 exactly).
IRON_LAT0 = -90.0
IRON_LAT_INC = 1.894736842105263
IRON_NLAT = 96


def pco2_latitudes():
    lat = np.empty(256)
    lat[0] = PCO2_LAT0
    lat[1:] = PCO2_LAT0 + np.cumsum(PCO2_LAT_INC)
    return lat


def iron_latitudes():
    return IRON_LAT0 + IRON_LAT_INC * np.arange(IRON_NLAT)


def days_in_year(year):
    return 366 if (year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)) else 365


def load_apco2_year(year, weights):
    path = os.path.join(APCO2_DIR, f"apCO2_{year}")
    ndays = days_in_year(year)
    expected_bytes = ndays * 256 * 2 * 4
    actual_bytes = os.path.getsize(path)
    assert actual_bytes == expected_bytes, (
        f"apCO2_{year}: {actual_bytes} bytes, expected {expected_bytes} "
        f"({ndays} days)"
    )
    raw = np.fromfile(path, dtype=">f4").reshape(ndays, 256, 2)
    # lon0/lon1 confirmed bit-identical -- use lon0 only.
    zonal = raw[:, :, 0]  # (ndays, 256), mole fraction (e.g. 0.00035 = 350 ppm)
    daily_mean_ppm = (zonal * weights[None, :]).sum(axis=1) / weights.sum() * 1e6
    return daily_mean_ppm


def load_iron_climatology(weights):
    expected_bytes = 144 * 96 * 12 * 4
    actual_bytes = os.path.getsize(IRON_PATH)
    assert actual_bytes == expected_bytes, (
        f"iron file: {actual_bytes} bytes, expected {expected_bytes}"
    )
    raw = np.fromfile(IRON_PATH, dtype=">f4").reshape(12, 96, 144)
    zonal_mean = raw.mean(axis=2)  # (12, 96) -- average over 144 lon points first
    global_mean = (zonal_mean * weights[None, :]).sum(axis=1) / weights.sum()  # (12,)
    return global_mean  # kg Fe/m^2/s, per data.darwin's darwin_inscal_iron comment


def plot_apco2(pdf):
    lat = pco2_latitudes()
    weights = np.cos(np.radians(lat))

    years = sorted(
        int(re.search(r"apCO2_(\d+)$", f).group(1))
        for f in glob.glob(os.path.join(APCO2_DIR, "apCO2_*"))
    )
    print(f"apCO2: found {len(years)} years: {years[0]}-{years[-1]}")

    all_ppm, all_year_labels = [], []
    for year in years:
        daily = load_apco2_year(year, weights)
        all_ppm.append(daily)
        all_year_labels.append(year)
        print(f"  {year}: {daily.size} days, mean {daily.mean():.2f} ppm, "
              f"range {daily.min():.2f}-{daily.max():.2f}")

    ppm = np.concatenate(all_ppm)
    time_axis = []
    for year, daily in zip(all_year_labels, all_ppm):
        n = daily.size
        time_axis.append(year + (np.arange(n) + 0.5) / n)
    time_axis = np.concatenate(time_axis)

    fig, ax = plt.subplots(figsize=(11, 6))
    ax.plot(time_axis, ppm, lw=0.8, color="tab:red")
    ax.set_xlabel("Year")
    ax.set_ylabel("Global-mean atmospheric CO2 (ppm)")
    ax.set_title(
        "apCO2 forcing (data.darwin pCO2File) -- global (cos-lat-weighted) "
        f"mean, {all_year_labels[0]}-{all_year_labels[-1]}\n"
        "Source: NOAA Marine Boundary Layer, steph_input/v06/apCO2/"
    )
    ax.grid(alpha=0.3)
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)

    mask = (time_axis >= 1992) & (time_axis <= 2026)
    fig, ax = plt.subplots(figsize=(11, 6))
    ax.plot(time_axis[mask], ppm[mask], lw=0.8, color="tab:red")
    ax.set_xlabel("Year")
    ax.set_ylabel("Global-mean atmospheric CO2 (ppm)")
    ax.set_title("apCO2 forcing -- zoom on this run's actual period (1992-2025)")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def plot_iron(pdf):
    lat = iron_latitudes()
    weights = np.cos(np.radians(lat))
    monthly = load_iron_climatology(weights)
    print("iron dust: global-mean monthly climatology (kg Fe/m^2/s):")
    for m, v in enumerate(monthly, start=1):
        print(f"  month {m:2d}: {v:.4e}")

    months = np.arange(1, 13)
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(months, monthly, marker="o", color="tab:brown")
    ax.set_xlabel("Month")
    ax.set_ylabel("Global-mean iron dust deposition (kg Fe/m^2/s)")
    ax.set_xticks(months)
    ax.set_title(
        "Iron dust forcing (data.darwin ironFile) -- global (cos-lat-weighted) "
        "mean monthly climatology\n"
        "Source: ecco_darwin CESM-MIMI Hamilton et al. 2020, "
        "steph_input/v06/forcing/ -- repeating cycle (ironPeriod=-12), not a "
        "real time series"
    )
    ax.grid(alpha=0.3)
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def main():
    os.makedirs(os.path.dirname(OUT_PDF), exist_ok=True)
    with PdfPages(OUT_PDF) as pdf:
        plot_apco2(pdf)
        plot_iron(pdf)
    print(f"\nDone: {OUT_PDF}")


if __name__ == "__main__":
    main()
