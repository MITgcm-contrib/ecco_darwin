"""
Build monthly SOCAT surface-pCO2 climatology files for the pkg/profiles
package, on the global_oce_biogeo_bling grid (128x64, 2.8125 deg, 2x2 tile
decomposition, sNx=64 sNy=32).

Source data: SOCATv2026_tracks_gridded_monthly.nc, a 1-degree gridded
monthly SOCAT product with a 'tmnth' time axis (days since 1970-01-01,
mid-month timestamps, e.g. 15.5 for Jan 1970), 'fco2_ave_weighted'
(per-cruise-weighted mean fCO2, uatm) and 'fco2_std_weighted' fields on
(tmnth, ylat, xlon), xlon in [-180, 180).

Output: socat_pco2_clim_month01.nc .. month12.nc, one per calendar month,
each a multi-year climatological average of that calendar month (across
all years present in the source file) binned onto the model's 128x64
grid and masked to ocean cells only (bathy.bin < 0).

PROVENANCE: this is a corrected reimplementation of the script originally
used to generate the socat_pco2_clim_month*.nc files shipped in
input_ad/, reconstructed and validated against the raw source file
(available at raw_data/SOCATv2026_tracks_gridded_monthly.nc) after the
original working copy was lost (it lived only in a scratch directory,
never committed). VERIFIED: running this script reproduces the exact
observation count AND exact per-tile maximum for all 12 months against
the values recorded from the original campaign (see results/output.txt):

    month  n     max_tile      month  n     max_tile
    01     2759  867           07     2050  901
    02     2733  883           08     1978  873
    03     2695  884           09     2012  871
    04     2378  891           10     2334  919
    05     2101  883           11     2344  926
    06     1931  874           12     2564  871

Two things had to be fixed relative to a naive reconstruction, both
confirmed against the real source data by exactly reproducing the counts
above:

  1. Remapping method: the model grid (2.8125 deg) is coarser than the
     source (1 deg), so each model cell spans ~2-3 source cells along
     each axis. BIN-AVERAGING all source cells whose center falls within
     a model cell (not nearest-neighbor sampling) is required to
     reproduce the original observation counts -- nearest-neighbor
     undercounts by ~40% (1661 vs. 2759 for January) because it discards
     valid source cells whenever the single "nearest" one happens to be
     masked.
  2. Mid-month observation dates: the original computed each month's
     nominal date by stepping 30 days at a time from Jan 15
     (`date(2007,1,15) + timedelta(days=30*(n-1))`), which drifts from
     true calendar mid-month by 2 days in July growing to 4-5 days by
     December. This version uses true calendar mid-month dates directly.
     This only affects the `prof_date` metadata written into each file,
     not the averaged pCO2 values or the bin-averaging above.

One thing that is NOT independently verified: the exact prof_PCOweight
formula (1/sigma**2 using fco2_std_weighted, floored at MIN_SIGMA_UATM)
matches the original by recollection, not by reproducing shipped bytes --
unlike the counts above, the input_ad/*.nc files don't carry the original
per-obs weight values in a form this script's output was diffed against.
If you regenerate the shipped .nc files with this script, prof_date will
change (this is the fix) and prof_PCOweight may differ slightly from the
original for low-observation-count grid cells (near-zero sigma), even
though prof_PCO/prof_lon/prof_lat will not. See README.md "Known
limitations."
"""
import datetime
import numpy as np
import netCDF4 as nc

SOCAT_SOURCE = "../../../../raw_data/SOCATv2026_tracks_gridded_monthly.nc"
BATHY_FILE = "bathy.bin"
OUT_PREFIX = "socat_pco2_clim_month"

NX, NY = 128, 64
DELX = DELY = 2.8125
SNX, SNY = 64, 32          # tile size
NOBSGLOB = 1200            # pkg/profiles cap, per tile per file (PROFILES_SIZE.h)
MIN_SIGMA_UATM = 1.0       # floor on fco2_std_weighted before inverting to a weight,
                           # avoids divide-by-zero / blown-up weights for single-obs cells


def month_mid_date(month, base_year=2007):
    """True calendar mid-month date, not a compounded 30-day-step approximation."""
    return datetime.date(base_year, month, 15)


def load_bathy_mask():
    bathy = np.fromfile(BATHY_FILE, dtype=">f4").reshape(NY, NX)
    return bathy < 0  # ocean = negative depth


def bin_average(src_field, xlon_src, ylat_src, lon_edges, lat_edges):
    """Average all finite source cells whose center falls in each model bin."""
    ix_bin = np.searchsorted(lon_edges, xlon_src, side="right") - 1
    iy_bin = np.searchsorted(lat_edges, ylat_src, side="right") - 1
    sum_grid = np.zeros((NY, NX))
    cnt_grid = np.zeros((NY, NX))
    for j in range(len(ylat_src)):
        if not (0 <= iy_bin[j] < NY):
            continue
        for i in range(len(xlon_src)):
            if not (0 <= ix_bin[i] < NX):
                continue
            v = src_field[j, i]
            if np.isfinite(v):
                sum_grid[iy_bin[j], ix_bin[i]] += v
                cnt_grid[iy_bin[j], ix_bin[i]] += 1
    return np.where(cnt_grid > 0, sum_grid / np.maximum(cnt_grid, 1), np.nan)


def build_month(month, src, ocean_mask, xlon_src, ylat_src, lon_edges, lat_edges):
    tmnth = src.variables["tmnth"][:]
    month_idx = np.arange(month - 1, len(tmnth), 12)  # this calendar month, all years

    fco2 = np.ma.filled(src.variables["fco2_ave_weighted"][month_idx, :, :], np.nan)
    clim = np.nanmean(fco2, axis=0)
    remapped = bin_average(clim, xlon_src, ylat_src, lon_edges, lat_edges)
    remapped[~ocean_mask] = np.nan

    fco2_std = np.ma.filled(src.variables["fco2_std_weighted"][month_idx, :, :], np.nan)
    sig_clim = np.nanmean(fco2_std, axis=0)
    sig_remapped = bin_average(sig_clim, xlon_src, ylat_src, lon_edges, lat_edges)
    sig_remapped = np.maximum(sig_remapped, MIN_SIGMA_UATM)

    valid = np.argwhere(np.isfinite(remapped))
    tile_counts = {}
    for (jy, ix_) in valid:
        tile = (ix_ // SNX, jy // SNY)
        tile_counts[tile] = tile_counts.get(tile, 0) + 1
    max_tile = max(tile_counts.values()) if tile_counts else 0
    assert max_tile <= NOBSGLOB, (
        f"month {month:02d}: tile obs count {max_tile} exceeds NOBSGLOB={NOBSGLOB} "
        "-- pkg/profiles will silently truncate or fail at runtime"
    )
    print(f"month{month:02d}: n={len(valid)} obs, max_tile={max_tile} (NOBSGLOB={NOBSGLOB})")

    mlon = (np.arange(NX) + 0.5) * DELX
    mlat = -90.0 + (np.arange(NY) + 0.5) * DELY
    prof_date = int(month_mid_date(month).strftime("%Y%m%d"))

    outname = f"{OUT_PREFIX}{month:02d}.nc"
    with nc.Dataset(outname, "w", format="NETCDF4_CLASSIC") as ds:
        ds.createDimension("iOBS", len(valid))
        ds.createDimension("iDEPTH", 1)
        lon_v = ds.createVariable("prof_lon", "f8", ("iOBS",))
        lat_v = ds.createVariable("prof_lat", "f8", ("iOBS",))
        date_v = ds.createVariable("prof_date", "i4", ("iOBS",))
        pco2_v = ds.createVariable("prof_PCO", "f8", ("iOBS", "iDEPTH"))
        wt_v = ds.createVariable("prof_PCOweight", "f8", ("iOBS", "iDEPTH"))
        for k, (jy, ix_) in enumerate(valid):
            lon_v[k] = mlon[ix_]
            lat_v[k] = mlat[jy]
            date_v[k] = prof_date
            pco2_v[k, 0] = remapped[jy, ix_]
            wt_v[k, 0] = 1.0 / (sig_remapped[jy, ix_] ** 2)
    return outname


def main():
    ocean_mask = load_bathy_mask()
    lon_edges = np.arange(NX + 1) * DELX
    lat_edges = -90.0 + np.arange(NY + 1) * DELY
    with nc.Dataset(SOCAT_SOURCE) as src:
        xlon_src = np.mod(src.variables["xlon"][:], 360.0)
        ylat_src = src.variables["ylat"][:]
        for month in range(1, 13):
            build_month(month, src, ocean_mask, xlon_src, ylat_src, lon_edges, lat_edges)


if __name__ == "__main__":
    main()
