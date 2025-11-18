#!/usr/bin/env python3
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import savefig
from scipy.ndimage import binary_dilation
import dask
import argparse
from tqdm import tqdm

def clean_duplicate_rivers(coastal_runoff, radius=9, rel_tol=0.10):
    """
    coastal_runoff : 2D array (coastal runoff, nan = land/no runoff)
    radius        : half-width of window (9 → 19x19 window)
    rel_tol       : relative tolerance defining "very similar" runoff

    Returns:
        cleaned_runoff : 2D array with duplicates removed
    """

    ny, nx = coastal_runoff.shape
    cleaned = coastal_runoff.copy()

    # Track processed cells so we don't re-process clusters
    visited = np.zeros_like(coastal_runoff, dtype=bool)

    for iy in range(ny):
        for ix in range(nx):

            # Skip zero or nan or already processed
            current_val = coastal_runoff[iy, ix]
            if (np.isnan(current_val)) or (current_val <= 0) or visited[iy, ix]:
                continue

            # Define search window
            y0 = max(0, iy - radius)
            y1 = min(ny, iy + radius + 1)
            x0 = max(0, ix - radius)
            x1 = min(nx, ix + radius + 1)

            # Extract the patch
            patch = coastal_runoff[y0:y1, x0:x1]

            # Boolean mask of positive runoff
            positive_mask = patch > 0

            if not np.any(positive_mask):
                continue

            # Compute similarity mask
            similar_mask = np.zeros_like(patch, dtype=bool)
            for py, px in zip(*np.where(positive_mask)):
                val = patch[py, px]
                if abs(val - current_val) <= rel_tol * current_val:
                    similar_mask[py, px] = True

            # If only one similar value → skip
            if np.sum(similar_mask) <= 1:
                visited[iy, ix] = True
                continue

            # Coordinates of similar cells in global grid
            sy, sx = np.where(similar_mask)
            gy = sy + y0
            gx = sx + x0
            vals = coastal_runoff[gy, gx]

            # Find max runoff cell
            idx_max = np.argmax(vals)
            gy_max, gx_max = gy[idx_max], gx[idx_max]

            # Set all others to zero
            for j in range(len(vals)):
                if j != idx_max:
                    cleaned[gy[j], gx[j]] = 0

            # Mark all cells in this cluster as visited
            visited[gy, gx] = True

    return cleaned

def main(glofas_input_NCDFfile,arg_plot):
    # glofas_input_NCDFfile : input netcdf file. Must be one full year. m3/sec
    # arg_plot : plot map of annual discharge

    # Convert arg_plot to boolean if passed as string
    if isinstance(arg_plot, str):
        arg_plot = arg_plot.lower() in ["true", "1", "yes"]

    # Open the NetCDF file
    ds = xr.open_dataset(glofas_input_NCDFfile, engine="netcdf4")
    # load year
    year = (np.floor(np.mean(ds["valid_time"].data.astype('datetime64[Y]').astype(int) + 1970))).astype(int)

    #load pixel area m2
    area_ds = xr.open_dataset('GLOFAS_pixarea_Global_03min.nc', engine="netcdf4")
    area = area_ds["Band1"]  # shape: (lat, lon)
    area = area.isel(lat=slice(0, 3000))

    # load water mask [0-1]
    water_ds = xr.open_dataset("GLOFAS_fracwater_Global_03min.nc", engine="netcdf4")
    water = water_ds["Band1"]  # shape: (lat, lon)
    ocean = water.isnull().astype(int)
    ocean = ocean.isel(lat=slice(0, 3000))

    # Compute annual discharge at each grid cell (km3/yr) for sanity check
    annual_dis = ds["dis24"].sum(dim="valid_time")*1e-9*60*60*24
    annual_dis = annual_dis.where(annual_dis > 0)

    # Convert discharge m3/s to m/s
    area = area.rename({'lat': 'latitude', 'lon': 'longitude'})
    runoff = ds["dis24"].data / area.data  # pure numpy, ignores coords
    runoff = xr.DataArray(runoff, coords=ds["dis24"].coords, dims=ds["dis24"].dims)

    # Ensure consistent coordinates ---
    # (Rename to match)
    ocean = ocean.rename({'lat': 'latitude', 'lon': 'longitude'})
    # Create boolean masks ---
    land_mask = ~ocean  # land = True, ocean = False
    # Identify "coastal" land cells ---
    # Dilate ocean mask (grow it by 2 pixel in all directions)
    ocean_dilated_2 = binary_dilation(ocean, iterations = 2)
    # Coastal land = land cells that touch ocean
    coastal_land_mask_2 = land_mask & ocean_dilated_2
    # Now masking annual discharge
    coastal_annual_dis = annual_dis * coastal_land_mask_2.data
    cleaned_coastal = clean_duplicate_rivers(coastal_annual_dis.data, radius=9, rel_tol=0.1)
    print(f"Sum Annual discharge = {np.nansum(cleaned_coastal)} km3/yr")
    print(f"Max discharge (Amazon) = {np.nanmax(cleaned_coastal)} km3/yr")
    # Get a boolean mask of all cells in cleaned_coastal that are zero
    zero_mask = cleaned_coastal == 0
    # Set the corresponding cells in coastal_land_mask_2 to zero
    coastal_land_mask_2.data[zero_mask] = 0
    # Now masking daily runoff
    coastal_runoff = runoff * coastal_land_mask_2.data

    # # --- (optional) Visualize annual discharge ---
    if arg_plot:
        vmax = float(np.nanmax(cleaned_coastal))  # cast to plain Python float
        coastal_annual_dis.plot(
            cmap="turbo",
            norm=plt.matplotlib.colors.LogNorm(vmin = 1, vmax = vmax),  # log scale for river discharge
            cbar_kwargs={'label': "Annual discharge (km3/yr)"}
        )
        savefig('Annual_discharge_GLOFAS')

    # Write to binary file (little-endian real*4 format)
    dtype = '>f4'  # big-endian float32
    out_file = f"GloFas_{year}"

    with open(out_file, "wb") as f:
        for t in tqdm(range(coastal_runoff.sizes['valid_time']), desc="Writing binary"):
            # Explicitly compute one 2D slice into memory
            arr_2d = coastal_runoff.isel(valid_time=t).compute().data.astype(np.float32)

            # Replace NaNs or fill values
            arr_2d = np.nan_to_num(arr_2d, nan=0.0)

            # Convert to big-endian and write
            arr_2d.astype(dtype).tofile(f)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process GloFAS runoff NetCDF and extract coastal runoff.")
    parser.add_argument("ncfile", help="Path to GloFAS NetCDF input file (one year).")
    parser.add_argument("plot", help="Whether to plot annual discharge (True/False).")
    args = parser.parse_args()

    main(args.ncfile, args.plot)