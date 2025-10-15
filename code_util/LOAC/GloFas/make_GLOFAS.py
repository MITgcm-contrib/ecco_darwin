#!/usr/bin/env python3
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import savefig
from scipy.ndimage import binary_dilation
import dask
import argparse

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

    # apply coastal mask
    # --- 1. Ensure consistent coordinates ---
    # (Rename to match)
    ocean = ocean.rename({'lat': 'latitude', 'lon': 'longitude'})
    ocean_aligned = ocean.interp_like(annual_dis, method='nearest')
    # --- 2. Create boolean masks ---
    ocean_mask = ocean_aligned == 1   # ocean = True, land = False
    land_mask = ~ocean_mask  # land = True, ocean = False
    # --- 3. Identify "coastal" land cells ---
    # Dilate ocean mask (grow it by 1 pixel in all directions)
    ocean_dilated = binary_dilation(ocean_mask)
    # Coastal land = land cells that touch ocean
    coastal_land_mask = land_mask & ocean_dilated
    # --- 4. Apply to runoff product ---
    coastal_runoff = runoff * coastal_land_mask
    # Chunk datasets
    coastal_runoff = coastal_runoff.chunk({"valid_time": 1, "latitude": 1000, "longitude": 1000})
    # Now masking
    coastal_runoff = coastal_runoff.where(coastal_runoff > 0)
    # --- 5. (optional) Visualize annual discharge ---
    coastal_annual_dis = annual_dis * coastal_land_mask
    coastal_annual_dis = coastal_annual_dis.chunk({"latitude": 1000, "longitude": 1000})
    coastal_annual_dis = coastal_annual_dis.where(coastal_annual_dis > 0)

    if arg_plot:
        vmax = float(coastal_annual_dis.max().compute())  # cast to plain Python float
        coastal_annual_dis.plot(
            cmap="turbo",  # colorful colormap
            norm=plt.matplotlib.colors.LogNorm(vmin = 1e-3, vmax = vmax),  # log scale for river discharge
            cbar_kwargs={'label': "Annual discharge (km3/yr)"}
        )
        savefig('Annual_discharge_GLOFAS')

    # Write to binary file (little-endian real*4 format)
    # Convert xarray DataArray to NumPy array
    with open(f"GloFas_{year}", "wb") as f:
        for t in range(coastal_runoff.shape[0]):
            for j in range(coastal_runoff.shape[1]):
                row = coastal_runoff.isel(valid_time=t, latitude=j).values.astype(np.float32)
                row.tofile(f)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process GloFAS runoff NetCDF and extract coastal runoff.")
    parser.add_argument("ncfile", help="Path to GloFAS NetCDF input file (one year).")
    parser.add_argument("plot", help="Whether to plot annual discharge (True/False).")
    args = parser.parse_args()

    main(args.ncfile, args.plot)