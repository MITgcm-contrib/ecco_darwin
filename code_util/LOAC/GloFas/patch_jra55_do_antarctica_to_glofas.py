#!/usr/bin/env python3

import numpy as np
import os
from glob import glob
import MITgcmutils
import sys
import argparse
# ================================================
# PARSE ARGUMENTS
# ================================================
parser = argparse.ArgumentParser(description="Patch GloFAS–ECCO runoff with JRA55 southern runoff.")
parser.add_argument(
    "--glofas_prefix",
    type=str,
    default=None,
    help="Prefix of previously regridded GloFAS→ECCO runoff binary files."
)
args = parser.parse_args()

# ================================================
# USER INPUTS
# ================================================

# Folder containing GloFAS→ECCO binary files
glofas_dir = '/Users/rsavelli/Documents/CMS_LOAC/'
# JRA55-do runoff file for year 2000 gridded on ECCO grid (any leap year works as it repeats over year)
jra_path = "/Volumes/G-DRIVEArmorATD/Documents/CMS_LOAC/jra55/river_runoff_ECCOv4r5/jra55_do_runoff_ECCO_V4r5_2000"
gridDir = '/Users/rsavelli/Documents/Models/grid/ECCO_V4r5/'   # contains XC.data, YC.data, RAC.data, hFacC.data

# dtype of all files
dtype = '>f4'

# ECCO grid
numFaces = 13
nx = 90
ny = numFaces * nx
nz = 50

ntime = 366       # 2000 = leap year

# ---- load ECCO grid (XC, YC, RAC, hFacC) ----
def load_grid(grid_dir):
    path_grid = grid_dir
    XC = MITgcmutils.mds.rdmds((f"{path_grid}/XC"))
    YC = MITgcmutils.mds.rdmds((f"{path_grid}/YC"))
    hfacc = MITgcmutils.mds.rdmds((f"{path_grid}/hFacC"))
    hfacc = hfacc[0, :, :].squeeze()
    RAC = MITgcmutils.mds.rdmds((f"{path_grid}/RAC"))
    return XC, YC, RAC, hfacc

print("Loading ECCO grid...")
XC, YC, RAC, hFacC = load_grid(gridDir)

# ================================================
# LOAD JRA55 RUNOFF (YEAR 2000)
# ================================================

print("Loading JRA55-do runoff gridded on ECCO...")
jra = np.memmap(
    jra_path,
    dtype=dtype,
    mode="r",
    shape=(ntime, ny, nx)
)

jra_flat = MITgcmutils.llc.flat(np.sum(jra, axis=0))

# Extract southern latitudes (< -59.975)
mask_jra_south = YC < -59.975

# ================================================
# PATCH EACH Glofas FILE
# ================================================

eccofiles = sorted(glob(os.path.join(glofas_dir, f"{args.glofas_prefix}*")))

print(f"\nUsing prefix: {args.glofas_prefix}")
print(f"\nFound {len(eccofiles)} ECCO files to patch.\n")

for fpath in eccofiles:

    fname = os.path.basename(fpath)
    print(f"--- Patching {fname} ---")

    # Determine time dimension from file size
    filesize = os.path.getsize(fpath)
    ntime = filesize // (nx * ny * np.dtype(dtype).itemsize)

    # memmap ECCO file in read/write mode
    ecco = np.memmap(
        fpath,
        dtype=dtype,
        mode="r+",
        shape=(ntime, ny, nx)
    )

    ecco_flat1 = MITgcmutils.llc.flat(np.sum(ecco, axis=0))

    # Compute annual discharge
    # --- compute annual mean discharge ---
    annual_dis = np.nansum(ecco, axis=0) * RAC * 1e-9 * 60 * 60 * 24

    # Patch per-day
    for day in range(ntime):
        jra_day = jra[day, :, :].copy()
        # Flatten JRA southern region for that day
        jra_day[~mask_jra_south] = 0  # (time, south_lat, lon)
        # Extract current ECCO day
        control = ecco[day, :, :].copy()
        control[~mask_jra_south] = 0
        # If ANY non-zero exists in the southern region -> already patched
        if np.sum(jra_day) == np.sum(control):
            print(f"Already patched — skipping.")
            sys.exit()
        else:
            ecco[day, :, :] = ecco[day, :, :] + jra_day

    ecco_flat = MITgcmutils.llc.flat(np.sum(ecco, axis=0))

    del ecco  # flush memmap

    print(f"Patched file: {fname}\n")

print("All files successfully patched.")