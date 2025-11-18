#!/usr/bin/env python3

import os
import glob
import numpy as np
import h5py
from tqdm import tqdm
import xarray as xr
import MITgcmutils
import pickle
import argparse

def parse_args():
    parser = argparse.ArgumentParser(
        description="Remap GloFAS runoff to ECCO using precomputed or new mapping."
    )

    parser.add_argument(
        "--compute-grid",
        type=str,
        default="true",
        help="Whether to compute the GloFAS to ECCO mapping (true/false)"
    )

    parser.add_argument(
        "--coast-mask",
        type=str,
        required=True,
        help="Path to the ECCO coastal mask file (e.g. ECCO_V4r5_coastMask_orig.mat)"
    )

    parser.add_argument(
        "--mapping-file",
        type=str,
        required=True,
        help="Mapping pickle file name (e.g. glofas_ECCO_V4r5_grid_orig.pkl)"
    )

    parser.add_argument(
        "--out-prefix",
        type=str,
        required=True,
        help="Output binary file prefix (e.g. GloFAS_runoff_ECCO_V4r5_)"
    )

    args = parser.parse_args()

    # Convert compute-grid to boolean
    val = args.compute_grid.lower().strip()
    if val in ["true", "1", "yes", "y"]:
        compute_grid = True
    elif val in ["false", "0", "no", "n"]:
        compute_grid = False
    else:
        raise ValueError("--compute-grid must be true or false")

    return compute_grid, args.coast_mask, args.mapping_file, args.out_prefix

# ----------------------
# User-editable settings
# ----------------------
gridDir = '/Users/rsavelli/Documents/Models/grid/ECCO_V4r5/'   # contains XC.data, YC.data, RAC.data, hFacC.data
dataDir = '/Users/rsavelli/Documents/CMS_LOAC/'  # contains coastal mask
binDir1 = '/Users/rsavelli/Documents/CMS_LOAC/'  # directory with GloFAS runoff binary files
cellarea_file = os.path.join(binDir1, 'GloFas/GLOFAS_pixarea_Global_03min.nc')  # cellarea per source grid cell
saveDir = '/Users/rsavelli/Documents/CMS_LOAC/'
writeDir = '/Users/rsavelli/Documents/CMS_LOAC/'

filename_pattern = 'GloFas_*'  # adjust to your filenames

# GloFAS grid specs
nLon = 7200
nLat = 3000
lon0 = -179.975
lat0 = -59.975
lonInc = 0.05
latInc = 0.05

# ECCO stitched grid settings
numFaces = 13
nx = 90
ny = numFaces * nx
nz = 50

dtype = '>f4'        # big-endian float32

compute_grid, coast_mask_file, mapping_file, out_prefix = parse_args()

# ----------------------
# Helper functions
# ----------------------
# ---- load ECCO grid (XC, YC, RAC, hFacC) ----
def load_grid(grid_dir):
    path_grid = grid_dir
    XC = MITgcmutils.mds.rdmds((f"{path_grid}/XC"))
    YC = MITgcmutils.mds.rdmds((f"{path_grid}/YC"))
    hfacc = MITgcmutils.mds.rdmds((f"{path_grid}/hFacC"))
    hfacc = hfacc[0, :, :].squeeze()
    RAC = MITgcmutils.mds.rdmds((f"{path_grid}/RAC"))
    print(f"ECCO min lon {np.min(XC)}")
    print(f"ECCO max lon {np.max(XC)}")
    print(f"ECCO min lat {np.min(YC)}")
    print(f"ECCO max lat {np.max(YC)}")
    return XC, YC, RAC, hfacc

def write_bin(fname, arr, dtype=np.float32):
    os.makedirs(os.path.dirname(fname), exist_ok=True)
    arr.astype(dtype).tofile(fname)

def build_glofas_grid(nLon, nLat, lon0, lat0, lonInc=0.05, latInc=0.05):
    lon = np.concatenate(([lon0], lon0 + np.cumsum(np.ones(nLon-1) * lonInc)))
    lat = np.concatenate(([lat0], lat0 + np.cumsum(np.ones(nLat-1) * latInc)))
    xx, yy = np.meshgrid(lon, lat)  # xx, yy shapes: (nLat, nLon)
    print(f"Glofas min lon {np.min(lon)}")
    print(f"Glofas max lon {np.max(lon)}")
    print(f"Glofas min lat {np.min(lat)}")
    print(f"Glofas max lat {np.max(lat)}")
    return lon, lat, xx, yy


def aggregate_day(src_values, ecco_idx, weights, n_ecco):
    """
    Aggregate source runoff to ECCO grid for one day,
    summing runoff and averaging the weights per unique ECCO cell.
    """
    unique_idx = np.unique(ecco_idx)

    sum_runoff = np.zeros_like(unique_idx, dtype=np.float64)
    mean_weight = np.zeros_like(unique_idx, dtype=np.float64)

    for i, uidx in enumerate(unique_idx):
        mask = (ecco_idx == uidx)
        # Sum runoff for this ECCO cell
        sum_runoff[i] = np.nansum(src_values[mask])
        # Average weight for this ECCO cell
        mean_weight[i] = np.nanmean(weights[mask])

    # Compute final aggregated runoff
    aggregated = np.zeros(n_ecco, dtype=np.float64)
    aggregated[unique_idx] = sum_runoff / mean_weight

    return aggregated

# ----------------------
# Start
# ----------------------
if __name__ == '__main__':
    os.makedirs(writeDir, exist_ok=True)
    os.makedirs(saveDir, exist_ok=True)

    # ---- load coast mask (MAT) ----
    matfile = os.path.join(dataDir, coast_mask_file)
    if not os.path.exists(matfile):
        raise FileNotFoundError(matfile + " not found")

    with h5py.File(matfile, 'r') as f:
            if 'coastMask' in f:
                ds = f['coastMask']
                coastMask = ds[()]

    print("Loading ECCO grid...")
    XC, YC, RAC, hFacC = load_grid(gridDir)

    # Build GloFAS lon/lat mesh
    print("Building GloFAS grid...")
    lon, lat, xx, yy = build_glofas_grid(nLon, nLat, lon0, lat0, lonInc, latInc)
    # xx, yy shapes: (nLat, nLon)

    # List GloFAS runoff files
    files = sorted(glob.glob(os.path.join(binDir1, filename_pattern)))
    if len(files) == 0:
        raise FileNotFoundError("No runoff files found with pattern: " + os.path.join(binDir1, filename_pattern))
    print(f"Found {len(files)} runoff files.")

    # load GloFAS pixel area
    cellarea = None
    if os.path.exists(cellarea_file):
        print("Loading GloFAS pixel area netCDF")
        area_ds = xr.open_dataset(cellarea_file, engine="netcdf4")
        area = area_ds["Band1"]  # shape: (lat, lon)
        area = area.isel(lat=slice(0, 3000))

        # Verify and possibly transpose
        if area.shape[-2:] == (nLat, nLon):
            cellarea = area
        elif area.shape[-2:] == (nLon, nLat):
            cellarea = area.T
        else:
            raise ValueError(f"Unexpected area shape {area.shape}, expected ({nLat},{nLon}) or ({nLon},{nLat})")

        # ensure it's float and aligned
        cellarea = np.array(cellarea.astype(np.float64))
        cellarea = np.flipud(cellarea)

    # ----------------------
    # compute grid mapping (GloFAS non-zero runoff -> nearest ECCO wet cell)
    # ----------------------
    mapping_file = os.path.join(saveDir, mapping_file)
    if compute_grid:
        print("Computing mapping...")

        # Accumulate total runoff across files (memmap per file) to find non-zero source cells
        total_runoff = np.zeros((nLat, nLon), dtype=np.float64)  # shape (nLat,nLon) like runoff transposed
        for fpath in tqdm(files, desc="Summing runoff"):
            filesize = os.path.getsize(fpath)
            floats = filesize // np.dtype(dtype).itemsize
            if floats % (nLon * nLat) != 0:
                raise RuntimeError(f"{fpath} size not multiple of nLon*nLat")
            numDays = floats // (nLon * nLat)
            data_mm = np.memmap(fpath, dtype=dtype, mode='r').reshape((numDays, nLat, nLon))
            s = np.nansum(data_mm, axis=0)  # shape (nLat,nLon)
            total_runoff += s
            del data_mm

        total_runoff = np.flipud(total_runoff)
        total_runoff_flat = total_runoff.ravel()
        ir_flat = np.where(total_runoff_flat != 0)[0]  # 1D array of positions
        runoffLon = xx.ravel()[ir_flat]   # 1D array length = number of non-zero cells
        runoffLat = yy.ravel()[ir_flat]

        # prepare ECCO wet cells: wet where coastMask is NaN
        wet_mask = np.isnan(coastMask)  # True where wet

        # flatten ECCO arrays in C-order and select wet ones
        XC_flat = XC.ravel(order='C')
        YC_flat = YC.ravel(order='C')
        RAC_flat = RAC.ravel(order='C')
        wet_flat_mask = wet_mask.ravel(order='C')

        wet_indices_fullflat = np.flatnonzero(wet_flat_mask)        # indices in full flattened ECCO grid
        wetX_flat = XC_flat[wet_flat_mask]
        wetY_flat = YC_flat[wet_flat_mask]
        wetRAC_flat = RAC_flat[wet_flat_mask]

        # src_area per non-zero source cell
        src_area = cellarea.ravel()[ir_flat]

        # flatten runoff source arrays
        n_src = ir_flat.size
        print(f"\nNumber of non-zero GloFAS source cells: {n_src}")

        # filter out any NaN wet points
        valid_mask = ~np.isnan(wetX_flat) & ~np.isnan(wetY_flat)
        wetX_valid = wetX_flat[valid_mask]
        wetY_valid = wetY_flat[valid_mask]
        wetRAC_valid = wetRAC_flat[valid_mask]
        wet_indices_valid = wet_indices_fullflat[valid_mask]

        # prepare mapping arrays
        llc_indices = np.empty(n_src, dtype=np.int64)   # stores full-flat ECCO index for each src cell
        llc_areas = np.empty(n_src, dtype=np.float64)
        llc_weights = np.empty(n_src, dtype=np.float64)

        print("Mapping non-zero GloFAS cells to nearest ECCO wet cell...")
        for i in range(len(runoffLon)):
            dlon2 = (runoffLon[i] - wetX_valid) ** 2 * np.cos(np.deg2rad(runoffLat[i]))
            dlat2 = (runoffLat[i] - wetY_valid) ** 2
            distWet = dlon2 + dlat2

            ix = np.argmin(distWet)

            # store mapped ECCO full-flat indices, area and weight
            llc_indices[i] = wet_indices_valid[ix]
            llc_areas[i] = wetRAC_valid[ix]
            # avoid division by zero
            with np.errstate(divide='ignore', invalid='ignore'):
                llc_weights[i] = wetRAC_valid[ix] / src_area[i]

        # build glofas mapping dict
        glofas = {
            'index': ir_flat,  # flat indices of source in C-order
            'area': src_area,
            'ECCOIndex': llc_indices,
            'ECCOArea': llc_areas,
            'ECCOWeight': llc_weights,
            'years': np.array([int(os.path.basename(fp)[-4:]) if os.path.basename(fp)[-4:].isdigit() else 0 for fp in files])
        }

        # Save the dictionary to a binary file
        print("Saving mapping to pickle file ...")
        with open(mapping_file, 'wb') as f:
            pickle.dump(glofas, f)
        print("Mapping saved.")

    else:
        print("Loading existing mapping ...")
        with open(mapping_file, 'rb') as f:
            glofas = pickle.load(f)

    # ----------------------
    # Second stage: loop files, per day aggregate -> ECCO grid and write per-year binary
    # ----------------------
    print("Starting per-file, per-day remapping and writing ...")
    ny_llc, nx_llc = XC.shape  # ECCO grid size

    for fpath in files:
        fname = os.path.basename(fpath)
        import re

        m = re.search(r'(\d{4})', fname[::-1])
        yr = int(fname[::-1][m.start():m.end()][::-1]) if m else None

        print(f"\nProcessing file {fname}, year {yr}")

        # Determine number of days
        filesize = os.path.getsize(fpath)
        floats = filesize // np.dtype(dtype).itemsize
        if floats % (nLon * nLat) != 0:
            raise RuntimeError(f"{fpath} size not integer multiple of nLon*nLat")
        numDays = floats // (nLon * nLat)

        # Memory-map the source GloFAS file (read-only)
        data_mm = np.memmap(fpath, dtype=dtype, mode='r').reshape((numDays, nLat, nLon))

        # Memory-map the output ECCO file as 3D (nx_llc, ny_llc, numDays)
        out_fname = f"{out_prefix}{yr:04d}"
        out_path = os.path.join(writeDir, out_fname)
        ECCO_year = np.memmap(out_path, dtype=dtype, mode='w+', shape=(numDays, ny_llc, nx_llc))

        # Precompute mapping
        src_idx = glofas['index'].astype(np.int64)
        ecco_idx = glofas['ECCOIndex'].astype(np.int64)
        weights = glofas['ECCOWeight'].astype(np.float64)

        unique_idx, inverse = np.unique(ecco_idx, return_inverse=True)

        # --- Loop over days with tqdm ---
        for j in tqdm(range(numDays), desc=f"Year {yr}", unit="day"):
            # Flatten daily runoff
            data_mm_flipped = np.flipud(data_mm[j, :, :].squeeze())
            runoff_flat = data_mm_flipped.ravel(order='C')

            # Aggregate
            src_values_day = runoff_flat[src_idx]  # flatten daily runoff and select non-zero sources
            aggregated_runoff = aggregate_day(src_values_day, ecco_idx, weights, n_ecco=XC.size)

            # Assign directly to 3D memory-mapped array
            # Work directly on memmap (no intermediate arrays)
            sl = ECCO_year[j, :, :]
            flat = sl.reshape(-1)
            flat[unique_idx] = aggregated_runoff[unique_idx].astype(dtype)

        # Flush to disk
        del data_mm, ECCO_year
        print(f"Year {yr} written to {out_path}")

    print("All done.")