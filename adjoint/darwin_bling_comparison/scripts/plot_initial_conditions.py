#!/usr/bin/env python3
"""Compare Darwin's and BLING's initial conditions for the 6 matched BGC
tracers (DIC, ALK, O2, NO3, PO4, FeT): surface maps and global-mean
vertical profiles, side by side.

Useful both as a one-off sanity check after any IC change (e.g. after
running match_ics_to_bling.py) and as a general diagnostic for spotting
IC mismatches like the ~7% ALK / ~2% DIC provenance mismatch documented
in readme.txt -- run this BEFORE and AFTER an IC fix to see the effect.

Land masking: MITgcm's ptracers_init_varia.F unconditionally zeros every
tracer where maskC==0 at cold start, so a cell is treated as dry wherever
the field is exactly 0 -- no separate bathymetry file needed.

Vertical profiles are cos(latitude)-weighted horizontal means at each
depth level (a reasonable area-weighting approximation on this lat-lon
grid; not volume-weighted, since per-cell dz doesn't vary horizontally
here).

Requires MITgcmutils is NOT needed here -- these are plain (15,64,128)
big-endian float32 binaries, read directly with numpy.

Configure via environment variables:
  DARWIN_INPUT_AD_DIR -- Darwin's input_ad/ (source of *_darwin_init.bin)
  BLING_INPUT_AD_DIR  -- BLING_SOCAT's input_ad/ (source of *_init.bin)
  IC_OUT_DIR          -- where to write the two output PNGs
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

DARWIN_DIR = os.environ.get("DARWIN_INPUT_AD_DIR", "../../global_oce_biogeo_darwin/input_ad")
BLING_DIR = os.environ.get("BLING_INPUT_AD_DIR", "../../global_oce_biogeo_bling_SOCAT/input_ad")
OUT_DIR = os.environ.get("IC_OUT_DIR", "../figures")

NZ, NY, NX = 15, 64, 128
DELZ = np.array([50., 70., 100., 140., 190., 240., 290., 340., 390., 440.,
                  490., 540., 590., 640., 690.])
DEPTH_CENTER = np.cumsum(DELZ) - DELZ / 2.0

# tracer -> (Darwin file, BLING file, Darwin unit->mol/m3 factor, display units)
TRACERS = {
    'DIC': ('dic_darwin_init.bin', 'dic_init.bin', 1e-3, 'mol m$^{-3}$'),
    'ALK': ('alk_darwin_init.bin', 'alk_init.bin', 1e-3, 'eq m$^{-3}$'),
    'O2':  ('o2_darwin_init.bin', 'o2_init.bin', 1e-3, 'mol m$^{-3}$'),
    'NO3': ('no3_darwin_init.bin', 'no3_init.bin', 1e-3, 'mol m$^{-3}$'),
    'PO4': ('po4_darwin_init.bin', 'po4_init.bin', 1e-3, 'mol m$^{-3}$'),
    'FeT': ('fet_darwin_init.bin', 'fe_init.bin', 1e-3, 'mol m$^{-3}$'),
}


def load(path):
    return np.fromfile(path, dtype='>f4').reshape(NZ, NY, NX)


def lat_axis():
    return -90.0 + (np.arange(NY) + 0.5) * (180.0 / NY)


def zonal_mean_profile(field):
    """cos(lat)-weighted horizontal mean at each depth level, dry (==0)
    cells excluded."""
    lat = lat_axis()
    w2d = np.cos(np.radians(lat))[:, None] * np.ones((1, NX))
    profile = np.full(NZ, np.nan)
    for k in range(NZ):
        wet = field[k] != 0
        if wet.any():
            profile[k] = np.average(field[k][wet], weights=w2d[wet])
    return profile


def plot_maps():
    fig, axes = plt.subplots(6, 2, figsize=(9, 18))
    for row, (name, (dfile, bfile, to_mol, units)) in enumerate(TRACERS.items()):
        d = load(f"{DARWIN_DIR}/{dfile}")[0] * to_mol
        b = load(f"{BLING_DIR}/{bfile}")[0]
        d_masked = np.ma.masked_where(d == 0, d)
        b_masked = np.ma.masked_where(b == 0, b)
        vmin = min(np.nanmin(d_masked), np.nanmin(b_masked))
        vmax = max(np.nanmax(d_masked), np.nanmax(b_masked))
        for col, (label, field) in enumerate([('Darwin', d_masked), ('BLING', b_masked)]):
            ax = axes[row, col]
            im = ax.imshow(field, origin='lower', extent=[0, 360, -90, 90],
                            vmin=vmin, vmax=vmax, cmap='viridis', aspect='auto')
            ax.set_title(f'{name} -- {label}', fontsize=10)
            if col == 0:
                ax.set_ylabel('latitude')
            if row == 5:
                ax.set_xlabel('longitude')
            fig.colorbar(im, ax=ax, label=units, shrink=0.8)
    fig.suptitle('Surface initial conditions: Darwin vs BLING', fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    out = f"{OUT_DIR}/initial_conditions_maps.png"
    fig.savefig(out, dpi=130)
    print(f"saved to {out}")


def plot_profiles():
    fig, axes = plt.subplots(2, 3, figsize=(13, 7))
    for ax, (name, (dfile, bfile, to_mol, units)) in zip(axes.flat, TRACERS.items()):
        d = load(f"{DARWIN_DIR}/{dfile}") * to_mol
        b = load(f"{BLING_DIR}/{bfile}")
        d_prof = zonal_mean_profile(d)
        b_prof = zonal_mean_profile(b)
        ax.plot(b_prof, DEPTH_CENTER, 's-', color='#D55E00', linewidth=3.5, label='BLING')
        ax.plot(d_prof, DEPTH_CENTER, 'o--', color='#0072B2', linewidth=1.5, label='Darwin')
        ax.invert_yaxis()
        ax.set_title(name, fontsize=11)
        ax.set_xlabel(units)
        ax.set_ylabel('depth (m)')
        ax.legend(fontsize=8, frameon=False)
    fig.suptitle('Global-mean vertical profiles (cos-lat weighted, wet cells only):\nDarwin vs BLING initial conditions', fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    out = f"{OUT_DIR}/initial_conditions_profiles.png"
    fig.savefig(out, dpi=130)
    print(f"saved to {out}")


if __name__ == '__main__':
    os.makedirs(OUT_DIR, exist_ok=True)
    plot_maps()
    plot_profiles()
