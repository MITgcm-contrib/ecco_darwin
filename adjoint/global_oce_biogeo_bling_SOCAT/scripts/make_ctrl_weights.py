"""
Generate physically-scaled per-field prior-uncertainty (sigma) weight
files for the 8 control-vector fields, replacing the uniform ones_32b.bin
weighting that leaves m1qn3 optimizing in raw physical units simultaneously
across fields spanning ~7 orders of magnitude (K vs psu vs mol/m3 vs nM Fe)
-- see README.md "Known limitations" for the full diagnosis.

Run this from input_ad/ (it reads the *_init.bin / lev_clim_*.bin fields
already there to sanity-check the chosen values, and writes wt_*.bin next
to them).

Each output file matches the exact convention of the existing
ones_32b.bin: a single global (not per-tile) big-endian float32 field of
shape (nz=15, ny=64, nx=128), spatially uniform (no land/ocean masking --
land values are ignored by the ctrl package at runtime, same as
ones_32b.bin). This is a spatially-uniform scalar per field, not a
spatially-varying uncertainty map -- the goal here is fixing the
cross-field scale mismatch, not adding untested spatial structure.

Sigma values are prior IC uncertainties (physical units), chosen from:
  - standard ECCO-style prior IC uncertainty conventions for theta/salt
    (~1 degC, ~0.2 psu are commonly used starting points)
  - ~1-10% of each BLING tracer's own field range for the 6 BLING
    tracers, cross-checked against this experiment's own initial
    condition files below

These are defensible, physically-grounded starting values, not a
definitive calibration -- if you have a better-informed estimate of
actual IC uncertainty for this domain/config, use it instead.
"""
import numpy as np

NZ, NY, NX = 15, 64, 128
N = NZ * NY * NX

SIGMA = {
    "wt_theta.bin": 1.0,      # degC    (ECCO-style convention)
    "wt_salt.bin":  0.2,      # psu     (ECCO-style convention)
    "wt_DIC.bin":   0.02,     # mol/m3  (20 mmol/m3, ~1-2% of ~2 mol/m3 field mean)
    "wt_ALK.bin":   0.02,     # mol/m3  (20 mmol/m3, same basis as DIC)
    "wt_O2.bin":    0.015,    # mol/m3  (15 mmol/m3; field max ~0.48 mol/m3)
    "wt_NO3.bin":   0.003,    # mol/m3  (3 mmol/m3; field max ~0.049 mol/m3)
    "wt_PO4.bin":   0.0002,   # mol/m3  (0.2 mmol/m3; field max ~0.0036 mol/m3, ~1/16 of NO3, Redfield-consistent)
    "wt_FE.bin":    5.0e-7,   # mol/m3  (0.5 nM; field max ~4.3 nM)
}

# sanity-check reference: (field, its own init/climatology file)
REFERENCE_FIELDS = {
    "wt_theta.bin": "lev_clim_temp.bin",
    "wt_salt.bin":  "lev_clim_salt.bin",
    "wt_DIC.bin":   "dic_init.bin",
    "wt_ALK.bin":   "alk_init.bin",
    "wt_O2.bin":    "o2_init.bin",
    "wt_NO3.bin":   "no3_init.bin",
    "wt_PO4.bin":   "po4_init.bin",
    "wt_FE.bin":    "fe_init.bin",
}


def main():
    for fname, sigma in SIGMA.items():
        ref = REFERENCE_FIELDS[fname]
        ref_field = np.fromfile(ref, dtype=">f4")
        assert ref_field.size == N, f"{ref} has unexpected size {ref_field.size}, expected {N}"
        print(f"{fname:14s} sigma={sigma:<10.4g} "
              f"(reference {ref}: mean={ref_field.mean():.4g}, max={ref_field.max():.4g})")

        arr = np.full(N, sigma, dtype=">f4")
        arr.tofile(fname)


if __name__ == "__main__":
    main()
