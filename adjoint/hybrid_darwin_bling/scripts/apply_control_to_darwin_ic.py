#!/usr/bin/env python3
"""Apply the current M1QN3 control-vector guess (xx_ptr1/xx_ptr2, BLING's
DIC/ALK controls) to Darwin's own physical initial-condition files.

Mirrors MITgcm's ctrl_map_genarr3d.F preconditioning transform exactly:
    physical = prior + xx_gen / sqrt(ctrl_weight)
using the SAME prior IC (BLING's own dic_init.bin/alk_init.bin, since
Darwin's ICs were already matched to BLING's via
../../global_oce_biogeo_darwin/scripts/match_ics_to_bling.py) and the
SAME weight files (wt_DIC.bin/wt_ALK.bin) BLING's own adjoint run used --
so this reproduces exactly what BLING did internally to its own physical
field, applied externally to Darwin's copy.

Usage: apply_control_to_darwin_ic.py <cycle>
Reads xx_ptr1.<cycle>.*.data / xx_ptr2.<cycle>.*.data (BLING's 4-tile
control files, in BLING_RUN_DIR) and writes dic_darwin_init.bin/
alk_darwin_init.bin (mmol/m3) into DARWIN_RUN_DIR for the next Darwin
forward run.

Requires MITgcmutils (ships with MITgcm under utils/python/MITgcmutils --
add that directory to PYTHONPATH before running this script).

Configure via environment variables:
  BLING_RUN_DIR       -- this pilot's run_bling/ (has xx_ptr1/xx_ptr2.<cycle>)
  DARWIN_RUN_DIR       -- this pilot's run_darwin/ (destination for the IC files)
  BLING_PRIOR_DIR      -- global_oce_biogeo_bling_SOCAT/input_ad/ (prior IC + weight files)
"""
import os
import sys
import numpy as np
from MITgcmutils import rdmds

NZ, NY, NX = 15, 64, 128

BLING_RUN = os.environ.get("BLING_RUN_DIR", "../run_bling")
DARWIN_RUN = os.environ.get("DARWIN_RUN_DIR", "../run_darwin")
BLING_PRIOR_DIR = os.environ.get("BLING_PRIOR_DIR", "../../global_oce_biogeo_bling_SOCAT/input_ad")

FIELDS = {
    "xx_ptr1": dict(prior="dic_init.bin", weight="wt_DIC.bin", out="dic_darwin_init.bin"),
    "xx_ptr2": dict(prior="alk_init.bin", weight="wt_ALK.bin", out="alk_darwin_init.bin"),
}


def main(cycle):
    for xxname, cfg in FIELDS.items():
        prior = np.fromfile(f"{BLING_PRIOR_DIR}/{cfg['prior']}", dtype=">f4").reshape(NZ, NY, NX)
        weight = np.fromfile(f"{BLING_PRIOR_DIR}/{cfg['weight']}", dtype=">f4").reshape(NZ, NY, NX)

        try:
            xx = rdmds(f"{BLING_RUN}/{xxname}", cycle)
        except IOError:
            # optimcycle 0's xx file may not exist yet (first guess = all
            # zero perturbation) -- fall back to zero, i.e. physical = prior.
            print(f"{xxname}.{cycle:010d} not found, treating as zero perturbation (first guess)")
            xx = np.zeros((NZ, NY, NX))

        sqrtw = np.sqrt(np.where(weight > 0, weight, np.nan))
        delta = np.where(weight > 0, xx / sqrtw, 0.0)
        physical_bling_units = prior + delta  # mol/m3 or eq/m3, BLING's own units

        physical_darwin_units = (physical_bling_units * 1000.0).astype(">f4")  # -> mmol/m3
        out_path = f"{DARWIN_RUN}/{cfg['out']}"
        physical_darwin_units.tofile(out_path)
        print(f"{xxname}.{cycle:010d} -> {out_path}  "
              f"surface min/mean/max = {physical_darwin_units[0].min():.6g}/"
              f"{physical_darwin_units[0].mean():.6g}/{physical_darwin_units[0].max():.6g}")


if __name__ == "__main__":
    main(int(sys.argv[1]))
