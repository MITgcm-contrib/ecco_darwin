#!/usr/bin/env python3
"""Same transform as apply_control_to_darwin_ic.py, defaulting to the
full-year run directories (run_bling_1yr / run_darwin_1yr) instead of
the 1-month pilot's."""
import os
import sys
import numpy as np
from MITgcmutils import rdmds

NZ, NY, NX = 15, 64, 128

BLING_RUN = os.environ.get("BLING_RUN_DIR", "../run_bling_1yr")
DARWIN_RUN = os.environ.get("DARWIN_RUN_DIR", "../run_darwin_1yr")
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
            print(f"{xxname}.{cycle:010d} not found, treating as zero perturbation (first guess)")
            xx = np.zeros((NZ, NY, NX))

        sqrtw = np.sqrt(np.where(weight > 0, weight, np.nan))
        delta = np.where(weight > 0, xx / sqrtw, 0.0)
        physical_bling_units = prior + delta

        physical_darwin_units = (physical_bling_units * 1000.0).astype(">f4")
        out_path = f"{DARWIN_RUN}/{cfg['out']}"
        physical_darwin_units.tofile(out_path)
        print(f"{xxname}.{cycle:010d} -> {out_path}  "
              f"surface min/mean/max = {physical_darwin_units[0].min():.6g}/"
              f"{physical_darwin_units[0].mean():.6g}/{physical_darwin_units[0].max():.6g}")


if __name__ == "__main__":
    main(int(sys.argv[1]))
