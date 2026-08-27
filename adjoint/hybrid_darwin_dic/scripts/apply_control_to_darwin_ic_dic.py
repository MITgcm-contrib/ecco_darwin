#!/usr/bin/env python3
"""Apply the current M1QN3 control-vector guess (xx_ptr1/xx_ptr2, pkg/dic's
DIC/ALK controls) to Darwin's own physical initial-condition files.

RESTRICTED TO DIC+ALK ONLY on Darwin's side, matching hybrid_darwin_bling's
own pilot scope -- but enforced here, not by de-registering controls in
run_dic/data.ctrl. All 7 fields (theta, salt, DIC, ALK, PO4, DOP, O2) stay
registered there: de-registering PO4/DOP/O2/theta/salt was tried first and
found to corrupt pkg/dic's own ALK adjoint gradient (confirmed via a
controlled A/B test -- see data.ctrl's own comment for the full
diagnosis). pkg/dic is therefore left free to optimize all 7 internally;
this script simply never reads xx_theta/xx_salt/xx_ptr3/xx_ptr4/xx_ptr5,
so Darwin's ICs for temperature, salinity, PO4, DOP, and O2 stay frozen at
their prior values regardless of what pkg/dic does with them.

Mirrors MITgcm's ctrl_map_genarr3d.F preconditioning transform exactly:
    physical = prior + xx_gen / sqrt(ctrl_weight)
using the SAME prior ICs (global_oce_biogeo_dic_SOCAT's own dic_init.bin/
alk_init.bin) and the SAME weight files (wt_DIC.bin/wt_ALK.bin) pkg/dic's
own adjoint run used.

Usage: apply_control_to_darwin_ic_dic.py <cycle>
Reads xx_ptr1/xx_ptr2.<cycle>.*.data (pkg/dic's 4-tile control files,
run_dic/) and writes dic_darwin_init.bin/alk_darwin_init.bin (mmol/m3,
run_darwin/) for the NEXT Darwin forward run.
"""
import sys
import numpy as np

sys.path.insert(0, "/Users/carrolld/Documents/research/ECCO/offline/scripts")
from MITgcmutils import rdmds

NZ, NY, NX = 15, 64, 128

DIC_RUN = "/Users/carrolld/Documents/research/adjoint/MITgcm/verification/hybrid_darwin_dic/run_dic"
DARWIN_RUN = "/Users/carrolld/Documents/research/adjoint/MITgcm/verification/hybrid_darwin_dic/run_darwin"
DIC_PRIOR_DIR = "/Users/carrolld/Documents/research/adjoint/MITgcm/verification/global_oce_biogeo_dic_SOCAT/input_ad"

# pkg/dic's own tracer order (dic_tr_register.F): DIC=1, Alk=2 (PO4/DOP/O2
# excluded from this pilot's control scope -- see module docstring).
FIELDS = {
    "xx_ptr1": dict(prior="dic_init.bin", weight="wt_DIC.bin", out="dic_darwin_init.bin"),
    "xx_ptr2": dict(prior="alk_init.bin", weight="wt_ALK.bin", out="alk_darwin_init.bin"),
}


def main(cycle):
    for xxname, cfg in FIELDS.items():
        prior = np.fromfile(f"{DIC_PRIOR_DIR}/{cfg['prior']}", dtype=">f4").reshape(NZ, NY, NX)
        weight = np.fromfile(f"{DIC_PRIOR_DIR}/{cfg['weight']}", dtype=">f4").reshape(NZ, NY, NX)

        try:
            xx = rdmds(f"{DIC_RUN}/{xxname}", cycle)
        except IOError:
            print(f"{xxname}.{cycle:010d} not found, treating as zero perturbation (first guess)")
            xx = np.zeros((NZ, NY, NX))

        sqrtw = np.sqrt(np.where(weight > 0, weight, np.nan))
        delta = np.where(weight > 0, xx / sqrtw, 0.0)
        physical_dic_units = prior + delta  # mol/m3 (or mol eq/m3 for ALK), pkg/dic's own units

        physical_darwin_units = (physical_dic_units * 1000.0).astype(">f4")  # -> mmol/m3
        out_path = f"{DARWIN_RUN}/{cfg['out']}"
        physical_darwin_units.tofile(out_path)
        print(f"{xxname}.{cycle:010d} -> {out_path}  "
              f"surface min/mean/max = {physical_darwin_units[0].min():.6g}/"
              f"{physical_darwin_units[0].mean():.6g}/{physical_darwin_units[0].max():.6g}")


if __name__ == "__main__":
    main(int(sys.argv[1]))
