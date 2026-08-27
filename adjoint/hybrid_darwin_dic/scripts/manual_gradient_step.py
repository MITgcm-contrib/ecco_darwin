#!/usr/bin/env python3
"""Hand-rolled gradient-descent step for the DIC-pkg hybrid pilot,
bypassing M1QN3/optim.x entirely.

WHY: the gradient returned for genarr3d slot 2 is not usable in this
configuration -- whichever field occupies that slot comes back with
surface values ~1e12 against a normal ~1e2-1e5. Because M1QN3 packs every
registered control into one vector, that one bad field dominates the
two-norm and destroys the step sizing even for fields that are themselves
clean, which is what stalled the M1QN3-driven loop (see
run_hybrid_loop_dic.sh and results/cost_trajectory_dic.csv: the control
never moved and J was identical across cycles). This script therefore
bypasses optim.x and reads the raw adxx_ptr* output directly.

Slots 1, 3 and 4 are verified clean, agreeing bit-for-bit across
independent runs. NOTE: registering DIC and ALK as separate single-control
runs (run_dic_A / run_dic_B, as this script does) is NOT necessary -- the
7-control run_dic config puts them in slots 3 and 4 and returns both
gradients, bit-identical, from ONE adjoint run at about half the cost.
Prefer that; this two-run path is kept because it is what produced the
committed trajectory.

MECHANICS:
  1. run_dic_A (DIC only) and run_dic_B (ALK only) have each just run
     mitgcmuv_ad against the CURRENT physical dic_init.bin/alk_init.bin.
  2. Convert each adxx field to a physical-space gradient:
         dJ/dphysical = adxx * sqrt(ctrl_weight)
     (same conversion used throughout this project, e.g.
     darwin_bling_comparison/scripts/compute_multipoint_comparison.py).
  3. Take a plain steepest-descent step, clipped per-point to +/-1 prior
     sigma (the same wt_*.bin sigma used for BLING/pkg-dic's own
     preconditioning) so no single grid point -- outlier or not -- can
     move further than one prior standard deviation in a single cycle.
  4. Write the updated physical field back into dic_init.bin/alk_init.bin
     in BOTH run_dic_A and run_dic_B (so next cycle's pkg/dic forward
     pass starts from the new state), and into Darwin's own
     dic_darwin_init.bin/alk_darwin_init.bin (mmol/m3 unit conversion).

run_dic_A/run_dic_B always run at optimcycle=0 (never incremented -- see
run_hybrid_loop_dic_manual.sh's own comment), so MITgcm always writes
adxx_ptr1/2.0000000000 regardless of which hybrid cycle this is. The
<cycle> argument here is used only for logging context, not file lookup.

Usage: manual_gradient_step.py <cycle>
"""
import sys
import numpy as np

sys.path.insert(0, "/Users/carrolld/Documents/research/ECCO/offline/scripts")
from MITgcmutils import rdmds

NZ, NY, NX = 15, 64, 128
HD = "/Users/carrolld/Documents/research/adjoint/MITgcm/verification/hybrid_darwin_dic"
RUN_A = f"{HD}/run_dic_A"
RUN_B = f"{HD}/run_dic_B"
DARWIN_RUN = f"{HD}/run_darwin"

# target: a "typical" (median-gradient) point moves by this fraction of
# its prior sigma per cycle; clip_frac bounds the MAXIMUM per-point move
# regardless of gradient magnitude, as a safety net against any
# undetected residual corruption.
TARGET_MEDIAN_STEP_FRAC = 0.3
CLIP_FRAC = 1.0


def load_field(path):
    return np.fromfile(path, dtype=">f4").reshape(NZ, NY, NX)


def main(cycle):
    hfac = rdmds(f"{RUN_A}/hFacC")
    wet = hfac[0] > 0

    wt_dic = load_field(f"{RUN_A}/wt_DIC.bin")
    wt_alk = load_field(f"{RUN_B}/wt_ALK.bin")

    adxx_dic = rdmds(f"{RUN_A}/adxx_ptr1", 0)
    adxx_alk = rdmds(f"{RUN_B}/adxx_ptr2", 0)

    g_dic = adxx_dic[0] * np.sqrt(wt_dic[0])
    g_alk = adxx_alk[0] * np.sqrt(wt_alk[0])

    for name, g, wt, path_a, path_b in [
        ("DIC", g_dic, wt_dic, f"{RUN_A}/dic_init.bin", f"{RUN_B}/dic_init.bin"),
        ("ALK", g_alk, wt_alk, f"{RUN_A}/alk_init.bin", f"{RUN_B}/alk_init.bin"),
    ]:
        median_abs_g = np.median(np.abs(g[wet]))
        sigma = wt[0][wet].mean()
        alpha = (TARGET_MEDIAN_STEP_FRAC * sigma) / median_abs_g if median_abs_g > 0 else 0.0

        raw_step = -alpha * g
        clip = CLIP_FRAC * wt[0]
        step = np.clip(raw_step, -clip, clip)
        step = np.where(wet, step, 0.0)

        current = load_field(path_a)  # dic_init.bin/alk_init.bin identical in run_dic_A/B
        updated = current.copy()
        updated[0] = current[0] + step

        updated.astype(">f4").tofile(path_a)
        updated.astype(">f4").tofile(path_b)

        darwin_out = f"{DARWIN_RUN}/{'dic' if name == 'DIC' else 'alk'}_darwin_init.bin"
        (updated * 1000.0).astype(">f4").tofile(darwin_out)

        print(f"{name}: alpha={alpha:.4g}  median|step|={np.median(np.abs(step[wet])):.4g}  "
              f"max|step|={np.abs(step[wet]).max():.4g} (sigma={sigma:.4g})  "
              f"-> {path_a}, {path_b}, {darwin_out}")


if __name__ == "__main__":
    main(int(sys.argv[1]))
