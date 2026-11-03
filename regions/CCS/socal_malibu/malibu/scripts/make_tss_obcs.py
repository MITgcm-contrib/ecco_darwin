#!/usr/bin/env python3
"""Zero TSS (pTracer 36) open-boundary files for a Malibu work dir (2026-10-14).

Without an OB?ptrFile, pkg/obcs sets a passive tracer on the open boundary to the adjacent interior value
every step (obcs_calc.F, darwin3 24885b71, lines 314-325: OBWptr = pTracer(I_obc+1)), so inflow carries
the interior concentration back in. For TSS next to an eroding bed this is a sediment source: the
pkg/sediment monitor of run_2025h1_v3 (tides) shows a net open-boundary import of 5.9e14 g by day 180
(S5, no tides: 5.5e13 g) vs an initial domain total of 2.3e11 g. The TSS tracer holds only river and
seabed sediment and the SoCal/V5r1 parents carry none, so inflow gets TSS = 0.

Writes <work dir>/obcs/<THETA file name with TSS>, e.g. malibu_128x72x50_TSS_{West,South,East}, float32 big-endian, with the same number of records
as the THETA file of each side (obcs_files.json) (e.g. 182 daily records with make_ic_obcs.py --pad-end).

Usage: python3 make_tss_obcs.py <work dir>     (then OB?ptrFile(36) in data.obcs)
"""
import json
import os
import sys

import numpy as np


def main():
    wk = sys.argv[1]
    ob = json.load(open(os.path.join(wk, 'obcs_files.json')))
    for side in ('West', 'South', 'East'):
        ref = os.path.join(wk, 'obcs', ob['files']['THETA'][side.lower()])
        out = ref.replace('_THETA_', '_TSS_')
        if os.path.exists(out):
            sys.exit(f'{out} exists')
        np.zeros(os.path.getsize(ref) // 4, '>f4').tofile(out)
        print(out, os.path.getsize(out), 'bytes (= THETA)')


if __name__ == '__main__':
    main()
