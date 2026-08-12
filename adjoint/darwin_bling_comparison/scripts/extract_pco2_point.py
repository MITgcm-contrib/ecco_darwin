#!/usr/bin/env python3
"""Read the darwinPCO2Snap diagnostic output from a run dir and print the
surface pCO2 value at one grid point (j,i).

Requires MITgcmutils (ships with MITgcm under utils/python/MITgcmutils --
add that directory to PYTHONPATH before running this script).
"""
import sys
from MITgcmutils import rdmds  # noqa: E402

run_dir, j, i = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
# snapshot triggers at t=1,252,800s (Jan 15 12:00) = timestep 1252800/900=1392
d = rdmds(f'{run_dir}/darwinPCO2Snap', 1392)
# d shape: (ny, nx) since only level 1 was requested
val = float(d[j, i])
print(val)
