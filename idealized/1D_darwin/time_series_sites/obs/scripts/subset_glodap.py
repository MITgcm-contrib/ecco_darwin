#!/usr/bin/env python3
"""Subset GLODAPv2.2023 per-ocean CSV zips to within RADIUS deg of each site.

Usage: python3 -I subset_glodap.py OBS_ROOT ZIP [ZIP ...]
Writes OBS_ROOT/<SITE>/raw/glodap/GLODAPv2.2023_<ocean>_within1deg.csv
(all original columns, unmodified rows).
"""
import os
import site as _site
import sys
import zipfile

sys.path.append(_site.getusersitepackages())
import pandas as pd

RADIUS = 1.0
SITES = {"HOT": (22.75, -158.0), "BATS": (31.67, -64.17), "HydroS": (32.17, -64.5),
         "Papa": (50.1, -144.9), "PAP": (49.0, -16.5)}

root = sys.argv[1]
for zp in sys.argv[2:]:
    ocean = os.path.basename(zp).replace("GLODAPv2.2023_", "").replace(".csv.zip", "")
    with zipfile.ZipFile(zp) as z:
        name = [n for n in z.namelist() if n.endswith(".csv")][0]
        hits = {s: [] for s in SITES}
        with z.open(name) as f:
            for chunk in pd.read_csv(f, chunksize=200000, low_memory=False):
                for s, (la, lo) in SITES.items():
                    m = ((chunk["G2latitude"] - la).abs() <= RADIUS) & \
                        ((chunk["G2longitude"] - lo).abs() <= RADIUS)
                    if m.any():
                        hits[s].append(chunk[m])
    for s, parts in hits.items():
        if not parts:
            continue
        d = os.path.join(root, s, "raw", "glodap")
        os.makedirs(d, exist_ok=True)
        out = os.path.join(d, f"GLODAPv2.2023_{ocean}_within1deg.csv")
        df = pd.concat(parts)
        df.to_csv(out, index=False)
        print(s, ocean, len(df), "rows ->", out)
