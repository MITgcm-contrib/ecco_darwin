#!/usr/bin/env python3
"""obs_summary.csv (site x variable x source) -> obs_summary_core.csv (site x core variable).

Usage: python3 -I obs/scripts/summarize_core.py obs/obs_summary.csv obs/obs_summary_core.csv
Columns: site,variable,n,t0,t1 (years),z0,z1 (m),src (comma-joined sources).
"""
import os
import site as _site
import sys

_usp = _site.getusersitepackages()
if os.path.isdir(_usp) and _usp not in sys.path:
    sys.path.append(_usp)
import pandas as pd

CORE = ["ALK", "bbp", "bSi_flux", "Chl", "Chl_HPLC", "Chl_sat", "DIC", "DOC", "fCO2", "Fe", "NO3", "O2",
        "pCO2", "pH", "PIC_flux", "PO4", "POC", "POC_flux", "PON", "PP", "salt", "SiO2", "temp", "theta"]


def main(src, dst):
    s = pd.read_csv(src)
    s = s[s["variable"].isin(CORE)]
    g = s.groupby(["site", "variable"])
    out = g.agg(n=("n", "sum"), t0=("t0", "min"), t1=("t1", "max"), z0=("zmin", "min"), z1=("zmax", "max"),
                src=("source", lambda x: ",".join(sorted(set(x)))))
    out["t0"] = out["t0"].str[:4].astype(int)
    out["t1"] = out["t1"].str[:4].astype(int)
    out = out.reset_index()
    out = out.sort_values(["site", "variable"])
    out.to_csv(dst, index=False)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
