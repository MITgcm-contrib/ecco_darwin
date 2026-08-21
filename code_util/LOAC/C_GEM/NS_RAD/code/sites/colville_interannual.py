"""
Colville, INTERANNUAL discharge + riverine TOC forcing variant. Not a different river --
`import *` inherits every parameter from sites/colville.py unchanged (geometry, marine
boundary, tides, all other BOUNDARIES) and overrides only the two things this variant
exists for. `colville.py` itself, runs/definitive/colville, and every report/validation
tool are completely unaffected by this file's existence.

DISCHARGE_FILE and BOUNDARY_FORCING["TOC"]["cub"] now point at multi-year (1980-2023)
daily series built from PWBM by tools/build_interannual_forcings.py, instead of the
2022 climatology colville.py uses. Everything else -- met forcing, tides, marine
boundary, ice -- is UNCHANGED and still repeats a single 2022-typical year: interannual
variability here is deliberately isolated to river discharge + upstream DOC loading.
See CLAUDE.md -> "Interannual forcing" for the full derivation, unit-conversion
methodology, and the DOC-only/no-POC caveat.

    CGEM_SITE=colville_interannual CGEM_MAXT_DAYS=<N> CGEM_WARMUP_DAYS=365 \\
        PYTHONPATH=code python code/main.py

or tools/run_interannual.sh colville. Writes to runs/interannual/colville/ by convention
(not enforced by this file -- see the run wrapper).
"""
from .colville import *  # noqa: F401,F403

DISCHARGE_FILE = "colville_river_discharge_interannual_1980-2023_m3sec.csv"
BOUNDARY_FORCING = {"TOC": {"cub": "colville_toc_interannual_1980-2023_mmolC_m3.csv"}}

LABEL = "Colville-interannual"
