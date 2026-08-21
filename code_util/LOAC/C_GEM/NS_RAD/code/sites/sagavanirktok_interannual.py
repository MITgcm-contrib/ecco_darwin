"""
Sagavanirktok, INTERANNUAL discharge + riverine TOC forcing variant. See
sites/colville_interannual.py's module docstring for the pattern and caveats -- this
file is the same thing for Sagavanirktok. `import *` inherits every parameter from
sites/sagavanirktok.py unchanged; sagavanirktok.py itself and
runs/definitive/sagavanirktok are unaffected.

    CGEM_SITE=sagavanirktok_interannual CGEM_MAXT_DAYS=<N> CGEM_WARMUP_DAYS=365 \\
        PYTHONPATH=code python code/main.py
"""
from .sagavanirktok import *  # noqa: F401,F403

DISCHARGE_FILE = "sagavanirktok_river_discharge_interannual_1980-2023_m3sec.csv"
BOUNDARY_FORCING = {"TOC": {"cub": "sagavanirktok_toc_interannual_1980-2023_mmolC_m3.csv"}}

LABEL = "Sagavanirktok-interannual"
