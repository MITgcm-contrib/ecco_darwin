"""
Kuparuk, INTERANNUAL discharge + riverine TOC forcing variant. See
sites/colville_interannual.py's module docstring for the pattern and caveats -- this
file is the same thing for Kuparuk. `import *` inherits every parameter from
sites/kuparuk.py unchanged; kuparuk.py itself and runs/definitive/kuparuk are
unaffected.

    CGEM_SITE=kuparuk_interannual CGEM_MAXT_DAYS=<N> CGEM_WARMUP_DAYS=365 \\
        PYTHONPATH=code python code/main.py
"""
from .kuparuk import *  # noqa: F401,F403

DISCHARGE_FILE = "kuparuk_river_discharge_interannual_1980-2023_m3sec.csv"
BOUNDARY_FORCING = {"TOC": {"cub": "kuparuk_toc_interannual_1980-2023_mmolC_m3.csv"}}

LABEL = "Kuparuk-interannual"
