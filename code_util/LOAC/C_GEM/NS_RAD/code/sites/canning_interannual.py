"""
Canning, INTERANNUAL discharge + riverine TOC forcing variant. See
sites/colville_interannual.py's module docstring for the pattern and caveats -- this
file is the same thing for Canning, with one more override: canning.py's
DISCHARGE_FILE is DISCHARGE_IS_RECONSTRUCTED (an ad-hoc Hulahula-proxy reconstruction,
since Canning has no 2022 gauge record at all). PWBM's discharge is real whole-basin
hydrological modelling of the Canning watershed itself, not a proxy scaled from a
neighboring river, so that flag is cleared here.

Still inherits EL=0 from canning.py (a deliberate placeholder pending an observed
estuary length -- see canning.py's module docstring), so this variant is NOT runnable
until that is resolved, same as canning.py today. Its forcing files are still built by
tools/build_interannual_forcings.py for consistency ("one recipe for all four rivers").

    CGEM_SITE=canning_interannual CGEM_MAXT_DAYS=<N> CGEM_WARMUP_DAYS=365 \\
        PYTHONPATH=code python code/main.py     # once EL is real
"""
from .canning import *  # noqa: F401,F403

DISCHARGE_FILE = "canning_river_discharge_interannual_1980-2023_m3sec.csv"
BOUNDARY_FORCING = {"TOC": {"cub": "canning_toc_interannual_1980-2023_mmolC_m3.csv"}}
DISCHARGE_IS_RECONSTRUCTED = False

LABEL = "Canning-interannual"
