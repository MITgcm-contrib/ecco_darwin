"""
Per-river site definitions for NS-RAD.

Each module in this package supplies only the parameters that genuinely differ
between rivers -- channel geometry, tidal amplitude, the discharge forcing filename,
and the upstream/downstream boundary concentrations. Everything else (Chezy,
sediment, biogeochemistry, numerics) stays shared in config.py, so a change to the
reaction network applies to all four sites at once.

config.py selects one of these at import time via the CGEM_SITE environment
variable. See SITES below for valid keys.
"""

SITES = ("colville", "kuparuk", "sagavanirktok", "canning", "idealized",
         # EXPERIMENTAL: the two Colville delta distributaries, run as separate
         # processes and aggregated afterwards -- the multi-SITE machinery doubling as
         # multi-CHANNEL machinery. Not rivers; see sites/_colville_delta.py.
         "colville_main", "colville_minor",
         # INTERANNUAL variants: same river, but discharge + riverine TOC follow a
         # real multi-year (1980-2023) PWBM-derived daily series instead of the 2022
         # climatology -- everything else (geometry, marine boundary, met, tides) is
         # inherited unchanged. See sites/colville_interannual.py and CLAUDE.md ->
         # "Interannual forcing". canning_interannual still can't run (EL=0).
         "colville_interannual", "kuparuk_interannual", "sagavanirktok_interannual",
         "canning_interannual")
