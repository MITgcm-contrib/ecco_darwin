"""
Colville delta -- MAIN (eastern, 1310 m) distributary. EXPERIMENTAL.

One of the two channels of the two-outlet Colville delta, run as its own process.
Carries 96.7% of the discharge. See sites/_colville_delta.py for the derivation, the
approximations, and why this experiment exists at all.

Not a river: aggregate with colville_minor (weighted by surface area) before comparing
to anything. Compare against runs/definitive/colville for the lumped equivalent.
"""
from ._colville_delta import (  # noqa: F401
    WATERTEMP_FILE, AMPL, pfun, BOUNDARIES, BOUNDARY_CHEM_SOURCE,
    DISCHARGE_FILE, SURGE_FILE, GAUGE, DISCHARGE_IS_UPSTREAM_PROXY,
    EL, L_FLARE, distance, channel,
)

LABEL = "Colville-main"

Q_FRACTION, B_lb, B_ub, DEPTH_lb, DEPTH_ub = channel("main")

# Single-thread channel: the distributaries are resolved explicitly here, so there is
# no thread count left to lump. These make the site inert under CGEM_MULTICHANNEL=on.
N_CHAN_LB = 1.0
N_CHAN_UP = 1.0

GEOMETRY_IS_PLACEHOLDER = False
IS_DISTRIBUTARY = True
