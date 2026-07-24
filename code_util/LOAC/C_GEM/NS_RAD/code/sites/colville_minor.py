"""
Colville delta -- MINOR (western, 240 m) distributary. EXPERIMENTAL.

The other channel of the two-outlet Colville delta. Carries only 3.3% of the discharge
(width partitions as W^1.98, so an 18%-of-width channel takes 3% of the flow), and its
depth follows from that small share -- 0.82 m against the main channel's 2.24 m. The
result is a mouth velocity of 0.082 m/s, HALF the main channel's, which is the whole
point of the experiment: this is the channel that should hold salt.

See sites/_colville_delta.py for the derivation and the approximations. Not a river:
aggregate with colville_main before comparing to anything.
"""
from ._colville_delta import (  # noqa: F401
    WATERTEMP_FILE, AMPL, pfun, BOUNDARIES, BOUNDARY_CHEM_SOURCE,
    DISCHARGE_FILE, SURGE_FILE, GAUGE, DISCHARGE_IS_UPSTREAM_PROXY,
    EL, L_FLARE, distance, channel,
)

LABEL = "Colville-minor"

Q_FRACTION, B_lb, B_ub, DEPTH_lb, DEPTH_ub = channel("minor")

N_CHAN_LB = 1.0
N_CHAN_UP = 1.0

GEOMETRY_IS_PLACEHOLDER = False
IS_DISTRIBUTARY = True
