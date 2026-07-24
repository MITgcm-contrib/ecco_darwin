"""
Shared derivation for the two Colville delta distributary runs (EXPERIMENTAL).

WHY THIS EXISTS. C-GEM is a single-channel solver, so a two-outlet delta can only be
represented as one lumped channel of the summed width. That is fine for bulk conveyance
but it cannot produce a RANGE: salt intrusion is nonlinear in velocity, so the intrusion
of the mean channel is not the mean intrusion across channels. Observed mouth salinity on
these rivers is 8-32 PSU in individual grabs while the lumped model sits near 0 in summer
-- and width, storm surge and the Chezy/dispersion `distance` gradient have all been
tested and all fail to close it (see CLAUDE.md). Resolving the distributaries is the one
remaining lever with the right physics, and it needs NO solver change: the model already
runs one process per site, so each distributary is simply its own site.

THE PARTITION IS DERIVED, NOT ASSUMED. Manning with a shared slope and roughness gives
Q_i ~ W_i * D_i^(5/3), and Colville's own at-a-station hydraulic geometry gives
D = 0.360*Q^0.297. Eliminating D:

    Q_i^(1 - 5f/3) ~ W_i    ->    Q_i ~ W_i^(1/(1-5f/3)) = W_i^1.98

so discharge partitions as very nearly the SQUARE of width, not linearly with it. The
minor channel is 18% of the mouth width but takes only 3.3% of the flow -- and because
its depth then follows from its own (small) Q, it ends up shallow. The result is the
mechanism this experiment is testing:

    channel   W(m)   Qfrac    Q      D(m)    U_mouth
    main      1310   0.967   468.2   2.236   0.160 m/s
    minor      240   0.034    16.3   0.824   0.082 m/s   <- HALF the main-channel velocity

A lumped 1550 m channel cannot express that spread.

UPSTREAM = STREAM TUBE. Above the delta apex the two distributaries are one trunk, so
each run's prismatic width is its conveyance SHARE of the trunk (B_ub_total * Q_FRACTION)
and its upstream depth is the trunk depth. Both channels then carry the same upstream
velocity as each other and as the multichannel-corrected lumped run (0.205 m/s at the
2022 open-water mean), which is what makes the three runs comparable: they differ ONLY
in the delta reach. That velocity follows the trunk width, so it moves if B_UB_TOTAL is
revised -- it is a consistency property of the construction, not a tuned number.

APPROXIMATION TO KEEP IN VIEW. Each distributary is extended over the full 27 km domain
rather than stopping at the delta apex, because C-GEM needs an upstream boundary. The
stream-tube construction makes that reach mass- and velocity-consistent with the trunk,
but it is not a real bifurcating network -- the two channels cannot exchange water. This
is a diagnostic experiment, not a delta model.
"""
from .colville import (  # noqa: F401
    WATERTEMP_FILE, AMPL, pfun, BOUNDARIES, BOUNDARY_CHEM_SOURCE,
    DISCHARGE_FILE, SURGE_FILE, GAUGE, DISCHARGE_IS_UPSTREAM_PROXY,
    EL, L_FLARE, distance,
)

# --- the trunk this delta is fed by (sites/colville.py) ----------------------
_B_UB_PERCHAN = 423.0     # SWORD per-channel prismatic width [m]
# OBSERVED braided total beyond the flare (SWORD v17b raw median, 102 nodes), NOT
# B_ub x a braid count -- see sites/colville.py and config.MULTICHANNEL for why the
# reconstruction overstates it (1052 observed vs 1269 from 423 x median n_chan_mod 3).
_B_UB_TOTAL = 1052.0
_DEPTH_TRUNK = 2.25       # trunk depth [m], = colville DEPTH_lb/ub
_QBAR = 484.5             # 2022 open-water mean discharge [m3/s] (the depth-relation argument)

# --- at-a-station hydraulic geometry, D = c * Q^f (208 USGS ADCP surveys) ----
_C_D, _F_D = 0.360, 0.297

# --- SWORD v17b seaward distributary widths at Harrison Bay ------------------
# These are the two entries of chan_widths in docs/sword_widths.json; they sum to
# colville.B_lb = 1550 m, which is exactly what the lumped run uses as its mouth width.
CHAN_WIDTHS = {"main": 1310.0, "minor": 240.0}

_P = 1.0 / (1.0 - 5.0 * _F_D / 3.0)          # = 1.980, the width exponent for Q
_WP = {k: w ** _P for k, w in CHAN_WIDTHS.items()}
_WPSUM = sum(_WP.values())
Q_SHARE = {k: v / _WPSUM for k, v in _WP.items()}


def channel(name):
    """Geometry for one distributary: (Q_FRACTION, B_lb, B_ub, DEPTH_lb, DEPTH_ub)."""
    f = Q_SHARE[name]
    return (f,
            CHAN_WIDTHS[name],          # B_lb  = its own SWORD mouth width
            _B_UB_TOTAL * f,            # B_ub  = its conveyance share of the trunk
            _C_D * (f * _QBAR) ** _F_D,  # DEPTH_lb from its OWN discharge share
            _DEPTH_TRUNK)               # DEPTH_ub = the shared trunk depth
