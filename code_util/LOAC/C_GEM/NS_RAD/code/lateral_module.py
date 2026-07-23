"""
Distributed lateral loading.

Tundra / thermokarst / tributary inputs that enter ALONG the channel rather than only
at the upstream boundary -- a first-order representation of the terrestrial-aquatic
coupling that dominates Arctic land-to-ocean carbon delivery. The active site declares
a total lateral inflow (config.LATERAL_INFLOW, spread uniformly over the domain into
variables.q_lat during init) and the chemistry of that water (config.LATERAL_CONC).

Each timestep, after transport, every declared species is nudged toward the lateral
concentration in proportion to the local volumetric inflow fraction:

    dc_i = (q_lat_i / V_i) * (c_lat - c_i) * dt ,   V_i = D_i * DELXI

This is the minimal (solute-source) form: it adds mass but does not yet grow the
downstream discharge (a hydrodynamic change left for a later pass). It is a no-op when
LATERAL_CONC is empty -- which it is for the four real rivers, so their results are
unchanged.
"""
from numba import njit

from config import M, DELXI, DELTI, LATERAL_CONC
from variables import v, D, q_lat


@njit(cache=True)
def _lateral_kernel(c, clat, q_lat, D, DELXI, DELTI, M):
    for i in range(1, M + 1):
        vol = D[i] * DELXI                       # cell water volume [m^3]
        if vol > 0.0:
            c[i] += (q_lat[i] / vol) * (clat - c[i]) * DELTI


def lateral():
    """Inject the per-cell lateral inflow as a mixing source into each declared species."""
    for name, clat in LATERAL_CONC.items():
        _lateral_kernel(v[name]['c'], float(clat), q_lat, D, DELXI, DELTI, M)
