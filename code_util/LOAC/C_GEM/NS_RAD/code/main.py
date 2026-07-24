"""
Main model routine (translated from main.c)
"""
import os
import math
from os import listdir
from typing import Union, Any

from file_module import exfread, close_output, icewrite
from fun_module import river_dispersion, set_surge
from init_module import init
from hyd_module import hyd
from transport_module import transport
from biogeo_module import biogeo
from sed_module import sed
from heat_module import heat_budget
from ice_module import ice_step
from lateral_module import lateral
from config import (M, MAXT, DELTI, WARMUP, TS, DEPTH_lb, B_lb, PI, G, EL, B_ub,
                    DELXI, SITE, SITE_LABEL, DISCHARGE_FILE, WATERTEMP_FILE, ICE_MODEL,
                    WIND_FILE, SOLAR_FILE, AIRTEMP_FILE, RELHUM_FILE, PCO2_FILE,
                    SEATEMP_FILE, BOUNDARY_FORCING, SURGE_FILE, Q_FRACTION)
from variables import dispersion, v

# Forcings live alongside the code tree, not on the original author's machine.
# Resolved relative to THIS file so the model can be run from any output directory
# (which matters -- main() deletes *.dat in the cwd and writes results there, so
# each river is run from its own directory).
FORCINGS = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "..", "forcing")


def forcing(fname):
    """Absolute path to a forcing file, with a clear error if it is missing."""
    path = os.path.normpath(os.path.join(FORCINGS, fname))
    if not os.path.exists(path):
        raise SystemExit(f"Missing forcing file: {path}")
    return path

def main():
    """Main model routine."""

    # Remove existing output files
    for fileName in listdir("."):
        # Check file extension
        if fileName.endswith('.dat'):
            # Remove File
            os.remove(fileName)

    init()  # Initialize the model

    # Reference discharge for the ice break-up trigger [m^3 s^-1]: the ANNUAL MEAN.
    # The spring freshet runs many times this, so it mechanically clears the cover;
    # winter flows stay far below it, so there is no false winter trigger. (The winter
    # low-flow percentile is ~0 on these ice-affected gauges and cannot be used.)
    import numpy as _np
    _q = _np.abs(_np.loadtxt(forcing(DISCHARGE_FILE)))
    q_ref = float(_np.mean(_q)) if _q.size else 0.0

    # Resolve forcing paths once -- they are constant for the whole run. This removes
    # eight os.path.exists() stat() calls per timestep (millions over a full run) and
    # is bit-identical: exfread receives the same absolute path string it did before
    # (and caches its parsed series on exactly that string).
    # Filenames come from config (per-site overridable). For the four real rivers the
    # met defaults are the shared regional records, so these resolve to exactly the
    # paths that were previously hardcoded here; the idealized site points them all at
    # its own analytic series.
    P_DISCHARGE = forcing(DISCHARGE_FILE)
    P_WIND = forcing(WIND_FILE)
    P_SOLAR = forcing(SOLAR_FILE)
    P_AIRTEMP = forcing(AIRTEMP_FILE)
    P_RELHUM = forcing(RELHUM_FILE)
    P_PCO2 = forcing(PCO2_FILE)
    P_WATERTEMP = forcing(WATERTEMP_FILE)
    P_SEATEMP = forcing(SEATEMP_FILE)
    # Optional wind-driven storm-surge series for the marine boundary (daily m).
    P_SURGE = forcing(SURGE_FILE) if SURGE_FILE else None

    # Time-varying solute boundaries (per-site; empty for the four real rivers, which
    # hold every solute boundary constant). Resolve each (species, end) -> absolute
    # path once here; the loop refreshes v[species][end] from it every timestep, the
    # same way the two temperature boundaries are driven. See config.BOUNDARY_FORCING.
    BC_FORCING = [(sp, end, forcing(fn))
                  for sp, ends in BOUNDARY_FORCING.items()
                  for end, fn in ends.items()]

    # Main simulation loop
    for t in range(0, MAXT + 1, DELTI):
        # Print the number of simulated days
        print(f"t: {float(t) / (24.0 * 60.0 * 60.0):.2f} days")

        # Read forcing files
        # River discharge (per-site; selected by CGEM_SITE)
        Qr, previousdays = exfread(P_DISCHARGE, t)
        # Q_FRACTION is 1.0 for a whole river and is the conveyance share for a delta
        # distributary run (sites/colville_main.py etc.) -- see config.Q_FRACTION.
        Qr = - Qr * Q_FRACTION

        # Dispersion coefficient [m^2/s], from LOCAL width/depth/velocity.
        # Replaces the Van der Burgh / Savenije estuary form, which under
        # observation-based geometry clamped to zero after a few grid points and
        # left transport purely advective. See config.DISPERSION_MODEL.
        river_dispersion(Qr)

        # Wind. Both the saline-zone and tidal-river wind read the same series --
        # there is only one regional record (NDBC PRDA2), so the split in
        # fun_module.piston_velocity is currently fed identical values.
        Uw_sal, previousdays = exfread(P_WIND, t)
        Uw_tid, previousdays = exfread(P_WIND, t)

        # Solar radiation
        I0, previousdays = exfread(P_SOLAR, t)
        # Air Temperature -- drives the surface heat budget (heat_module)
        air_temp, previousdays = exfread(P_AIRTEMP, t)
        # Relative humidity for the latent-heat term, OBSERVED at Deadhorse Airport
        # (colocated with PRDA2), replacing the assumed constant. See tools/build_humidity.py.
        rel_hum, previousdays = exfread(P_RELHUM, t)
        # pCO2  (note the capital B -- the file is pCO2_Barrow_2022.csv, which the
        # original lowercase spelling only resolved on case-insensitive filesystems)
        pCO2, previousdays = exfread(P_PCO2, t)  #microatm
        pCO2 = pCO2 * 1e-6  # atm
        # Water Temperature -- MUST stay last: the ice gate downstream uses the
        # `previousdays` returned here, and every call overwrites it.
        #
        # This is MODELLED RIVER temperature (tools/build_river_temp.py), not the
        # shipped watertemp.csv, which was PRDA2 sea-water temperature -- it biased
        # every rate cold by 6-8 C and gated the model off coastal sea ice, opening
        # on day 173 when Colville and Kuparuk peak on day 153. See CLAUDE.md.
        water_temp, previousdays = exfread(P_WATERTEMP, t)

        # Water temperature is now a TRANSPORTED field (variables.v['T']), not a
        # scalar. Refresh its two boundary values each timestep -- openbound() reads
        # them per call, so this is all that is needed to drive it:
        #   downstream = SEA temperature (PRDA2 buoy, watertemp.csv)
        #   upstream   = RIVER temperature (modelled, WATERTEMP_FILE)
        sea_temp, _ = exfread(P_SEATEMP, t)
        v['T']['clb'] = sea_temp     # marine end-member
        v['T']['cub'] = water_temp   # riverine end-member

        # Time-varying solute boundaries (per-site; the list is empty for the four
        # real rivers, so this is a no-op for them). Refresh each declared downstream
        # (clb) / upstream (cub) boundary from its forcing series; openbound reads
        # these fresh on the next transport call, so this is all that is required.
        for _sp, _end, _path in BC_FORCING:
            _bv, _ = exfread(_path, t)
            v[_sp][_end] = _bv

        # Wind-driven storm surge at the marine boundary: refresh the surge offset that
        # fun_module.Tide adds to the harmonic elevation (0 if the site has no surge
        # file). MUST precede hyd, which reads Tide(t) for the downstream water level.
        if P_SURGE is not None:
            _surge, _ = exfread(P_SURGE, t)
            set_surge(_surge)

        # Run hydrodynamics and transport processes
        hyd(t, Qr)
        transport(t, previousdays)

        # Distributed lateral loading: inject tundra/thermokarst inputs ALONG the
        # channel as a mixing source, AFTER transport has moved the upstream-boundary
        # water. A no-op unless the active site sets LATERAL_INFLOW/LATERAL_CONC
        # (the four real rivers leave it off, so their results are unchanged).
        lateral()

        # Surface heat budget: warm/cool the transported temperature field from the
        # atmosphere, AFTER advection+dispersion have moved it. Uses the saline-zone
        # wind (Uw_sal); at these sites both wind series are identical anyway. This
        # is what makes the interior temperature physical rather than a pure mixing
        # tracer between the two boundaries. See heat_module for its assumptions.
        heat_budget(t, air_temp, Uw_sal, I0, previousdays, rel_hum)

        # Advance the prognostic ice cover: grow from the freeze energy heat_budget just
        # booked, thicken by conduction, melt in spring, and break up on the freshet.
        # Must follow heat_budget (consumes its per-step deficit) and precede biogeo
        # (which reads ice_frac to shut off gas exchange and dim under-ice light).
        ice_step(t, air_temp, Uw_sal, I0, Qr, q_ref)

        # Write the ice cover at the same cadence as the other state fields.
        if ICE_MODEL and (float(t) / float(TS * DELTI)) % 1 == 0:
            icewrite(t)

        # Run biogeochemical and sediment processes after warmup
        if t > WARMUP:
            biogeo(t, Uw_sal, Uw_tid, v['T']['c'], pCO2, I0, previousdays)
            sed(t, previousdays)

    # Flush/close the NetCDF output (no-op for the .dat path)
    close_output()

if __name__ == "__main__":
    main()
