"""discharge2grid -- river discharge datasets onto ocean-model grids (MITgcm, ECCO LLC, any
curvilinear grid), conserving volume; optionally with river water temperature
(flux-weighted, heat-conserving).

    from discharge2grid import ModelGrid, Mapper, SourceSet
    grid  = ModelGrid.from_netcdf("grid.nc")
    src   = ...                               # a reader in discharge2grid.readers
    deliv = Mapper(grid, sigma_km=18, max_snap_km=95, exclude_outside=True).map_sources(src)
    deliv.report()
    write_mitgcm_runoff("runoff_X_2015", deliv) / write_sparse_netcdf("X_2015.nc", deliv)

See README.md (quick start) and docs/METHODS.md (what it does and why).
"""
from .grid import ModelGrid
from .sources import SourceSet, daily_axis
from .mapping import Mapper, Delivery, Snap, KEPT, TOO_FAR, OUTSIDE
from .writers import (write_mitgcm_runoff, write_mitgcm_temperature, write_sparse_netcdf,
                      read_sparse_netcdf)

__version__ = "0.2.0"
__all__ = ["ModelGrid", "SourceSet", "daily_axis", "Mapper", "Delivery", "Snap", "KEPT",
           "TOO_FAR", "OUTSIDE", "write_mitgcm_runoff", "write_mitgcm_temperature", "write_sparse_netcdf",
           "read_sparse_netcdf"]
