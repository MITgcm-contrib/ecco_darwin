"""Dataset readers. Each returns a SourceSet (point sources with discharge in m3 s-1).

glofas   GloFAS v4 raw files: outlet finding, dead ends, fixed-outlet extraction;
         plus the legacy "runoff map" reader used for the pan_AO/WestAO forcing
points   river/station tables (CSV or netCDF)
gridded  coastal runoff fields (e.g. JRA55-do friver)
temperature  river water temperature: daily climatology / dated series from gridded
         products, sampled at the sources; filling of unknown temperatures

To add a dataset, write a function that returns a SourceSet (docs/INSTRUCTIONS.md, step 4).
"""
from . import glofas, points, gridded, temperature  # noqa: F401
