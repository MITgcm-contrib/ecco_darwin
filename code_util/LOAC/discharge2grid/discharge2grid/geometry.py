"""Spherical geometry helpers. Everything is done with 3-D unit vectors so the code works
across the dateline, at the poles and on any grid topology (lat-lon, curvilinear, LLC)."""
import math

import numpy as np

EARTH_R_M = 6371.0e3


def to_xyz(lon_deg, lat_deg):
    """Unit vectors on the sphere, shape (..., 3)."""
    la = np.radians(np.asarray(lat_deg, dtype=float))
    lo = np.radians(np.asarray(lon_deg, dtype=float))
    c = np.cos(la)
    return np.stack([c * np.cos(lo), c * np.sin(lo), np.sin(la)], axis=-1)


def chord_to_km(d):
    """Chord length on the unit sphere -> great-circle distance in km."""
    return 2.0 * EARTH_R_M * np.arcsin(np.clip(np.asarray(d) / 2.0, 0.0, 1.0)) / 1e3


def km_to_chord(km):
    """Great-circle distance in km -> chord length on the unit sphere."""
    th = km * 1e3 / EARTH_R_M
    return 2.0 * math.sin(min(th, math.pi) / 2.0)


def gc_km(lon1, lat1, lon2, lat2):
    """Great-circle distance (km) between points, broadcasting."""
    return chord_to_km(np.linalg.norm(to_xyz(lon1, lat1) - to_xyz(lon2, lat2), axis=-1))


def wrap_lon(lon):
    """Longitudes to [-180, 180)."""
    return (np.asarray(lon, dtype=float) + 180.0) % 360.0 - 180.0


def regular_cell_area(lon1d, lat1d):
    """Exact spherical area (m2) of the cells of a regular lat-lon grid given its cell
    CENTRES (1-D, uniformly spaced). Returns (nlat, nlon)."""
    lon1d = np.asarray(lon1d, float)
    lat1d = np.asarray(lat1d, float)
    dlon = abs(float(np.median(np.diff(lon1d)))) if lon1d.size > 1 else 360.0
    dlat = abs(float(np.median(np.diff(lat1d)))) if lat1d.size > 1 else 180.0
    s = np.radians(np.clip(lat1d - dlat / 2, -90, 90))
    n = np.radians(np.clip(lat1d + dlat / 2, -90, 90))
    band = EARTH_R_M ** 2 * np.radians(dlon) * np.abs(np.sin(n) - np.sin(s))
    return np.repeat(band[:, None], lon1d.size, axis=1)
