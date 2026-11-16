"""Checks on a multi-year run.

scan_dropouts: find rivers that vanish from the sources in some years.

This is how the intermittent-outlet problem of the pan_AO / WestAO input was found
(Olfusa, Iceland, ~260 m3/s, missing in 11 of 36 years; Pur, ~1,000 m3/s, missing in
1999). Sources from all years are merged into rivers (within merge_km, largest first);
a GAP is a year whose mean discharge is below gap_frac x the river's median yearly mean.
For each gap, `lost_frac` compares the drop with the change of ALL water within
deficit_km: ~1 means the water is really missing; ~0 means it only moved to a nearby
cell (e.g. the Gaula, Norway: an outlet that alternates between cells -- filling it would
count it twice).
"""
import numpy as np

from .geometry import gc_km


def yearly_table(sources):
    """(lon, lat, mean_q) of a SourceSet: what scan_dropouts needs per year."""
    return np.asarray(sources.lon, float), np.asarray(sources.lat, float), sources.mean_q()


def scan_dropouts(tables, merge_km=8.0, min_q=20.0, gap_frac=0.3, deficit_km=30.0,
                  keep_frac=0.5):
    """tables: {year: (lon, lat, mean_q)}. Returns a list of dicts, largest river first:
    lon, lat, median_q, gap_years, lost_frac (median over gaps), keep (lost_frac >= keep_frac)."""
    years = sorted(tables)
    lo = np.concatenate([np.asarray(tables[y][0], float) for y in years])
    la = np.concatenate([np.asarray(tables[y][1], float) for y in years])
    qq = np.concatenate([np.asarray(tables[y][2], float) for y in years])
    yi = np.concatenate([np.full(len(tables[y][0]), k) for k, y in enumerate(years)])
    key = np.round(lo, 3) * 1e4 + np.round(la, 3)
    ukey, inv = np.unique(key, return_inverse=True)
    ulo = np.zeros(ukey.size)
    ula = np.zeros(ukey.size)
    ulo[inv], ula[inv] = lo, la
    Qs = np.zeros((ukey.size, len(years)))              # per position, per year
    np.add.at(Qs, (inv, yi), qq)
    tot = Qs.sum(axis=1)
    centre = -np.ones(ukey.size, dtype=np.int64)
    for i in np.argsort(-tot):
        if centre[i] >= 0:
            continue
        free = centre < 0
        centre[free & (gc_km(ulo[i], ula[i], ulo, ula) <= merge_km)] = i
    rid = np.unique(centre)
    pos = {c: k for k, c in enumerate(rid)}
    Qr = np.zeros((rid.size, len(years)))
    for s in range(ukey.size):
        Qr[pos[centre[s]]] += Qs[s]
    out = []
    for k, c in enumerate(rid):
        med = float(np.median(Qr[k]))
        if med < min_q:
            continue
        gaps = np.flatnonzero(Qr[k] < gap_frac * med)
        if gaps.size == 0:
            continue
        near = gc_km(ulo[c], ula[c], ulo, ula) <= deficit_km
        tsum = Qs[near].sum(axis=0)
        tmed = float(np.median(tsum))
        lost = [(tmed - tsum[g]) / max(med - Qr[k, g], 1e-9) for g in gaps]
        lf = float(np.median(lost))
        out.append(dict(lon=float(ulo[c]), lat=float(ula[c]), median_q=med,
                        gap_years=[years[g] for g in gaps], lost_frac=lf, keep=lf >= keep_frac))
    out.sort(key=lambda r: -r["median_q"])
    return out
