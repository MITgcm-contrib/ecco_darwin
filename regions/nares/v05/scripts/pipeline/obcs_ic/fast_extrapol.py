"""Vectorized drop-in replacement for gen_obcs.Hextrapol.

WHY THIS EXISTS
---------------
gen_obcs.py's Hextrapol() fills zero-valued wet cells with their nearest
already-nonzero wet neighbor via a legacy "one ring per wave" loop: each
wave recomputes the distance from every still-zero wet cell to every
currently-nonzero wet cell from scratch (np.where + full pairwise
distance), and only fills cells within sqrt(2) (immediate 8-neighbors).
That's O(n_waves x n_remaining x n_wet) per call -- fine when few cells
need filling (e.g. fiji's west boundary, where nearly every level has a
native parent point nearby), catastrophic when a whole level has zero
valid parent points near the boundary line (n_valid==0 in
gen_obcs._gen_obcs_batched's per-level loop) and Hextrapol must fill the
ENTIRE wet line from its own internal seeds. gen_obcs_fast.py's own
_gen_obcs_batched calls this once per timestep per level (it deliberately
kept Hextrapol/Vextrapol untouched, on the documented assumption they
"never called griddata themselves, so they were never the bottleneck") --
but the fill pattern (which cell copies from which neighbor) only depends
on the wet-mask geometry, which is IDENTICAL across all timesteps; only
the copied values differ. Recomputing it 8784 times (once per timestep) is
pure waste, and is the confirmed root cause of fiji's benchmark job
hanging for 3+ hours on UVEL/south and TRAC01/south after finishing
THETA/west in under a minute.

EXACTNESS
---------
Every boundary mask this is applied to (mask1[k, :, :] in
gen_obcs._gen_obcs_batched / gen_obcs_fast._gen_obcs_batched) has one
spatial dimension of size 1 -- it's a 1D line, not a genuine 2D field (see
gen_obcs.get_regbnd_info: west/east slice XC1/mask1 down to width 1 in one
axis, north/south the other). Hextrapol's wave process is exactly a
multi-source BFS flood-fill through 8-connected wet cells, copying each
cell's value from the wet neighbor it's reached from. In 1D there is only
one possible path between any two cells in the same connected wet run --
no detours, no shortcuts -- so "nearest by raw Euclidean distance" and
"nearest by BFS graph distance through wet cells" are identical, and a
connected-component-restricted distance_transform_edt reproduces
Hextrapol's converged output exactly, not approximately. (This does not
claim to be exact for genuinely 2D wet masks with concave/non-convex
shapes, where BFS graph distance and raw Euclidean distance can diverge --
not needed here since every real call site uses a size-1-dimension mask.)
"""

import numpy as np
from scipy import ndimage


def hextrapol_fast(var_grid, wet_grid):
    """Same signature/semantics as gen_obcs.Hextrapol(var_grid, wet_grid):
    returns (filled_grid, n_remaining_unfilled_wet_cells)."""
    wet = wet_grid == 1

    # Original's `for _ in range(n_remaining)` never executes at all when
    # nothing needs filling -- effectively free. label()/distance_transform_edt
    # below are not: measured 8x SLOWER than the original on fiji's THETA/west
    # (the common, nothing-to-fill case) without this early exit, since their
    # fixed per-call overhead dominates when there's no real work. Mirror the
    # original's zero-cost path explicitly instead of paying that overhead
    # unconditionally on every timestep x level call.
    if not np.any(wet & (var_grid == 0)):
        return var_grid, 0

    seed = wet & (var_grid != 0)
    out = var_grid.copy()

    if not seed.any():
        n_remaining = int(np.sum(wet & (out == 0)))
        return out, n_remaining

    labels, nlab = ndimage.label(wet, structure=np.ones((3, 3), dtype=int))
    for lab in range(1, nlab + 1):
        comp = labels == lab
        comp_seed = comp & seed
        if not comp_seed.any():
            continue  # this connected wet run has no data at all -- matches
            # the original's "len(Wet_Vals)==0 -> continue_iter=False" break:
            # stays zero, contributes to n_remaining below.
        _, (iy, ix) = ndimage.distance_transform_edt(~comp_seed, return_indices=True)
        targets = comp & (out == 0)
        out[targets] = out[iy[targets], ix[targets]]

    n_remaining = int(np.sum(wet & (out == 0)))
    return out, n_remaining
