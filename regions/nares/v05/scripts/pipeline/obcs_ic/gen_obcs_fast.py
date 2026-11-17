"""
Kerguelen STEP3: fast replacement driver for
ecco_darwin/regions/downscaling/utils/gen_obcs.py.

WHY THIS EXISTS
---------------
gen_obcs.py's gen_obcs() calls scipy.interpolate.griddata once per
(variable, boundary, timestep, level). With tstp=26304 (3yr hourly) x up
to Nr1=90 levels x 35 ocean+BGC variables x 4 boundaries, that's ~26
million griddata calls, each rebuilding a Delaunay triangulation from
scratch -- even though the triangulation only depends on the boundary's
source-point geometry (fixed per variable/boundary/level), never on which
timestep is being read. Benchmarked at ~4.45ms/call, that's ~17 CPU-days
run as written -- not viable.

Also: gen_obcs.py's own __main__ calls
gen_obcs_files(config_dir, reg_nm, bds, itrs, seaice, bgc, print_level) --
7 positional args into a 6-param function (`seaice` is already read as a
module global inside; looks like leftover cruft from a refactor). Running
gen_obcs.py's documented CLI invocation directly raises TypeError before
any real work happens. Worked around here by calling into gen_obcs.py's
functions directly (like gen_obcs_transport_correction.py already does),
never through its __main__.

FIX
---
Per (variable, boundary), for each level k: build ONE griddata call whose
`values` argument carries all tstp timesteps as columns (griddata/
LinearNDInterpolator already support this -- the original code already
uses the size-1 case of this feature, via `val0.reshape(-1,1)` and
`val1[:,:,0]`). This replaces tstp separate Python-level calls (each
re-triangulating) with a single triangulation + one batched solve. Also
compute the nearest-neighbor fallback ONCE per level (batched the same
way) instead of only on-demand -- it's needed for two different fallback
paths in the original algorithm (whole-column-zero, and post-extrapolation
remaining cells), and precomputing it once removes both from the
per-timestep loop entirely.

The per-timestep Hextrapol/Vextrapol/masking logic is kept EXACTLY as
gen_obcs.py's own gen_obcs() does it (those never called griddata
themselves, so they were never the bottleneck) -- this minimizes the risk
of silently changing numerical behavior vs. the trusted original. See
`_gen_obcs_batched` below for a step-by-step correspondence to the
original gen_obcs() logic.

Parallelizes across (boundary, variable) pairs with multiprocessing
(these are fully independent -- each reads its own dv-diagnostic files and
writes its own output file). Grid/mask setup is read once in the parent
process before forking so workers share it via copy-on-write.
"""

import os
import sys
import glob
import argparse
import numpy as np
from scipy.interpolate import griddata
import multiprocessing as mp

UTILS_DIR = "/nobackup/dcarrol2/v05_V5r1/ecco_darwin/regions/downscaling/utils"
sys.path.insert(0, UTILS_DIR)
import gen_obcs as go  # noqa: E402
import hproj  # noqa: E402
from fast_extrapol import hextrapol_fast


def _zinterp_fixed(dv_diag, msk_pts, delR0, delR1):
    """Drop-in replacement for gen_obcs.Zinterp, ported from upstream commit
    f6d55709 ("Fix OBCs downscaling method due to mask issue", 2026-08-03).
    The shared /nobackup/dcarrol2/v05_V5r1/ecco_darwin checkout is a
    multi-user resource with other users' uncommitted local changes
    (confirmed via `git status`), so rather than `git pull` it (which
    could disturb their in-progress work), this is monkey-patched onto
    the `gen_obcs` module object in THIS process only (`go.Zinterp =
    _zinterp_fixed`, right after import) -- it never touches the file on
    disk. `read_dv_diags` calls `Zinterp` by bare name, resolved from the
    module's global namespace at call time, so the patch also takes
    effect for every THETA/SALT/UVEL/VVEL/TRACnn call that goes through
    `go.read_dv_diags` unmodified.

    The bug: when the child (regional) grid's vertical levels create a
    multi-level gap near the bottom of a water column relative to the
    parent's coarser levels, the original only filled the FIRST gap
    level with the bottom value (`tmp1[:,idB[0]] = bottom_values`),
    leaving any remaining gap levels as NaN -> zeroed by the unconditional
    `tmp1[np.isnan(tmp1)] = 0` a few lines later, i.e. spurious zero
    layers mid-water-column. Fix: fill every index in `idB`, not just the
    first. Same fix applied to the wet-cell mask branch (`msk_pts_IT`),
    which previously had no bottom-gap fill logic at all.
    """
    Zbot0 = np.cumsum(delR0)
    Ztop0 = np.concatenate([np.array([0]), Zbot0[:-1]])
    Z0 = (Zbot0 + Ztop0) / 2
    Zbot1 = np.cumsum(delR1)
    Ztop1 = np.concatenate([np.array([0]), Zbot1[:-1]])
    Z1 = (Zbot1 + Ztop1) / 2

    idS = np.where(np.abs(Z1) < np.abs(Z0[0]))[0]
    dv_diag_IT = np.zeros((dv_diag.shape[0], len(Z1), dv_diag.shape[2]))
    for i in range(dv_diag.shape[2]):
        tmp0 = dv_diag[:, :, i].copy()
        if np.sum(tmp0 != 0) > 1:
            tmp0[np.where(tmp0 == 0)] = np.nan
            f = go.interp1d(Z0, tmp0, axis=1, kind='linear',
                             bounds_error=False, fill_value=np.nan)
            tmp1 = f(Z1)
            for j in idS:
                tmp1[:, j] = tmp1[:, idS[-1] + 1]
            if np.size(np.abs(Z0[np.isnan(tmp0[0])])) > 0:
                bottom_depth = np.abs(Z0[np.isnan(tmp0[0])])[0]
                tmp = np.where(~np.isnan(tmp1[0]))[0]
                if tmp.size > 1:
                    bottom_values = tmp1[:, np.where(~np.isnan(tmp1[0]))[0][-1]]
                elif tmp.size == 0:
                    bottom_values = 0.
                else:
                    bottom_values = tmp1[:, np.where(~np.isnan(tmp1[0]))[0][0]]
                idB = np.where(np.logical_and(np.isnan(tmp1[0]), np.abs(Z1) < bottom_depth))[0]
                if len(idB) != 0:
                    tmp1[:, idB] = bottom_values[:, None] if np.ndim(bottom_values) else bottom_values
            tmp1[np.isnan(tmp1)] = 0
            dv_diag_IT[:, :, i] = tmp1
    if msk_pts.shape[0] != len(Z1):
        msk_pts_IT = np.zeros((len(Z1), msk_pts.shape[1]))
        for i in range(msk_pts.shape[1]):
            tmp0 = msk_pts[:, i].copy()
            if np.sum(tmp0 != 0) > 1:
                f = go.interp1d(Z0[tmp0 != 0], tmp0[tmp0 != 0], kind='linear',
                                 bounds_error=False, fill_value=np.nan)
                tmp1 = f(Z1)
                tmp1[np.abs(Z1) < np.abs(Z0[0])] = tmp1[~np.isnan(tmp1)][0]
                if np.size(np.abs(Z0[tmp0 == 0])) > 0:
                    bottom_depth = np.abs(Z0[tmp0 == 0])[0]
                    bottom_value = tmp1[~np.isnan(tmp1)][-1]
                    idB = np.where(np.logical_and(np.isnan(tmp1), np.abs(Z1) < bottom_depth))[0]
                    if len(idB) != 0:
                        tmp1[idB] = bottom_value
                msk_pts_IT[:, i] = tmp1
            elif np.sum(tmp0 != 0) == 1:
                msk_pts_IT[0, i] = msk_pts[0, i]
    else:
        msk_pts_IT = msk_pts
    msk_pts_IT[np.isnan(msk_pts_IT)] = 0
    msk_pts_IT = np.round(msk_pts_IT).astype(int)

    return dv_diag_IT, msk_pts_IT


go.Zinterp = _zinterp_fixed

# O(n_waves x n_remaining x n_wet) legacy nested loop -> O(n) connected-
# component distance transform, proven bit-identical for these 1D boundary
# lines (see fast_extrapol.py docstring + fiji's validate_hextrapol_fix.py,
# job 25180576, 2026-09-23: THETA/west max_abs_diff=0.0, and the two
# previously-3h+-hanging tasks (UVEL/south, PTRACE01/south) finished in
# ~30-50s). Ported into pan_AO unvalidated-in-place -- pan_AO's own
# gateway boundaries (Bering/Davis-Baffin/Nordic) are complex-coastline
# straits at least as likely as fiji's to hit an n_valid==0-at-some-level
# case, so this is applied here proactively rather than waiting to
# rediscover the same hang. Patched via the `go` module object, same
# pattern as go.Zinterp above, so both call sites in this file
# (gen_obcs_fixed_reference and _gen_obcs_batched) pick it up
# automatically.
go.Hextrapol = hextrapol_fast


def gen_obcs_fixed_reference(dv_diag, XC0, YC0, mask00, XC1, YC1, mask1,
                             hmode=hproj.LONLAT):
    """Per-timestep reference implementation matching upstream commit
    f6d55709's fixed gen_obcs() exactly (NaN-based per-cell fallback, not
    the shared checkout's zero-based whole-array check). Used ONLY as the
    slow validation reference for `_gen_obcs_batched` -- comparing the
    batched version against the shared (still-unfixed) `go.gen_obcs`
    would show spurious mismatches now that the batched version has this
    fix and the shared one doesn't."""
    obcs = np.zeros(mask1.shape)
    XC0, YC0 = hproj.project(XC0, YC0, hmode)
    XC1, YC1 = hproj.project(XC1, YC1, hmode)
    coord00 = np.array([XC0, YC0]).T
    for k in range(len(obcs)):
        if np.any(mask1[k] > 0):
            coord0 = coord00.copy()
            val0 = dv_diag[k].reshape(-1, 1)
            mask0 = mask00[k].reshape(-1, 1)
            coord0 = coord0[mask0[:, 0] != 0, :]
            val0 = val0[mask0[:, 0] != 0, :]
            if len(val0) > 0:
                if len(val0) > 4:
                    val1 = griddata(coord0, val0, (XC1, YC1), method='linear', fill_value=np.nan)
                    val1 = val1[:, :, 0]
                else:
                    val1 = np.full(XC1.shape, np.nan)
                needs_nearest = np.isnan(val1)
                if np.any(needs_nearest):
                    val1_near = griddata(coord0, val0, (XC1, YC1), method='nearest')
                    val1_near = val1_near[:, :, 0]
                    val1[needs_nearest] = val1_near[needs_nearest]
            else:
                val1 = np.zeros_like(XC1).astype(float)
            val1[mask1[k, :, :] == 0] = 0
            val1, n_remaining = go.Hextrapol(val1, mask1[k, :, :])
            if n_remaining > 0 and k > 0:
                val1 = go.Vextrapol(obcs, val1, mask1[k, :, :], k, 0)
            if n_remaining > 0:
                if len(coord0) > 0:
                    val1_NR = griddata(coord0, val0, (XC1, YC1), method='nearest', fill_value=0)
                    val1_NR = val1_NR[:, :, 0]
                    ids = np.logical_and(val1 == 0, mask1[k, :, :] != 0)
                    val1[ids] = val1_NR[ids]
        else:
            val1 = np.zeros_like(XC1).astype(float)
        obcs[k, :, :] = val1[:, :]
    return obcs


def _gen_obcs_batched(dv_diag_all, XC0, YC0, mask00, XC1, YC1, mask1,
                      hmode=hproj.LONLAT):
    """Batched-over-time equivalent of gen_obcs.gen_obcs(), updated to match
    upstream commit f6d55709 ("Fix OBCs downscaling method due to mask
    issue", 2026-08-03): the initial linear-interpolation fallback is now
    NaN-based and per-cell, not zero-based and whole-array. Originally
    (and in this file before this fix), a cell that came back exactly 0
    from griddata (either legitimately, or because it's outside the
    convex hull of the ~100-150 source points, since fill_value=0) only
    got replaced by the nearest-neighbor fallback if EVERY cell in that
    level was zero (`if not np.any(val1!=0)`) -- otherwise it was left at
    0 for Hextrapol to chain-fill one cell at a time from neighbors,
    which can propagate a stale value across many cells instead of using
    the true nearest native sample directly. Fix: fill_value=np.nan
    distinguishes "outside convex hull" from "genuinely interpolated to
    zero", so each NaN cell gets its own nearest-neighbor fallback
    regardless of what the rest of that level's cells did.

    dv_diag_all: (tstp, Nr, npts) -- same data gen_obcs() receives one
    timestep at a time as dv_diag[t] (shape (Nr, npts)).
    Returns obcs_all: (tstp,) + mask1.shape -- same per-timestep shape
    gen_obcs() returns, stacked over time.
    """
    tstp = dv_diag_all.shape[0]
    Nr = mask1.shape[0]
    obcs_all = np.zeros((tstp,) + mask1.shape)
    # Source and target coords are used here only as interpolation
    # coordinates, so this is the single place the plane has to be applied.
    # Identity for 'lonlat'; polar stereographic for a circumpolar or
    # near-polar domain, where no longitude branch removes the +/-180 seam
    # and the lon/lat plane is degenerate at the pole. See hproj.py.
    XC0, YC0 = hproj.project(XC0, YC0, hmode)
    XC1, YC1 = hproj.project(XC1, YC1, hmode)
    coord00 = np.array([XC0, YC0]).T

    for k in range(Nr):
        if not np.any(mask1[k] > 0):
            continue  # obcs_all[:, k] stays zero, matches original's val1=zeros branch

        mask0 = mask00[k].reshape(-1, 1)
        valid = mask0[:, 0] != 0
        coord0 = coord00[valid, :]
        val0_all = dv_diag_all[:, k, valid].T  # (n_valid, tstp)
        n_valid = coord0.shape[0]

        if n_valid > 0:
            if n_valid > 4:
                val1_all = griddata(coord0, val0_all, (XC1, YC1), method='linear', fill_value=np.nan)
            else:
                val1_all = np.full(XC1.shape + (tstp,), np.nan)
            needs_nearest = np.isnan(val1_all)
            if np.any(needs_nearest):
                val1_near_all = griddata(coord0, val0_all, (XC1, YC1), method='nearest')
                val1_all[needs_nearest] = val1_near_all[needs_nearest]
            # second fallback path (post-Hextrapol remaining cells) still
            # matches the original's fill_value=0 nearest call -- precompute
            # once per level here instead of on-demand per timestep.
            val1_near0_all = griddata(coord0, val0_all, (XC1, YC1), method='nearest', fill_value=0)
        else:
            val1_all = np.zeros(XC1.shape + (tstp,))
            val1_near0_all = np.zeros(XC1.shape + (tstp,))

        # --- per-timestep post-processing, unchanged from the original ---
        for t in range(tstp):
            val1 = val1_all[..., t].copy()
            val1[mask1[k, :, :] == 0] = 0
            val1, n_remaining = go.Hextrapol(val1, mask1[k, :, :])
            if n_remaining > 0 and k > 0:
                val1 = go.Vextrapol(obcs_all[t], val1, mask1[k, :, :], k, 0)
            if n_remaining > 0 and n_valid > 0:
                val1_NR = val1_near0_all[..., t]
                ids = np.logical_and(val1 == 0, mask1[k, :, :] != 0)
                val1[ids] = val1_NR[ids]
            obcs_all[t, k, :, :] = val1

    return obcs_all


def _squeeze_bnd(obcs_all, bnd):
    """Same squeeze gen_obcs_files applies to gen_obcs()'s per-timestep
    output -- explicit east=-1/west=0/north=-1/south=0 indexing, matching
    upstream commit f6d55709 (was `east/west: [...,0]`, `north/south:
    [...,0,:]` for all four, a latent bug that happened to be a no-op
    under this file's own pre-sliced-to-width-1 XC1/mask1 convention, but
    matched here anyway to track upstream exactly and remove any doubt)."""
    if bnd == 'east':
        return obcs_all[:, :, :, -1]
    elif bnd == 'west':
        return obcs_all[:, :, :, 0]
    elif bnd == 'north':
        return obcs_all[:, :, -1, :]
    else:  # south
        return obcs_all[:, :, 0, :]


# Bug D: gen_obcs.py's hardcoded seaice vnms=['AREA','HEFF','HSNOW','UICE','VICE']
# don't match the real diagnostic field names this project's STEP2 config
# actually extracted (confirmed against data.diagnostics_vec and the real
# dv/ filenames): 'SIarea','SIheff','SIhsnow','SIuice','SIvice' (MITgcm's own
# SEAICE-package names, 8-char-padded). Passing the un-prefixed names into
# go.read_dv_diags makes its glob match zero files (itrs=[]), so the read
# loop body never runs and `dv_diag` is referenced undefined -> NameError.
# Simply substituting the correct SI-prefixed names isn't enough either:
# go.read_dv_diags's rotation-branch detection (`'ICE' in vnm`, case-sensitive)
# and its Nr/mask-shape membership checks are also keyed off the wrong
# (unprefixed) names, and its U/V filename reconstruction
# (`sfx[:-4]+"U"+vnm[-3:]`) only happens to work for 'UVEL'/'VVEL' because
# those are exactly 4 chars with the direction letter first -- it silently
# rebuilds the WRONG filename for 6-char 'SIuice'/'SIvice'. Rather than
# patch go.read_dv_diags's string-slicing to handle a second naming
# convention, these two functions replace it for seaice fields only,
# writing out the exact real filenames and doing the same
# AngleCS/AngleSN rotation math UVEL/VVEL already use (validated correct).
SEAICE_SCALAR = ('SIarea', 'SIheff', 'SIhsnow')
SEAICE_VECTOR = ('SIuice', 'SIvice')

# Map from this script's own output-side field names to the field names the
# parent run's diagnostics_vec package actually writes into dv/. Fixed
# 2026-08-10: data.diagnostics_vec used to request 'SIarea'/'SIuice'/.../
# 'TRAC01'..31, but diagnostics_vec_prepare_subfield.F's SET_SUBFIELDS only
# recognizes 'AREA'/'UICE'/'VICE'/'HEFF'/'HSNOW' and 'PTRACE01'..31 -- those
# requested names never matched anything in that IF/ELSEIF chain, so the
# fields were silently written as all-zero for the whole 3-year run.
# data.diagnostics_vec now requests the recognized names; this map lets the
# OBCS output filenames (SIarea_*.bin, TRAC01_*.bin, ...) stay unchanged
# while reading the newly-correct dv/ files.
_DV_READ_NAME = {'SIarea': 'AREA', 'SIheff': 'HEFF', 'SIhsnow': 'HSNOW',
                  'SIuice': 'UICE', 'SIvice': 'VICE'}


def _dv_read_name(vnm):
    if vnm.startswith('TRAC') and vnm[4:].isdigit():
        return 'PTRACE' + vnm[4:]
    return _DV_READ_NAME.get(vnm, vnm)


def _read_dv_diags_seaice(vnm, bnd, bnd_domain):
    dv_diags_dir = os.path.join(go.config_dir, 'parent/outputs/OBCS/')
    nm_pt = len(bnd_domain[bnd]['XC'])

    def _read_field(fld):
        sfx = dv_diags_dir + f"{bnd}_BC_mask_{fld}"
        files = sorted(glob.glob(sfx + ".*"))
        if not files:
            raise FileNotFoundError(f"no dv output found for {sfx}.*")
        chunks = []
        for f in files:
            itr = int(f[len(sfx) + 1:-4])
            d = np.fromfile(sfx + f".{itr:010d}.bin", '>f4')
            tstp = int(np.size(d) / nm_pt)
            chunks.append(d.reshape(tstp, 1, nm_pt))
        return np.concatenate(chunks, axis=0)

    if vnm in SEAICE_VECTOR:
        dv_diagU = _read_field(_dv_read_name('SIuice'))
        dv_diagV = _read_field(_dv_read_name('SIvice'))
        dv_diag = np.zeros_like(dv_diagU)
        for m in range(nm_pt):
            if vnm == 'SIuice':
                dv_diag[:, :, m] = (bnd_domain[bnd]['AngleCS'][m] * dv_diagU[:, :, m] -
                                    bnd_domain[bnd]['AngleSN'][m] * dv_diagV[:, :, m])
            else:
                dv_diag[:, :, m] = (bnd_domain[bnd]['AngleSN'][m] * dv_diagU[:, :, m] +
                                    bnd_domain[bnd]['AngleCS'][m] * dv_diagV[:, :, m])
    else:
        dv_diag = _read_field(_dv_read_name(vnm))

    msk_pts = bnd_domain[bnd]['maskC'][:1, :]  # surface-only, no Zinterp needed
    return dv_diag, msk_pts


def _get_regbnd_info_seaice(vnm, bnd, tstp, grid1):
    if vnm == 'SIuice':
        mask = grid1['maskW']
    elif vnm == 'SIvice':
        mask = grid1['maskS']
    else:
        mask = grid1['maskC']

    if bnd == 'west':
        XC1 = grid1['XC'][:, :1]; YC1 = grid1['YC'][:, :1]; mask1 = mask[:, :, :1]
    elif bnd == 'east':
        XC1 = grid1['XC'][:, -1:]; YC1 = grid1['YC'][:, -1:]; mask1 = mask[:, :, -1:]
    elif bnd == 'north':
        XC1 = grid1['XC'][-1:, :]; YC1 = grid1['YC'][-1:, :]; mask1 = mask[:, -1:, :]
    elif bnd == 'south':
        XC1 = grid1['XC'][:1, :]; YC1 = grid1['YC'][:1, :]; mask1 = mask[:, :1, :]
    else:
        raise ValueError('Boundary ' + bnd + ' not recognized')

    mask1 = mask1[:1, :, :]
    obcs_init = np.zeros((tstp, 1, np.size(XC1)))
    return XC1, YC1, mask1, obcs_init


# ---------------------------------------------------------------------------
# Disc (non-rectangular) open boundary -- pan_AO only
# ---------------------------------------------------------------------------
# Every other region's open boundary IS the edge of the model array, and both
# go.get_regbnd_info (gen_obcs.py:68-83) and _get_regbnd_info_seaice above
# hard-code exactly that: XC1 = grid1['XC'][:, :1] for west, and so on.
#
# pan_AO's domain is a disc inscribed in a 1680x1680 square, so the square's
# perimeter is not the boundary. Measured on the model-facing bathymetry:
#
#     square array perimeter : 6,716 cells, only    174 wet
#     disc rim ring          : 4,748 cells,       1,700 wet
#
# Using the edge slices would build boundary forcing along a line that is
# 97% land in this model, and -- the dangerous part -- it would NOT fail.
# Every array still has the right shape, so the run would start and quietly
# integrate with no Atlantic inflow at the Greenland-Scotland ridge, which is
# the sector that sets the Arctic heat budget.
#
# MITgcm supports this natively. pkg/obcs stores the boundary as per-row /
# per-column INDEX ARRAYS OB_Iwest(j), OB_Ieast(j), OB_Jsouth(i), OB_Jnorth(i)
# with OB_indexUnset = 0 meaning "no open boundary in this row/column"
# (obcs_readparms.F:237-238). Crucially the FORCING FILE LAYOUT IS UNCHANGED:
# E/W files are still (Nr, Ny) indexed by j, N/S still (Nr, Nx) indexed by i.
# So only the GATHER changes -- row j of a west file must come from model cell
# (j, OB_Iwest(j)-1) instead of (j, 0). Rows whose index is 0 are written as
# zeros and never read by pkg/obcs.
#
# WHY THE GATHER IS JUST mask[:, j, i], with no face offset. gen_obcs_files_fast
# trims HFacS to [:, 1:, :] and HFacW to [:, :, :-1] (see the comment there).
# After those trims, index m of HFacS is the NORTH face of cell m and index m
# of HFacW is the WEST face of cell m. Check the trimmed index the box path
# actually uses against the box boundary CELL:
#
#     west  cell i=0     -> mask[:, :, :1]  index 0     == cell i
#     east  cell i=nx-1  -> mask[:, :, -1:] index nx-1  == cell i
#     south cell j=0     -> mask[:, :1, :]  index 0     == cell j
#     north cell j=ny-1  -> mask[:, -1:, :] index ny-1  == cell j
#
# In all four the trimmed mask index EQUALS the boundary cell index, so
# gathering at the OB cell reproduces the box convention exactly on a box and
# extends it to a staircase without inventing a new one. This matters: the
# box convention is the validated one, and the disc must not quietly differ.
#
# The indices come from grid/disc_ob_index.npz, written by
# grid/gen_disc_ob_index.py -- the SAME file data.obcs and gen_dvmasks_disc.py
# are built from. One definition of the boundary, three consumers; that is
# what keeps the dv masks, the forcing files and the namelist in agreement.

_DISC_AXIS = {'west': ('OB_Iwest', 'row'), 'east': ('OB_Ieast', 'row'),
              'south': ('OB_Jsouth', 'col'), 'north': ('OB_Jnorth', 'col')}

# Which C-grid FACE each velocity component is masked on, per boundary.
#
# Transcribed directly from darwin3_ckpt69k/pkg/obcs/obcs_apply_uv.F:145-213,
# which is the only authority on this -- the model multiplies the prescribed
# value by exactly these masks, so anything written on a different face is
# either silently discarded or silently missing:
#
#     West:  uFld(Iobc+1,j,k) = OBWu(j,k) * _maskW(Iobc+1,j,k)
#     East:  uFld(Iobc,  j,k) = OBEu(j,k) * _maskW(Iobc,  j,k)
#     North: vFld(i,Jobc,  k) = OBNv(i,k) * _maskS(i,Jobc,  k)
#     South: vFld(i,Jobc+1,k) = OBSv(i,k) * _maskS(i,Jobc+1,k)
#
# The normal component is always masked on the INTERIOR face -- but
# "interior" is +1 at west/south and +0 at east/north. The offset is
# therefore DIRECTION-DEPENDENT, and that is why this cannot be expressed as
# the single axis-wide slice it used to be
# (`HFacW[:, :, :-1]` / `HFacS[:, 1:, :]`, applied once in
# gen_obcs_files_fast). One global trim is necessarily right at one end of an
# axis and wrong at the other: it gave west-normal 0 wet faces of 737 and
# north-normal 58 of 522, while east and south came out correct. See
# DECISIONS.md 2026-10-01 (15:50).
#
# This was a no-op for every BOX domain (fiji, kerguelen) -- there the cell
# outside the boundary is outside the ARRAY, so the outer face is wet and
# both choices agree (verified: fiji west outer 1677 == interior 1677). The
# disc is the first geometry where the cell outside the boundary is LAND
# inside the array, which is what made a long-latent bug reachable.
#
# Offsets apply to x for 'W' masks and to y for 'S' masks. Only west-U and
# south-V are offset; the other six are read at the boundary cell itself.
_MASK_FACE = {
    ('west',  'U'): ('W', 1), ('west',  'V'): ('S', 0),
    ('east',  'U'): ('W', 0), ('east',  'V'): ('S', 0),
    ('north', 'U'): ('W', 0), ('north', 'V'): ('S', 0),
    ('south', 'U'): ('W', 0), ('south', 'V'): ('S', 1),
}

# Velocity variables, by the component they carry. Everything else is a
# scalar on maskC, which sits at the cell centre and has no face ambiguity --
# which is why the other 160 OBCS files were never affected.
_UCOMP = ('UVEL', 'UICE', 'SIuice')
_VCOMP = ('VVEL', 'VICE', 'SIvice')


def _disc_mask_for(vnm, bnd, grid1):
    """(mask array, x-offset, y-offset) for this variable on this boundary.

    Reads the UNTRIMMED HFac-derived masks and applies the offset at gather
    time. Pre-trimming cannot work here: the same array is needed with two
    different offsets depending on which boundary is being written.
    """
    if vnm in _UCOMP:
        comp = 'U'
    elif vnm in _VCOMP:
        comp = 'V'
    else:
        return grid1['maskC'], 0, 0
    which, off = _MASK_FACE[(bnd, comp)]
    if which == 'W':
        return grid1['maskW_full'], off, 0
    return grid1['maskS_full'], 0, off

# Mirrors go.get_regbnd_info:88-99 exactly. Kept as named tuples rather than
# re-written as a condition so a future change upstream is easy to diff.
_SURF_ONLY_SEAICE = ('AREA', 'UICE', 'VICE', 'HEFF', 'HSNOW', 'ETAN')
_SURF_ONLY_PLAIN = ('AREA', 'ETAN')


def load_disc_ob(config_dir, mode='auto'):
    """Load the disc OB index arrays, or None for an ordinary box region.

    mode 'on' refuses to fall back. That asymmetry is deliberate: a missing
    npz under 'auto' is the normal case for fiji/kerguelen, but for pan_AO it
    would silently restore the square-perimeter boundary described above, so
    its job script passes 'on' and gets a hard failure instead.
    """
    if mode == 'off':
        return None
    npz = os.path.join(config_dir, 'grid', 'disc_ob_index.npz')
    if not os.path.exists(npz):
        if mode == 'on':
            raise SystemExit(
                f"FAIL: -disc on, but {npz} is missing. Run "
                f"grid/gen_disc_ob_index.py first. Falling back to the box "
                f"edge slices here would build the open boundary along the "
                f"square perimeter without any error.")
        return None
    z = np.load(npz)
    return {k: z[k].astype(int) for k in
            ('OB_Iwest', 'OB_Ieast', 'OB_Jsouth', 'OB_Jnorth')}


def _disc_fill_coords(c, has):
    """Nearest-valid fill for the rows/columns that have no open boundary.

    Those entries are masked to zero and their results are discarded, but
    they are still handed to griddata as TARGET coordinates, so they must be
    finite. NaN would survive the linear pass and then crash the cKDTree in
    the nearest-neighbour fallback; a wild value (0,0) would stretch the
    target cloud across the globe. Copying the nearest real boundary cell's
    coordinate keeps the cloud compact and is provably inert, because
    mask1 == 0 at exactly these entries and _gen_obcs_batched zeroes them
    (`val1[mask1[k] == 0] = 0`) before Hextrapol, which also never reads a
    mask1 == 0 cell as a source.
    """
    n = has.size
    out = np.asarray(c, dtype=float).copy()
    if not has.any():
        return np.zeros(n)
    k = np.arange(n)
    prev = np.maximum.accumulate(np.where(has, k, -1))
    nxt = np.minimum.accumulate(np.where(has, k, n)[::-1])[::-1]
    for m in np.flatnonzero(~has):
        a, b = prev[m], nxt[m]
        if a < 0:
            out[m] = out[b]
        elif b >= n:
            out[m] = out[a]
        else:
            out[m] = out[a] if (m - a) <= (b - m) else out[b]
    return out


def _get_regbnd_info_disc(mask, bnd, grid1, ob, surface_only, tstp,
                          off_x=0, off_y=0):
    """Box-equivalent (XC1, YC1, mask1, obcs_init) gathered along the OB
    index arrays instead of sliced off an array edge. Shapes are identical to
    the box path, so _gen_obcs_batched and _squeeze_bnd are untouched.

    off_x/off_y shift the gather onto the C-grid FACE the model masks this
    component on -- see _MASK_FACE. They are 0 for every scalar; only
    west-U (off_x=1) and south-V (off_y=1) are shifted.
    """
    key, axis = _DISC_AXIS[bnd]
    idx = ob[key]
    XC, YC = grid1['XC'], grid1['YC']
    ny, nx = XC.shape
    n = ny if axis == 'row' else nx
    if idx.size != n:
        raise SystemExit(
            f"FAIL: {key} has {idx.size} entries but the regional grid is "
            f"{ny}x{nx} -- disc_ob_index.npz and the stitched ncgrid "
            f"describe different domains")
    # 2026-10-01: negatives are rejected BEFORE `has = idx > 0`. Previously
    # the `idx[has].min() < 1` test could never fire -- `has` had already
    # excluded every negative -- so a malformed array would have been
    # gathered as "no open boundary in this row" with no complaint. Found by
    # test_disc_transport.py while porting the same guard into the transport
    # correction. pkg/obcs uses OB_indexUnset = 0 and nothing else.
    if (idx < 0).any():
        raise SystemExit(
            f"FAIL: {key} holds {int((idx < 0).sum())} negative entries "
            f"(min {idx.min()}); 0 is the only valid 'no boundary' value")
    has = idx > 0
    if has.any() and idx[has].max() > (nx if axis == 'row' else ny):
        raise SystemExit(
            f"FAIL: {key} holds an index outside 1..{nx if axis == 'row' else ny}")

    k = np.arange(n)
    pos = np.where(has, idx - 1, 0)          # placeholder index where unset
    rows, cols = (k, pos) if axis == 'row' else (pos, k)

    xs = _disc_fill_coords(XC[rows, cols], has)
    ys = _disc_fill_coords(YC[rows, cols], has)
    # Coordinates always come from the cell centre (XC/YC are (ny,nx)); only
    # the MASK is read on the shifted face. The interpolation target location
    # is the boundary cell -- what the offset selects is which face the model
    # will apply the result to, not where the value lives.
    mrows = rows + off_y
    mcols = cols + off_x
    if mrows.max() >= mask.shape[1] or mcols.max() >= mask.shape[2]:
        raise SystemExit(
            f"FAIL: {bnd} mask gather at offset (x={off_x}, y={off_y}) runs "
            f"past the mask array {mask.shape} -- the untrimmed HFac is "
            f"required here, not the trimmed one")
    m = np.where(has[None, :], mask[:, mrows, mcols], 0.0)   # (Nr, n)

    if axis == 'row':
        XC1, YC1, mask1 = xs[:, None], ys[:, None], m[:, :, None]
    else:
        XC1, YC1, mask1 = xs[None, :], ys[None, :], m[:, None, :]
    if surface_only:
        mask1 = mask1[:1, :, :]
    return XC1, YC1, mask1, np.zeros((tstp, mask1.shape[0], n))


def _regbnd_disc_plain(vnm, bnd, tstp, grid1, ob):
    """Disc counterpart of go.get_regbnd_info.

    NOTE: unlike the box path (gen_obcs.py:60-65) this does NOT pick between
    two globally pre-trimmed masks. The face a velocity component is masked
    on depends on the BOUNDARY as well as the component -- see _MASK_FACE.
    """
    mask, off_x, off_y = _disc_mask_for(vnm, bnd, grid1)
    surf = vnm in (_SURF_ONLY_SEAICE if go.seaice else _SURF_ONLY_PLAIN)
    return _get_regbnd_info_disc(mask, bnd, grid1, ob, surf, tstp,
                                 off_x=off_x, off_y=off_y)


def _regbnd_disc_seaice(vnm, bnd, tstp, grid1, ob):
    """Disc counterpart of _get_regbnd_info_seaice -- always surface-only."""
    mask, off_x, off_y = _disc_mask_for(vnm, bnd, grid1)
    return _get_regbnd_info_disc(mask, bnd, grid1, ob, True, tstp,
                                 off_x=off_x, off_y=off_y)


def report_disc_ob(ob, grid1, bnd_ls):
    """Loud, independent re-verification that the OB indices and the STITCHED
    ncgrid still agree, run once before any forcing is written.

    gen_disc_ob_index.py checked the indices against the model-facing .bin it
    built them from; this checks them against the ncgrid the interpolation
    actually reads. If those two ever diverge the forcing would be written
    onto land and nothing downstream would notice."""
    maskC = grid1['maskC']
    ny, nx = grid1['XC'].shape
    total, dry = 0, 0
    for bnd in bnd_ls:
        key, axis = _DISC_AXIS[bnd]
        idx = ob[key]
        has = idx > 0
        k = np.arange(idx.size)[has]
        pos = idx[has] - 1
        rows, cols = (k, pos) if axis == 'row' else (pos, k)
        nd = int((maskC[0, rows, cols] == 0).sum()) if k.size else 0
        total += int(has.sum())
        dry += nd
        print(f'    {bnd:5s}: {int(has.sum()):5d} open-boundary entries'
              + (f'   *** {nd} ON LAND ***' if nd else ''))
    print(f'  {total} open-boundary entries total on a {ny}x{nx} disc')
    if dry:
        raise SystemExit(
            f"FAIL: {dry} open-boundary cells are DRY in the stitched ncgrid "
            f"-- disc_ob_index.npz and the ncgrid disagree, so the forcing "
            f"would be written onto land.")
    if total == 0:
        raise SystemExit("FAIL: every OB index array is empty")


_WORKER_STATE = {}


def _init_worker(grid0, grid1, bnd_domain, out_dir, config_dir, seaice, hmode,
                 disc_ob):
    _WORKER_STATE['grid0'] = grid0
    _WORKER_STATE['grid1'] = grid1
    _WORKER_STATE['bnd_domain'] = bnd_domain
    _WORKER_STATE['out_dir'] = out_dir
    _WORKER_STATE['hmode'] = hmode
    _WORKER_STATE['disc_ob'] = disc_ob
    go.config_dir = config_dir
    go.seaice = seaice


def _run_one(task):
    vnm, bnd = task
    grid0 = _WORKER_STATE['grid0']
    grid1 = _WORKER_STATE['grid1']
    bnd_domain = _WORKER_STATE['bnd_domain']
    out_dir = _WORKER_STATE['out_dir']
    Nr0 = grid0['HFacC'].shape[0]
    Nr1 = grid1['HFacC'].shape[0]

    # The PARENT side needs no change for a disc: bnd_domain comes from the
    # *_BC_mask.bin files, which gen_dvmasks_disc.py already built along the
    # disc rim. Only the REGIONAL target geometry differs -- see the disc
    # section above.
    ob = _WORKER_STATE['disc_ob']

    if vnm in SEAICE_SCALAR or vnm in SEAICE_VECTOR:
        dv_diag, msk_pts = _read_dv_diags_seaice(vnm, bnd, bnd_domain)
        tstp = len(dv_diag)
        if ob is None:
            XC1, YC1, mask1, _ = _get_regbnd_info_seaice(vnm, bnd, tstp, grid1)
        else:
            XC1, YC1, mask1, _ = _regbnd_disc_seaice(vnm, bnd, tstp, grid1, ob)
    else:
        dv_diag, msk_pts = go.read_dv_diags([0], _dv_read_name(vnm), bnd, bnd_domain, Nr0, Nr1, grid0, grid1)
        tstp = len(dv_diag)
        if ob is None:
            XC1, YC1, mask1, _ = go.get_regbnd_info(vnm, bnd, Nr1, tstp, grid1)
        else:
            XC1, YC1, mask1, _ = _regbnd_disc_plain(vnm, bnd, tstp, grid1, ob)

    obcs_all = _gen_obcs_batched(dv_diag, bnd_domain[bnd]['XC'], bnd_domain[bnd]['YC'],
                                  msk_pts, XC1, YC1, mask1,
                                  _WORKER_STATE['hmode'])
    obcs = _squeeze_bnd(obcs_all, bnd)
    output_file = os.path.join(out_dir, f'{vnm}_{bnd}.bin')
    obcs.ravel('C').astype('>f4').tofile(output_file)
    return (vnm, bnd, tstp, output_file)


def gen_obcs_files_fast(config_dir, reg_nm, boundaries, seaice, bgc, nproc, print_level,
                        hmode=None, disc='auto'):
    grd_ls = ['XC', 'YC', 'AngleCS', 'AngleSN', 'HFacC', 'HFacS', 'HFacW', 'drF']

    if print_level >= 1:
        print('> Reading in the parent model tile geometry')
    tmp = go.read_eccogrid(config_dir, reg_nm)
    grid0 = dict(zip(grd_ls + ['rF'], tmp))
    Nr0 = tmp[grd_ls.index("HFacC")].shape[0]
    llc = tmp[grd_ls.index("XC")].shape[1]
    del tmp

    if print_level >= 1:
        print('> Reading in the regional model tile geometry')
    tmp = go.read_ncgrid(config_dir, reg_nm, grd_ls)

    # Keep the UNTRIMMED face arrays before anything is sliced. HFacW is
    # (Nr, ny, nx+1) and HFacS is (Nr, ny+1, nx) -- they carry BOTH faces of
    # every cell, and the disc path needs that, because which face a velocity
    # component is masked on depends on the boundary (see _MASK_FACE). Any
    # single global trim throws away exactly the half that two of the four
    # boundaries need.
    # bool, not float: these are 0/1 and are only ever indexed down to a
    # (Nr, n) sliver before use, so the full array never needs to be float.
    # As float64 the pair would add ~4 GB to every one of the 14 forked
    # workers; as bool it is ~254 MB each.
    maskW_full = tmp[grd_ls.index("HFacW")] > 0
    maskS_full = tmp[grd_ls.index("HFacS")] > 0

    # The trims below remain for the BOX path only (go.get_regbnd_info).
    #
    # HFacS has one extra row (ny+1) on the C-grid -- the open north/south
    # boundary lives on the INTERIOR face (row 1..ny), not the outer domain
    # edge (row 0..ny-1). Trimming [:, :-1, :] (dropping the last row)
    # instead of [:, 1:, :] (dropping the first) reads the south boundary
    # off a degenerate non-existent outer-wall row and the north boundary
    # one row short of the true edge -- upstream commit f6d55709 ("Fix
    # OBCs downscaling method due to mask issue", 2026-08-03) fixed this
    # exact bug in gen_obcs.py; ported here since this file duplicated the
    # buggy trim (copied before that fix existed upstream).
    #
    # 2026-10-01: that fix states the right principle ("INTERIOR face") but
    # encodes it as one axis-wide slice, which can only be correct at ONE end
    # of each axis -- "interior" is +1 at west/south and +0 at east/north.
    # Harmless on a box domain (outer face is wet, both choices agree);
    # it zeroed west-normal and 89% of north-normal on the disc. The disc
    # path no longer uses these arrays. See DECISIONS.md 2026-10-01 (15:50).
    tmp[grd_ls.index("HFacS")] = tmp[grd_ls.index("HFacS")][:, 1:, :]
    tmp[grd_ls.index("HFacW")] = tmp[grd_ls.index("HFacW")][:, :, :-1]
    grid1 = dict(zip(grd_ls, tmp))
    del tmp

    bnd_ls, dv_masks_ls = go.read_dv_masks(config_dir, boundaries, llc)
    for nm in ['C', 'S', 'W']:
        tmp = np.copy(grid1["HFac" + nm])
        tmp[tmp > 0] = 1
        grid1.update({"mask" + nm: tmp})
        del tmp
    grid1['maskW_full'] = maskW_full
    grid1['maskS_full'] = maskS_full

    # BOX path, west-U and north-V (2026-10-09, found on WestAO; Nares has it too). The trims above leave the
    # west-U / north-V gathers on the OUTER box faces (0 / ny), which in an ncgrid from gen_ncgrid wrap periodically
    # to the opposite edge, so their mask is the INTERSECTION of both edges' wetness. Nares: west 10678 vs 15508
    # wet (31% of Baffin-edge U zeroed), north 2260 vs 10587 (~80% of north-edge V, ocean and ice, zeroed).
    # pkg/obcs prescribes the normal velocity on the INTERIOR face (west Iobc+1 -> face 1, north Jobc -> face ny-1).
    grid1['maskW'][:, :, 0] = maskW_full[:, :, 1]
    grid1['maskS'][:, -1, :] = maskS_full[:, -2, :]

    out_dir = os.path.join(config_dir, 'forcings/OBCS/')
    os.makedirs(out_dir, exist_ok=True)

    # Choose the interpolation plane from the regional domain. pan_AO is
    # circumpolar and reaches 89.978N, so the lon/lat plane has the +/-180
    # seam running through its interior and is degenerate at the pole; the
    # 2026-09-23 pickups generated that way came out 12.1% NaN. hproj stays
    # on 'lonlat' for every non-polar region, so this changes nothing for
    # fiji or kerguelen. See hproj.py.
    if hmode is None:
        hmode = hproj.choose_mode(grid1['XC'], grid1['YC'])
    if print_level >= 1:
        print('> Boundary interpolation plane: %s' % hmode)

    disc_ob = load_disc_ob(config_dir, disc)
    if disc_ob is None:
        if print_level >= 1:
            print('> Open boundary: the four array edges (box region)')
    else:
        print('> Open boundary: DISC -- gathering along OB_Iwest/OB_Ieast/'
              'OB_Jsouth/OB_Jnorth')
        report_disc_ob(disc_ob, grid1, bnd_ls)

    if print_level >= 1:
        print('> Getting boundary domain')
    bnd_domain = go.gen_bnd_domain(bnd_ls, dv_masks_ls, grid0, Nr0)

    tasks = []
    for vnm in ['THETA', 'SALT', 'UVEL', 'VVEL', 'ETAN']:
        for bnd in bnd_ls:
            tasks.append((vnm, bnd))
    if seaice:
        # Real diagnostic/output names (SIarea/SIheff/SIhsnow/SIuice/SIvice),
        # NOT gen_obcs.py's hardcoded AREA/HEFF/HSNOW/UICE/VICE -- see Bug D
        # note above _read_dv_diags_seaice. Output filenames follow the same
        # names; regional_run/input*/data.obcs's OBWaFile/OBWhFile/OBWsnFile/
        # OBWuiceFile/OBWviceFile (currently unset -- separate STEP4 gap) must
        # be pointed at these exact filenames.
        for vnm in SEAICE_SCALAR + SEAICE_VECTOR:
            for bnd in bnd_ls:
                tasks.append((vnm, bnd))
    if bgc:
        sfx = os.path.join(config_dir, 'parent/outputs/OBCS/') + f"{bnd_ls[0]}_BC_mask_PTRACE*"
        ptr_files = glob.glob(sfx)
        ptr_ls = np.unique([int(f[len(sfx) - 1:-15]) for f in ptr_files])
        for i in ptr_ls:
            vnm = f'TRAC{i:02d}'
            for bnd in bnd_ls:
                tasks.append((vnm, bnd))

    if print_level >= 1:
        print(f'> {len(tasks)} (variable, boundary) tasks, {nproc} parallel workers')

    ctx = mp.get_context('fork')
    with ctx.Pool(processes=nproc, initializer=_init_worker,
                  initargs=(grid0, grid1, bnd_domain, out_dir, config_dir, seaice,
                            hmode, disc_ob)) as pool:
        for vnm, bnd, tstp, output_file in pool.imap_unordered(_run_one, tasks):
            if print_level >= 1:
                print(f'    - done: {vnm} {bnd} (tstp={tstp}) -> {output_file}')

    print("OBCS files successfully generated!")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--config_dir", required=True)
    parser.add_argument("-n", "--reg_nm", required=True)
    parser.add_argument("-bnd", "--bnd_nm", required=True, help="e.g. EWNS")
    parser.add_argument("-seaice", "--si", action="store_true")
    parser.add_argument("-bgc", "--darwin", action="store_true")
    parser.add_argument("-p", "--nproc", type=int, default=os.cpu_count())
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-hp", "--hproj", dest="hproj", type=str, default=None,
                        choices=list(hproj.MODES),
                        help="plane to interpolate in; default is chosen from "
                             "the domain. See hproj.py")
    parser.add_argument("-disc", "--disc", default="auto",
                        choices=("auto", "on", "off"),
                        help="gather the open boundary along "
                             "grid/disc_ob_index.npz instead of the four array "
                             "edges. 'auto' uses it if present; 'on' FAILS if "
                             "it is missing (use this for pan_AO -- a silent "
                             "fallback would put the boundary on the square "
                             "perimeter, which is 97%% land)")
    args = parser.parse_args()

    print_level = 1 if args.verbose else 0
    gen_obcs_files_fast(args.config_dir, args.reg_nm, args.bnd_nm, args.si, args.darwin,
                         args.nproc, print_level, hmode=args.hproj, disc=args.disc)
