"""
plot_boundary_velocity_profile.py

Plots the TIME-MEAN, depth-resolved boundary-normal velocity at each open
boundary, side by side for the parent (native LLC270) and the downscaled
child model -- a depth x along-boundary "section" plot, as opposed to
check_obcs_transport.py's depth-INTEGRATED transport-per-unit-length
transect. This shows the vertical (baroclinic) structure directly, so it's
possible to see exactly where (both along the boundary AND at depth) the
child's velocity departs from the parent's, rather than only in the
depth-integrated total.

Reuses check_obcs_transport.py's native/child velocity reading (same
rotation convention, same multi-file discovery/concatenation) and
gen_obcs.py's grid-reading helpers, so results stay consistent with -- and
safe to run alongside -- check_obcs_transport.py. Read-only: never writes
any OBCS or grid file, only PNG plots.

Usage:
    python3 plot_boundary_velocity_profile.py -d /path/to/grid -n GoM_1km \\
                                               -bnd ENS -i 26281 -o /path/to/plots -v
"""
import os
import argparse
import warnings
import numpy as np

from gen_obcs import transp_tiles, read_eccogrid, read_ncgrid, read_dv_masks, gen_bnd_domain
from check_obcs_transport import NORMAL_VNM, read_native_normal_velocity, read_child_velocity

############################################################
#                    DEPTH HELPERS                          #
############################################################

def cell_center_depths(drF):
    """Cell-center depth (negative downward, m) from cell thickness drF (Nr,)."""
    drF = np.asarray(drF)
    return -(np.cumsum(drF) - drF / 2)


def masked_time_mean(vel, mask):
    """
    vel: (tstp, Nr, n) velocity; mask: (Nr, n) wet fraction (0/1 or hFac).
    Returns (Nr, n) time-mean velocity, NaN wherever never wet, so dry
    cells plot as blank rather than a misleading 0 m/s.
    """
    wet = mask > 0
    vel_masked = np.where(wet[None, :, :], vel, np.nan)
    with warnings.catch_warnings(), np.errstate(invalid='ignore'):
        warnings.simplefilter('ignore', category=RuntimeWarning)   # all-dry columns -> empty-slice mean
        mean = np.nanmean(vel_masked, axis=0)
    mean = np.where(wet, mean, np.nan)
    return mean

############################################################
#                        PLOTTING                            #
############################################################

def shared_xlim(coord0, mean0, coord1, mean1, pad_frac=0.1, min_pad=0.25):
    """
    pcolormesh auto-scales each subplot to its OWN x-data -- since coord1
    (child) always spans the full domain width while coord0 (native) only
    spans wherever native samples exist, the two panels end up on wildly
    different x-axes, making it impossible to visually compare them and
    making any genuinely narrow "before" wet strip on the child panel
    nearly invisible against its own much wider axis. This computes one
    x-range -- the union of wherever EITHER model actually has finite
    (wet) data -- with a little padding, to use on both panels instead.
    Falls back to each array's own full range if nothing is finite.
    """
    wet0 = coord0[np.any(np.isfinite(mean0), axis=0)] if mean0.size else np.array([])
    wet1 = coord1[np.any(np.isfinite(mean1), axis=0)] if mean1.size else np.array([])
    combined = np.concatenate([wet0, wet1])
    if combined.size == 0:
        combined = np.concatenate([coord0, coord1])
    lo, hi = float(np.min(combined)), float(np.max(combined))
    pad = max((hi - lo) * pad_frac, min_pad)
    return lo - pad, hi + pad


def plot_velocity_profile(bnd, vnm, reg_nm, coord0, Z0, mean0, coord1, Z1, mean1, plot_dir):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    finite_vals = np.concatenate([mean0[np.isfinite(mean0)].ravel(), mean1[np.isfinite(mean1)].ravel()])
    vmax = np.max(np.abs(finite_vals)) if finite_vals.size > 0 else 1.0
    vmax = vmax if vmax > 0 else 1.0

    order0 = np.argsort(coord0)
    order1 = np.argsort(coord1)
    xlim = shared_xlim(coord0, mean0, coord1, mean1)

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(14, 5), sharey=True)

    pc0 = ax0.pcolormesh(coord0[order0], Z0, mean0[:, order0], shading='nearest',
                          cmap='RdBu_r', vmin=-vmax, vmax=vmax)
    ax0.set_title(f'parent (native)\n{reg_nm} {bnd} {vnm}')
    ax0.set_xlabel('Latitude' if bnd in ('east', 'west') else 'Longitude')
    ax0.set_ylabel('Depth (m)')
    ax0.set_xlim(xlim)

    ax1.pcolormesh(coord1[order1], Z1, mean1[:, order1], shading='nearest',
                    cmap='RdBu_r', vmin=-vmax, vmax=vmax)
    ax1.set_title(f'downscaled (child)\n{reg_nm} {bnd} {vnm}')
    ax1.set_xlabel('Latitude' if bnd in ('east', 'west') else 'Longitude')
    ax1.set_xlim(xlim)

    fig.colorbar(pc0, ax=[ax0, ax1], label='time-mean velocity (m/s)', shrink=0.85)

    os.makedirs(plot_dir, exist_ok=True)
    fpath = os.path.join(plot_dir, f'{reg_nm}_{bnd}_{vnm}_velocity_profile.png')
    fig.savefig(fpath, dpi=150, bbox_inches='tight')
    plt.close(fig)
    return fpath

############################################################
#                          MAIN                              #
############################################################

def main(config_dir, reg_nm, boundaries, itrs, output_dir, verbose, legacy_bug_mode=False):
    grd_ls = ['XC', 'YC', 'AngleCS', 'AngleSN', 'HFacC', 'HFacS', 'HFacW', 'drF']
    obcs_dir = os.path.join(config_dir, 'forcings/OBCS/')

    if verbose:
        print('> Reading parent and downscaled grid geometry')
    tmp = read_eccogrid(config_dir, reg_nm)
    grid0 = dict(zip(grd_ls + ['rF'], tmp))
    Nr0 = tmp[grd_ls.index("HFacC")].shape[0]
    llc = tmp[grd_ls.index("XC")].shape[1]

    tmp = read_ncgrid(config_dir, reg_nm, grd_ls)
    tmp[grd_ls.index("HFacS")] = tmp[grd_ls.index("HFacS")][:, 1:, :]
    tmp[grd_ls.index("HFacW")] = tmp[grd_ls.index("HFacW")][:, :, :-1]
    grid1 = dict(zip(grd_ls, tmp))
    Nr1 = tmp[grd_ls.index("HFacC")].shape[0]
    for nm in ['S', 'W']:
        m = np.copy(grid1["HFac" + nm]); m[m > 0] = 1
        grid1["mask" + nm] = m

    bnd_ls, dv_masks_ls = read_dv_masks(config_dir, boundaries, llc)
    bnd_domain = gen_bnd_domain(bnd_ls, dv_masks_ls, grid0, Nr0)

    Z0 = cell_center_depths(grid0['drF'])
    Z1 = cell_center_depths(grid1['drF'])

    for bnd in bnd_ls:
        vnm = NORMAL_VNM[bnd]
        if bnd in ('east', 'west'):
            mask1 = grid1['maskW'][:, :, -1:] if bnd == 'east' else grid1['maskW'][:, :, :1]
            coord1 = (grid1['YC'][:, -1:] if bnd == 'east' else grid1['YC'][:, :1]).ravel()
            mask0 = bnd_domain[bnd]['maskW']
            coord0 = bnd_domain[bnd]['YC']
        else:
            mask1 = grid1['maskS'][:, -1:, :] if bnd == 'north' else grid1['maskS'][:, :1, :]
            coord1 = (grid1['XC'][-1:, :] if bnd == 'north' else grid1['XC'][:1, :]).ravel()
            mask0 = bnd_domain[bnd]['maskS']
            coord0 = bnd_domain[bnd]['XC']
        mask1_2d = mask1[:, :, 0] if bnd in ('east', 'west') else mask1[:, 0, :]
        n1 = mask1_2d.shape[1]

        if verbose:
            print(f'> [{bnd}] reading native parent velocity and downscaled {vnm}')
        native_vel = read_native_normal_velocity(config_dir, bnd, itrs, bnd_domain, Nr0, verbose=verbose,
                                                  legacy_bug_mode=legacy_bug_mode)
        child_vel = read_child_velocity(obcs_dir, bnd, vnm, Nr1, n1)
        if native_vel.shape[0] != child_vel.shape[0]:
            raise SystemExit(
                f"[{bnd}] timestep mismatch: native has {native_vel.shape[0]} timestep(s), child "
                f"has {child_vel.shape[0]}. Pass the exact iteration(s) with -i to bypass auto-"
                f"discovery, same as with check_obcs_transport.py.")

        mean0 = masked_time_mean(native_vel, mask0)
        mean1 = masked_time_mean(child_vel, mask1_2d)

        fpath = plot_velocity_profile(bnd, vnm, reg_nm, coord0, Z0, mean0, coord1, Z1, mean1,
                                       os.path.join(output_dir, 'diagnostics'))
        print(f'    [{bnd}] saved velocity profile plot: {fpath}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-d", "--config_dir", required=True, type=str,
                        help="Directory where gen_obcs.py's inputs/outputs live (same -d you gave it)")
    parser.add_argument("-n", "--reg_nm", required=True, type=str, help="Name of the regional cutout")
    parser.add_argument("-bnd", "--bnd_nm", required=True, type=str,
                        help="Open boundaries to plot, e.g. ENS")
    parser.add_argument("-i", "--itr", nargs='+', type=int, required=True,
                        help="Same iteration(s) used when running gen_obcs.py/gen_obcs2.py")
    parser.add_argument("-o", "--output_dir", type=str, default=None,
                        help="Where to save plots (default: <config_dir>; saved under "
                             "<output_dir>/diagnostics/)")
    parser.add_argument("-legacy-bug-mode", "--legacy_bug_mode", action="store_true",
                        help="Replicate the original gen_obcs.py's multi-file merge bug when "
                             "reading native files -- only needed if the OBCS files you're "
                             "comparing against were generated by the unfixed gen_obcs.py")
    parser.add_argument("-v", "--verbose", action="store_true")

    args = parser.parse_args()
    output_dir = args.output_dir or args.config_dir
    main(args.config_dir, args.reg_nm, args.bnd_nm, args.itr, output_dir, args.verbose,
         args.legacy_bug_mode)
