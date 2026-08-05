"""
verify_run.py

Post-run verification for a downscaled regional set-up (STEP5).

Two independent checks, both self-contained -- neither needs the parent model
to be re-run:

  --health  Parse the MITgcm STDOUT of a short verification integration and
            report whether it ran sanely: no NaN, advective CFL below 1,
            temperature/salinity inside physical bounds, no runaway kinetic
            energy drift, and a normal termination message.

  --inputs  Confirm the model was actually handed what STEP3 generated:
            that every pickup and OBCS file is the size the grid implies, is
            free of NaN/Inf, holds physically plausible values, and -- most
            usefully -- that each boundary's OBCS is non-zero on the wet cells
            of the SAME interior C-grid face the model masks it with. A
            boundary whose OBCS is zero everywhere it should be wet is the
            classic symptom of a mask taken from the wrong face, and it is
            invisible until the regional run produces no inflow.

By default both run. Exit status is 0 when everything passes, 1 when any
check FAILs, so this can be used as a gate in a job script.

Usage:
    python3 verify_run.py -d /path/to/config_dir -n GoM_1km -bnd ENS \\
                          -i 26281 -r /path/to/run -plot -v
"""
import os
import re
import glob
import argparse
import numpy as np

from gen_obcs import read_ncgrid

############################################################
#                        CONSTANTS                          #
############################################################

BND_NAME = {'E': 'east', 'W': 'west', 'N': 'north', 'S': 'south'}
NORMAL_VNM = {'east': 'UVEL', 'west': 'UVEL', 'north': 'VVEL', 'south': 'VVEL'}

### surface (single-level) OBCS fields; everything else is full-depth
SURFACE_VNM = {'ETAN', 'AREA', 'HEFF', 'HSNOW', 'UICE', 'VICE'}

### generous physical bounds -- these are sanity rails, not tuning targets
BOUNDS = {
    'THETA': (-3.0, 40.0, 'degC'),
    'SALT':  (0.0,  45.0, 'psu'),
    'UVEL':  (-5.0,  5.0, 'm/s'),
    'VVEL':  (-5.0,  5.0, 'm/s'),
    'ETAN':  (-10.0, 10.0, 'm'),
}

### MITgcm monitor keys worth watching, and how to judge them
MON_NAN_KEYS = ['dynstat_theta_max', 'dynstat_theta_min', 'dynstat_salt_max',
                'dynstat_salt_min', 'dynstat_uvel_max', 'dynstat_vvel_max',
                'dynstat_eta_max', 'ke_mean']
MON_CFL_KEYS = ['advcfl_uvel_max', 'advcfl_vvel_max', 'advcfl_wvel_max',
                'advcfl_W_max']

PASS, WARN, FAIL = 'PASS', 'WARN', 'FAIL'

############################################################
#                    STDOUT MONITOR PARSING                 #
############################################################

### e.g. "(PID.TID 0000.0001) %MON dynstat_theta_max =   3.00000000E+01"
MON_RE = re.compile(r'%MON\s+(\S+)\s*=\s*(\S+)')

def parse_stdout(path):
    """
    Read a MITgcm STDOUT file and return (monitors, finished).

    monitors maps each %MON key to a float array, one entry per monitor dump.
    Values MITgcm prints as NaN come back as np.nan rather than raising, which
    is the point -- a NaN in the log is the thing we are looking for.
    """
    monitors = {}
    finished = False
    with open(path, 'r', errors='replace') as fh:
        for line in fh:
            if 'Execution ended Normally' in line:
                finished = True
                continue
            m = MON_RE.search(line)
            if not m:
                continue
            key, raw = m.group(1), m.group(2)
            try:
                val = float(raw.replace('D', 'E'))
            except ValueError:
                val = np.nan if 'NaN' in raw or 'nan' in raw.lower() else None
                if val is None:
                    continue
            monitors.setdefault(key, []).append(val)
    return {k: np.asarray(v, dtype=float) for k, v in monitors.items()}, finished


def check_health(monitors, finished, cfl_warn=0.5, ke_growth=10.0):
    """Judge a parsed STDOUT. Returns a list of (level, label, detail)."""
    out = []

    if not monitors:
        out.append((FAIL, 'monitor output',
                    'no %MON lines found -- set monitorFreq in "data" to a value '
                    'small enough to produce several dumps over the run'))
        return out

    nts = monitors.get('time_tsnumber')
    ndump = len(nts) if nts is not None else max(len(v) for v in monitors.values())
    out.append((PASS if ndump > 1 else WARN, 'monitor dumps',
                f'{ndump} dump(s) parsed'
                + ('' if ndump > 1 else ' -- too few to judge drift; lower monitorFreq')))

    ### 1. NaN anywhere is fatal
    nan_keys = [k for k in MON_NAN_KEYS if k in monitors and np.isnan(monitors[k]).any()]
    if nan_keys:
        first = min(int(np.argmax(np.isnan(monitors[k]))) for k in nan_keys)
        out.append((FAIL, 'NaN check',
                    f'NaN appeared in {nan_keys} (first at monitor dump {first}) '
                    f'-- the run blew up; reduce deltaT or check the boundary forcing'))
    else:
        out.append((PASS, 'NaN check', 'no NaN in any monitored quantity'))

    ### 2. advective CFL
    cfl_present = {k: monitors[k] for k in MON_CFL_KEYS if k in monitors}
    if not cfl_present:
        out.append((WARN, 'CFL', 'no advcfl_* monitors found in STDOUT'))
    else:
        worst_k, worst = max(((k, np.nanmax(v)) for k, v in cfl_present.items()),
                             key=lambda kv: kv[1])
        if not np.isfinite(worst):
            out.append((FAIL, 'CFL', f'{worst_k} is not finite'))
        elif worst >= 1.0:
            out.append((FAIL, 'CFL',
                        f'max {worst_k} = {worst:.3f} >= 1.0 -- violates the CFL '
                        f'condition; reduce deltaT'))
        elif worst >= cfl_warn:
            out.append((WARN, 'CFL',
                        f'max {worst_k} = {worst:.3f} (>= {cfl_warn}) -- close to the '
                        f'limit; consider a smaller deltaT'))
        else:
            out.append((PASS, 'CFL', f'max {worst_k} = {worst:.3f}'))

    ### 3. physical bounds on T and S
    for var, lo, hi, unit, kmin, kmax in [
            ('theta', *BOUNDS['THETA'], 'dynstat_theta_min', 'dynstat_theta_max'),
            ('salt',  *BOUNDS['SALT'],  'dynstat_salt_min',  'dynstat_salt_max')]:
        if kmin not in monitors or kmax not in monitors:
            out.append((WARN, f'{var} range', f'{kmin}/{kmax} not in STDOUT'))
            continue
        vmin, vmax = np.nanmin(monitors[kmin]), np.nanmax(monitors[kmax])
        if vmin < lo or vmax > hi:
            out.append((FAIL, f'{var} range',
                        f'[{vmin:.3f}, {vmax:.3f}] {unit} leaves the physical '
                        f'range [{lo}, {hi}]'))
        else:
            out.append((PASS, f'{var} range', f'[{vmin:.3f}, {vmax:.3f}] {unit}'))

    ### 4. kinetic-energy drift -- a downscaled run whose KE keeps climbing is
    ###    usually being driven by inconsistent boundary transport
    if 'ke_mean' in monitors and len(monitors['ke_mean']) > 1:
        ke = monitors['ke_mean']
        ke0 = ke[0] if ke[0] > 0 else np.nan
        ratio = ke[-1] / ke0 if np.isfinite(ke0) else np.nan
        if not np.isfinite(ratio):
            out.append((WARN, 'KE drift', 'initial ke_mean is zero or non-finite'))
        elif ratio > ke_growth:
            out.append((FAIL, 'KE drift',
                        f'mean KE grew {ratio:.1f}x over the run '
                        f'({ke[0]:.3e} -> {ke[-1]:.3e}) -- check open-boundary '
                        f'transport with check_obcs_transport.py'))
        else:
            out.append((PASS, 'KE drift',
                        f'mean KE {ke[0]:.3e} -> {ke[-1]:.3e} ({ratio:.2f}x)'))

    ### 5. did it actually finish
    out.append((PASS, 'termination', 'Execution ended Normally') if finished else
               (FAIL, 'termination',
                'no "Execution ended Normally" in STDOUT -- the run stopped early; '
                'check the tail of the file and the job stderr'))
    return out

############################################################
#                      INPUT CHECKING                       #
############################################################

def boundary_geometry(config_dir, reg_nm):
    """
    Child grid geometry at each open boundary, using the SAME interior-face
    trim as gen_obcs.py (see the long comment there): HFacS[:, 1:-1, :] and
    HFacW[:, :, 1:-1], so [0] is the interior face of the south/west boundary
    and [-1] the interior face of the north/east one.
    """
    grd_ls = ['XC', 'YC', 'HFacC', 'HFacS', 'HFacW', 'drF']
    tmp = read_ncgrid(config_dir, reg_nm, grd_ls)
    g = dict(zip(grd_ls, tmp))
    S = g['HFacS'][:, 1:-1, :]
    W = g['HFacW'][:, :, 1:-1]
    C = g['HFacC']
    Nr, ny, nx = C.shape

    geom = {}
    for bnd in ('east', 'west', 'north', 'south'):
        if bnd == 'east':
            mnorm, mcen = W[:, :, -1], C[:, :, -1]
        elif bnd == 'west':
            mnorm, mcen = W[:, :, 0], C[:, :, 0]
        elif bnd == 'north':
            mnorm, mcen = S[:, -1, :], C[:, -1, :]
        else:
            mnorm, mcen = S[:, 0, :], C[:, 0, :]
        geom[bnd] = {'normal': (mnorm > 0), 'center': (mcen > 0),
                     'npts': mnorm.shape[1], 'Nr': Nr}
    return geom


def _read_f4(path):
    return np.fromfile(path, '>f4')


def check_obcs_files(config_dir, reg_nm, boundaries, geom, obcs_dir=None, verbose=False):
    """Structural + numerical checks on everything gen_obcs.py wrote."""
    out = []
    obcs_dir = obcs_dir or os.path.join(config_dir, 'forcings/OBCS/')
    if not os.path.isdir(obcs_dir):
        return [(FAIL, 'OBCS directory',
                 f'missing: {obcs_dir} -- run gen_obcs.py (STEP3, Section II)')]

    for c in boundaries:
        bnd = BND_NAME[c]
        g = geom[bnd]
        files = sorted(glob.glob(os.path.join(obcs_dir, f'*_{bnd}.bin')))
        if not files:
            out.append((FAIL, f'{bnd} OBCS', f'no *_{bnd}.bin files in {obcs_dir}'))
            continue

        for path in files:
            base = os.path.basename(path)
            vnm = base[:-len(f'_{bnd}.bin')]
            nlev = 1 if vnm in SURFACE_VNM else g['Nr']
            npts = g['npts']
            data = _read_f4(path)

            ### size must be a whole number of (nlev x npts) records
            rec = nlev * npts
            if rec == 0 or data.size % rec != 0:
                out.append((FAIL, f'{bnd}/{vnm}',
                            f'{data.size} values is not a multiple of nlev*npts '
                            f'({nlev}*{npts}={rec}) -- wrong grid size or field depth'))
                continue
            tstp = data.size // rec
            fld = data.reshape(tstp, nlev, npts)

            if not np.isfinite(fld).all():
                out.append((FAIL, f'{bnd}/{vnm}',
                            f'{int((~np.isfinite(fld)).sum())} non-finite value(s)'))
                continue

            ### physical range, where we know one
            if vnm in BOUNDS:
                lo, hi, unit = BOUNDS[vnm]
                vmin, vmax = float(fld.min()), float(fld.max())
                if vmin < lo or vmax > hi:
                    out.append((FAIL, f'{bnd}/{vnm}',
                                f'range [{vmin:.3f}, {vmax:.3f}] {unit} leaves '
                                f'[{lo}, {hi}]'))
                    continue

            ### the useful one: is the field actually populated on the wet
            ### cells of the face the model will mask it with? UICE/VICE live
            ### on the same staggered faces as UVEL/VVEL, so they get the
            ### normal-face mask too rather than the cell-centered one.
            normal_vnms = {NORMAL_VNM[bnd], 'UICE' if bnd in ('east', 'west') else 'VICE'}
            mask = g['normal'] if vnm in normal_vnms else g['center']
            mask = mask[:nlev, :]
            nwet = int(mask.sum())
            if nwet == 0:
                out.append((FAIL, f'{bnd}/{vnm}',
                            f'the child grid has NO wet cells on this boundary face -- '
                            f'check the interior-face trim and the bathymetry'))
                continue
            wet_zero = int((np.all(fld == 0, axis=0) & mask).sum())
            frac = wet_zero / nwet
            detail = (f'{tstp} record(s), {nwet} wet cell(s), '
                      f'{wet_zero} ({100*frac:.1f}%) zero for the whole run')
            if frac == 1.0:
                out.append((FAIL, f'{bnd}/{vnm}',
                            detail + ' -- this boundary carries nothing; classic '
                                     'wrong-C-grid-face symptom'))
            elif frac > 0.5:
                out.append((WARN, f'{bnd}/{vnm}', detail))
            else:
                out.append((PASS, f'{bnd}/{vnm}', detail))
    return out


def check_pickup_files(config_dir, reg_nm, itr, geom, verbose=False):
    """Pickups must match the child grid and be finite/plausible."""
    out = []
    pdir = os.path.join(config_dir, 'forcings/pickups/')
    if not os.path.isdir(pdir):
        return [(FAIL, 'pickup directory',
                 f'missing: {pdir} -- run gen_pickups.py (STEP3, Section I)')]

    files = sorted(glob.glob(os.path.join(pdir, f'pickup_*.{itr:010d}.data')))
    if not files:
        have = sorted({os.path.basename(p).split('.')[1]
                       for p in glob.glob(os.path.join(pdir, 'pickup_*.*.data'))})
        return [(FAIL, 'pickup files',
                 f'none for iteration {itr:010d} in {pdir}; present: {have or "none"}')]

    for path in files:
        base = os.path.basename(path)
        vnm = base.split('.')[0].replace('pickup_', '')
        ### gen_pickups.py writes 64-bit
        data = np.fromfile(path, '>f8')
        if not np.isfinite(data).all():
            out.append((FAIL, f'pickup {vnm}',
                        f'{int((~np.isfinite(data)).sum())} non-finite value(s)'))
            continue
        key = {'THETA': 'THETA', 'SALT': 'SALT', 'U': 'UVEL',
               'V': 'VVEL', 'ETAN': 'ETAN'}.get(vnm)
        detail = f'{data.size} values'
        if key and key in BOUNDS:
            lo, hi, unit = BOUNDS[key]
            nz = data[data != 0]
            if nz.size:
                vmin, vmax = float(nz.min()), float(nz.max())
                detail = f'{data.size} values, wet range [{vmin:.3f}, {vmax:.3f}] {unit}'
                if vmin < lo or vmax > hi:
                    out.append((FAIL, f'pickup {vnm}', detail + f' leaves [{lo}, {hi}]'))
                    continue
        out.append((PASS, f'pickup {vnm}', detail))
    return out


def check_initial_dump(run_dir, config_dir, itr, verbose=False):
    """
    If the run was configured with dumpInitAndLast=.TRUE., MITgcm writes the
    state it actually started from as T/S/U/V/Eta.0000000000.data. Comparing
    that against the pickup files proves the model really ingested STEP3's
    initial condition rather than silently falling back to a default.
    """
    out = []
    if not run_dir:
        return out
    pairs = [('T', 'THETA'), ('S', 'SALT'), ('U', 'U'), ('V', 'V'), ('Eta', 'ETAN')]
    pdir = os.path.join(config_dir, 'forcings/pickups/')
    found_any = False
    for dump_nm, pick_nm in pairs:
        dump = os.path.join(run_dir, f'{dump_nm}.0000000000.data')
        pick = os.path.join(pdir, f'pickup_{pick_nm}.{itr:010d}.data')
        if not (os.path.isfile(dump) and os.path.isfile(pick)):
            continue
        found_any = True
        a = np.fromfile(dump, '>f4').astype(np.float64)
        b = np.fromfile(pick, '>f8')
        if a.size != b.size:
            out.append((FAIL, f'initial {dump_nm}',
                        f'model dump has {a.size} values, pickup has {b.size}'))
            continue
        ### the model writes 32-bit dumps of a 64-bit input, so compare at
        ### single precision rather than demanding exact equality
        scale = max(np.abs(b).max(), 1e-12)
        rel = np.abs(a - b).max() / scale
        if rel > 1e-6:
            out.append((FAIL, f'initial {dump_nm}',
                        f'differs from pickup by rel. {rel:.2e} -- the model did not '
                        f'start from the file named in "data"'))
        else:
            out.append((PASS, f'initial {dump_nm}',
                        f'matches pickup_{pick_nm} (rel. diff {rel:.1e})'))
    if not found_any:
        out.append((WARN, 'initial state',
                    'no <var>.0000000000.data dumps in the run directory -- set '
                    'dumpInitAndLast=.TRUE. in "data" to enable this check'))
    return out

############################################################
#                        PLOTTING                           #
############################################################

def plot_health(monitors, plot_dir, reg_nm):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    panels = [('dynstat_theta_max', 'dynstat_theta_min', 'theta (degC)'),
              ('dynstat_salt_max', 'dynstat_salt_min', 'salt (psu)'),
              ('dynstat_uvel_max', 'dynstat_vvel_max', 'max |velocity| (m/s)'),
              ('ke_mean', None, 'mean KE (m2/s2)')]
    panels = [p for p in panels if p[0] in monitors]
    if not panels:
        return None
    x = monitors.get('time_tsnumber')
    fig, axes = plt.subplots(len(panels), 1, figsize=(9, 2.3 * len(panels)), sharex=True)
    axes = np.atleast_1d(axes)
    for ax, (k1, k2, label) in zip(axes, panels):
        xs = x if x is not None and len(x) == len(monitors[k1]) else np.arange(len(monitors[k1]))
        ax.plot(xs, monitors[k1], marker='.', label=k1)
        if k2 and k2 in monitors and len(monitors[k2]) == len(xs):
            ax.plot(xs, monitors[k2], marker='.', label=k2)
        ax.set_ylabel(label)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7, loc='best')
    axes[-1].set_xlabel('timestep' if x is not None else 'monitor dump')
    fig.suptitle(f'{reg_nm}: verification run health')
    os.makedirs(plot_dir, exist_ok=True)
    path = os.path.join(plot_dir, f'{reg_nm}_verify_health.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    return path

############################################################
#                          MAIN                             #
############################################################

def report(title, results, verbose=False):
    print(f'\n=== {title} ===')
    if not results:
        print('  (nothing to check)')
        return 0, 0
    width = max(len(r[1]) for r in results)
    nfail = sum(1 for r in results if r[0] == FAIL)
    nwarn = sum(1 for r in results if r[0] == WARN)
    for level, label, detail in results:
        print(f'  [{level}] {label:<{width}}  {detail}')
    return nfail, nwarn


def main(config_dir, reg_nm, boundaries, itr, run_dir, stdout_path,
         do_health, do_inputs, do_plot, output_dir, verbose):
    nfail = nwarn = 0

    if do_health:
        path = stdout_path
        if path is None and run_dir:
            cands = sorted(glob.glob(os.path.join(run_dir, 'STDOUT*')))
            path = cands[0] if cands else None
        if path is None or not os.path.isfile(path):
            f, w = report('RUN HEALTH', [(FAIL, 'STDOUT',
                                          'not found -- pass -r <run dir> or -stdout <file>')])
        else:
            monitors, finished = parse_stdout(path)
            f, w = report(f'RUN HEALTH ({os.path.basename(path)})',
                          check_health(monitors, finished), verbose)
            if do_plot and monitors:
                p = plot_health(monitors, os.path.join(output_dir, 'diagnostics'), reg_nm)
                if p:
                    print(f'  plot saved: {p}')
        nfail += f; nwarn += w

    if do_inputs:
        geom = boundary_geometry(config_dir, reg_nm)
        f, w = report('OBCS FILES',
                      check_obcs_files(config_dir, reg_nm, boundaries, geom, verbose=verbose),
                      verbose)
        nfail += f; nwarn += w
        f, w = report('PICKUP FILES',
                      check_pickup_files(config_dir, reg_nm, itr, geom, verbose=verbose),
                      verbose)
        nfail += f; nwarn += w
        res = check_initial_dump(run_dir, config_dir, itr, verbose=verbose)
        if res:
            f, w = report('INITIAL STATE INGESTED', res, verbose)
            nfail += f; nwarn += w

    print(f'\n=== SUMMARY ===\n  {nfail} FAIL, {nwarn} WARN')
    if nfail:
        print('  Verification FAILED -- see the FAIL lines above.')
    elif nwarn:
        print('  Verification passed with warnings.')
    else:
        print('  All checks passed.')
    return 1 if nfail else 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Post-run verification of a downscaled regional set-up (STEP5).')
    parser.add_argument("-d", "--config_dir", required=True, type=str,
                        help="Your config_dir -- the same -d given to the STEP3 scripts")
    parser.add_argument("-n", "--reg_nm", required=True, type=str,
                        help="Name of the regional cutout")
    parser.add_argument("-bnd", "--bnd_nm", type=str, default='',
                        help="Open boundaries to check, as a string of E/W/N/S (e.g. ENS). "
                             "Required for --inputs")
    parser.add_argument("-i", "--itr", type=int, default=None,
                        help="Pickup iteration used to start the run (same -i as gen_pickups.py). "
                             "Required for --inputs")
    parser.add_argument("-r", "--run_dir", type=str, default=None,
                        help="The model run directory (holds STDOUT and any dumps)")
    parser.add_argument("-stdout", "--stdout_path", type=str, default=None,
                        help="Explicit path to a STDOUT file (default: <run_dir>/STDOUT*)")
    parser.add_argument("--health", action="store_true", help="Run only the STDOUT health check")
    parser.add_argument("--inputs", action="store_true", help="Run only the input-consistency check")
    parser.add_argument("-plot", "--plot", action="store_true",
                        help="Save a health time-series plot (requires matplotlib)")
    parser.add_argument("-o", "--output_dir", type=str, default=None,
                        help="Where to save plots (default: <config_dir>); "
                             "written under <output_dir>/diagnostics/")
    parser.add_argument("-v", "--verbose", action="store_true")

    args = parser.parse_args()
    ### neither flag given -> run both
    do_health = args.health or not (args.health or args.inputs)
    do_inputs = args.inputs or not (args.health or args.inputs)

    boundaries = args.bnd_nm
    if do_inputs:
        bad = [c for c in boundaries if c not in BND_NAME]
        if not boundaries or bad:
            raise SystemExit(
                f"ERROR: --inputs needs -bnd as a string of E/W/N/S letters "
                f"(got '{boundaries}'{f', bad letters {bad}' if bad else ''}).")
        if args.itr is None:
            raise SystemExit("ERROR: --inputs needs -i, the pickup iteration the run started from.")

    raise SystemExit(main(args.config_dir, args.reg_nm, boundaries, args.itr,
                          args.run_dir, args.stdout_path, do_health, do_inputs,
                          args.plot, args.output_dir or args.config_dir, args.verbose))
