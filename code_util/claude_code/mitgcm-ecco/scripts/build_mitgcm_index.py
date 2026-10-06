#!/usr/bin/env python3
"""Build a searchable Markdown index of an MITgcm source tree (and optionally darwin3).

Usage:
  build_mitgcm_index.py <mitgcm_root> <out_dir> [--darwin3 <darwin3_root>] [--label TEXT]

Writes into <out_dir>:
  README.md            what was indexed (commit/label) and how to regenerate
  packages.md          catalogue of every package
  pkg/<name>.md        per package: description, docs, deps, data.* file, namelist params
                       (with PARAMS-header descriptions + CPP guards), CPP options, headers,
                       routines, call sites outside the package, verification usage
  verification.md      catalogue of every verification experiment
  verification/<exp>.md per experiment: README, build variants, packages, SIZE.h,
                       code overrides, input variants (enabled packages, key params), results
  core-params.md       model/src namelists (data PARM01-05, data.pkg, eedata)
  core-calltree.md     static call tree from THE_MODEL_MAIN through model/src
  docs.md              section map of the Sphinx manual (every heading + line, params each section cites)
  bibliography.md      manual bib entries (MITgcm + darwin3) and every place the manual cites each one
  docs-ext.md          section map of ECCO / ecco_darwin / darwin3 docs (with --docroot)
  tools.md             genmake2 / testreport help text
The parser is regex-based and heuristic; always confirm details in the source.
"""
import os, re, sys, subprocess, collections, argparse

ap = argparse.ArgumentParser()
ap.add_argument('root'); ap.add_argument('out')
ap.add_argument('--darwin3'); ap.add_argument('--label', default='')
ap.add_argument('--extra', action='append', default=[],
                help='LABEL=ROOT[:pkg1,pkg2] index packages/experiments of a branch clone that are new or differ from upstream (or only the listed pkgs)')
ap.add_argument('--docroot', action='append', default=[],
                help='LABEL=ROOT[:substr1,substr2] map the .rst/.md/readme docs of another repo (ECCO, ecco_darwin, darwin3 manual) into docs-ext.md; optional path substrings filter files')
A = ap.parse_args()
ROOT, OUT = os.path.abspath(A.root), os.path.abspath(A.out)

def rd(p):
    try:
        with open(p, errors='replace') as f: return f.read()
    except OSError: return ''

def wr(rel, text):
    p = os.path.join(OUT, rel); os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, 'w') as f: f.write(text.rstrip() + '\n')

def is_comment(line):
    return bool(line) and line[0] in 'Cc*!' and not line.startswith('#')

def strip_inline(line):
    return line.split('!')[0]

SRC_EXT = ('.F', '.F90', '.h', '.FOR')

def fortran_files(d):
    for dp, dn, fn in os.walk(d):
        for f in sorted(fn):
            if f.endswith(SRC_EXT): yield os.path.join(dp, f)

# ---------------------------------------------------------------- descriptions of variables
DESC_RE = re.compile(r'^[Cc!]\s+(\w+(?:\s*,\s*\w+)*)\s*::\s*(.*)$')
def var_descriptions(files):
    d = {}
    for p in files:
        if not p.endswith('.h'): continue
        lines = rd(p).splitlines()
        for i, l in enumerate(lines):
            m = DESC_RE.match(l)
            if not m: continue
            txt = m.group(2).strip()
            # absorb continuation comment lines (deeply indented, no ::)
            j = i + 1
            while j < len(lines) and re.match(r'^[Cc!]\s{8,}\S', lines[j]) and '::' not in lines[j] and len(txt) < 220:
                txt += ' ' + lines[j][1:].strip(); j += 1
            for n in m.group(1).split(','):
                d.setdefault(n.strip().lower(), txt)
    return d

# ---------------------------------------------------------------- namelists
def namelists(path):
    """Return {group: [(param, cppguard)]} and data file names opened."""
    lines = rd(path).splitlines()
    groups, files, stack = collections.OrderedDict(), [], []
    i = 0
    while i < len(lines):
        l = lines[i]
        for m in re.finditer(r"OPEN_COPY_DATA_FILE\s*\(\s*$", l, re.I):
            pass
        m = re.search(r"'(data(?:\.[\w.]+)?)'", l)
        if m and ('OPEN_COPY_DATA_FILE' in l.upper() or 'OPEN_COPY_DATA_FILE' in lines[i-1].upper()):
            files.append(m.group(1))
        if l.startswith('#'):
            d = l[1:].strip()
            if re.match(r'if(n?def)?\b', d): stack.append(d)
            elif d.startswith('else') and stack: stack[-1] = 'NOT(' + stack[-1] + ')'
            elif d.startswith('endif') and stack: stack.pop()
        m = re.match(r'^\s+NAMELIST\s*/\s*(\w+)\s*/(.*)$', l, re.I) if not is_comment(l) else None
        if m:
            g = m.group(1).upper(); rest = m.group(2)
            params = groups.setdefault(g, [])
            def add(txt):
                guard = ' & '.join(s for s in stack if s)
                for p in strip_inline(txt).split(','):
                    p = p.strip().strip('&').strip()
                    if re.match(r'^\w+$', p): params.append((p, guard))
            add(rest)
            j = i + 1
            while j < len(lines):
                n = lines[j]
                if n.startswith('#'):
                    d = n[1:].strip()
                    if re.match(r'if(n?def)?\b', d): stack.append(d)
                    elif d.startswith('else') and stack: stack[-1] = 'NOT(' + stack[-1] + ')'
                    elif d.startswith('endif') and stack: stack.pop()
                    j += 1; continue
                if is_comment(n) or not n.strip(): j += 1; continue
                if len(n) > 5 and n[5] not in ' 0' and not n[:5].strip():
                    add(n[6:]); j += 1; continue
                break
            i = j; continue
        i += 1
    return groups, files

# ---------------------------------------------------------------- CPP options
def cpp_options(path, pkg=None):
    out, cmt = [], []
    for l in rd(path).splitlines():
        if is_comment(l):
            t = l[1:].strip().strip('|').strip()
            if t and not set(t) <= set('=-*o '): cmt.append(t)
            continue
        m = re.match(r'^#\s*(define|undef)\s+(\w+)', l)
        if m:
            name = m.group(2)
            skip = name.endswith('_H') or (pkg and name == 'ALLOW_' + pkg.upper())
            if not skip:
                out.append((name, m.group(1), ' '.join(cmt[-3:])[:220]))
            cmt = []
        elif not l.strip():
            cmt = []
    return out

def header_blurb(path):
    """First meaningful descriptive comment lines of a file."""
    txt = []
    lines = rd(path).splitlines()
    for idx, l in enumerate(lines[:80]):
        if not is_comment(l):
            if txt and not l.startswith('#'): break
            continue
        t = l[1:].strip().strip('|').strip()
        t = re.sub(r'^o\s+', '', t)
        if not t or set(t) <= set('=-*o ') or t.startswith('!') and len(t) < 3: continue
        if re.match(r'^(!?(ROUTINE|INTERFACE|DESCRIPTION|USES|INPUT|OUTPUT|LOCAL|REVISION|CALLING|CBOP|CEOP)\b)', t.upper().lstrip('!')):
            if 'DESCRIPTION' in t.upper(): continue
            if txt: break
            continue
        if re.match(r'^#?include|^(SUBROUTINE|FUNCTION)\b', t, re.I): continue
        if t.upper().startswith(os.path.basename(path).upper()): continue
        txt.append(t)
        if len(' '.join(txt)) > 200: break
    return ' '.join(txt)[:240]

# ---------------------------------------------------------------- routine map + calls
SUB_RE = re.compile(r'^\s+(?:\w+\s+)*?(SUBROUTINE|FUNCTION)\s+(\w+)', re.I)
CALL_RE = re.compile(r'\bCALL\s+(\w+)', re.I)

def scan_routines(root, subdirs):
    defs = {}   # NAME -> (relfile, owner)
    calls = []  # (caller_relfile, lineno, CALLEE, guards)
    for sd in subdirs:
        base = os.path.join(root, sd)
        if not os.path.isdir(base): continue
        for p in fortran_files(base):
            rel = os.path.relpath(p, root)
            parts = rel.split(os.sep)
            owner = parts[1] if parts[0] == 'pkg' else parts[0] + '/' + parts[1] if len(parts) > 2 else parts[0]
            stack = []
            for n, l in enumerate(rd(p).splitlines(), 1):
                if l.startswith('#'):
                    d = l[1:].strip()
                    if re.match(r'if(n?def)?\b', d): stack.append(d)
                    elif d.startswith('else') and stack: stack[-1] = 'NOT(' + stack[-1] + ')'
                    elif d.startswith('endif') and stack: stack.pop()
                    continue
                if is_comment(l): continue
                code = strip_inline(l)
                m = SUB_RE.match(code)
                if m and not code.strip().upper().startswith('END'):
                    defs.setdefault(m.group(2).upper(), (rel, owner))
                for c in CALL_RE.findall(code):
                    calls.append((rel, n, c.upper(), tuple(stack)))
    return defs, calls

# ---------------------------------------------------------------- verification
def parse_size(path):
    t = rd(path); d = {}
    for k in ('sNx', 'sNy', 'OLx', 'OLy', 'nSx', 'nSy', 'nPx', 'nPy', 'Nr'):
        m = re.search(r'\b' + k + r'\s*=\s*(\d+)', t)
        if m: d[k] = int(m.group(1))
    return d

def grid_str(d):
    try:
        return '%dx%dx%d' % (d['sNx']*d['nSx']*d['nPx'], d['sNy']*d['nSy']*d['nPy'], d['Nr'])
    except KeyError: return '?'

def nml_values(path, keys):
    t = rd(path); out = {}
    for k in keys:
        m = re.search(r'^\s*' + k + r'\s*=\s*([^,\n]+)', t, re.I | re.M)
        if m: out[k] = m.group(1).strip()
    return out

def pkgs_conf(path):
    return [l.split('#')[0].strip() for l in rd(path).splitlines() if l.split('#')[0].strip()]

def load_groups(root):
    g = {}
    for l in rd(os.path.join(root, 'pkg/pkg_groups')).splitlines():
        l = l.split('#')[0]
        if ':' in l:
            k, v = l.split(':', 1); g[k.strip()] = v.split()
    return g

def expand(pkgs, groups, seen=None):
    out = []
    for p in pkgs:
        neg = p.startswith('-')
        if neg: continue
        if p in groups: out += expand(groups[p], groups)
        else: out.append(p)
    negs = {p[1:] for p in pkgs if p.startswith('-')}
    return [p for p in dict.fromkeys(out) if p not in negs]

def git_label(root):
    try:
        return subprocess.run(['git', '-C', root, 'log', '-1', '--format=%h %ad %s', '--date=short'],
                              capture_output=True, text=True).stdout.strip()
    except Exception: return ''


PKG_DESC = {
 'admtlm': 'Adjoint/tangent-linear combined (ADM-TLM) driver for singular-vector / Hessian-type calculations (legacy).',
 'aim_v23': 'AIM intermediate-complexity atmospheric physics (SPEEDY v23: convection, clouds, radiation, surface fluxes) for atmosphere set-ups.',
 'atm2d': '2-D (zonally averaged) atmosphere coupled to an MITgcm ocean (IGSM-style climate coupling).',
 'atm_common': 'Shared diagnostics/variables for atmospheric physics packages.',
 'atm_compon_interf': 'Atmosphere-side interface of the coupler (cpl_aim+ocn style atmosphere-ocean coupling).',
 'atm_ocn_coupler': 'Stand-alone coupler component exchanging fields between atmosphere and ocean MITgcm components.',
 'atm_phys': 'Grey-radiation idealized atmospheric physics (Frierson/O\'Gorman-type moist physics).',
 'autodiff': 'Automatic differentiation support (TAF/Tapenade/OpenAD): checkpointing, tape storage directives, adjoint I/O and dumps.',
 'bbl': 'Bottom boundary layer: a thin bottom layer with its own T/S exchanging with the interior and downslope (dense-overflow representation; mutually exclusive with down_slope).',
 'bling': 'BLING biogeochemistry (Biogeochemistry with Light, Iron, Nutrients and Gases; Galbraith et al.) via gchem/ptracers.',
 'bulk_force': 'Simple bulk-formula surface forcing (older alternative to exf).',
 'cal': 'Calendar tool (Gregorian/360-day/model calendar) used by exf, ecco, ctrl, diagnostics with calendarDumps.',
 'cd_code': 'C-D scheme: carries D-grid velocities to stabilise Coriolis on the C-grid at coarse resolution (useCDscheme).',
 'cfc': 'CFC-11/CFC-12 ocean tracers with air-sea gas exchange (OCMIP protocol) via gchem/ptracers.',
 'cheapaml': 'Cheap atmospheric mixed layer: simple prognostic atmospheric boundary layer above the ocean for surface fluxes.',
 'chronos': 'Clock/date utilities (legacy).',
 'compon_communic': 'MPI component communication layer for coupled (multi-executable) runs.',
 'cost': 'Cost-function framework for adjoint/state estimation (accumulates and writes the cost).',
 'ctrl': 'Control-vector handling for adjoint/optimisation: generic 2-D/3-D/time-varying controls, packing/unpacking, preconditioning, smoothing hooks.',
 'debug': 'Debug utilities (debugMode, debugLevel, field stats printing, DEBUG_ENTER/LEAVE).',
 'diagnostics': 'Diagnostics framework: data.diagnostics output streams (time averages/snapshots, levels, frequencies) and statistics-diagnostics; packages register fields via DIAGNOSTICS_ADDTOLIST and fill with DIAGNOSTICS_FILL.',
 'dic': 'Simple dissolved inorganic carbon biogeochemistry (DIC, ALK, PO4, DOP, O2, Fe) with carbonate chemistry and air-sea CO2 flux (OCMIP-style).',
 'down_slope': 'Down-slope flow parameterisation (Campin & Goosse 1999) moving dense bottom water downslope (excludes bbl).',
 'ebm': 'Energy-balance atmosphere model coupled to the ocean.',
 'ecco': 'ECCO state-estimation package: model-data misfit cost terms, generic cost (gencost), averaging, ECCO-specific I/O.',
 'embed_files': 'Embeds input files into the executable (e.g. for testing/portability).',
 'exch2': 'Generalised tile exchanges for cubed-sphere and LLC grids (facets, blank tiles via data.exch2, W2_mapIO).',
 'exf': 'External forcing: reads/interpolates/time-interpolates atmospheric state or fluxes, bulk formulae (Large & Yeager), runoff, open-boundary-friendly interpolation; standard ECCO forcing path.',
 'fizhi': 'Fizhi atmospheric physics (NASA GEOS-like) for atmosphere configurations.',
 'flt': 'Lagrangian floats/particles advected online (2-D/3-D, profiling floats).',
 'frazil': 'Frazil ice formation: removes supercooling and adjusts heat/salt (used with shelfice/seaice set-ups).',
 'gchem': 'Geochemistry driver: interface between ptracers and BGC packages (dic, bling, cfc, darwin); separate forcing step, surface forcing, chemistry calls.',
 'generic_advdiff': 'Generic advection-diffusion operators (2nd/3rd/4th order, DST, flux-limited, OS7MP, Prather SOM, multi-dim) used for T, S and ptracers.',
 'ggl90': 'Gaspar-Gregoris-Lefevre (1990) TKE vertical mixing scheme (ECCO v4 default), with IDEMIX and Langmuir options.',
 'gmredi': 'Gent-McWilliams / Redi isopycnal eddy parameterisation (GM, Redi, tapering, Visbeck, GEOMETRIC, advective/skew forms).',
 'grdchk': 'Gradient check: compares adjoint gradient against finite differences (data.grdchk).',
 'gridalt': 'Alternative vertical grid support for atmospheric physics (fizhi).',
 'icefront': 'Vertical ice-front (tidewater glacier face) melt parameterisation at side walls.',
 'kl10': 'Klymak & Legg (2010) internal-wave-breaking vertical mixing scheme.',
 'kpp': 'K-Profile Parameterization (Large et al. 1994) vertical mixing with nonlocal transport.',
 'land': 'Simple land model (bucket hydrology, soil temperature) for atmospheric configurations.',
 'layers': 'Diagnostics of transport in isopycnal/temperature/salinity layers (residual overturning).',
 'longstep': 'Longer time step for passive tracers than for dynamics (ptracers stepped every N steps with averaged velocities).',
 'matrix': 'Transport-matrix method: extracts explicit/implicit tracer transport matrices.',
 'mdsio': 'MDS (.data/.meta) binary I/O routines; always compiled.',
 'mnc': 'NetCDF I/O (per-tile files) for state, pickups, diagnostics when diag_mnc is on.',
 'mom_common': 'Momentum code shared by flux-form and vector-invariant: viscosity, bottom/side drag, Smagorinsky/Leith, metric terms.',
 'mom_fluxform': 'Flux-form momentum equations (default unless vectorInvariantMomentum=.TRUE.).',
 'mom_vecinv': 'Vector-invariant momentum equations (used by ECCO/cubed-sphere/LLC configs; enstrophy/energy-conserving Coriolis).',
 'monitor': 'Monitor statistics (%MON lines in STDOUT: CFL, KE, field min/max/mean/sd) used by testreport.',
 'my82': 'Mellor-Yamada (1982) level-2 vertical mixing.',
 'mypackage': 'TEMPLATE package: skeleton showing how to write a new package (readparms, check, init, diagnostics, pickup, tendencies). Start new packages from it.',
 'obcs': 'Open boundary conditions: prescribed/relaxed boundary values (T,S,U,V,eta,seaice, ptracers), sponges, Orlanski radiation, Stevens and Flather schemes, tidal forcing, balancing.',
 'obsfit': 'Model-data comparison for generic (unstructured) observations in the cost function.',
 'ocn_compon_interf': 'Ocean-side interface of the coupler for coupled atmosphere-ocean runs.',
 'offline': 'Offline mode: reads precomputed velocities/diffusivities/forcing (e.g. ECCO archive) to drive passive tracers/BGC without dynamics.',
 'openad': 'OpenAD automatic-differentiation support (legacy).',
 'opps': 'OPPS convection (Paluszkiewicz & Romea penetrative plume scheme).',
 'pp81': 'Pacanowski & Philander (1981) Richardson-number vertical mixing.',
 'profiles': 'Model-data comparison for in situ profiles (Argo, CTD, XBT...) for ECCO cost.',
 'ptracers': 'Passive tracers: arbitrary number of tracers with their own advection schemes, diffusivities, initial files, surface/relaxation; carrier for BGC.',
 'rbcs': 'Relaxation boundary conditions: 3-D relaxation of T, S, ptracers (and U/V) toward prescribed fields with masks/timescales.',
 'regrid': 'Regridding utilities for diagnostics output.',
 'runclock': 'Wall-clock based run termination (stop gracefully before walltime).',
 'rw': 'Read/write utilities on top of mdsio (READ_FLD_XY_RL etc.); always compiled.',
 'salt_plume': 'Brine rejection from sea-ice growth distributed vertically (salt plume parameterisation, ECCO v4).',
 'sbo': 'Solid-body ocean diagnostics: angular momentum, centre of mass (Earth rotation studies).',
 'seaice': 'Dynamic-thermodynamic sea ice (viscous-plastic/EVP/JFNK/Krylov dynamics, zero-layer/ITD thermodynamics, snow, advection).',
 'shap_filt': 'Shapiro filter for grid-scale noise.',
 'shelfice': 'Ice-shelf cavities: pressure loading, three-equation melt, boundary-layer options, conservative remeshing.',
 'showflops': 'FLOP counting/timing (PAPI).',
 'smooth': 'Diffusion-operator smoothing for control/covariance (2-D/3-D correlation operators).',
 'sphere': 'Spherical-harmonic utilities.',
 'steep_icecavity': 'Ice cavities with steep ice draft (alternative shelfice treatment).',
 'streamice': 'Ice-stream/shelf dynamics model (shallow-shelf / hybrid), optionally coupled to shelfice.',
 'tapenade': 'Tapenade automatic-differentiation support.',
 'thsice': 'Thermodynamic sea ice (Winton 2000 two-layer) usable alone or with seaice dynamics.',
 'zonal_filt': 'Zonal (polar) Fourier filter for lat-lon grids.',
 'darwin': 'darwin3 Darwin ecosystem model: many plankton types/groups (allometric or random traits), nutrients, carbon chemistry, iron, radtrans light; called via gchem.',
 'radtrans': 'darwin3 spectral radiative transfer (direct/diffuse irradiance in nlam wavebands, OASIM forcing) used by Darwin.',
}

# ================================================================= main
os.makedirs(OUT, exist_ok=True)
groups = load_groups(ROOT)
depend = collections.defaultdict(list)
for l in rd(os.path.join(ROOT, 'pkg/pkg_depend')).splitlines():
    l = l.split('#')[0].split()
    if len(l) > 1: depend[l[0]] += l[1:]

import filecmp
def dir_differs(a, b):
    """True if directory trees a and b differ in file names or contents (ignores build junk)."""
    if not os.path.isdir(b): return True
    c = filecmp.dircmp(a, b, ignore=['__pycache__', '.DS_Store'])
    if c.left_only or c.right_only: return True
    _, mism, err = filecmp.cmpfiles(a, b, c.common_files, shallow=False)
    if mism or err: return True
    return any(dir_differs(os.path.join(a, d), os.path.join(b, d)) for d in c.common_dirs)

def diff_summary(a, b):
    if not os.path.isdir(b): return 'new (not in its upstream base)'
    added, changed = [], []
    for dp, dn, fn in os.walk(a):
        for f in fn:
            p = os.path.join(dp, f); rel = os.path.relpath(p, a); q = os.path.join(b, rel)
            if not os.path.exists(q): added.append(rel)
            elif not filecmp.cmp(p, q, shallow=False): changed.append(rel)
    removed = [os.path.relpath(os.path.join(dp, f), b) for dp, dn, fn in os.walk(b) for f in fn
               if not os.path.exists(os.path.join(a, os.path.relpath(os.path.join(dp, f), b)))]
    out = []
    if added: out.append('added: ' + ', '.join(sorted(added)[:40]))
    if changed: out.append('changed: ' + ', '.join(sorted(changed)[:40]))
    if removed: out.append('removed: ' + ', '.join(sorted(removed)[:40]))
    return '; '.join(out) or 'identical'

# ---- jobs: (root, pkg, label, page_name, diffnote)
jobs = [(ROOT, p, 'MITgcm', p, '', ROOT) for p in sorted(d for d in os.listdir(os.path.join(ROOT, 'pkg'))
                                                   if os.path.isdir(os.path.join(ROOT, 'pkg', d)))]
scans = {ROOT: scan_routines(ROOT, ['model', 'pkg', 'eesupp'])}
defs, calls = scans[ROOT]
if A.darwin3:
    D3 = os.path.abspath(A.darwin3)
    scans[D3] = scan_routines(D3, ['model', 'pkg', 'eesupp'])
    for d in ('darwin', 'radtrans'):
        if os.path.isdir(os.path.join(D3, 'pkg', d)): jobs.append((D3, d, 'darwin3', d, '', D3))
    for k, v in scans[D3][0].items(): defs.setdefault(k, v)

EXTRAS = []   # (label, root, restrict-list or None)
for spec in A.extra:
    label, rest = spec.split('=', 1)
    root, _, only = rest.partition(':')
    EXTRAS.append((label, os.path.abspath(os.path.expanduser(root)), only.split(',') if only else None))

import tempfile
def branch_base(root):
    """Unpack the clone's merge-base with origin/master so we diff the user's own changes only."""
    try:
        sha = subprocess.run(['git', '-C', root, 'merge-base', 'HEAD', 'origin/master'],
                             capture_output=True, text=True, check=True).stdout.strip()
        d = tempfile.mkdtemp(prefix='mitgcm_base_')
        arc = subprocess.run(['git', '-C', root, 'archive', sha, 'pkg', 'verification', 'model'],
                             capture_output=True, check=True).stdout
        subprocess.run(['tar', '-x', '-C', d], input=arc, check=True)
        return d, sha[:9]
    except Exception:
        return ROOT, 'upstream index root'

def snapshot(root):
    """Copy of the clone's tracked + untracked-not-ignored files (drops build-generated junk)."""
    try:
        names = subprocess.run(['git', '-C', root, 'ls-files', '-co', '--exclude-standard', '-z', '--',
                                'pkg', 'verification', 'model', 'eesupp'],
                               capture_output=True, check=True).stdout
        if not names: return root
        d = tempfile.mkdtemp(prefix='mitgcm_base_snap_')
        tar = subprocess.run(['tar', '-c', '--null', '-T', '-', '-C', root], input=names,
                             capture_output=True, check=True, cwd=root).stdout
        subprocess.run(['tar', '-x', '-C', d], input=tar, check=True)
        return d
    except Exception:
        return root

# ---- verification experiments
KEYS = ['deltaT', 'deltaTmom', 'deltaTtracer', 'nTimeSteps', 'nIter0', 'endTime', 'startTime',
        'nonlinFreeSurf', 'select_rStar', 'eosType', 'implicitFreeSurface', 'nonHydrostatic',
        'usingCurvilinearGrid', 'usingSphericalPolarGrid', 'usingCartesianGrid', 'usingCylindricalGrid',
        'buoyancyRelation', 'tempAdvScheme', 'saltAdvScheme', 'momStepping', 'useRealFreshWaterFlux']
exp_info, pkg_used_by = {}, collections.defaultdict(set)

def parse_exp(ed, name, pkg_key=lambda p: p):
    entries = sorted(os.listdir(ed))
    readme = ''
    for r in ('README', 'README.md', 'README.txt', 'readme'):
        if os.path.isfile(os.path.join(ed, r)): readme = rd(os.path.join(ed, r)); break
    codes = [x for x in entries if x.startswith('code') and os.path.isdir(os.path.join(ed, x))]
    inputs = [x for x in entries if x.startswith('input') and os.path.isdir(os.path.join(ed, x))]
    info = {'readme': readme, 'codes': {}, 'inputs': {}, 'diff': '',
            'results': sorted(os.listdir(os.path.join(ed, 'results'))) if 'results' in entries else []}
    for c in codes:
        cd = os.path.join(ed, c)
        files = sorted(os.listdir(cd))
        pc = pkgs_conf(os.path.join(cd, 'packages.conf')) if 'packages.conf' in files else []
        exp_p = expand(pc, groups) if pc else []
        for p in exp_p: pkg_used_by[pkg_key(p)].add(name)
        std = {'packages.conf', 'SIZE.h', 'SIZE.h_mpi', 'CPP_OPTIONS.h', 'CPP_EEOPTIONS.h', 'genmake_local', 'README'}
        overrides = [f for f in files if f not in std and not f.endswith('_OPTIONS.h') and not f.endswith('_SIZE.h')]
        opts = [f for f in files if f.endswith('_OPTIONS.h') or f.endswith('_SIZE.h')]
        info['codes'][c] = {'pkgs': pc, 'expanded': exp_p, 'size': parse_size(os.path.join(cd, 'SIZE.h')),
                            'size_mpi': os.path.isfile(os.path.join(cd, 'SIZE.h_mpi')),
                            'overrides': overrides, 'opts': opts}
    for i in inputs:
        idir = os.path.join(ed, i)
        files = sorted(os.listdir(idir))
        dp = rd(os.path.join(idir, 'data.pkg'))
        on = re.findall(r'^\s*(use\w+)\s*=\s*\.?(?:true|t)\.?', dp, re.I | re.M)
        info['inputs'][i] = {'use': on, 'vals': nml_values(os.path.join(idir, 'data'), KEYS),
                             'files': [f for f in files if f.startswith('data') or f in ('eedata', 'prepare_run')]}
    return info

def exp_dirs(root):
    V = os.path.join(root, 'verification')
    if not os.path.isdir(V): return []
    return sorted(d for d in os.listdir(V) if os.path.isdir(os.path.join(V, d)) and
                  any(x.startswith('code') for x in os.listdir(os.path.join(V, d))))

VER = os.path.join(ROOT, 'verification')
exps = exp_dirs(ROOT)
for e in exps: exp_info[e] = parse_exp(os.path.join(VER, e), e)

# extras: packages and experiments that are new or differ from upstream
extra_rows, extra_exps = [], []
BASES = {}
for label, root, only in EXTRAS:
    BASE, bsha = branch_base(root); BASES[label] = bsha
    shown_root, root = root, snapshot(root)
    scans[root] = scan_routines(root, ['model', 'pkg', 'eesupp'])
    for k, v in scans[root][0].items(): defs.setdefault(k, v)
    plist = only or sorted(d for d in os.listdir(os.path.join(root, 'pkg')) if os.path.isdir(os.path.join(root, 'pkg', d)))
    for p in plist:
        a, b = os.path.join(root, 'pkg', p), os.path.join(BASE, 'pkg', p)
        if os.path.isdir(a) and (only or dir_differs(a, b)):
            jobs.append((root, p, label, f'{p}@{label}', diff_summary(a, b), shown_root))
    if only is None:
        for e in exp_dirs(root):
            a, b = os.path.join(root, 'verification', e), os.path.join(BASE, 'verification', e)
            if dir_differs(a, b):
                nm = f'{e}@{label}'
                exp_info[nm] = parse_exp(a, nm, pkg_key=lambda p, l=label, r=root: p)
                exp_info[nm]['diff'] = diff_summary(a, b)
                extra_exps.append(nm)
    # core files the branch changes (model/src, model/inc); skipped for package-restricted extras
    for sub in (('model/src', 'model/inc') if only is None else ()):
        a, b = os.path.join(root, sub), os.path.join(BASE, sub)
        if os.path.isdir(a) and dir_differs(a, b):
            extra_rows.append((label, root, sub, diff_summary(a, b)))

def exp_line(e):
    info = exp_info[e]
    c = info['codes'].get('code') or next(iter(info['codes'].values()), {})
    first = next((l.strip() for l in info['readme'].splitlines() if l.strip() and not set(l.strip()) <= set('=-#*')), '')
    tests = [r for r in info['results'] if r.startswith('output')]
    return grid_str(c.get('size', {})), ' '.join(c.get('pkgs', [])), first[:110], len(tests)

vlines = ['# Verification experiments', '',
          'Source: `verification/` in the indexed tree. Grid = global nx×ny×Nr from `code/SIZE.h` '
          '(sNx·nSx·nPx etc.). "tests" counts reference outputs in `results/` (incl. adjoint/TLM/secondary).',
          'Per-experiment details: `verification/<exp>.md`. Packages listed as written in packages.conf '
          '(groups like `gfd`/`oceanic` unexpanded; `-pkg` = excluded).', '',
          '| experiment | grid | packages (code/packages.conf) | tests | description |', '|---|---|---|---|---|']
for e in exps:
    g, p, f, n = exp_line(e)
    vlines.append(f'| `{e}` | {g} | {p} | {n} | {f.replace("|", "/")} |')
if extra_exps:
    vlines += ['', "## Experiments in the user's branches (new or modified vs upstream)", '',
               '`name@label` — label = branch clone (see README.md).', '',
               '| experiment | grid | packages | tests | description |', '|---|---|---|---|---|']
    for e in extra_exps:
        g, p, f, n = exp_line(e)
        vlines.append(f'| `{e}` | {g} | {p} | {n} | {f.replace("|", "/")} |')
wr('verification.md', '\n'.join(vlines))

for e in exps + extra_exps:
    info = exp_info[e]
    L = [f'# verification/{e}', '']
    if info['diff']: L += [f"**vs its upstream base (merge-base):** {info['diff']}", '']
    if info['readme']:
        rl = info['readme'].splitlines()
        L += ['## README (first 40 lines)', '```'] + rl[:40] + (['... (truncated)'] if len(rl) > 40 else []) + ['```', '']
    L.append('## Build variants (code*/)')
    for c, ci in info['codes'].items():
        s = ci['size']
        L.append(f"- **{c}**: packages.conf = `{' '.join(ci['pkgs']) or '(inherits/none)'}`")
        if ci['expanded']: L.append(f"  - expanded: {' '.join(ci['expanded'])}")
        if s: L.append(f"  - SIZE.h: grid {grid_str(s)}; " + ', '.join(f'{k}={v}' for k, v in s.items())
                       + ('; has SIZE.h_mpi' if ci['size_mpi'] else ''))
        if ci['opts']: L.append(f"  - option/size headers: {', '.join(ci['opts'])}")
        if ci['overrides']: L.append(f"  - modified/extra source: {', '.join(ci['overrides'])}")
    L += ['', '## Input variants (input*/)']
    for i, ii in info['inputs'].items():
        L.append(f"- **{i}**: data.pkg on: {', '.join(ii['use']) or '-'}")
        if ii['vals']: L.append('  - data: ' + ', '.join(f'{k}={v}' for k, v in ii['vals'].items()))
        if ii['files']: L.append('  - namelist files: ' + ' '.join(ii['files']))
    L += ['', '## Reference results', ' '.join(f'`{r}`' for r in info['results']) or '-', '',
          f"Run: `cd verification; ./testreport -of <optfile> -t {e.split('@')[0]}` (add `-adm`/`-tlm`/`-tap` for adjoint variants)."]
    wr(f'verification/{e}.md', '\n'.join(L))

# ---- docs map
docs = []
for dp, dn, fn in os.walk(os.path.join(ROOT, 'doc')):
    for f in sorted(fn):
        if not f.endswith('.rst'): continue
        p = os.path.join(dp, f); lines = rd(p).splitlines(); title = ''
        for k in range(len(lines) - 1):
            if lines[k].strip() and re.match(r'^[=\-~^#*]{3,}\s*$', lines[k + 1]) and not re.match(r'^[=\-~^#*]{3,}$', lines[k].strip()):
                title = lines[k].strip(); break
        docs.append((os.path.relpath(p, ROOT), title, rd(p)))

def rst_sections(text):
    """[(line, level, title, names)] for every heading; names = parameters/CPP flags/files the section cites."""
    lines, heads, order = text.splitlines(), [], []
    for k in range(len(lines) - 1):
        t, u = lines[k].rstrip(), lines[k + 1].rstrip()
        if (t.strip() and not t.startswith(' ') and re.match(r'^([=\-~^#*"+`\'])\1{2,}$', u)
                and len(u) >= len(t.strip()) - 1 and not re.match(r'^([=\-~^#*"+`\'])\1+$', t.strip())
                and not t.strip().startswith(('..', ':'))):
            over = k > 0 and lines[k - 1].rstrip() == u
            key = ('o' if over else 'u') + u[0]
            if key not in order: order.append(key)
            heads.append((k + 1, order.index(key), t.strip()))
    out = []
    for n, (ln, lev, title) in enumerate(heads):
        body = '\n'.join(lines[ln:heads[n + 1][0] - 1] if n + 1 < len(heads) else lines[ln:])
        names = re.findall(r':varlink:`~?([A-Za-z_][\w]*)', body)
        names += [w for w in re.findall(r'``([A-Za-z_][\w]{3,})``', body) if '_' in w or re.search(r'[a-z][A-Z]', w)]
        names += [m.split('/')[-1] for m in re.findall(r':filelink:`(?:[^`<]*<)?([\w./-]+\.(?:F|h|F90))>?`', body)]
        seen = []
        for w in names:
            if w not in seen: seen.append(w)
        out.append((ln, lev, title, seen))
    return out

DL = ['# MITgcm manual section map (doc/**/*.rst)', '',
      'The online manual (https://mitgcm.readthedocs.io) is built from these files. Every heading is listed with',
      'its line number (`L123`) and the parameters / CPP flags / source files the section cites, so',
      '`grep -n <name> docs.md` finds the section that explains a parameter. Read that line range from the',
      'source tree (`git -C <clone> show origin/master:doc/... | sed -n <a>,<b>p`) for equations and intent.', '']
for r, t, text in sorted(docs):
    DL += ['', f'## `{r}` — {t}']
    for ln, lev, title, names in rst_sections(text):
        if lev > 3: continue
        s = f"{'  ' * lev}- L{ln} {title}"
        if names: s += ' — ' + ', '.join(f'`{w}`' for w in names[:40]) + (f' (+{len(names) - 40})' if len(names) > 40 else '')
        DL.append(s)
wr('docs.md', '\n'.join(DL))

# ---- bibliography: manual_references.bib entries + where the manual cites each one
def parse_bib(text):
    macros, entries = {}, {}
    for m in re.finditer(r'@string\s*\{\s*(\w+)\s*=\s*[{"](.*?)[}"]\s*\}', text, re.I): macros[m.group(1).lower()] = m.group(2)
    for chunk in re.split(r'\n(?=\s*@)', '\n' + text):
        m = re.match(r'\s*@(\w+)\s*\{\s*([^,\s]+)\s*,(.*)', chunk, re.S)
        if not m or m.group(1).lower() in ('string', 'comment', 'preamble'): continue
        f = {}
        for fm in re.finditer(r'(\w+)\s*=\s*(\{(?:[^{}]|\{(?:[^{}]|\{[^{}]*\})*\})*\}|"[^"]*"|\w+)', m.group(3)):
            v = fm.group(2)
            v = v[1:-1] if v[0] in '{"' else macros.get(v.lower(), v)
            f[fm.group(1).lower()] = v
        entries[m.group(2)] = f
    return entries

def tex_clean(s):
    s = re.sub(r'\\[\'"`^~=.]\{?(\w)\}?', r'\1', s)  # accents
    s = re.sub(r'\\(?:emph|textit|textbf|it|bf|em)\s*', '', s).replace('~', ' ').replace('\\ ', ' ').replace('--', '–')
    s = re.sub(r'\\&', '&', s); s = re.sub(r'\$([^$]*)\$', r'\1', s)
    return re.sub(r'\s+', ' ', s.replace('{', '').replace('}', '').replace('\\', '')).strip()

def fmt_authors(a):
    names = [tex_clean(x) for x in re.split(r'\s+and\s+', a) if x.strip()]
    def last(n): return n.split(',')[0].strip() if ',' in n else n.split()[-1]
    if not names: return '?'
    return last(names[0]) if len(names) == 1 else (f'{last(names[0])} & {last(names[1])}' if len(names) == 2 else f'{last(names[0])} et al.')

def fmt_ref(f):
    s = f"{fmt_authors(f.get('author', f.get('editor', '?')))} ({f.get('year', '?')}). {tex_clean(f.get('title', '?'))}."
    venue = f.get('journal') or f.get('booktitle') or f.get('publisher') or f.get('institution') or f.get('school') or ''
    if venue: s += f" {tex_clean(venue)}"
    if f.get('volume'): s += f" {f['volume']}"
    if f.get('pages'): s += f", {tex_clean(f['pages'])}"
    if f.get('doi'): s += f". doi:{re.sub(r'^https?://(dx\.)?doi\.org/', '', f['doi'].strip())}"
    elif f.get('url'): s += f". {f['url'].strip()}"
    return s

def bib_cites(docroot, files, label):
    """{key: [(label, relpath, line, section, varlink_on_line)]}"""
    out = collections.defaultdict(list)
    for rel in files:
        text = rd(os.path.join(docroot, rel)); lines = text.splitlines()
        secs = rst_sections(text)
        for i, l in enumerate(lines, 1):
            for keys in re.findall(r':cite(?:[tp]|alp|alt)?:`([^`]+)`', l):
                sec = ''
                for ln, lev, title, _ in secs:
                    if ln <= i: sec = title
                    else: break
                v = re.findall(r':varlink:`~?(\w+)`', l)
                for k in keys.split(','): out[k.strip()].append((label, rel, i, sec, v[0] if v else ''))
    return out

bibs, cites = {}, collections.defaultdict(list)
rst_rel = lambda root: sorted(os.path.relpath(os.path.join(dp, f), root) for dp, dn, fn in os.walk(root) for f in fn
                              if f.endswith('.rst') and '/old_doc' not in dp)
mdoc = os.path.join(ROOT, 'doc')
if os.path.exists(os.path.join(mdoc, 'manual_references.bib')):
    bibs.update(parse_bib(rd(os.path.join(mdoc, 'manual_references.bib'))))
    for k, v in bib_cites(mdoc, rst_rel(mdoc), 'mitgcm').items(): cites[k] += v
if A.darwin3 and os.path.exists(os.path.join(A.darwin3, 'doc', 'manual_references.bib')):
    d3doc = os.path.join(os.path.abspath(A.darwin3), 'doc')
    for k, v in parse_bib(rd(os.path.join(d3doc, 'manual_references.bib'))).items(): bibs.setdefault(k, v)
    d3files = [r for r in rst_rel(d3doc) if re.search(r'darwin|radtrans', r)]
    for k, v in bib_cites(d3doc, d3files, 'darwin3').items(): cites[k] += v
if bibs:
    B = ['# Bibliography: papers behind MITgcm / darwin3 schemes', '',
         'Generated from the manuals\' `manual_references.bib` (MITgcm, plus darwin3\'s extra ecosystem entries) and',
         'every `:cite:` in the manual. **Part 1** lists, per manual page, which papers it cites, where (`L<line>`), under which',
         'section, and the parameter on that line. **Part 2** gives the full references. `grep -n <param or author> bibliography.md`.',
         'Papers the manuals don\'t cite (ECCO, ECCO-Darwin) are in `references/literature.md`.', '', '## Part 1: citations by manual page', '']
    byfile = collections.defaultdict(list)
    for k, occ in cites.items():
        for lab_, rel, i, sec, v in occ: byfile[(lab_, rel)].append((i, k, sec, v))
    for (lab_, rel) in sorted(byfile):
        B.append(f"### {'darwin3 ' if lab_ == 'darwin3' else ''}`doc/{rel}`")
        for i, k, sec, v in sorted(byfile[(lab_, rel)]):
            f = bibs.get(k)
            who = f"{fmt_authors(f.get('author', f.get('editor', '?')))} {f.get('year', '?')}" if f else '(not in bib)'
            B.append(f"- L{i} `{k}` {who} — {sec}" + (f" — `{v}`" if v else ''))
        B.append('')
    B += ['## Part 2: full references (`key`: citation; cited N times)', '']
    for k in sorted(bibs, key=str.lower):
        B.append(f"- `{k}`: {fmt_ref(bibs[k])}" + (f" ({len(cites[k])}×)" if cites.get(k) else ' (not cited in manual text)'))
    wr('bibliography.md', '\n'.join(B))

# ---- packages
catalog = []
for root, pkg, label, page, diffnote, shown in jobs:
    rdefs, rcalls = scans[root]
    pd = os.path.join(root, 'pkg', pkg)
    files = sorted(os.listdir(pd))
    srcs = [os.path.join(pd, f) for f in files if f.endswith(SRC_EXT)]
    if not srcs: continue
    vdesc = var_descriptions(srcs)
    optf = [f for f in files if f.upper() == pkg.upper() + '_OPTIONS.H'] or [f for f in files if f.endswith('_OPTIONS.h')]
    blurb = ''
    for cand in optf + [f for f in files if f.lower() == pkg + '.h'] + [f for f in files if 'readparms' in f.lower()] + files:
        if cand.endswith(SRC_EXT):
            blurb = header_blurb(os.path.join(pd, cand))
            if blurb and 'CPP options file' not in blurb and len(blurb) > 25: break
    readme = next((rd(os.path.join(pd, f)) for f in files if f.lower().startswith('readme')), '')
    nls, dfiles = collections.OrderedDict(), []
    own = 'ifdef ALLOW_' + pkg.upper()
    for s in srcs:
        g, df = namelists(s)
        for k, v in g.items(): nls.setdefault(k, []).extend(v)
        dfiles += df
    opts = []
    for f in optf: opts += [(f,) + o for o in cpp_options(os.path.join(pd, f), pkg)]
    mydefs = {n: f for n, (f, o) in rdefs.items() if o == pkg and f.startswith('pkg/' + pkg + '/')}
    ext_calls = [(f, n, c, g) for (f, n, c, g) in rcalls if c in mydefs and not f.startswith('pkg/' + pkg + '/')]
    if label == 'MITgcm':
        used = sorted(e for e in pkg_used_by.get(pkg, []) if '@' not in e)
    else:
        used = sorted(e for e in pkg_used_by.get(pkg, []) if e.endswith('@' + label))
    docrefs = [r for r, t, txt in docs if re.search(r'pkg/' + re.escape(pkg) + r'\b', txt)]
    docrefs.sort(key=lambda r: (pkg not in os.path.basename(r), r))
    has_ad = [f for f in files if f.endswith('_diff.list') or 'ad_check' in f]
    blurb = PKG_DESC.get(pkg, blurb) if label in ('MITgcm', 'darwin3') else (blurb or PKG_DESC.get(pkg, ''))
    catalog.append((page, label, blurb, sorted(set(dfiles)), sum(len(v) for v in nls.values()), len(opts), bool(has_ad), used))

    L = [f'# pkg/{pkg}' + (f'  ({label}: {shown})' if label != 'MITgcm' else ''), '']
    if blurb: L += [blurb, '']
    if diffnote: L += [f'**vs its upstream base (merge-base, see README):** {diffnote}', '']
    if depend.get(pkg): L.append('**pkg_depend:** ' + ' '.join(depend[pkg]) + '  (`+` requires, `-` excludes)')
    ingroups = [g for g, v in groups.items() if pkg in v]
    if ingroups: L.append('**in groups:** ' + ', '.join(ingroups))
    L.append('**runtime switch:** `use' + pkg.upper() + '`-style flag in `data.pkg` (check exact name in packages_boot.F)' if pkg not in ('mdsio', 'rw', 'mom_common', 'generic_advdiff', 'debug', 'monitor') else '**always-on utility package** (no data.pkg switch)')
    if dfiles: L.append('**reads:** ' + ', '.join(f'`{d}`' for d in sorted(set(dfiles))))
    if docrefs: L.append('**manual:** ' + ', '.join(f'`{d}`' for d in docrefs[:5]))
    if has_ad: L.append('**adjoint support files:** ' + ', '.join(has_ad))
    if readme: L += ['', '## README', '```'] + readme.splitlines()[:30] + ['```']
    if nls:
        L += ['', '## Namelist parameters']
        for g, ps in nls.items():
            L.append(f'### {g}')
            seen = set()
            for p, guard in ps:
                if p.lower() in seen: continue
                seen.add(p.lower())
                d = vdesc.get(p.lower(), '')
                guard = ' & '.join(x for x in guard.split(' & ') if x != own)
                L.append(f'- `{p}`' + (f' — {d}' if d else '') + (f'  _[{guard}]_' if guard else ''))
    if opts:
        L += ['', '## CPP options (defaults as shipped)']
        for f, n, kind, c in opts:
            L.append(f'- `{n}` ({kind}, {f})' + (f' — {c}' if c else ''))
    hdrs = [f for f in files if f.endswith('.h')]
    if hdrs:
        L += ['', '## Headers'] + [f'- `{h}` — {header_blurb(os.path.join(pd, h))[:150]}' for h in hdrs]
    if mydefs:
        byfile = collections.defaultdict(list)
        for n, f in sorted(mydefs.items()): byfile[os.path.basename(f)].append(n)
        L += ['', f'## Routines ({len(mydefs)})',
              ', '.join(f'`{os.path.basename(f)}`' for f in sorted(byfile))]
    if ext_calls:
        L += ['', '## Called from outside the package']
        agg = collections.defaultdict(list)
        for f, n, c, g in ext_calls: agg[(c, f)].append(n)
        for (c, f), ns in sorted(agg.items(), key=lambda x: (x[0][1], x[0][0])):
            L.append(f'- `{c}` ← `{f}:' + ','.join(map(str, ns[:4])) + '`')
    if used:
        L += ['', f'## Verification experiments compiling it ({len(used)})', ' '.join(f'`{u}`' for u in used)]
    wr(f'pkg/{page}.md', '\n'.join(L))

C = ['# MITgcm packages', '',
     'One row per `pkg/` directory. Details (namelists with descriptions, CPP options, call sites, '
     'experiments): `pkg/<name>.md`. "#nml" = namelist parameters, "#cpp" = CPP options in its *_OPTIONS.h, '
     'AD = ships adjoint diff lists. Experiments = verification experiments whose packages.conf (groups '
     'expanded) includes it. Rows named `pkg@label` come from the user\'s branch clones and are new or '
     'differ from upstream (see README.md for label → path).', '',
     '| pkg | data file | #nml | #cpp | AD | experiments | description |', '|---|---|---|---|---|---|---|']
for page, label, blurb, dfs, nn, no, ad, used in sorted(catalog, key=lambda r: (r[1] not in ('MITgcm', 'darwin3'), r[0])):
    nm = page + (' (darwin3)' if label == 'darwin3' else '')
    C.append(f"| `{nm}` | {' '.join(dfs)} | {nn} | {no} | {'y' if ad else ''} | {len(used)}"
             f"{': ' + ', '.join(used[:4]) + (', …' if len(used) > 4 else '') if used else ''} | {blurb[:130].replace('|', '/')} |")
if extra_rows:
    C += ['', '## Core (model/) changes in the user\'s branches', '']
    for label, root, sub, d in extra_rows:
        C.append(f'- **{label}** `{sub}`: {d}')
wr('packages.md', '\n'.join(C))


# ---- core params
core = collections.OrderedDict()
cdesc = var_descriptions(list(fortran_files(os.path.join(ROOT, 'model/inc'))) + list(fortran_files(os.path.join(ROOT, 'eesupp/inc'))))
for f in ('model/src/ini_parms.F', 'model/src/packages_boot.F', 'eesupp/src/eeset_parms.F'):
    g, df = namelists(os.path.join(ROOT, f))
    core[f] = (g, df)
L = ['# Core namelists', '', 'Descriptions come from `::` comments in model/inc and eesupp/inc headers '
     '(PARAMS.h etc.). Guards in _[italics]_ are the CPP conditions around the namelist entry.', '']
for f, (g, df) in core.items():
    L.append(f'## {f}  (reads {", ".join(sorted(set(df))) or "?"})')
    for grp, ps in g.items():
        L.append(f'### {grp}'); seen = set()
        for p, guard in ps:
            if p.lower() in seen: continue
            seen.add(p.lower()); d = cdesc.get(p.lower(), '')
            L.append(f'- `{p}`' + (f' — {d}' if d else '') + (f'  _[{guard}]_' if guard else ''))
wr('core-params.md', '\n'.join(L))

# ---- call tree
by_caller = collections.defaultdict(list)
caller_routine = {}
# map (file) -> routine name(s): use defs reverse map
file_routines = collections.defaultdict(list)
for n, (f, o) in defs.items(): file_routines[f].append(n)
# assign calls to the routine they're in by re-scanning model/src files in order
def calls_in_routine(name):
    if name not in defs: return []
    f = defs[name][0]; out = []; cur = None; stack = []
    for l in rd(os.path.join(ROOT, f)).splitlines():
        if l.startswith('#'):
            d = l[1:].strip()
            if re.match(r'if(n?def)?\b', d): stack.append(d)
            elif d.startswith('else') and stack: stack[-1] = 'NOT(' + stack[-1] + ')'
            elif d.startswith('endif') and stack: stack.pop()
            continue
        if is_comment(l): continue
        code = strip_inline(l); m = SUB_RE.match(code)
        if m and not code.strip().upper().startswith('END'): cur = m.group(2).upper(); continue
        if cur == name:
            for c in CALL_RE.findall(code):
                g = [s.replace('ifdef ', '').replace('ifndef ', '!') for s in stack if 'ALLOW_' in s or 'ifn' in s]
                out.append((c.upper(), g))
    return out

SKIP = {'TIMER_START', 'TIMER_STOP', 'PRINT_MESSAGE', 'PRINT_ERROR', 'DEBUG_ENTER', 'DEBUG_LEAVE',
        'DEBUG_CALL', 'DEBUG_STATS_RL', 'BARRIER', 'MONITOR', 'DEBUG_MSG', 'WRITE_FULLARRAY_RL',
        'EXCH_UV_XYZ_RL', 'EXCH_XYZ_RL', 'EXCH_XY_RL', 'EXCH_UV_XY_RL', 'COMP_FLDS', 'DIAGNOSTICS_FILL',
        'WRITE_FLD_XYZ_RL', 'WRITE_FLD_XY_RL', 'WRITE_FLD_XY_RS', 'WRITE_FLD_XYZ_RS', 'GLOBAL_SUM_TILE_RL',
        'GLOBAL_MAX_R8', 'GLOBAL_SUM_R8', 'PLOT_FIELD_XYRL', 'PLOT_FIELD_XYZRL', 'PLOT_FIELD_XYRS',
        'DEBUG_FLD_STATS_RL', 'DIAGNOSTICS_SCALE_FILL', 'DIAGNOSTICS_FRACT_FILL', 'MDS_WRITELOCAL',
        'ALL_PROC_DIE', 'DIAGNOSTICS_IS_ON', 'STOP', 'PACKAGES_PRINT_MSG', 'PACKAGES_ERROR_MSG', 'LCASE',
        'OPEN_COPY_DATA_FILE', 'MDSFINDUNIT', 'WRITE_0D_RL', 'WRITE_0D_I', 'WRITE_0D_L', 'WRITE_0D_C',
        'WRITE_1D_RL', 'WRITE_1D_I', 'WRITE_XY_XLINE_RS', 'WRITE_XY_YLINE_RS', 'WRITE_COPY1D_RS',
        'EXCH_XY_RS', 'EXCH_XYZ_RS', 'EXCH_UV_AGRID_3D_RL', 'EXCH_UV_DGRID_3D_RL', 'EXCH_UV_XY_RS',
        'EXCH_SM_3D_RL', 'EXCH_3D_RL', 'EXCH_3D_RS', 'FLUSH', 'TIMER_PRINTALL', 'GLOBAL_SUM_TILE_RS',
        'GLOBAL_MAX_R4', 'GLOBAL_MIN_R8', 'MDS_FLUSH'}
tree = ['# Core call tree (static)', '',
        'Static CALL graph from `THE_MODEL_MAIN`, expanded through routines in model/src (depth ≤ 7). '
        'Package routines are leaves tagged `[pkg/x]`. `{ALLOW_X}` = enclosing #ifdef. Timers, '
        'exchanges, I/O and debug calls are omitted. A routine expanded earlier is marked "(↑)". '
        'Runtime IF tests (e.g. `IF (useKPP)`) are not shown — read the source for those.', '', '```']
expanded = set()
def walk(name, depth, guards):
    shown = set()
    for c, g in calls_in_routine(name):
        if c in SKIP or c.startswith('DEBUG_') or c.startswith('PLOT_') or c in shown: continue
        shown.add(c)
        owner = defs.get(c, ('?', '?'))
        tag = (' [pkg/' + owner[1] + ']') if owner[0].startswith('pkg/') else ('' if owner[0] != '?' else ' [?]')
        gs = (' {' + ','.join(g) + '}') if g else ''
        is_core = owner[0].startswith('model/')
        mark = ' (↑)' if (is_core and c in expanded) else ''
        tree.append('  ' * depth + c + tag + gs + mark)
        if is_core and c not in expanded and depth < 7:
            expanded.add(c); walk(c, depth + 1, g)
expanded.add('THE_MODEL_MAIN'); tree.append('THE_MODEL_MAIN')
walk('THE_MODEL_MAIN', 1, [])
tree.append('```')
wr('core-calltree.md', '\n'.join(tree))

# ---- tools help
def helptext(cmd, cwd):
    try:
        r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=60)
        return (r.stdout + r.stderr).strip()
    except Exception as ex: return f'(could not run: {ex})'
T = ['# tools: genmake2 and testreport', '', '## genmake2 -help', '```',
     helptext(['sh', 'tools/genmake2', '-help'], ROOT)[:12000], '```', '',
     '## testreport -help', '```', helptext(['sh', './testreport', '-help'], VER)[:8000], '```', '',
     '## Optfiles shipped (tools/build_options)', ' '.join(f'`{f}`' for f in sorted(os.listdir(os.path.join(ROOT, 'tools/build_options'))))]
wr('tools.md', '\n'.join(T))

# ---- external docs (ECCO v4 docs, ecco_darwin readmes, darwin3 manual pages)
def md_sections(text):
    out, lines, fence = [], text.splitlines(), False
    heads = []
    for k, l in enumerate(lines):
        if l.startswith('```'): fence = not fence
        m = None if fence else re.match(r'^(#{1,4})\s+(.+?)\s*#*\s*$', l)
        if m: heads.append((k + 1, len(m.group(1)) - 1, m.group(2)))
    for n, (ln, lev, title) in enumerate(heads):
        body = '\n'.join(lines[ln:heads[n + 1][0] - 1] if n + 1 < len(heads) else lines[ln:])
        names = [w for w in re.findall(r'`([A-Za-z_][\w.]{3,})`', body) if '_' in w or re.search(r'[a-z][A-Z]', w) or w.endswith(('.F', '.h'))]
        out.append((ln, lev, title, list(dict.fromkeys(names))))
    return out

def readme_facts(text):
    """first informative line + MITgcm version / optfile / key settings a plain-text readme mentions."""
    lines = [l.strip(' #*=-~\t') for l in text.splitlines()]
    lines = [l for l in lines if len(re.findall(r'[A-Za-z]', l)) >= 3]
    first = re.sub(r'\s+', ' ', lines[0])[:160] if lines else ''
    facts = []
    for pat in (r'checkpoint\d+[a-z]?', r'darwin3[\w./-]*', r'\b[\w-]+_(?:gfortran|ifort|intel|mpi)[\w+.-]*', r'\bv4r\d\b', r'\b(?:llc|cs)\d+\b', r'nTimeSteps\s*=\s*\d+',
                r'git (?:checkout|clone)\s+\S+(?:\s+\S+)?', r'cvs (?:co|checkout)\s+[^\n]{0,60}'):
        for m in re.findall(pat, text, re.I):
            m = re.sub(r'\s+', ' ', m.strip())
            if m not in facts: facts.append(m)
    return first, facts[:14]

if A.docroot:
    X = ['# External docs map (ECCO v4, ecco_darwin, darwin3 manual)', '',
         'Section headings with line numbers (`L123`) and the names each section cites, for docs outside MITgcm',
         'proper. Plain-text readmes show their first line and the MITgcm checkpoint / optfile / grid / git',
         'commands they mention. Read the real file at the path shown (roots listed below).', '']
    for spec in A.docroot:
        lab_, _, rest = spec.partition('=')
        root_, _, filt = rest.partition(':')
        root_ = os.path.expanduser(root_); subs = [s for s in filt.split(',') if s]
        X += ['', f'# @{lab_} = `{root_}` ({git_label(root_)})']
        for dp, dn, fn in os.walk(root_):
            # skip hidden/build dirs and vendored MITgcm/darwin3 trees (their manual is mapped in docs.md already)
            dn[:] = sorted(d for d in dn if not d.startswith(('.', '_build')) and d not in ('figs', '_static', 'images', 'node_modules')
                           and not os.path.isdir(os.path.join(dp, d, 'eesupp')) and not os.path.isdir(os.path.join(dp, d, 'phys_pkgs')))
            for f in sorted(fn):
                p = os.path.join(dp, f); rel = os.path.relpath(p, root_)
                if subs and not any(s in rel for s in subs): continue
                low = f.lower(); text = None
                if low.endswith('.rst'): text = rd(p); secs = rst_sections(text)
                elif low.endswith('.md'): text = rd(p); secs = md_sections(text)
                elif low.startswith('readme') or (low.endswith('.txt') and 'readme' in low): text = rd(p); secs = None
                if text is None or not text.strip(): continue
                if secs is None:
                    first, facts = readme_facts(text)
                    X.append(f'- `{rel}` ({len(text.splitlines())} lines) — {first}' + (' — ' + ', '.join(f'`{x}`' for x in facts) if facts else ''))
                    continue
                X.append(f'\n## `{rel}`')
                for ln, lev, title, names in secs:
                    if lev > 3: continue
                    s = f"{'  ' * lev}- L{ln} {title}"
                    if names: s += ' — ' + ', '.join(f'`{w}`' for w in names[:30]) + (f' (+{len(names) - 30})' if len(names) > 30 else '')
                    X.append(s)
    wr('docs-ext.md', '\n'.join(X))

# ---- README
lab = A.label or git_label(ROOT)
d3lab = (' + darwin3 pkg/darwin, pkg/radtrans: ' + git_label(A.darwin3)) if A.darwin3 else ''
def git_branch(r):
    try: return subprocess.run(['git', '-C', r, 'branch', '--show-current'], capture_output=True, text=True).stdout.strip()
    except Exception: return ''
xlab = ''.join(f"\n- `@{l}` = `{r}` (branch `{git_branch(r)}`, {git_label(r)}; working tree incl. uncommitted edits; diffed against merge-base {BASES.get(l, '')})"
               + (f" — only {', '.join(o)}" if o else '') for l, r, o in EXTRAS)
wr('README.md', f"""# MITgcm source index

Indexed: MITgcm {lab}{d3lab}
Counts: {len(catalog)} package pages, {len(exps)} upstream + {len(extra_exps)} branch verification experiments, {len(docs)} manual pages.\n\nUser branch clones indexed (only packages/experiments that are new or differ from upstream):{xlab or " none"}

Files: `packages.md`, `pkg/<name>.md`, `verification.md`, `verification/<exp>.md`, `core-params.md`,
`core-calltree.md`, `docs.md`, `bibliography.md`, `tools.md`{', `docs-ext.md` (' + ', '.join('@' + s.split('=')[0] for s in A.docroot) + ')' if A.docroot else ''}.

Generated by `scripts/build_mitgcm_index.py`; heuristic parsing, so treat it as a map and confirm
in the source (MITgcm: `git -C ~/Documents/research/ECCO/BBL/MITgcm show origin/master:<path>`,
darwin3: `~/Documents/GitHub/darwin3`).

Regenerate:
  git -C <MITgcm clone> fetch origin && mkdir -p /tmp/m && git -C <clone> archive origin/master | tar -x -C /tmp/m
  python3 ~/.claude/skills/mitgcm-ecco/scripts/build_mitgcm_index.py /tmp/m ~/.claude/skills/mitgcm-ecco/references/mitgcm-index \\
      --darwin3 ~/Documents/GitHub/darwin3 --label "<commit / tag>" [--extra LABEL=ROOT[:pkgs] ...]
  (or just run scripts/refresh_index.sh, which does all of the above for the standard clones)
""")
print(f'indexed {len(catalog)} packages, {len(exps)} experiments, {len(docs)} docs -> {OUT}')
