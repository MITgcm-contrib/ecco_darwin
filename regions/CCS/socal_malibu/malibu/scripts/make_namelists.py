#!/usr/bin/env python3
"""Write the grid-dependent malibu_v01 namelists into ../input/.

Reads ../work/grid_info.json and ../work/obcs_files.json (from make_grid_bathy.py
and make_ic_obcs.py), the parent delYFile, and the CCS v05/v05_kelp inputs as
templates. Writes: data, data.obcs, data.exf, data.ptracers, data.darwin.
Static files (data.pkg, data.radtrans, data.sediment, data.diagnostics, ...)
live in ../input/ already.

All atmospheric, wave, river and iron forcing is read from the CCS parent's
80x100 files and interpolated in the model (exf / darwin interp, method 1);
nothing is regridded offline.

Usage: python3 make_namelists.py <ecco_darwin>/regions/CCS <parent run_template dir>
       [--case v5r1 --box <extract_v5r1_box.py output> --start YYYY-MM-DD
        --days N --work <make_ic_obcs.py --out dir> --out <input dir>]
--case v5r1 (PACE years): ICs/OBCS from the ECCO-Darwin v05 LLC270 V5r1 run
(make_ic_obcs.py --mode v5r1), atmosphere and iron from the same run's box
(v5r1/ in the run dir, in-model interpolation from the 1/3-deg box), OASIM
monthly light (oasim/, Oliver Jahn's C61extrap below-surface files), no waves.
"""
import argparse
import json
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, '..', 'work')
INP = os.path.join(HERE, '..', 'input')

DELTAT = 90.0
NTIMESTEPS = 40          # short test (1 h); set endTime for real runs
START = 19920116         # data.cal start date (parent IC date)

# parent CCS v05 grid, forcing at cell centres
PNX, PNY, PXG0, PDX, PYG0 = 80, 100, -131.66667, 1.0 / 3.0, 19.8998


def wrap(prefix, vals, fmt='{:.6f}', per=6):
    out, line = [], []
    for v in vals:
        line.append(fmt.format(v))
        if len(line) == per:
            out.append(', '.join(line))
            line = []
    if line:
        out.append(', '.join(line))
    return prefix + (',\n' + ' ' * len(prefix)).join(out) + ','


PGRID = {}   # forcing source grid: lon0, dlon, nlon, lat (cell centres)


def parent_interp(prefix, dely=None):
    if not PGRID:     # CCS v05 80x100 parent
        yc = PYG0 + np.r_[0.0, np.cumsum(dely)][:-1] + 0.5 * dely
        PGRID.update(lon0=PXG0 + 0.5 * PDX, dlon=PDX, nlon=PNX, lat=yc)
    yc = np.asarray(PGRID['lat'])
    return '\n'.join([
        f'  {prefix}_lon0 = {PGRID["lon0"]:.6f},',
        f'  {prefix}_lon_inc = {PGRID["dlon"]:.8f},',
        f'  {prefix}_lat0 = {yc[0]:.6f},',
        wrap(f'  {prefix}_lat_inc = ', np.diff(yc)),
        f'  {prefix}_nlon = {PGRID["nlon"]},',
        f'  {prefix}_nlat = {yc.size},',
        f'  {prefix}_interpMethod = 1,'])


def sub_once(s, old, new):
    if s.count(old) != 1:
        sys.exit(f'expected one match of:\n{old}')
    return s.replace(old, new)


def main():
    global INP, START, NTIMESTEPS
    ap = argparse.ArgumentParser()
    ap.add_argument('ccs')
    ap.add_argument('ptmpl')
    ap.add_argument('--case', choices=('ccs1992', 'v5r1'), default='ccs1992')
    ap.add_argument('--box', default=None)
    ap.add_argument('--start', default=None)
    ap.add_argument('--days', type=float, default=None)
    ap.add_argument('--work', default=WORK)
    ap.add_argument('--out', default=INP)
    ap.add_argument('--waves', default=None,
                    help='v5r1: dir from make_waves.py (ERA5): pkg/sediment sedWaveMode = 3 and '
                         'EXF WvHeight/WvPeriod for kelp wave mortality; linked as waves/ in the run dir')
    ap.add_argument('--sed-bed', action='store_true',
                    help='v5r1: initial erodible bed sediment/sed_bed_<nx>x<ny> (make_sed_bed.py, '
                         'linked from <work>/sediment) w/ literature erosion rate and bulk density')
    ap.add_argument('--sed-layered', action='store_true',
                    help='with --sed-bed: consolidating layered bed (SEDIMENT_PARM02)')
    ap.add_argument('--diag-pertype', action='store_true',
                    help='v5r1: daily per-plankton PP/GR + c01-c07 diagnostics, top 20 levels '
                         '(needs a build with DARWIN_DIAG_PERTYPE, BUILD=build_diag)')
    ap.add_argument('--kelp-death', type=float, default=None,
                    help='data.darwin mp_spp_death for kelp [1/d] (default template 1.15e-7 1/s = 0.0099/d)')
    ap.add_argument('--diag-budget', action='store_true',
                    help='v5r1: daily advective/diffusive fluxes of c03, c05 (TRAC23, TRAC25), top 20 levels')
    ap.add_argument('--chl2cmax-small', type=float, default=None,
                    help='data.traits CHL2CMAX for phytoplankton classes 3-5 [mg Chl/mmol C] (default 0.2)')
    ap.add_argument('--phyto-mort', type=float, default=None,
                    help='data.traits MORT for phytoplankton classes 1-5 [1/d] (default 0.1)')
    ap.add_argument('--no-waves', action='store_true',
                    help='drop the CCS_wave WvHeight/WvPeriod forcing (files not on the Mac)')
    a = ap.parse_args()
    ccs, ptmpl = a.ccs, a.ptmpl
    v5 = a.case == 'v5r1'
    INP = a.out
    os.makedirs(INP, exist_ok=True)
    if v5:
        box = json.load(open(os.path.join(a.box, 'box.json')))
        lon = np.array(box['lon'])
        PGRID.update(lon0=lon[0], dlon=float(np.mean(np.diff(lon))), nlon=lon.size, lat=box['lat'])
        START = int(a.start.replace('-', ''))
        NTIMESTEPS = int(round(a.days * 86400 / DELTAT))
    g = json.load(open(os.path.join(WORK, 'grid_info.json')))
    ob = json.load(open(os.path.join(a.work, 'obcs_files.json')))
    dely = np.fromfile(os.path.join(ptmpl, 'delYFile'), '>f4').astype(float)
    nx, ny, nr, tag = g['nx'], g['ny'], g['nr'], ob['tag']
    delr = np.array(g['delr'])

    # tRef/sRef: horizontal mean of the initial state on each level
    T = np.fromfile(os.path.join(a.work, 'init', f'THETA_{tag}'), '>f4').reshape(nr, ny, nx)
    S = np.fromfile(os.path.join(a.work, 'init', f'SALT_{tag}'), '>f4').reshape(nr, ny, nx)
    tref, sref = T.mean(axis=(1, 2)), S.mean(axis=(1, 2))

    # ---------------- data
    d = f"""# ====================
# | Model parameters |  malibu_v01: Malibu / Pt. Dume kelp nest, ~460 m
# ====================
 &PARM01
{wrap(' tRef = ', tref, '{:.2f}', 10)}
{wrap(' sRef = ', sref, '{:.3f}', 10)}
 no_slip_sides  = .TRUE.,
 no_slip_bottom = .TRUE.,
 viscAr=0.5E-4,
 viscAh=1.E0,
 viscAhGrid=1.E-2,
 viscC4Leith=1.5,
 viscC4Leithd=1.5,
 viscA4GridMax=0.5,
 useAreaViscLength=.TRUE.,
 diffKhT=1.E0,
 diffKhS=1.E0,
 diffKrNrT={nr}*1.E-5,
 diffKrNrS={nr}*1.E-5,
 bottomDragQuadratic = 0.0025,
 ivdc_kappa=10.,
 implicitDiffusion=.TRUE.,
 implicitViscosity=.TRUE.,
 useRealFreshWaterFlux=.TRUE.,
# 2-m levels near the surface: hFacMinDr in m, not the 50 m of CCS v05
 hFacMinDr=1.,
 hFacMin=0.3,
 hFacInf=0.1,
 hFacSup=5.,
 select_rStar=2,
 nonlinFreeSurf=4,
 gravity=9.81,
 rhonil=1029.,
 rhoConst=1029.,
 rhoConstFresh=1000.,
 convertFW2Salt=-1.,
 eosType='JMD95Z',
 implicitFreeSurface=.TRUE.,
 exactConserv=.TRUE.,
 useSingleCpuIO=.TRUE.,
 tempAdvScheme=30,
 saltAdvScheme=30,
 tempVertAdvScheme=3,
 saltVertAdvScheme=3,
 tempImplVertAdv=.TRUE.,
 saltImplVertAdv=.TRUE.,
 staggerTimeStep=.TRUE.,
 vectorInvariantMomentum=.TRUE.,
 useJamartWetPoints=.TRUE.,
 readBinaryPrec=32,
 writeBinaryPrec=32,
 debugLevel=-1,
 /

 &PARM02
 cg2dMaxIters=300,
 cg2dTargetResWunit=1.E-12,
 /

 &PARM03
 nIter0=0,
 nTimeSteps={NTIMESTEPS},
# endTime=31536000.,
 deltaT={DELTAT},
 forcing_In_AB=.FALSE.,
 momDissip_In_AB=.FALSE.,
 doAB_onGtGs=.FALSE.,
 alph_AB=0.5,
 beta_AB=0.281105,
 pChkptFreq=2592000.,
 monitorFreq=900.,
 dumpInitAndLast=.TRUE.,
 pickupStrictlyMatch=.FALSE.,
 /

 &PARM04
 usingSphericalPolarGrid=.TRUE.,
 xgOrigin = {g['lon0']},
 ygOrigin = {g['lat0']},
 delX = {nx}*{g['dlon']},
 delYfile = '{g['delyfile']}',
{wrap(' delR = ', delr, '{:.2f}', 10)}
 /

 &PARM05
 bathyFile       = '{g['bathy']}',
 hydrogThetaFile = 'init/THETA_{tag}',
 hydrogSaltFile  = 'init/SALT_{tag}',
 uVelInitFile    = 'init/U_{tag}',
 vVelInitFile    = 'init/V_{tag}',
 pSurfInitFile   = 'init/ETAN_{nx}x{ny}',
 /
"""
    if v5:
        # flux-limited DST3 (33) for T and S, as for the ptracers: unlimited DST3 (30, from CCS v05)
        # overshot S by ~1 psu at the Malibu Creek plume front in the PACE-year runs (2026-10-07)
        assert d.count('AdvScheme=30,') == 2
        d = d.replace(' tempAdvScheme=30,', '# flux-limited DST3 for T, S (unlimited 30 overshot at the creek plume front)\n'
                      ' tempAdvScheme=33,').replace(' saltAdvScheme=30,', ' saltAdvScheme=33,')
    open(os.path.join(INP, 'data'), 'w').write(d)

    # ---------------- data.obcs
    sides = ob['sides']
    cap = {'west': 'W', 'south': 'S', 'east': 'E'}
    lines = [' &OBCS_PARM01']
    if 'west' in sides:
        lines.append(f' OB_Iwest = {ny}*1,')
    if 'east' in sides:
        lines.append(f' OB_Ieast = {ny}*-1,')
    if 'south' in sides:
        lines.append(f' OB_Jsouth = {nx}*1,')
    lines += [' useOBCSprescribe = .TRUE.,', ' useOBCSsponge = .TRUE.,',
              ' useOBCSbalance = .TRUE.,', '#']
    fmap = {'THETA': 't', 'SALT': 's', 'U': 'u', 'V': 'v'}
    for side in sides:
        c = cap[side]
        for f, k in fmap.items():
            lines.append(f" OB{c}{k}File = 'obcs/{ob['files'][f][side]}',")
        for n, tr in enumerate(ob['tracers'], 1):
            if tr in ob['files']:
                lines.append(f" OB{c}ptrFile({n}) = 'obcs/{ob['files'][tr][side]}',")
        lines.append('#')
    lines += [' /', '', '# sponge: 8 cells (~3.7 km), 1 h at the boundary, 1 d inside',
              ' &OBCS_PARM03', ' spongeThickness = 8,',
              ' Urelaxobcsinner = 86400.,', ' Urelaxobcsbound = 3600.,',
              ' Vrelaxobcsinner = 86400.,', ' Vrelaxobcsbound = 3600.,', ' /', '']
    open(os.path.join(INP, 'data.obcs'), 'w').write('\n'.join(lines))

    # ---------------- data.exf (from v05_kelp; parent-grid interpolation)
    exf = open(os.path.join(ccs, 'v05_kelp', 'input', 'data.exf')).read()
    exf = exf.replace("  runofffile        = 'era_xx_it42_v2/CCS_kelp_80x100_Fekete_runoff',\n", '')
    exf = exf.replace('runoffperiod      = 86400.0,,', 'runoffperiod      = 86400.0,')
    fields = ['atemp', 'aqh', 'precip', 'uwind', 'vwind', 'swdown', 'lwdown',
              'WvHeight', 'WvPeriod']
    if a.no_waves and not v5:
        exf = re.sub(r"  WvHeight.*\n|  WvPeriod.*\n", '', exf)
        fields = fields[:-2]
    if v5:
        # V5r1 atmosphere (6-hourly, ECCO V5r1 iter70) cut to the box; no
        # wave forcing for this case (sedWaveMode = 0, no wave mortality)
        exf = exf.replace("'era_xx_it42_v2/CCS_kelp_80x100_EXF", "'v5r1/EXF")
        exf = re.sub(r"  WvHeight.*\n|  WvPeriod.*\n", '', exf)
        fields = fields[:-2]
    # rivers (make_rivers.py): point sources on the nest grid replace the
    # parent's JRA55-do runoff, which smears river input over 1/3-deg cells
    exf = sub_once(exf, "  runofffile        = './JRA55-do/jra55_do_runoff_CCS_cutout',\n",
                   "# Malibu Creek point source on the nest grid (python/make_rivers.py)\n"
                   "  runofffile        = 'rivers/rivers_runoff_128x72',\n")
    nml04 = (' &EXF_NML_04\n' + '\n'.join(parent_interp(f, dely) for f in fields)
             + '\n  runoff_interpMethod = 0,\n /')
    exf = re.sub(r' &EXF_NML_04.*?\n /', lambda m: nml04, exf, flags=re.S)
    obcs = ' &EXF_NML_OBCS\n' + ''.join(
        f'  obcs{cap[s]}startdate1 = {START},\n  obcs{cap[s]}startdate2 = 000000,\n'
        f'  obcs{cap[s]}period = 2592000.0,\n' for s in sides) + ' /'
    exf = re.sub(r' &EXF_NML_OBCS.*?\n/', lambda m: obcs, exf, flags=re.S)
    if v5 and a.waves:
        wv = json.load(open(os.path.join(a.waves, 'waves.json')))
        lat_w = np.array(wv['lat'])
        add = ''
        for f in ('WvHeight', 'WvPeriod'):
            add += (f"  {f}file = 'waves/{wv.get('exf_prefix', 'era5')}_{f}',\n  {f}startdate1 = 19920101,\n"
                    f"  {f}startdate2 = 000000,\n  {f}period = 3600.0,\n")
        exf = exf.replace(' &EXF_NML_02\n', ' &EXF_NML_02\n# hourly waves (python/make_waves.py, ' + wv.get('source', 'era5') + ') for kelp wave mortality\n' + add, 1)
        interp = ''.join(f"  {f}_lon0 = {wv['lon0']:.4f},\n  {f}_lon_inc = {wv['dlon']:.4f},\n"
                         f"  {f}_lat0 = {lat_w[0]:.4f},\n  {f}_lat_inc = {lat_w.size - 1}*{np.diff(lat_w).mean():.4f},\n"
                         f"  {f}_nlon = {wv['nx']},\n  {f}_nlat = {wv['ny']},\n  {f}_interpMethod = 1,\n"
                         for f in ('WvHeight', 'WvPeriod'))
        exf = exf.replace('\n  runoff_interpMethod = 0,\n /', '\n' + interp + '  runoff_interpMethod = 0,\n /', 1)
        assert exf.count('WvHeightfile') == 1 and exf.count('WvHeight_lon0') == 1
    if v5:
        assert 'era_xx' not in exf, 'unreplaced CCS forcing in data.exf'
    open(os.path.join(INP, 'data.exf'), 'w').write(exf)

    # ---------------- data.ptracers (v05_kelp settings, nest tracer list)
    # Tracer order follows darwin DARWIN_INDICES.h with DARWIN_ALLOW_CDOM:
    # 1-19 DIC..O2, 20 CDOM, 21-27 c01-c07, 28-32 Chl01-05, 33-35 kelp, 36 TSS
    # (obcs_files.json from make_ic_obcs.py gives the same order)
    names = (['DIC', 'NO3', 'NO2', 'NH4', 'PO4', 'FeT', 'SiO2', 'DOC', 'DON', 'DOP', 'DOFe',
              'POC', 'PON', 'POP', 'POFe', 'POSi', 'PIC', 'ALK', 'O2', 'CDOM']
             + [f'c{n:02d}' for n in range(1, 8)] + [f'Chl{n:02d}' for n in range(1, 6)]
             + ['Bmag', 'Cmag', 'Qmag', 'TSS'])
    units = (['uM C', 'uM N', 'uM N', 'uM N', 'uM P', 'uM Fe', 'uM Si', 'uM C', 'uM N', 'uM P',
              'uM Fe', 'uM C', 'uM N', 'uM P', 'uM Fe', 'uM Si', 'uM C', 'meq/m^3', 'uM O', 'uM P']
             + ['uM C'] * 7 + ['mg/m^3'] * 5 + ['gDW/m^3', 'gC/m^3', 'mgN/gDryWeight', 'g/m^3'])
    nt = len(names)
    assert nt == len(ob['tracers']) == 36, (nt, len(ob['tracers']))
    ik = names.index('Bmag') + 1
    lines = [' &PTRACERS_PARM01', f' PTRACERS_numInUse= {nt},', '#',
             f' PTRACERS_advScheme(:)= {nt}*33,', f' PTRACERS_diffKh(:)= {nt}*0.E3,',
             f' PTRACERS_diffKr(:)= {nt}*1.E-5,', f' PTRACERS_useGMRedi(:)= {nt}*.FALSE.,',
             f'# A13: kelp is attached: no advection, diffusion or OBCS for tracers {ik}-{ik + 2}',
             f' PTRACERS_startStepFwd({ik}:{ik + 2}) = 3*1.E20,', '#',
             f' PTRACERS_EvPrRn(:)= {nt}*0.E3,', '#']
    lines += [f" PTRACERS_names({n:2d})= '{v}'," for n, v in enumerate(names, 1)]
    lines += ['#', '# uM is shorthand for mmol m^-3; CDOM is in P units (darwin)', '#']
    lines += [f" PTRACERS_units({n:2d})= '{v}'," for n, v in enumerate(units, 1)]
    lines += ['#'] + [f" PTRACERS_initialFile({n:2d})= 'init/{tr}_{tag}',"
                      for n, tr in enumerate(ob['tracers'], 1)]
    lines += [' &', '']
    open(os.path.join(INP, 'data.ptracers'), 'w').write('\n'.join(lines))

    # ---------------- data.darwin (from v05_kelp)
    dw = open(os.path.join(ccs, 'v05_kelp', 'input', 'data.darwin')).read()
    dw = dw.replace(' darwin_useSEAICE = .TRUE.,',
                    ' darwin_useSEAICE = .FALSE.,\n# no sea ice in the domain (radtrans: RT_iceconst = 0)\n iceconst = 0.,')
    dw = re.sub(r' DOPrunoff_interpMethod = 0,\n(.*?) DOPrunoff_interpMethod = 0,\n',
                r' DOPrunoff_interpMethod = 0,\n\1', dw, flags=re.S)
    dfields = ['iron']
    # river nutrients / carbon from the point sources (make_rivers.py, nest
    # grid, no interpolation); NO3 and NH4 given directly instead of DIN
    riv = {'DOC': 'DOC', 'DON': 'DON', 'DOP': 'DOP', 'IP': 'IP', 'DSi': 'DSi',
           'POC': 'POC', 'PON': 'PON', 'POP': 'POP', 'DIC': 'DIC'}
    for fld, tag in riv.items():
        dw = re.sub(rf" {fld}runofffile = '[^']*',",
                    f" {fld}runofffile = 'rivers/rivers_{tag}flux_128x72',", dw)
    dw = sub_once(dw, " DINrunofffile = './BGC_runoff/DIN_CCS_cutout',\n",
                  "# NO3 and NH4 from the river files (no DIN split)\n"
                  " NO3runofffile = 'rivers/rivers_NO3flux_128x72',\n"
                  " NO3runoffstartdate1 = 19920101,\n NO3runoffstartdate2 = 120000,\n"
                  " NO3runoffperiod     = 86400.0,\n"
                  " NH4runofffile = 'rivers/rivers_NH4flux_128x72',\n"
                  " NH4runoffstartdate1 = 19920101,\n NH4runoffstartdate2 = 120000,\n"
                  " NH4runoffperiod     = 86400.0,\n")
    dw = re.sub(r' R_DIP_IP_runoff\s*=\s*[0-9.Ee+-]+,', ' R_DIP_IP_runoff    = 1.0,', dw)
    # TSS from the river point sources (make_rivers.py), on the nest grid
    dw = sub_once(dw, " TSSrunofffile = './BGC_runoff/TSS_CCS_cutout',\n",
                  "# Malibu Creek sediment (Q x rating-curve SSC) on the nest grid\n"
                  " TSSrunofffile = 'rivers/rivers_TSSflux_128x72',\n")
    dw = re.sub(r' TSSrunoff_interpMethod = 0,\n', '', dw)
    for f in dfields + [x + 'runoff' for x in ['DOC', 'DON', 'DOP', 'DIN', 'IP', 'DSi',
                                                  'POC', 'POP', 'PON', 'DIC']]:
        dw = re.sub(rf' {f}_interpMethod = 0,\n', '', dw)
    interp = ('\n'.join(parent_interp(f, dely) for f in dfields) + '\n TSSrunoff_interpMethod = 0,\n'
              + '\n'.join(f' {x}runoff_interpMethod = 0,' for x in
                           ['DOC', 'DON', 'DOP', 'IP', 'DSi', 'POC', 'PON', 'POP', 'DIC', 'NO3', 'NH4']))
    dw = dw.replace(' &DARWIN_INTERP_PARAMS\n', ' &DARWIN_INTERP_PARAMS\n' + interp + '\n#\n')
    dw = dw.replace(' &DARWIN_TSS_PARAMS\n/',
                    ' &DARWIN_TSS_PARAMS\n'
                    '# settling, deposition and resuspension of TSS are done by pkg/sediment\n'
                    ' wTSS_sink = 0.,\n katten_TSS = 0.07,\n/\n'
                    ' &DARWIN_RADTRANS_PARAMS\n'
                    " darwin_waterAbsorbFile='OPTICS_COEFF2/optics_water.txt',\n"
                    " darwin_phytoAbsorbFile='OPTICS_COEFF2/optics_plankton.txt',\n"
                    " darwin_particleAbsorbFile='OPTICS_COEFF2/optics_detritus.txt',\n"
                    ' darwin_part_size_P = 1e-15,\n'
                    ' /')
    # kelp light saturation / compensation in umol photons m-2 s-1 (darwin PAR
    # units): G-MACMODS Macrocystis values 212.4 and 20.45 (mag_species.py,
    # stored there as /4.57 W m-2; Arzeno-Soltero et al. 2023 Suppl. Table 1).
    # v05_kelp has 213.7942 / 20.5841, the same values after a 4.57 -> 4.6
    # round trip.
    dw = sub_once(dw, ' mp_spp_PARs = 213.7942\n mp_spp_PARc = 20.5841\n',
                  '# malibu_v01: G-MACMODS Macrocystis (umol photons m-2 s-1)\n'
                  ' mp_spp_PARs = 212.4\n mp_spp_PARc = 20.45\n')
    # DARWIN_ALLOW_CDOM: group read right after DARWIN_PARAMS; darwin3 defaults
    dw = sub_once(dw, ' /\n &DARWIN_MACROALGAE_PARAMS\n',
                  ' /\n# CDOM tracer (DARWIN_ALLOW_CDOM), darwin3 24885b71 code defaults: fracCDOM = 0.02,\n'
                  '# CDOMdegrd = 1/200 d, CDOMbleach = 1/15 d, PARCDOM = 20, R_NP/R_CP/R_FeP_CDOM\n'
                  '# = 16/120/1e-3, CDOMcoeff = 100 m-1/(mmol P m-3) (a at 450 nm)\n'
                  ' &DARWIN_CDOM_PARAMS\n /\n &DARWIN_MACROALGAE_PARAMS\n')
    if v5:
        dw = dw.replace("'era_xx_it42_v2/CCS_kelp_80x100_iron_dust'", "'v5r1/iron_dust.bin'")
        assert 'era_xx' not in dw, 'unreplaced CCS forcing in data.darwin'
    open(os.path.join(INP, 'data.darwin'), 'w').write(dw)
    if v5:
        write_v5r1_extras(INP)
        if a.diag_pertype:
            dg = open(os.path.join(INP, 'data.diagnostics')).read()
            lv = ', '.join(f'{k}.' for k in range(1, 21))
            add = ('#--- [K8] per-plankton rates and carbon, top 20 levels (0-40 m), daily: PPnnnn primary\n'
                   '#    production, GRnnnn grazing loss (mmol C/m3/s), TRAC21-27 = c01-c07 (Chl:C with chl)\n'
                   '  frequency(15) = 86400.0,\n'
                   '   fields(1:10,15) = ' + ', '.join(f"'PP{j:04d}  '" for j in range(1, 6)) + ', '
                   + ', '.join(f"'GR{j:04d}  '" for j in range(1, 6)) + ',\n'
                   f'   levels(1:20,15) = {lv},\n'
                   "   filename(15) = 'diags/daily/rates',\n"
                   '  frequency(16) = 86400.0,\n'
                   '   fields(1:7,16) = ' + ', '.join(f"'TRAC{k}  '" for k in range(21, 28)) + ',\n'
                   f'   levels(1:20,16) = {lv},\n'
                   "   filename(16) = 'diags/daily/carbon',\n")
            assert dg.count(' /\n &DIAG_STATIS_PARMS') == 1
            dg = dg.replace(' /\n &DIAG_STATIS_PARMS', add + ' /\n &DIAG_STATIS_PARMS')
            open(os.path.join(INP, 'data.diagnostics'), 'w').write(dg)
        if a.diag_budget:
            dg = open(os.path.join(INP, 'data.diagnostics')).read()
            lv = ', '.join(f'{k}.' for k in range(1, 21))
            flds = [f'{q}Tr{t}' for t in (23, 25) for q in ('ADVx', 'ADVy', 'ADVr', 'DFrI')]
            add = ('#--- transport of the dominant phytoplankton C (c03 = TRAC23, c05 = TRAC25): advective\n'
                   '#    fluxes (x, y, r) and implicit vertical diffusion, top 20 levels, daily (budget closure)\n'
                   '  frequency(17) = 86400.0,\n'
                   f'   fields(1:{len(flds)},17) = ' + ', '.join(f"'{x:<8s}'" for x in flds) + ',\n'
                   f'   levels(1:20,17) = {lv},\n'
                   "   filename(17) = 'diags/daily/phyto_flux',\n")
            dg = dg.replace(' /\n &DIAG_STATIS_PARMS', add + ' /\n &DIAG_STATIS_PARMS')
            open(os.path.join(INP, 'data.diagnostics'), 'w').write(dg)
        if a.kelp_death is not None:
            dw = open(os.path.join(INP, 'data.darwin')).read()
            assert ' mp_spp_death = 1.15E-7\n' in dw, 'mp_spp_death line not found'
            dw = dw.replace(' mp_spp_death = 1.15E-7\n',
                            '# kelp background mortality %g /d: SBC LTER Macrocystis winter (Dec-Feb) frond loss 0.0151 /d\n'
                            '# (Rassweiler et al. 2008, Ecology 89:2068, Ecol. Archives E089-119 loss-rate data); within\n'
                            '# the G-MACMODS Monte Carlo range 0.003-0.017 /d (magpy/mag_montecarlo.py)\n'
                            ' mp_spp_death = %.4E\n' % (a.kelp_death, a.kelp_death / 86400.))
            open(os.path.join(INP, 'data.darwin'), 'w').write(dw)
        if a.chl2cmax_small is not None:
            tr = open(os.path.join(INP, 'data.traits')).read()
            old = ' CHL2CMAX        = 2*0.500000000000000       , 3*0.200000000000000       , 2*0.000000000000000E+000  ,\n'
            assert old in tr, 'CHL2CMAX line not found in data.traits'
            new = ('# CHL2CMAX classes 3-5 raised from 0.2 (C:Chl >= 60 g/g) to %g mg Chl/mmol C (C:Chl >= %.0f g/g):\n'
                   '# inshore southern CCE winter mixed-layer C:Chl 18-63 g/g, mean ~38 (Taylor, Landry, Selph &\n'
                   '# Wokuluk 2015, DSR II 112:117, Table 1)\n'
                   ' CHL2CMAX        = 2*0.500000000000000       , 3*%.15f       , 2*0.000000000000000E+000  ,\n'
                   % (a.chl2cmax_small, 12.011 / a.chl2cmax_small, a.chl2cmax_small))
            open(os.path.join(INP, 'data.traits'), 'w').write(tr.replace(old, new))
        if a.phyto_mort is not None:
            tr = open(os.path.join(INP, 'data.traits')).read()
            old = ' MORT    = 5*1.157407407407407E-006  , 2*3.858024691358024E-007  ,\n'
            assert old in tr, 'MORT line not found in data.traits'
            new = ('# MORT phytoplankton classes 1-5 cut from 0.1 to %g /d: S2 coastal losses (grazing 0.15-0.27\n'
                   '# + MORT 0.1 /d) exceed winter dilution-experiment mortality at Port of LA, 0-0.24 /d\n'
                   '# (Connell et al. 2017, DSR I 121:210, Table 1); sensitivity test\n'
                   ' MORT    = 5*%.15E  , 2*3.858024691358024E-007  ,\n' % (a.phyto_mort, a.phyto_mort / 86400.))
            open(os.path.join(INP, 'data.traits'), 'w').write(tr.replace(old, new))
        if a.waves:
            wv = json.load(open(os.path.join(a.waves, 'waves.json')))
            sd = open(os.path.join(INP, 'data.sediment')).read()
            sd = re.sub(r'# waves off.*?\n &SEDIMENT_PARM03\n sedWaveMode = 0,\n',
                        ('# ERA5 waves (python/make_waves.py; offshore 0.5-deg, nearest value at the\n'
                           '# coast, no nearshore transformation), hourly from\n'
                           if wv.get('source', 'era5') == 'era5' else
                           '# WAVEWATCH III wcoast 1/6-deg offshore + CDIP MOP 10-m points nearshore\n'
                           '# (python/make_waves.py --source mop_ww3), hourly from\n')
                        + '# the run start; orbital velocity from the local depth (linear theory)\n'
                        ' &SEDIMENT_PARM03\n sedWaveMode = 3,\n'
                        f" sedWaveHsFile = 'waves/sed_Hs_{nx}x{ny}',\n"
                        f" sedWavePerFile = 'waves/sed_Tp_{nx}x{ny}',\n"
                        f" sedWaveDirFile = 'waves/sed_Dir_{nx}x{ny}',\n"
                        f" sedWaveNrec = {wv['nrec_sed']},\n sedWaveDt = {wv['sedWaveDt']},\n", sd, flags=re.S)
            assert 'sedWaveMode = 3' in sd, 'data.sediment wave block not found'
            open(os.path.join(INP, 'data.sediment'), 'w').write(sd)
        if a.sed_bed:
            sd = open(os.path.join(INP, 'data.sediment')).read()
            new = ('# erodible bed (python/make_sed_bed.py): usSEABED mud % x Palos Verdes erodible\n'
                   '# mass 130 g/m2 at 65 % fines (Ferre et al. 2010, CSR, p18/Fig. 11d)\n'
                   f" sedBedFile = 'sediment/sed_bed_{nx}x{ny}',\n"
                   '# M = 0.01 kg/m2/s x (1 - porosity 0.9) (Moriarty et al. 2017, BG 14:1919, Table 3)\n'
                   ' sedEroM(1) = 1.0,\n'
                   '# dry bulk density ~1100 kg/m3, Palos Verdes shelf mud (Ferre et al. 2010, p20)\n'
                   ' sedBulkRho = 1.1E6,\n')
            for pat in (r" sedBedFile = ' ',\n", r' sedEroM\(1\) = 5\.E-2,\n', r' sedBulkRho = 1\.5E6,\n'):
                assert re.search(pat, sd), pat
                sd = re.sub(pat, '', sd)
            sd = sd.replace(' sedTauCD(1)', new + ' sedTauCD(1)', 1)
            if a.sed_layered:
                # layered bed: E = sedEroM (1 - phi) f (tau_b/tau_ce - 1); keep the effective rate of the
                # single-store bed: M = 0.01 kg/m2/s with phi = 0.9 (Moriarty et al. 2017, Table 3)
                sd = sd.replace('# M = 0.01 kg/m2/s x (1 - porosity 0.9) (Moriarty et al. 2017, BG 14:1919, Table 3)\n'
                                ' sedEroM(1) = 1.0,\n',
                                '# M = 0.01 kg/m2/s = 10 g/m2/s; layered bed applies (1 - sedPorosity) (Moriarty et al.\n'
                                '# 2017, BG 14:1919, Table 3, phi = 0.9)\n sedEroM(1) = 10.0,\n')
                assert 'sedEroM(1) = 10.0' in sd
                parm02 = (' &SEDIMENT_PARM02\n'
                          '# consolidating layered bed (Sherwood et al. 2018, GMD 11:1849; masses in g/m2 here,\n'
                          '# so ln-offsets in kg/m2 + ln(1000) = 6.908)\n'
                          ' sedNlay = 5,\n'
                          '# 5 layers ASSUMED\n'
                          ' sedCohesive(1) = .TRUE.,\n'
                          ' sedPorosity = 0.9,\n'
                          '# tau_cb,eq = exp((ln z_m - offset)/slope): slope 2 (Sherwood et al. 2018 p1861; range\n'
                          '# 1.7-2 in Suppl. Table S4); offset fitted to Palos Verdes shelf mud: 0.12-0.14 kg/m2\n'
                          '# remobilised at 0.32 Pa (Ferre et al. 2010, CSR, p18/Fig. 11d) -> ln(0.13) - 2 ln(0.32)\n'
                          '# = 0.24 (in Sherwood Table S4 range -0.47..0.3) + 6.908\n'
                          ' sedTcbOffset = 7.15,\n'
                          ' sedTcbSlope = 2.0,\n'
                          '# 0.05 / 1.5 Pa within 0.03-0.05 / 1.5-2.2 Pa (Sherwood et al. 2018 p1861-1862)\n'
                          ' sedTcbMin = 0.05,\n'
                          ' sedTcbMax = 1.5,\n'
                          '# consolidation 8 h, swelling 100 x (Sherwood et al. 2018 p1861, Suppl. Table S1)\n'
                          ' sedTcons = 28800.,\n'
                          ' sedTswell = 2880000.,\n'
                          '# active layer za = 0.007 m/Pa (tau_b - tau_ce) (Sherwood et al. 2018 Suppl. Table S1)\n'
                          ' sedZaCoef = 0.007,\n'
                          '# new-layer mass 10 g/m2 ASSUMED\n'
                          ' sedLayMass = 10.,\n'
                          ' /\n')
                assert ' &SEDIMENT_PARM02\n /\n' in sd
                sd = sd.replace(' &SEDIMENT_PARM02\n /\n', parm02, 1)
            open(os.path.join(INP, 'data.sediment'), 'w').write(sd)
    print('wrote data, data.obcs, data.exf, data.ptracers, data.darwin in', os.path.abspath(INP))


def write_v5r1_extras(inp):
    """data.cal, OASIM data.radtrans, daily diagnostics, and the static files."""
    import shutil
    src = os.path.join(HERE, '..', 'input')
    for f in ['data.pkg', 'data.sediment', 'data.kelpoptics', 'data.traits', 'data.gchem',
              'data.ggl90', 'eedata']:
        shutil.copy(os.path.join(src, f), os.path.join(inp, f))
    open(os.path.join(inp, 'data.cal'), 'w').write(f"""#
# Calendar Parameters (malibu_v01, PACE-year case: ICs from ECCO-Darwin v05 V5r1)
 &CAL_NML
 TheCalendar='gregorian',
 startDate_1={START},
 startDate_2=000000,
 calendarDumps = .TRUE.,
 &
""")
    files = '\n'.join([f" RT_EdFile({l:2d})='oasim/edsEd{l:02d}'," for l in range(1, 14)]
                      + [f" RT_EsFile({l:2d})='oasim/edsEs{l:02d}'," for l in range(1, 14)])
    rt = open(os.path.join(src, 'data.radtrans')).read()
    head = rt[:rt.index(' &RADTRANS_FORCING_PARAMS')]
    tail = rt[rt.index(' &RADTRANS_PARAMS'):]
    open(os.path.join(inp, 'data.radtrans'), 'w').write(
        "# radtrans: 13 wavebands, 400-700 nm (as in ECCO-Darwin v06)\n"
        "# OASIM monthly below-surface Ed/Es (Oliver Jahn, C61extrap/below, 1 deg global).\n"
        "# darwin3 24885b71 EXF reads RT files only as 12 monthly means (period -12, no\n"
        "# year suffix), so oasim/ holds run-window files from python/oasim_window.py\n"
        "# (Nov-Dec from 2024, Jan-Mar from 2025). Daily-mean light: no diel cycle.\n"
        " &RADTRANS_FORCING_PARAMS\n" + files + "\n RT_E_period = -12.,\n"
        " RT_iceconst = 0.,\n /\n"
        " &RADTRANS_INTERP_PARAMS\n RT_E_lon0 = -179.5,\n RT_E_lon_inc = 1.,\n RT_E_lat0 = -89.5,\n"
        " RT_E_lat_inc = 179*1.,\n RT_E_nlon = 360,\n RT_E_nlat = 180,\n RT_E_interpMethod = 1,\n /\n"
        + tail)

    def lst(n, names, fn, freq=86400., levels=None):
        f = ', '.join(f"'{x:<8s}'" for x in names)
        out = f"  frequency({n}) = {freq},\n   fields(1:{len(names)},{n}) = {f},\n"
        if levels:
            out += f"   levels(1:{len(levels)},{n}) = " + ', '.join(f'{v}.' for v in levels) + ',\n'
        return out + f"   filename({n}) = 'diags/daily/{fn}',\n"
    b13 = lambda p: [f'{p}{k:03d}' for k in range(1, 14)]
    d = ['# malibu_v01 PACE-year run (V5r1 BCs): daily means',
         ' &diagnostics_list', '  dumpatlast = .TRUE.,',
         lst(1, ['THETA', 'SALT', 'UVEL', 'VVEL', 'WVEL'], 'state_3d'),
         lst(2, ['ETAN', 'MXLDEPTH', 'oceQsw', 'TFLUX', 'SFLUX'], 'state_2d'),
         '#--- PAR, NO3, Chl01, CDOM (3-D)',
         lst(3, ['PAR', 'TRAC02', 'TRAC28', 'TRAC20'], 'light_nut'),
         '#--- kelp',
         lst(4, ['TRAC33', 'TRAC34', 'TRAC35', 'mu_mag', 'gT_mag', 'gQ_mag', 'gE_mag'], 'kelp'),
         '#--- TSS and sediment',
         lst(5, ['TRAC36'], 'TSS'),
         lst(6, ['SEDBED01', 'SEDDEP01', 'SEDERO01', 'SEDTAUB'], 'sediment'),
         '#--- optics for PACE (surface IOPs, Rirr, Ed/Es to 10 m)',
         lst(7, b13('Rirr'), 'Rirr'),
         lst(8, b13('a'), 'iop_a_sfc', levels=[1]),
         lst(9, b13('bb'), 'iop_bb_sfc', levels=[1]),
         lst(10, b13('Ed'), 'Ed_top', levels=[1, 2, 3, 4, 5, 6]),
         lst(11, b13('Es'), 'Es_top', levels=[1, 2, 3, 4, 5, 6]),
         lst(12, b13('aCDOM'), 'iop_acdom_sfc', levels=[1]),
         lst(13, b13('aprt'), 'iop_aprt_sfc', levels=[1]),
         '#--- all five Chl tracers (Chl01-05 = TRAC28-32; total Chl for PACE chlor_a)',
         lst(14, [f'TRAC{k}' for k in range(28, 33)], 'chl'),
         ' /', ' &DIAG_STATIS_PARMS', ' /', '']
    open(os.path.join(inp, 'data.diagnostics'), 'w').write('\n'.join(d))


if __name__ == '__main__':
    main()
