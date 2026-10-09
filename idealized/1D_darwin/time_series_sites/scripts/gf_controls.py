"""
Green's-function control parameters for the 1-D ECCO-Darwin v05 columns.

Each control is a multiplicative perturbation (factor = 1 + delta) applied to one
or more namelist entries.  Global: the same factor is applied at every site.

Plankton types (v05; Carroll et al. 2020, JAMES, sec. 2.3; mapping checked against data.traits):
  1 diatoms (HASSI)            2 other large eukaryotes (large PCMAX, sinks, PIC)
  3 Synechococcus (small, uses NO3/NO2, PIC)   4 high-light Prochlorococcus (NH4 only)
  5 low-light Prochlorococcus (NH4 + NO2)
  6 zooplankton grazing small picoplankton (PALAT 1 on 3-5)
  7 zooplankton grazing large eukaryotes (PALAT 0.85/0.9 on 1-2)
Trait entries are edited in data.traits (it overrides DARWIN_RANDOM_PARAMS, see
darwin_init_fixed.F: GENERATE_RANDOM then READ_TRAITS).  Scalar rates are edited
in data.darwin &DARWIN_PARAMS; entries absent from data.darwin are added with the
darwin_ckpt68g default (value taken from darwin_params.txt of the control run).

entry = (file, variable, indices or None for scalar, 1-based)
"""

CONTROLS = [
    # name,        delta, entries, description
    ('PCMAX_diat',  0.2, [('data.traits', 'PCMAX', [1])], 'max growth, diatom'),
    ('PCMAX_LgEukSyn', 0.2, [('data.traits', 'PCMAX', [2, 3])], 'max growth, other large eukaryotes + Synechococcus'),
    ('PCMAX_Pro',   0.2, [('data.traits', 'PCMAX', [4, 5])], 'max growth, Prochlorococcus (HL + LL)'),
    ('KSAT_phyto',  0.5, [('data.traits', v, [1, 2, 3, 4, 5]) for v in
                          ['KSATNO3', 'KSATNO2', 'KSATNH4', 'KSATPO4', 'KSATFET', 'KSATSIO2']],
                         'nutrient half-saturation, all phyto'),
    ('CHL2CMAX',    0.2, [('data.traits', 'CHL2CMAX', [1, 2, 3, 4, 5])], 'max Chl:C'),
    ('MORT_phyto',  0.5, [('data.traits', 'MORT', [1, 2, 3, 4, 5])], 'linear mortality, phyto'),
    ('GRAZEMAX',    0.2, [('data.traits', 'GRAZEMAX', [6, 7])], 'max grazing rate'),
    ('MORT_zoo',    0.5, [('data.traits', 'MORT', [6, 7])], 'linear mortality, zooplankton'),
    ('R_PICPOC',    0.5, [('data.traits', 'R_PICPOC', [2, 3])], 'PIC:POC of large eukaryotes + Synechococcus'),
    ('wPOM_sink',   0.5, [('data.darwin', v, None) for v in
                          ['wC_sink', 'wN_sink', 'wP_sink', 'wFe_sink', 'wSi_sink']],
                         'POM sinking speed'),
    ('KPOM',        1.0, [('data.darwin', v, None) for v in ['KPOC', 'KPON', 'KPOP', 'KPOFe']],
                         'POM remineralization rate'),
    ('KDOM',        0.5, [('data.darwin', v, None) for v in ['KDOC', 'KDON', 'KDOP', 'KDOFe']],
                         'DOM remineralization rate'),
    ('KDISSC',      0.5, [('data.darwin', 'Kdissc', None)], 'PIC dissolution rate'),
    ('SCAV_RAT',    0.5, [('data.darwin', 'SCAV_RAT', None)], 'Fe scavenging rate'),
    ('ALPFE',      -0.2, [('data.darwin', 'ALPFE', None)], 'soluble fraction of dust Fe'),
    # gyre-bias controls (added for ensemble 2): light, photoacclimation, grazing
    ('KATTEN_W',   -0.25, [('data.darwin', 'katten_w', None)], 'light attenuation by water'),
    ('MQYIELD',     0.3, [('data.traits', 'MQYIELD', [1, 2, 3, 4, 5])], 'max quantum yield (P-I initial slope), phyto'),
    ('KGRAZESAT',   0.5, [('data.traits', 'KGRAZESAT', [6, 7])], 'grazing half-saturation'),
    # palat(prey, predator) column-major: prey 3-5 (Syn, Pro HL, Pro LL) of zooplankton 6 = entries 38-40
    ('PALAT_pico', -0.3, [('data.traits', 'PALAT', [38, 39, 40])], 'palatability of Syn + Pro to small grazer'),
]

# Run-directory aliases. Ensemble 1 (PBS 25239309, runs_v1) used PCMAX_cocco/PCMAX_small for
# PCMAX_LgEukSyn/PCMAX_Pro; ensemble 2 (PBS 25241021) uses the control names directly.
RUNDIR_V1 = {'PCMAX_LgEukSyn': 'PCMAX_cocco', 'PCMAX_Pro': 'PCMAX_small'}
RUNDIR = {}

# Optional 20th control (added 2026-10-08, gf_solve.py --kz): background vertical diffusivity at 75-250 m,
# a stand-in for eddy/internal-wave nutrient supply. Not a namelist entry: edits diffkr_1x1x50, so its
# perturbation runs come from scripts/gf_kz_runs.py (Mac, vs runs/baseline) and its sensitivity column
# from scripts/gf_kz_cache.py. (name, delta, (z_top, z_bottom) in m, description)
KZ_CONTROL = ('KZ_BG', 1.0, (75., 250.), 'background diffusivity at 75-250 m')
