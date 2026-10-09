"""
Load the tidy observation files (obs/<SITE>/<SITE>_obs.csv) and convert to model units.

Model units: nutrients/DIC/ALK/O2/POC in mmol m-3, Chl mg m-3, PP mg C m-3 d-1,
pCO2 uatm, POC flux mg C m-2 d-1, temp degC (in-situ; model THETA is potential,
difference < 0.1 degC above 1000 m), salt psu.
"""
import os
import numpy as np
import pandas as pd

OBS_DIR = {'PAPA': 'Papa'}          # run-dir name -> obs dir name
RHO = 1.025                          # kg/L, per-kg -> per-m3 (umol/kg * 1.025 = mmol/m3)
CORE = ['NO3', 'PO4', 'SiO2', 'DIC', 'ALK', 'O2', 'Chl', 'PP', 'pCO2', 'POC_flux', 'temp', 'salt']


def _to_model_units(var, units, v):
    u = units.fillna('').str.lower()
    out = v.astype(float).copy()
    if var in ('NO3', 'PO4', 'SiO2', 'DIC', 'ALK', 'O2'):
        perkg = u.str.contains('kg')
        out[perkg] *= RHO
        if var == 'O2':
            ml = u.str.contains('ml/l')
            out[ml] *= 44.66                            # mL/L -> umol/L
    elif var == 'Chl':
        out[u.str.contains('kg')] *= RHO
    elif var == 'POC_flux':
        out[u.str.contains('mmol')] *= 12.011
        y = u.str.contains('mol/m2/yr') & ~u.str.contains('mmol')
        out[y] *= 12.011e3 / 365.25
    return out


def load_obs(site, obs_root, variables=CORE, t0='1992-01-01', t1='2026-01-01',
             drop_units=('fluorescence-derived',)):
    """Return tidy obs for one site in model units, with fCO2 folded into pCO2."""
    name = OBS_DIR.get(site, site)
    f = os.path.join(obs_root, name, '%s_obs.csv' % name)
    o = pd.read_csv(f, usecols=['time', 'depth_m', 'variable', 'value', 'units', 'source', 'qc'],
                    low_memory=False)
    o.loc[o['variable'] == 'theta', 'variable'] = 'theta'
    fc = o['variable'] == 'fCO2'
    o.loc[fc, 'value'] = o.loc[fc, 'value'].astype(float) * 1.0035   # fCO2 -> pCO2 (approx.)
    o.loc[fc, 'variable'] = 'pCO2'
    o = o[o['variable'].isin(variables)].copy()
    for d in drop_units:                     # e.g. uncalibrated fluorescence Chl
        o = o[~o['units'].fillna('').str.contains(d, case=False)]
    o['t'] = pd.to_datetime(o['time'], utc=True, errors='coerce').dt.tz_localize(None)
    o = o[(o['t'] >= t0) & (o['t'] < t1)].dropna(subset=['t', 'value'])
    o['depth_m'] = pd.to_numeric(o['depth_m'], errors='coerce').fillna(0.).clip(lower=0.)
    o['value_m'] = np.nan
    for v in o['variable'].unique():
        m = o['variable'] == v
        o.loc[m, 'value_m'] = _to_model_units(v, o.loc[m, 'units'], o.loc[m, 'value'])
    return o
