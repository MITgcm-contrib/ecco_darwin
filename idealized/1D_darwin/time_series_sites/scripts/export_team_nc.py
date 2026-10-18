#!/usr/bin/env python3
"""
Export 1-D run output as self-describing daily NetCDF time series for the team (2026-10-15; team meeting action item).

  python3 scripts/export_team_nc.py <out_dir> <label>=<run_dir>:<SITE> [...]

Writes <out_dir>/<SITE>_<label>_1992-2025.nc: daily means 1992-01-01..2025-12-31 on the wet model levels, with named
variables (no TRACnn), units, plankton-type labels and the run's provenance (experiment.txt, opt_manifest.txt) as
global attributes. Float32, zlib-compressed. Source: mnc_out/diags3D + diags2D + grid.
Unit conversions: PP mmol C m-3 s-1 -> mmol C m-3 d-1; fluxCO2 mmol C m-2 s-1 -> mmol C m-2 d-1, positive = ocean
uptake (FluxCO2 is added to the surface DIC tendency, code/darwin_surfforcing.F:287-324); pCO2 + apCO2 atm -> uatm.
"""
import os, sys
import numpy as np
import netCDF4 as nc

T0 = '1992-01-01 00:00:00'
# v05 plankton types (Carroll et al. 2020, JAMES); with the closure, type 2 is a nitrogen fixer
PFT = ['diatoms', 'other large eukaryotes', 'Synechococcus', 'high-light Prochlorococcus',
       'low-light Prochlorococcus']
ZOO = ['zooplankton grazing small picoplankton', 'zooplankton grazing large eukaryotes']
V3 = [  # name, source, scale, units, long name
    ('T', 'THETA', 1., 'degC', 'potential temperature'),
    ('S', 'SALT', 1., 'g kg-1', 'salinity'),
    ('Kz', 'GGL90Kr', 1., 'm2 s-1', 'GGL90 vertical diffusivity'),
    ('PAR', 'PAR', 1., 'uEin m-2 s-1', 'photosynthetically available radiation'),
    ('PP', 'PP', 86400., 'mmol C m-3 d-1', 'net primary production'),
    ('DIC', 'TRAC01', 1., 'mmol C m-3', 'dissolved inorganic carbon'),
    ('NO3', 'TRAC02', 1., 'mmol N m-3', 'nitrate'),
    ('NO2', 'TRAC03', 1., 'mmol N m-3', 'nitrite'),
    ('NH4', 'TRAC04', 1., 'mmol N m-3', 'ammonium'),
    ('PO4', 'TRAC05', 1., 'mmol P m-3', 'phosphate'),
    ('FeT', 'TRAC06', 1., 'mmol Fe m-3', 'total dissolved iron'),
    ('SiO2', 'TRAC07', 1., 'mmol Si m-3', 'silicic acid'),
    ('DOC', 'TRAC08', 1., 'mmol C m-3', 'dissolved organic carbon'),
    ('DON', 'TRAC09', 1., 'mmol N m-3', 'dissolved organic nitrogen'),
    ('DOP', 'TRAC10', 1., 'mmol P m-3', 'dissolved organic phosphorus'),
    ('DOFe', 'TRAC11', 1., 'mmol Fe m-3', 'dissolved organic iron'),
    ('POC', 'TRAC12', 1., 'mmol C m-3', 'detrital particulate organic carbon'),
    ('PON', 'TRAC13', 1., 'mmol N m-3', 'detrital particulate organic nitrogen'),
    ('POP', 'TRAC14', 1., 'mmol P m-3', 'detrital particulate organic phosphorus'),
    ('POFe', 'TRAC15', 1., 'mmol Fe m-3', 'detrital particulate organic iron'),
    ('POSi', 'TRAC16', 1., 'mmol Si m-3', 'detrital particulate silica'),
    ('PIC', 'TRAC17', 1., 'mmol C m-3', 'particulate inorganic carbon'),
    ('ALK', 'TRAC18', 1., 'meq m-3', 'alkalinity'),
    ('O2', 'TRAC19', 1., 'mmol O2 m-3', 'dissolved oxygen'),
] + [('c%02d' % (i + 1), 'TRAC%02d' % (20 + i), 1., 'mmol C m-3', 'plankton type %d carbon: %s' % (i + 1, n))
     for i, n in enumerate(PFT + ZOO)] \
  + [('Chl%02d' % (i + 1), 'TRAC%02d' % (27 + i), 1., 'mg Chl m-3', 'plankton type %d chlorophyll: %s' % (i + 1, n))
     for i, n in enumerate(PFT)]
V2 = [
    ('MLD', 'MXLDEPTH', 1., 'm', 'mixed-layer depth (MITgcm MXLDEPTH)'),
    ('pCO2', 'pCO2', 1e6, 'uatm', 'surface ocean pCO2'),
    ('apCO2', 'apCO2', 1e6, 'uatm', 'atmospheric pCO2'),
    ('fluxCO2', 'fluxCO2', 86400., 'mmol C m-2 d-1', 'air-sea CO2 flux, positive = ocean uptake'),
    ('pH', 'pH', 1., '1', 'surface pH'),
    ('Qnet', 'oceQnet', 1., 'W m-2', 'net surface heat flux into the ocean'),
    ('wspeed', 'EXFwspee', 1., 'm s-1', 'wind speed'),
    ('SIarea', 'SIarea', 1., '1', 'sea-ice fraction'),
]


def read(path):
    return open(path).read().strip() if os.path.exists(path) else ''


def export(out, lab, run, site):
    d3 = nc.Dataset(os.path.join(run, 'mnc_out', 'diags3D.0000000000.t001.nc'))
    d2 = nc.Dataset(os.path.join(run, 'mnc_out', 'diags2D.0000000000.t001.nc'))
    g = nc.Dataset(os.path.join(run, 'mnc_out', 'grid.t001.nc'))
    h = np.array(g['HFacC'][:]).ravel(); nz = int((h > 0).sum())
    zc = -np.array(g['RC'][:]).ravel()[:nz]; drf = np.array(g['drF'][:]).ravel()[:nz]
    te = np.array(d3['T'][:]) / 86400.                           # record = mean over the interval ending at T
    t = te - 0.5 * np.diff(np.r_[0., te])                         # interval centers (last record is a half day)
    assert np.array_equal(np.array(d2['T'][:]), np.array(d3['T'][:]))
    f = os.path.join(out, '%s_%s_1992-2025.nc' % (site, lab))
    o = nc.Dataset(f, 'w', format='NETCDF4')
    o.title = '1-D ECCO-Darwin v05 column at %s (%s set-up), daily means 1992-2025' % (site, lab)
    o.institution = 'Moss Landing Marine Laboratories, San Jose State University'
    o.contact = 'Dustin Carroll (dustin.carroll@sjsu.edu)'
    o.source = ('1-D ECCO-Darwin v05 (darwin3 darwin_ckpt68g) extracted from LLC270 V5r1; '
                'MITgcm-contrib/ecco_darwin idealized/1D_darwin/time_series_sites')
    o.site = site; o.setup = lab
    o.longitude = float(np.array(g['XC'][:]).ravel()[0]); o.latitude = float(np.array(g['YC'][:]).ravel()[0])
    o.bottom_depth_m = float(drf[:-1].sum() + drf[-1] * h[nz - 1])
    o.bgc_relaxation = 'NO3, PO4, FeT, SiO2, DIC, ALK, O2 relaxed (1 yr) to V5r1 below 200 m (linear ramp to 300 m)'
    o.experiment = read(os.path.join(run, 'experiment.txt'))
    if read(os.path.join(run, 'opt_manifest.txt')):
        o.parameter_changes = read(os.path.join(run, 'opt_manifest.txt'))
    o.history = 'created %s by scripts/export_team_nc.py from %s' % (np.datetime64('today'), os.path.basename(run))
    o.createDimension('time', len(t)); o.createDimension('depth', nz)
    v = o.createVariable('time', 'f8', ('time',)); v[:] = t
    v.units = 'days since ' + T0; v.calendar = 'standard'; v.long_name = 'center of the daily-mean interval'
    v = o.createVariable('depth', 'f4', ('depth',)); v[:] = zc
    v.units = 'm'; v.positive = 'down'; v.long_name = 'depth of model level center'
    v = o.createVariable('dz', 'f4', ('depth',)); v[:] = drf * h[:nz]
    v.units = 'm'; v.long_name = 'wet thickness of model level'
    for name, src, sc, u, ln in V3:
        a = np.ascontiguousarray(np.array(d3[src][:, :nz]).reshape(len(t), nz) * sc, 'f4')
        v = o.createVariable(name, 'f4', ('time', 'depth'), zlib=True, complevel=4)
        v[:] = a; v.units = u; v.long_name = ln
    a = sum(np.array(d3['TRAC%02d' % (27 + i)][:, :nz]).reshape(len(t), nz) for i in range(5))
    v = o.createVariable('Chl', 'f4', ('time', 'depth'), zlib=True, complevel=4)
    v[:] = np.ascontiguousarray(a, 'f4'); v.units = 'mg Chl m-3'; v.long_name = 'total chlorophyll (Chl01-Chl05)'
    for name, src, sc, u, ln in V2:
        v = o.createVariable(name, 'f4', ('time',), zlib=True, complevel=4)
        v[:] = np.ascontiguousarray(np.array(d2[src][:]).ravel() * sc, 'f4'); v.units = u; v.long_name = ln
    if 'DIAZO = 0, 1' in read(os.path.join(run, 'data.traits')).replace('  ', ' ') or 'diazotroph' in o.experiment:
        o['c02'].long_name = 'plankton type 2 carbon: nitrogen fixer (non-quota diazotroph, closure)'
        o['Chl02'].long_name = 'plankton type 2 chlorophyll: nitrogen fixer (non-quota diazotroph, closure)'
    o.close()
    print(f, '%.0f MB' % (os.path.getsize(f) / 1e6), flush=True)


if __name__ == '__main__':
    out = sys.argv[1]; os.makedirs(out, exist_ok=True)
    for spec in sys.argv[2:]:
        lab, rs = spec.split('='); run, s = rs.split(':')
        export(out, lab, run, s)
