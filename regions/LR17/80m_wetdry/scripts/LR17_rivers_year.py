#!/usr/bin/env python3
"""Finish the Hub'Eau river files for a long LR17 run (after LR17_rivers.py START END).

Usage:  python3 LR17_rivers_year.py START END

1. Boutonne: LR17_rivers.py keeps the first gauge with a complete series, which in 2025 is
   Saint-Severin (R608001003), upstream of Saint-Jean-d'Angely (R614292901), the gauge of the
   2022 runs. Rewrites LR17_rivers_START_END.json with Saint-Jean, its missing days filled
   from Saint-Severin times the ratio of their means over the common days.
2. River temperature: writes LR17_river_temp_START_END.csv (date,temp) from the Hub'Eau
   qualite_rivieres samples (parameter 1301) of the Charente at Taillebourg (05006900), the
   largest river. The samples are about monthly; gendata.py interpolates between them.
"""
import json
import os
import subprocess
import sys
import urllib.parse

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(os.environ.get('LR17_ROOT', os.path.dirname(HERE)), 'data', 'LR17', 'hubeau')
BASE = 'https://hubeau.eaufrance.fr/api/v2/'


def get(path, **params):
    url = BASE + path + '?' + urllib.parse.urlencode(params)
    return json.loads(subprocess.run(['curl', '-s', '-m', '90', url], capture_output=True, check=True).stdout)


def boutonne(start, end):
    fn = os.path.join(D, f'LR17_rivers_{start}_{end}.json')
    riv = json.load(open(fn))
    b = riv['Boutonne']
    if b['station'] == 'R614292901':
        print('Boutonne: already Saint-Jean')
        return
    days = [d[:10] for d in b['dates']]
    raw = get('hydrometrie/obs_elab', code_entite='R614292901', grandeur_hydro_elab='QmnJ',
              date_debut_obs_elab=start, date_fin_obs_elab=end, size=1000)
    sj = {o['date_obs_elab'][:10]: o['resultat_obs_elab'] / 1000.0 for o in raw['data']
          if o.get('resultat_obs_elab') is not None}
    up = np.array(b['q'])
    x = np.array([sj.get(d, np.nan) for d in days])
    g = np.isfinite(x)
    ratio = x[g].mean() / up[g].mean()
    print(f'Boutonne: Saint-Jean {g.sum()}/{len(days)} days, ratio {ratio:.3f}, r {np.corrcoef(x[g], up[g])[0, 1]:.3f}')
    riv['Boutonne'] = {'station': 'R614292901', 'name': "La Boutonne à Saint-Jean-d'Angély", 'lon': -0.527, 'lat': 45.937,
                       'dates': b['dates'], 'q': np.where(g, x, up * ratio).round(4).tolist(),
                       'filled_from': b['station'], 'filled_ratio': round(float(ratio), 4), 'filled_days': int((~g).sum())}
    json.dump(riv, open(fn, 'w'), indent=1)
    print('rewrote', fn)


def temperature(start, end):
    raw = get('qualite_rivieres/analyse_pc', code_parametre=1301, code_station='05006900', size=200,
              date_debut_prelevement=start, date_fin_prelevement=end, fields='date_prelevement,resultat')
    rows = {}
    for o in raw['data']:
        rows.setdefault(o['date_prelevement'][:10], o['resultat'])
    fn = os.path.join(D, f'LR17_river_temp_{start}_{end}.csv')
    with open(fn, 'w') as f:
        f.write('date,temp\n')
        for d in sorted(rows):
            f.write(f'{d},{rows[d]}\n')
    print('wrote', fn, len(rows), 'samples', min(rows.values()), '..', max(rows.values()), 'C')


if __name__ == '__main__':
    boutonne(*sys.argv[1:3])
    temperature(*sys.argv[1:3])
