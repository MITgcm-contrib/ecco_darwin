#!/usr/bin/env python3
"""Daily river discharge for LR17 from Hub'Eau (hydrometrie, obs_elab QmnJ).

Usage:  python3 LR17_rivers.py [START END]      (default 2022-02-01 2022-03-31)

For each river the gauges are tried from downstream to upstream (tidal
gauges near the mouths give levels only); the first with a complete daily
series over START..END is kept. Writes data/LR17/hubeau/<river>.json (raw
API answer) and data/LR17/hubeau/LR17_rivers_<START>_<END>.json:
{river: {station, name, lon, lat, dates, q (m3/s)}}.

The Boutonne joins the Charente at Carillon, downstream of the Charente
gauge (Beillant), so both are kept and the model adds them at the Charente
source. No correction is made for the catchment between gauge and mouth.
"""
import json
import os
import subprocess
import sys
import time
import urllib.parse
from datetime import date

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(os.path.dirname(HERE), 'data', 'LR17', 'hubeau')
API = 'https://hubeau.eaufrance.fr/api/v2/hydrometrie/obs_elab'

# river: [(station code, name, lon, lat), ...] downstream first
GAUGES = {
    'Charente': [('R423001001', 'La Charente à Chaniers - Beillant', -0.529, 45.706),
                 ('R520001001', 'La Charente à Saintes - Pont-Palissy', -0.630, 45.747)],
    'Boutonne': [('R618001001', 'La Boutonne à Tonnay-Boutonne', -0.703, 45.965),
                 ('R614292901', "La Boutonne à Saint-Jean-d'Angély", -0.527, 45.937),
                 ('R608001003', 'La Boutonne à Saint-Séverin-sur-Boutonne', -0.404, 46.080)],
    'Seudre': [('S012001001', 'La Seudre à Saujon', -0.925, 45.674),
               ('S011451001', 'La Seudre à Corme-Écluse', -0.852, 45.646),
               ('S011401003', 'La Seudre à Saint-André-de-Lidon', -0.748, 45.603)],
    'Sevre': [('N602000101', 'La Sèvre niortaise à Maillé [Bazoin]', -0.776, 46.319),
              ('N611061001', 'La Sèvre niortaise à Marans', -0.957, 46.321),
              ('N430062002', 'La Sèvre niortaise à Niort', -0.520, 46.328)],
    'Lay': [('N330161010', 'Le Lay à Mareuil-sur-Lay-Dissais', -1.220, 46.536)],
}


def get(params, tries=5):
    url = API + '?' + urllib.parse.urlencode(params)
    for k in range(tries):
        try:                        # curl: this Python lacks the CA chain
            r = subprocess.run(['curl', '-s', '-m', '90', url], capture_output=True, check=True)
            return json.loads(r.stdout)
        except Exception as e:                                  # flaky API: retry
            print(f'  retry {k + 1}: {e}')
            time.sleep(5 * (k + 1))
    raise RuntimeError('Hub\'Eau not answering: ' + url)


def main(start, end):
    os.makedirs(D, exist_ok=True)
    ndays = (date.fromisoformat(end) - date.fromisoformat(start)).days + 1
    out = {}
    for river, gauges in GAUGES.items():
        for code, name, lon, lat in gauges:
            raw = get({'code_entite': code, 'grandeur_hydro_elab': 'QmnJ',
                       'date_debut_obs_elab': start, 'date_fin_obs_elab': end, 'size': 1000})
            obs = sorted((o['date_obs_elab'], o['resultat_obs_elab']) for o in raw.get('data', [])
                         if o.get('resultat_obs_elab') is not None)
            print(f'{river:9s} {code} {len(obs):3d}/{int(ndays)} days  {name}')
            if len(obs) >= ndays:
                json.dump(raw, open(os.path.join(D, f'{river}.json'), 'w'))
                q = [v / 1000.0 for _, v in obs]                 # l/s -> m3/s
                out[river] = {'station': code, 'name': name, 'lon': lon, 'lat': lat,
                              'dates': [d for d, _ in obs], 'q': q}
                print(f'          kept: {min(q):.1f}..{max(q):.1f} m3/s, mean {sum(q) / len(q):.1f}')
                break
        else:
            print(f'{river}: no complete series')
    f = os.path.join(D, f'LR17_rivers_{start}_{end}.json')
    json.dump(out, open(f, 'w'), indent=1)
    print('wrote', f)


if __name__ == '__main__':
    main(*(sys.argv[1:3] if len(sys.argv) > 2 else ('2022-02-01', '2022-03-31')))
