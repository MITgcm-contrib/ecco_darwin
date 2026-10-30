#!/usr/bin/env python3
"""Download CMEMS IBI reanalysis fields (IBI_MULTIYEAR_PHY_005_002, 1/36 deg,
and IBI_MULTIYEAR_WAV_005_006) over the LR17 domain for the prototype run.

Usage:  python3 LR17_ibi.py START END [OUTDIR] [--only=ssh,cur]
        (dates YYYY-MM-DD; default OUTDIR data/LR17/ibi; --only downloads just the
        listed datasets, e.g. ssh,cur for the band-pass of gendata.py --fes-boundary)
Box: lon -1.80 .. -0.85, lat 45.65 .. 46.55 (LR17 80 m grid: lon -1.7001 ..
-0.9513, lat 45.7499 .. 46.4499, plus one IBI cell of margin).
Writes one NetCDF file per dataset:
  ssh    zos, hourly                     (tides + surge + mean)
  cur    ubar, vbar (depth mean) and surface uo, vo, hourly
  cur3d  uo, vo (3-D), daily
  temp3d thetao (3-D), daily
  sal    so (3-D), daily
  (the reanalysis has hourly fields at the surface only: tides and surge
  come from the hourly sea level and depth-mean current, the baroclinic
  structure from the daily 3-D fields)
  wav    VHM0, VTPK/VTM02, VMDR, hourly (waves for pkg/sediment)
  static deptho, mask
Needs Copernicus Marine credentials (~/.copernicusmarine).
"""
import os
import sys

import copernicusmarine as cm

TOP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BOX = dict(minimum_longitude=-1.80, maximum_longitude=-0.85,
           minimum_latitude=45.65, maximum_latitude=46.55)
SETS = {
    'ssh':  ('cmems_mod_ibi_phy-ssh_my_0.027deg_PT1H-m', ['zos']),
    'cur':  ('cmems_mod_ibi_phy-cur_my_0.027deg_PT1H-m', ['ubar', 'vbar', 'uo', 'vo']),
    'cur3d': ('cmems_mod_ibi_phy-cur_my_0.027deg_P1D-m', ['uo', 'vo']),
    'temp3d': ('cmems_mod_ibi_phy-temp_my_0.027deg_P1D-m', ['thetao']),
    'sal':  ('cmems_mod_ibi_phy-sal_my_0.027deg_P1D-m', ['so']),
    'wav':  ('cmems_mod_ibi_wav_my_0.027deg_PT1H-i', ['VHM0', 'VTPK', 'VTM02', 'VMDR']),
}


def main(start, end, out, only=None):
    os.makedirs(out, exist_ok=True)
    if not os.path.exists(os.path.join(out, 'ibi_static.nc')):
        st = cm.open_dataset(dataset_id='cmems_mod_ibi_phy_my_0.027deg-3D_static', **BOX)
        st.to_netcdf(os.path.join(out, 'ibi_static.nc'))
    for key, (dsid, var) in SETS.items():
        if only and key not in only:
            continue
        fn = os.path.join(out, f'ibi_{key}_{start}_{end}.nc')
        if os.path.exists(fn):
            print(key, 'exists')
            continue
        try:
            ds = cm.open_dataset(dataset_id=dsid, variables=var, **BOX,
                                 start_datetime=start + 'T00:00:00',
                                 end_datetime=end + 'T00:00:00')
        except Exception as e:                   # a variable name may differ
            print(key, 'failed:', e)
            continue
        ds.to_netcdf(fn)
        print(key, dict(ds.sizes), f'{os.path.getsize(fn) / 1e6:.0f} MB')


if __name__ == '__main__':
    only = next((a.split('=', 1)[1].split(',') for a in sys.argv if a.startswith('--only=')), None)
    args = [a for a in sys.argv[1:] if not a.startswith('--only=')]
    main(args[0], args[1], args[2] if len(args) > 2 else os.path.join(TOP, 'data', 'LR17', 'ibi'), only)
