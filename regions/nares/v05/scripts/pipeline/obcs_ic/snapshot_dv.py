#!/usr/bin/env python3
"""Copy the first NREC hourly records of each Nares boundary field from the (still growing)
shared parent extraction into a stable directory, renamed to the ecco_darwin convention
<bnd>_BC_mask_<FLD>.<iter>.bin. gen_obcs infers the record count from file size, so reading
live files would give fields of different lengths.
The shared extraction writes big-endian float64 (check_dv.py reads '>f8'); gen_obcs reads
'>f4' (pan_AO's extraction was float32), so records are converted to float32 here.
Usage: snapshot_dv.py <src dv dir> <mask dir> <dst dir> <nrec> [fields...]"""
import sys, os, numpy as np
src, mdir, dst, nrec = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
flds = sys.argv[5:] or ['THETA', 'SALT', 'UVEL', 'VVEL', 'ETAN', 'AREA', 'HEFF', 'HSNOW', 'UICE', 'VICE']
NR3, IT = 50, '0000604873'          # parent LLC270 levels; first dv record iteration
os.makedirs(dst, exist_ok=True)
for b in ('east', 'west', 'north', 'south'):
    npts = int(np.count_nonzero(np.fromfile(os.path.join(mdir, f'{b}_BC_mask.bin'), '>f4')))
    for f in flds:
        nz = NR3 if (f in ('THETA', 'SALT', 'UVEL', 'VVEL') or f.startswith('PTRACE')) else 1   # 2026-10-29: PTRACEnn are 3-D
        fn = os.path.join(src, f'nares_{b}_BC_mask_{f}.{IT}.bin')
        have = os.path.getsize(fn) // (8 * npts * nz)
        if have < nrec:
            sys.exit(f'FATAL {fn}: {have} complete records < {nrec} requested')
        d = np.fromfile(fn, '>f8', count=nrec * npts * nz)
        assert np.isfinite(d).all() and np.abs(d).max() < 1e6, (fn, np.abs(d).max())
        d.astype('>f4').tofile(os.path.join(dst, f'{b}_BC_mask_{f}.{IT}.bin'))
    print(f'{b}: {npts} points, {nrec} records x {len(flds)} fields')
print('snapshot done')
