"""Make Nares copies of the west_AO BGC tools (2026-10-29). Never edits the west_AO files."""
import re, shutil, os, sys
WP = '/nobackup/dcarrol2/downscaling/west_AO/regional_run/prep'
NP = '/nobackup/dcarrol2/downscaling/nares/prep'
def new(path):
    if os.path.exists(path): sys.exit(f'FATAL exists: {path}')
    return path
# 1. clip: Nares grid (Nr 90, Ny 360, Nx 1120; N/S edges 1120 points, E/W 360)
s = open(f'{WP}/clip_bgc_negatives.py').read()
assert s.count('NR, NY, NX = 90, 680, 860') == 1
s = s.replace('NR, NY, NX = 90, 680, 860', 'NR, NY, NX = 90, 360, 1120   # Nares copy 2026-10-29 (west_AO: 90, 680, 860)')
s = s.replace('"""Clip inherited', '"""NARES COPY (2026-10-29) of west_AO prep/clip_bgc_negatives.py: only NR/NY/NX changed.\nClip inherited', 1)
open(new(f'{NP}/clip_bgc_negatives_nares.py'), 'w').write(s)
# 2. seed: import the Nares clip module; OB mode seeds ALL FOUR edges (Dustin 2026-10-21 Nares note)
s = open(f'{WP}/seed_plankton_westao.py').read()
assert s.count('import clip_bgc_negatives as C') == 1
s = s.replace('import clip_bgc_negatives as C', 'import clip_bgc_negatives_nares as C')
s = s.replace('"""Seed plankton', '"""NARES COPY (2026-10-29) of west_AO prep/seed_plankton_westao.py: Nares grid via clip_bgc_negatives_nares,\nand OB mode seeds ALL FOUR edges (Dustin, 2026-10-21 Nares note: IC and all OBs, esp. the Lincoln Sea east edge).\nSeed plankton', 1)
old_start = s.index('    n = C.NX\n    names = {nm: f\'{C.OB[nm]}_north.bin\'')
old_end = s.index("    print(f'north: {nrec} records')\n") + len("    print(f'north: {nrec} records')\n")
newbody = '''    changed = set()
    for edge, n in [('north', C.NX), ('south', C.NX), ('east', C.NY), ('west', C.NY)]:
        names = {nm: f'{C.OB[nm]}_{edge}.bin' for nm in NAMES}
        paths = {nm: os.path.realpath(os.path.join(src, fn)) for nm, fn in names.items()}
        nrec = os.path.getsize(paths['DIC']) // (4 * C.NR * n)
        assert all(os.path.getsize(p) == nrec * 4 * C.NR * n for p in paths.values()), edge
        inp = {nm: np.memmap(p, '>f4', 'r', shape=(nrec, C.NR, n)) for nm, p in paths.items()}
        fh = {nm: open(os.path.join(dst, fn), 'wb') for nm, fn in names.items()}
        for r0 in range(0, nrec, 240):
            r1 = min(nrec, r0 + 240)
            f = {nm: np.array(inp[nm][r0:r1], 'f8') for nm in names}
            top = {nm: f[nm][:, :KTOP] for nm in names}            # views: seeding edits f in place
            seed(top, ratio)
            for nm in names:
                f[nm].astype('>f4').tofile(fh[nm])
        for h in fh.values():
            h.close()
        changed |= set(names.values())
        print(f'{edge}: {nrec} records', flush=True)
    C.link_rest(src, dst, changed)
'''
s = s[:old_start] + newbody + s[old_end:]
assert "C.link_rest(src, dst, set(names.values()))\n    print(f'north" not in s
s = s.replace("    C.link_rest(src, dst, set(names.values()))\n    for nm in PH + CHL", "    for nm in PH + CHL", 1) if "    C.link_rest(src, dst, changed)\n    C.link_rest" in s else s
open(new(f'{NP}/seed_plankton_nares.py'), 'w').write(s)
# 3. derived-tracer scripts: verbatim copies (grid sizes are arguments)
for f in ['gen_derived_tracer_pickup.py', 'gen_derived_tracer_obcs.py']:
    shutil.copy2(f'{WP}/{f}', new(f'{NP}/{f}'))
# 4. snapshot_dv.py: PTRACEnn are 50-level fields (backward compatible for physics)
p = f'{NP}/snapshot_dv.py'; s = open(p).read()
o = "nz = NR3 if f in ('THETA', 'SALT', 'UVEL', 'VVEL') else 1"
assert s.count(o) == 1
shutil.copy2(p, p + '.bak_pre_ptrace_20261029')
s = s.replace(o, "nz = NR3 if (f in ('THETA', 'SALT', 'UVEL', 'VVEL') or f.startswith('PTRACE')) else 1   # 2026-10-29: PTRACEnn are 3-D")
open(p, 'w').write(s)
print('TOOLS_DONE')
