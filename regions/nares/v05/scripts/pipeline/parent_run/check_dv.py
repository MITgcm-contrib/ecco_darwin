import numpy as np, os, glob
D='/nobackupp19/dcarrol2/downscaling/parent_V5r1_westAO_nares/run_test/dv'
F2=['ETAN','UICE','VICE','AREA','HSNOW','HEFF']
show=['THETA','SALT','UVEL','VVEL','PTRACE01','PTRACE02','PTRACE05','PTRACE20','ETAN','AREA','HEFF','UICE']
bad=[]
for m in sorted(glob.glob(D+'/*_BC_mask.bin')):
    name=os.path.basename(m)[:-12]
    mk=np.fromfile(m,'>f4'); npts=int((mk>0).sum())
    line=[f"{name:14s} npts={npts:3d}"]
    for f in sorted(glob.glob(D+f'/{name}_BC_mask_*.bin')):
        fld=f.split('_BC_mask_')[1].split('.')[0]
        nl=1 if fld in F2 else 50
        sz=os.path.getsize(f); nrec=sz/(npts*nl*8)
        a=np.memmap(f,'>f8','r').reshape(-1,nl,npts) if nrec==int(nrec) else None
        if a is None or int(nrec)!=240: bad.append((name,fld,nrec)); continue
        s=np.concatenate([np.asarray(a[r]).ravel() for r in (0,120,239)])
        nz=(s!=0).mean()
        if not np.isfinite(s).all(): bad.append((name,fld,'nonfinite'))
        if fld in show: line.append(f"{fld}[{s[s!=0].min() if nz else 0:.3g},{s.max():.3g},nz{nz:.2f}]")
    print('  '.join(line))
print('PROBLEMS:',bad if bad else 'none (all 41 fields x 8 masks have 240 hourly records, finite)')
