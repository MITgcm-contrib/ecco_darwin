#!/bin/bash
# Nares OBCS for Jan 1 - Sep 1 2015 (2026-10-23, arch evaluation idea 1: baseline AM98+ICT to Aug 31; Dustin "yes" to staging).
# Same chain as job_obcs_v2.sh (OBCS2_rot_tf), but 5856 hourly parent records (244 d) instead of 4344, in a NEW shadow config
# (nares_obcs_JanAug) so nothing that existing runs read is touched. The throughflow offset is NOT re-derived: Q = 0.331 Sv,
# the value job_obcs_v2 used for OBCS2_rot_tf, so Jan-Jun must be identical to OBCS2_rot_tf (checked at the end).
# Output: nares_obcs_JanAug/forcings/OBCS_rot_tf, linked as nares/forcings/OBCS2_rot_tf_JanAug. Ask Dustin before qsub.
#PBS -l select=1:ncpus=28:model=bro_ele
#PBS -q normal
#PBS -l walltime=2:00:00
#PBS -N nares_obcsJA
#PBS -W group_list=s1353
#PBS -j oe
#PBS -o /nobackup/dcarrol2/downscaling/nares/prep/job_obcs_2015JanAug.out
set -u
D=/nobackup/dcarrol2/downscaling; N=$D/nares; X=$D/nares_obcs_JanAug; Q=0.331; NREC=5856
PR=/nobackupp19/dcarrol2/downscaling/parent_V5r1_westAO_nares/run   # dv output is in run/dv (fixed after 25303875 failed)
export MAMBA_ROOT_PREFIX="$D/micromamba_root"
eval "$("$D/micromamba/bin/micromamba" shell hook --shell bash)"
micromamba activate downscaling
die() { echo "FATAL: $1"; exit 1; }
grep -q "offset Q = 0.331 Sv" $D/nares_obcsfix/prep/job_obcs_v2.out || die "OBCS2 offset is not 0.331 Sv"
[ -e $X ] && die "$X exists (remove its symlinks first)"
mkdir -p $X/forcings $X/parent/outputs
for f in $N/*; do b=$(basename $f); case $b in forcings|parent) ;; *) ln -sfn $f $X/$b;; esac; done
ln -sfn $N/forcings/pickups $X/forcings/pickups
ln -sfn $N/parent/inputs $X/parent/inputs
# mirror parent/outputs (grid, pickups, ...) except the OBCS snapshot (fixed after 25303895 failed on parent/outputs/grid)
for f in $N/parent/outputs/*; do b=$(basename $f); case $b in OBCS*) ;; *) ln -sfn $(readlink -f $f) $X/parent/outputs/$b;; esac; done
SNAP=$N/parent/outputs/OBCS_2015JanAug
echo "===== snapshot $NREC records $(date +%H:%M)"
python3 $N/prep/snapshot_dv.py $PR/dv $N/parent/inputs $SNAP $NREC || die snapshot
ln -sfn $SNAP $X/parent/outputs/OBCS
cd $N
echo "===== gen_obcs_fast $(date +%H:%M)"
python3 prep/gen_obcs_fast.py -d $X -n Nares -bnd EWNS -seaice -p 28 -v -disc off || die gen_obcs
echo "===== rotate $(date +%H:%M)"
python3 prep/rotate_obcs_velocity.py $X || die rotate
echo "===== offset Q = $Q Sv (as OBCS2_rot_tf)"
python3 prep/nares_throughflow_offset.py $X $Q OBCS_rot OBCS_rot_tf || die offset
python3 prep/obcs_boundary_budget.py $X OBCS_rot_tf
echo "===== check: records 1-4344 identical to OBCS2_rot_tf; record count $NREC"
python3 - <<P || die check
import os, glob, numpy as np
new = '$X/forcings/OBCS_rot_tf'; old = '$N/forcings/OBCS2_rot_tf'
fs = sorted(glob.glob(new + '/*.bin')); assert len(fs) == 40, len(fs)
worst = 0
for f in fs:
    b = os.path.basename(f); so = os.path.getsize(old + '/' + b); sn = os.path.getsize(f)
    assert so % 4344 == 0 and sn % $NREC == 0 and so // 4344 == sn // $NREC, (b, so, sn)
    a = np.fromfile(old + '/' + b, '>f4'); c = np.fromfile(f, '>f4', count=a.size)
    d = np.abs(a - c).max(); s = max(np.abs(a).max(), 1e-30); worst = max(worst, d / s)
    assert np.isfinite(np.fromfile(f, '>f4')).all(), b
    print(f'{b:22s} Jan-Jun max|diff| {d:.3e} (rel {d/s:.1e})')
print('worst relative difference', worst); assert worst < 1e-5
P
ln -sfn $X/forcings/OBCS_rot_tf $N/forcings/OBCS2_rot_tf_JanAug
echo "===== ALL DONE $(date +%H:%M)"
