#!/bin/bash
# Nares BGC inputs, Jan 1 - Sep 1 2015 (2026-10-29, BGC_PLAN.md, Dustin "Approve, stage all"). Ask Dustin before qsub.
# Same steps and rules as west_AO (job_westAO_inputs.sh + job_prod3_inputs.sh): gen_pickups -bgc, derived rDOC (0.50 DOC) /
# CDOM (0.02 DOC, >= 0), clip negatives with mass compensation, seed c_k < 1e-5 -> Xmin in the upper 100 m (IC + ALL 4 OBs).
# Nutrient fix (step 8): pan_AO factor v3, same rules as west_AO (nutfix_nares.py). Dustin "Check, then fix"; check in BGC_PLAN.md section 10.
# Works in a NEW shadow config (nares_bgc_inputs) so nothing that run_arch_Y16 or other runs read is touched.
# Results: nares/forcings/pickups_bgc -> shadow forcings/pickups_clip_seed_nutfix; nares/forcings/OBCS_bgc_JanAug -> shadow forcings/OBCSn.
#PBS -l select=1:ncpus=28:model=bro_ele
#PBS -q normal
#PBS -l walltime=6:00:00
#PBS -N nares_bgcin
#PBS -W group_list=s1353
#PBS -j oe
#PBS -o /nobackup/dcarrol2/downscaling/nares/prep/job_bgc_inputs.out
set -u
D=/nobackup/dcarrol2/downscaling; N=$D/nares; X=$D/nares_bgc_inputs; NREC=5856; IT=0000604872
PR=/nobackupp19/dcarrol2/downscaling/parent_V5r1_westAO_nares/run
UTILS=/nobackup/dcarrol2/v05_V5r1/ecco_darwin/regions/downscaling/utils
P=$N/prep
export MAMBA_ROOT_PREFIX="$D/micromamba_root"
eval "$("$D/micromamba/bin/micromamba" shell hook --shell bash)"
micromamba activate downscaling
S0=${S0:-0}                    # resume: qsub -v S0=<first step>
step() { echo; echo "===== $1  $(date +%H:%M:%S)"; }
run() { [ "$1" -ge "$S0" ]; }
die() { echo "FATAL: $1"; exit 1; }
F=$X/forcings; SNAP=$N/parent/outputs/OBCS_2015JanAug_bgc

if run 0; then
step "0 shadow config $X"
[ -e $X ] && die "$X exists (resume with S0>=1, or remove it)"
mkdir -p $F $X/parent/outputs
for f in $N/*; do b=$(basename $f); case $b in forcings|parent) ;; *) ln -sfn $f $X/$b;; esac; done
ln -sfn $N/parent/inputs $X/parent/inputs
for f in $N/parent/outputs/*; do b=$(basename $f); case $b in OBCS*) ;; *) ln -sfn $(readlink -f $f) $X/parent/outputs/$b;; esac; done
ls -la $X/parent/outputs/pickups/pickup_ptracers.$IT.data $X/parent/outputs/pickups/pickup_darwin.$IT.data || die "parent BGC pickups missing"
fi

if run 1; then
step "1 snapshot PTRACE01-31, $NREC records (physics snapshot files linked from OBCS_2015JanAug)"
[ -e $SNAP ] && die "$SNAP exists"
mkdir -p $SNAP
for f in $N/parent/outputs/OBCS_2015JanAug/*.bin; do ln -s $(readlink -f $f) $SNAP/$(basename $f); done
python3 $P/snapshot_dv.py $PR/dv $N/parent/inputs $SNAP $NREC $(for i in $(seq -w 1 31); do echo PTRACE$i; done) || die snapshot
ln -sfn $SNAP $X/parent/outputs/OBCS
ls $SNAP | wc -l
fi

if run 2; then
step "2 gen_pickups -bgc (physics + sea ice + ptracers + darwin) -> shadow forcings/pickups"
cd $N
python3 $UTILS/gen_pickups.py -d $X -n Nares -i ${IT#0000} -bgc -v -oh raise || die gen_pickups
ls -la $F/pickups | head -80
cmp $F/pickups/pickup_THETA.$IT.data $N/forcings/pickups/pickup_THETA.$IT.data && echo "THETA IC identical to the physics IC" || echo "NOTE: THETA IC differs from nares/forcings/pickups (only ptracer pickups are used)"
fi

if run 3; then
step "3 derived rDOC (0.50 DOC), CDOM (0.02 DOC, >= 0) -- IC"
python3 $P/gen_derived_tracer_pickup.py --pickup-dir $F/pickups --name rDOC --frac 0.50 --iter ${IT#0000} --force || die rDOC_ic
python3 $P/gen_derived_tracer_pickup.py --pickup-dir $F/pickups --name CDOM --frac 0.02 --clip-min 0 --iter ${IT#0000} --force || die CDOM_ic
fi

if run 4; then
step "4 clip + seed IC"
python3 $P/clip_bgc_negatives_nares.py ic $F/pickups $F/pickups_clip $IT || die clip_ic
python3 $P/seed_plankton_nares.py selftest || die selftest
python3 $P/seed_plankton_nares.py ic $F/pickups_clip $F/pickups_clip_seed $IT || die seed_ic
fi

if run 5; then
step "5 gen_obcs_fast -bgc (EWNS, $NREC records) -> shadow forcings/OBCS"
cd $N
python3 $P/gen_obcs_fast.py -d $X -n Nares -bnd EWNS -bgc -p 12 -v -disc off || die gen_obcs
ls $F/OBCS | wc -l
fi

if run 6; then
step "6 derived rDOC/CDOM -- OBCS"
python3 $P/gen_derived_tracer_obcs.py -d $X --parent-ptracers /nobackupp19/dcarrol2/v05_V5r1/darwin3/run3/data.ptracers --frac-rdoc 0.50 --frac-cdom 0.02 --clip-min-cdom 0 --nx-ew 360 --nx-ns 1120 --nrec $NREC --force || die derived_obcs
fi

if run 7; then
step "7 clip + seed OBCS (all 4 edges)"
python3 $P/clip_bgc_negatives_nares.py obcs $F/OBCS $F/OBCSc || die clip_obcs
python3 $P/seed_plankton_nares.py $F/OBCSc $F/OBCSs || die seed_obcs
fi

if run 8; then
step "8 nutrient fix (pan_AO factor v3): IC x min(FOB,3); OB x min(FOB, WOA23 winter / Jan-Apr max)"
FACN=/nobackup/dcarrol2/claude_scratch/panao_nutfix/nutfix_factor_v3.npz; NCG=$N/grid/Nares_ncgrid.nc
python3 $P/nutfix_nares.py $FACN $NCG limit $F/OBCSs || die nutfix_limit
python3 $P/nutfix_nares.py $FACN $NCG ic $F/pickups_clip_seed $F/pickups_clip_seed_nutfix || die nutfix_ic
python3 $P/nutfix_nares.py $FACN $NCG ob $F/OBCSs $F/OBCSn || die nutfix_ob
fi

if run 9; then
step "9 checks + links"
python3 - <<PY || die check
import os, numpy as np
F = '$F'; NREC = $NREC
stems = [f'TRAC{i:02d}' for i in range(1, 32)] + ['rDOC', 'CDOM']
for s in stems:
    for e, n in [('north', 1120), ('south', 1120), ('east', 360), ('west', 360)]:
        p = f'{F}/OBCSn/{s}_{e}.bin'
        assert os.path.getsize(p) == NREC * 90 * n * 4, p
        a = np.memmap(p, '>f4', 'r', shape=(NREC, 90, n))
        for r in (0, NREC // 2, NREC - 1):
            assert np.isfinite(a[r]).all(), (p, r)
print('OB files: 33 tracers x 4 edges, sizes and sampled records OK')
for i in list(range(1, 32)):
    nm = f'pickup_pTr{i:02d}.0000604872.data'
    a = np.fromfile(f'{F}/pickups_clip_seed_nutfix/{nm}', '>f8'); assert a.size == 90 * 360 * 1120 and np.isfinite(a).all(), nm
for nm in ('pickup_rDOC.0000604872.data', 'pickup_CDOM.0000604872.data'):
    a = np.fromfile(f'{F}/pickups_clip_seed_nutfix/{nm}', '>f8'); assert a.size == 90 * 360 * 1120 and np.isfinite(a).all() and a.min() >= 0, nm
print('IC pickups OK')
PY
ln -sfn $F/pickups_clip_seed_nutfix $N/forcings/pickups_bgc
ln -sfn $F/OBCSn $N/forcings/OBCS_bgc_JanAug
fi
step "ALL DONE"
