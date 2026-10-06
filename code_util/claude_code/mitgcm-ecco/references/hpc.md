# Running: NASA Pleiades and the local Mac

## Pleiades access (from Claude Code)

- Hosts: `pfe` (front end), `sfe` (2FA gateway), `lfe` (Lou archive). `~/.ssh/config`
  has a ControlMaster block (`ControlPersist 12h`) for `sfe pfe lfe`.
- The user opens the master connection themselves: suggest they type `! ssh -fN pfe`
  (passcode + PIN — never ask for or handle these). Then use non-interactive calls:
  `ssh -o BatchMode=yes pfe '...'`, `scp -o BatchMode=yes ...`.
  Check with `ssh -O check pfe`. If no master is up, give the user the commands to run.
- Use full PBS paths: `/PBS/bin/qstat -u <nas_user>`, `/PBS/bin/qsub`, `/PBS/bin/qalter`.
- **Ask before every qsub/qdel/qalter.** Work in copies of run dirs. Do analysis next to
  the data on a compute node (`qsub -I -q devel -l select=1:ncpus=40:model=sky_ele
  -l walltime=1:00:00`, or a small PBS script), never on the pfe login nodes: they are for
  editing, qsub/qstat and checks that take seconds. A ~20-min multi-core python run on
  pfe25 (Oct 2026) drew a NAS usage-policy violation and a CPU penalty. Vectorise first
  (e.g. KD-tree nearest-cell, not per-cell loops). Copy back only small results. Never work around NAS security
  controls. More detail: `ECCO/BBL/docs/Pleiades_from_Claude_Code.docx`,
  `research/HPC/NAS/notes/`.

## PBS job script conventions (csh)

```
#PBS -S /bin/csh
#PBS -l select=3:ncpus=40:model=sky_ele
#PBS -l walltime=08:00:00
#PBS -q long            # >8 h; use -q debug for short tests
#PBS -j oe
#PBS -m bea

module purge
module load comp-intel/2020.4.304 mpi-hpe/mpt hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt python3/3.9.12
limit stacksize unlimited
setenv FORT_BUFFERED 1
setenv MPI_BUFS_PER_PROC 128
unsetenv MPI_IB_RECV_MSGS
unsetenv MPI_UD_RECV_MSGS
cd $PBS_O_WORKDIR
mpiexec -np 113 /u/scicon/tools/bin/mbind.x ./mitgcmuv
```

- Node types: `sky_ele` (default, 40 cores), `bro` (28), `has` (24). `select × ncpus` must
  cover the rank count. A job stuck in queue can be moved, e.g.
  `qalter -l select=5:ncpus=28:model=bro` (ask first).
- Long runs are chained in segments (e.g. 5 yr) restarting from pickups; check
  `pickup*.data` exists and `nIter0` matches before resubmitting.
- Inode quota is tight: prefer monthly diagnostics and fewer streams. `useSingleCpuIO=.TRUE.` gives one
  global file per record (fewer inodes) but funnels I/O through rank 0; the offline LLC90 job uses
  `.FALSE.` deliberately for speed (~0.76 s/step). Create `diags/<stream>/` dirs before launch if used.
- The LatLon run_template on pfe has broken symlinks; copy `eedata` in rather than linking.
- Python on pfe: `module load python3/3.9.12; pip install --user MITgcmutils`.

## Remote paths

- User: `/nobackup/<nas_user>/` — e.g. `ECCO_V5r6_offline/{MITgcm,custom,archive}`,
  `v05_1deg_V4r5` (read-only ECCO-Darwin v05 LLC90; copied to `_bbl`),
  `sea_ice_BCs_latlon`, `downscaling/mac_delta/darwin3/run`; darwin3 under `/nobackupp19/<nas_user>/`.
- Colleagues (read-only): ECCO V4r6 ancillary `/nobackup/owang/runs/V4r6/PO.DAAC/ancillary_data/`;
  Oliver Jahn `/nobackup/ojahn/ecco_darwin/v06/llc270/data_darwin/`, `/nobackup/ojahn/forcing/oasim/`;
  Mike Wood `/nobackup/mwood7/Darwin/darwin3/`.
- Steph Dutkiewicz's MIT Engaging/ORCD paths (`/orcd/data/stephdut/001/`, `/home/stephdut/`)
  appear in `from_steph` namelists and must be repointed before running elsewhere.
- Syncing local edits: `rsync -avz <local code_*/run_*> pfe:/nobackup/<nas_user>/<run>/custom/`.
  The remote copy has drifted before — re-sync and diff before every submit.

## Local Mac runs

- Other Claude sessions run models here concurrently. Launch with
  `nohup ./mitgcmuv > run.log 2>&1 & echo $! > mitgcm.pid`; stop only that PID after
  checking `lsof -p <pid> | grep cwd`. Never `pkill`/`killall`/`pgrep -f mitgcmuv`.
- Wait on a run by watching for a marker (`grep -q "NORMAL END\|Execution ended" STDOUT.0000`)
  with Monitor/until-loops, not by polling process names. Never edit a script while it runs.
- MPI locally: `mpirun -np N ./mitgcmuv` with `nPx*nPy = N`. The Mac is shared: check `uptime` first
  and use ≤ 4 ranks.
