#!/bin/bash
# Stage run_bgc_test2d (2026-10-29, BGC_PLAN.md section 7; Dustin "Approve, stage all"). Ask Dustin before qsub.
# = run_arch_Y16 (P* 82500, eccen 1.6, scheme 77, capload, OBCS2_rot_tf_JanAug, Jan 1 cold start) with:
#   build_bgc binary; usePTRACERS/GCHEM/RADTRANS/OASIM on; input_darwin namelists (= pan_AO/west_AO ecosystem);
#   BGC IC from forcings/pickups_bgc and BGC OBCS from forcings/OBCS_bgc_JanAug (job_bgc_inputs.sh);
#   2 days (2880 steps); diagnostics = Y16 daily/3-hourly physics lists + BGC lists (no hourly, 6-hourly, monthly).
# Submit (after OK): qsub -q normal -l select=5:ncpus=40:mpiprocs=39:model=sky_ele -l walltime=4:00:00 -W group_list=s1353
#                    -N nares_bgct -v RUN=$R/run_bgc_test2d $N/prep/job_run_arch_long.sh
set -eu
N=/nobackup/dcarrol2/downscaling/nares; R=$N/regional_run; Y=$R/run_arch_Y16; T=$R/run_bgc_test2d
[ -e $T ] && { echo "FATAL: $T exists"; exit 1; }
mkdir $T; cd $T
# Y16 inputs: copy regular files, recreate symlinks (Y16 is not touched)
for f in $Y/*; do b=$(basename $f); case $b in diags|STD*|job_run.out) continue;; esac
  if [ -L $f ]; then ln -s $(readlink $f) $b; else cp -p $f $b; fi; done
mkdir -p diags/3hourly diags/daily
ln -sfn $R/build_bgc/mitgcmuv mitgcmuv
sed -i -E 's/^( *nTimeSteps=)[0-9]+,/\12880,/' data
sed -i -E 's/^( *(usePTRACERS|useGCHEM|useRADTRANS|useOASIM) *= *)\.FALSE\./\1.TRUE./' data.pkg
cp -p $R/input_darwin/data.darwin $R/input_darwin/data.gchem $R/input_darwin/data.oasim $R/input_darwin/data.traits $R/input_darwin/data.ptracers $R/input_darwin/data.radtrans .
# BGC IC (31 parent tracers + derived rDOC/CDOM) and OBCS (33 tracers x 4 edges)
for f in $N/forcings/pickups_bgc/pickup_pTr*.0000604872.* $N/forcings/pickups_bgc/pickup_rDOC.0000604872.* $N/forcings/pickups_bgc/pickup_CDOM.0000604872.* $N/forcings/pickups_bgc/pickup_darwin.0000604872.*; do ln -sfn $(readlink -f $f) $(basename $f); done
for f in $N/forcings/OBCS_bgc_JanAug/TRAC*_*.bin $N/forcings/OBCS_bgc_JanAug/rDOC_*.bin $N/forcings/OBCS_bgc_JanAug/CDOM_*.bin; do ln -sfn $(readlink -f $f) $(basename $f); done
# BGC forcing (same sources as west_AO stage_run_westAO.sh)
ln -sfn /nobackup/ojahn/forcing/oasim oasim
ln -sfn /nobackupp19/dcarrol2/v05_V5r1/ecco_darwin/regions/mac_delta/llc4320/AO_ecosystem/optics optics
ln -sfn /nobackupp19/dcarrol2/v05_V5r1/ecco_darwin/v05/3deg/data_darwin/3deg_Mahowald_2009_soluble_iron_dust.bin .
for f in /nobackupp19/dcarrol2/forcing/apCO2/NOAA_MBL/apCO2_*; do ln -sfn $f $(basename $f); done
cat > data.diagnostics <<'D'
# Nares BGC 2-day test (2026-10-29). Physics lists = run_arch_Y16 (1,2,3,4,6) for a day-1/2 comparison; BGC lists added.
# BGC names checked in code_bgc/darwin_diagnostics_init.F and radtrans_diagnostics_init.F.
 &DIAGNOSTICS_LIST
  frequency(1) = 10800.,
   fields(1:4,1) = 'SIuice  ','SIvice  ','SIarea  ','SIheff  ',
   filename(1) = 'diags/3hourly/SIdyn',
  frequency(2) = 86400.,
   fields(1:6,2) = 'SIshear ','SIdelta ','SItensil','SIsig1  ','SIsig2  ','SIhsnow ',
   filename(2) = 'diags/daily/SIdeform',
  frequency(3) = 86400.,
   fields(1:3,3) = 'ETAN    ','MXLDEPTH','oceQnet ',
   filename(3) = 'diags/daily/surf',
  frequency(4) = 86400.,
   fields(1:4,4) = 'SHIfwFlx','SHIhtFlx','SHIuStar','SHIgammT',
   filename(4) = 'diags/daily/SHI',
  frequency(6) = 86400.,
   fields(1:2,6) = 'THETA   ','SALT    ',
   levels(1,6) = 1.,
   filename(6) = 'diags/daily/surfTS',
# BGC surface (level 1): Chl, PP, plankton C, NO3 (TRAC02), PAR, pH, pCO2
  frequency(5) = 86400.,
   fields(1:7,5) = 'Chl     ','PP      ','plankC  ','TRAC02  ','PARF    ','pH      ','pCO2    ',
   levels(1,5) = 1.,
   filename(5) = 'diags/daily/bgcsurf',
# air-sea CO2
  frequency(7) = 86400.,
   fields(1:2,7) = 'fluxCO2 ','gDICsurf',
   filename(7) = 'diags/daily/co2',
# 3-D light and PP (check: PARF/RTtrans = 0 under the Petermann shelf)
  frequency(8) = 86400.,
   fields(1:3,8) = 'PARF    ','RTtrans ','PP      ',
   filename(8) = 'diags/daily/bgc3d',
 &
 &DIAG_STATIS_PARMS
 /
D
echo "--- diffs vs Y16"; for f in data data.pkg; do diff $Y/$f $f || true; done
echo "--- dangling links (BGC inputs not built yet?)"; find . -maxdepth 1 -xtype l | wc -l
echo STAGED $T
