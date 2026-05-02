#!/bin/bash

#SBATCH -A k10023
#SBATCH -p workq
#SBATCH -J ED_test
#SBATCH -t 23:59:00      # Specify Runtime in D-HH:MM  , for e.g 24 hrs.
#SBATCH --ntasks=1649   # Specify number of tasks, i.e. cores you need
#SBATCH -o %j.out       # Specify File to which standard out will be written
#SBATCH -e %j.err       # Specify File to which standard err will be written
#SBATCH --hint=nomultithread

date
echo "beginning: MITgcm id is: ${SLURM_JOB_ID}" 

 nsteps=58560
## nsteps=350400

start=`cat tmp.start`
if [ $start == 1 ]; then
nsteps=`expr $nsteps - 1`
fi
endstep=`expr $start + $nsteps`
cp data.temp2 data
mv data data.temp
cp data.temp data.temp2
mkdir -p diag3D
mkdir -p diag2D
###cat data.temp |sed 's/_stTIME_/'$start'/'|sed 's/_nsteps_/'$nsteps'/'|sed 's/_delt_/'$timestep'/'>data

cat data.temp |sed 's/_stTIME_/'$start'/'|sed 's/_nsteps_/'$nsteps'/'>data
echo $start $nsteps

ierr=0

sbatch --dependency=afterok:${SLURM_JOB_ID} run_sbatch.sh
#time srun -n 1656 --hint=nomultithread ./mitgcmuv65x  ||  ierr=1
#time srun -n 1656 ./mitgcmuv65x ||  ierr=1
/usr/bin/time srun -n 1649 --hint=nomultithread ./mitgcmuv ||  ierr=1

if [ $ierr == 1 ]; then
exit
fi

echo $endstep >tmp.start
npick1=`printf "%010i" $endstep`
cp data.temp2 data


## Rename pickup.ckptB with time for next run
##VAR=`grep timeStepNumber pickup.ckptA.meta | cut -d [ -f 2| cut -d ] -f 1`
##pickup_new_meta=$(printf "%s%010d%s" "pickup".$VAR."meta")
##pickup_new_data=$(printf "%s%010d%s" "pickup".$VAR."data")
##`cp pickup.ckptA.meta "$pickup_new_meta"`
##`cp pickup.ckptA.data "$pickup_new_data"`


##rm U.00* V.00* S.00* T.00* W.00* Eta.00*  ETAtave* wVeltave*

exit

