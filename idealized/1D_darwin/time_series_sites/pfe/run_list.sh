#!/bin/bash
# Run every run dir listed in $1 concurrently on this node (one serial mitgcmuv each).
source /usr/share/modules/init/bash
module purge
module load comp-intel/2020.4.304 hdf4/4.2.12 hdf5/1.8.18_serial netcdf/4.4.1.1_serial
ulimit -s unlimited
cd /nobackup/dcarrol2/1-D
while read d; do
  ( cd "$d" && rm -rf mnc_out/* STDOUT.* STDERR.* && ./mitgcmuv > output.txt 2>&1 ) &
done < "$1"
wait
while read d; do
  echo "$(hostname) $d : $(grep -c 'NORMAL END' $d/output.txt) NORMAL END"
done < "$1"
