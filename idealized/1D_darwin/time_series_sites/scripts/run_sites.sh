#!/bin/bash
# Run 1-D columns in parallel and wait: run_sites.sh <runs_dir> site1 site2 ...
# Writes <site>/<site>.pid; only these PIDs are ever touched.
cd "$1" || exit 1; shift
ulimit -s 65520
for s in "$@"; do
  ( cd "$s" && rm -rf mnc_out/* STDOUT.* STDERR.* output.txt pickup*ckpt* && exec ./mitgcmuv > output.txt 2>&1 ) &
  echo $! > "$s/$s.pid"
done
wait
for s in "$@"; do
  echo "$s: $(grep -c 'NORMAL END' $s/output.txt) NORMAL END, $(grep -c -i error $s/STDERR.0000) STDERR errors"
done
