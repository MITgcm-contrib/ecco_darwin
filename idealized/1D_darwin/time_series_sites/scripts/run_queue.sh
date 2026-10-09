#!/bin/bash
# Run 1-D columns, at most MAXJOBS (default 4) at once, and wait for all:
#   run_queue.sh <runs_dir> site1 site2 ...
# Each <site>/<site>.pid holds the mitgcmuv PID. The Mac is shared with other
# model runs, so keep MAXJOBS small.
MAXJOBS=${MAXJOBS:-4}
cd "$1" || exit 1; shift
ulimit -s 65520
for s in "$@"; do
  while [ "$(jobs -rp | wc -l)" -ge "$MAXJOBS" ]; do sleep 20; done
  ( cd "$s" && rm -rf mnc_out/* STDOUT.* STDERR.* output.txt pickup*ckpt* && exec ./mitgcmuv > output.txt 2>&1 ) &
  echo $! > "$s/$s.pid"
done
wait
for s in "$@"; do
  echo "$s: $(grep -c 'NORMAL END' $s/output.txt) NORMAL END, $(grep -c -i error $s/STDERR.0000) STDERR errors"
done
