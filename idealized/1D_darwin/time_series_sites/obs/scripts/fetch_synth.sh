#!/bin/bash
# Robust ERDDAP subset downloads for synthesis products (SOCAT, BGC-Argo).
# Usage: fetch_synth.sh OBS_ROOT [socat|bgcargo] [SITE ...]
# Each request is chunked in time and retried until curl succeeds and the file
# ends with a newline (this network truncates long transfers silently).
set -u
ROOT=$1; WHAT=$2; shift 2
SITES=${@:-HOT BATS HydroS Papa PAP}

latlon() { case $1 in
  HOT) echo "22.75 -158.0";; BATS) echo "31.67 -64.17";; HydroS) echo "32.17 -64.5";;
  Papa) echo "50.1 -144.9";; PAP) echo "49.0 -16.5";; esac; }

# fetch URL OUT -> 0 ok (data), 2 no-data, 1 failed
fetch() {
  local url=$1 out=$2 i code
  for i in 1 2 3 4 5 6 7 8; do
    code=$(curl -sg --max-time 1800 -o "$out.part" -w '%{http_code}' "$url")
    rc=$?
    if [ "$code" = "404" ] && grep -q 'no matching results\|nRows = 0\|Your query produced no matching' "$out.part" 2>/dev/null; then
      rm -f "$out.part"; return 2; fi
    if [ $rc -eq 0 ] && [ "$code" = "200" ] && [ "$(tail -c1 "$out.part" | od -An -c | tr -d ' ')" = '\n' ]; then
      mv "$out.part" "$out"; return 0; fi
    echo "  retry $i ($out rc=$rc http=$code)" >&2; sleep 5
  done
  rm -f "$out.part"; return 1
}

for site in $SITES; do
  read lat lon <<< "$(latlon $site)"
  d=$ROOT/$site/raw/$WHAT; mkdir -p "$d"
  if [ "$WHAT" = socat ]; then
    Q='expocode,platform_name,platform_type,qc_flag,time,latitude,longitude,depth,sal,temp,fCO2_recommended,WOCE_CO2_water'
    box=$(python3 -I -c "print(f'latitude%3E={$lat-1}&latitude%3C={$lat+1}&longitude%3E={$lon-1}&longitude%3C={$lon+1}')")
    base="https://data.pmel.noaa.gov/socat/erddap/tabledap/socat_v2026_fulldata.csv?$Q&$box&WOCE_CO2_water=%222%22"
    prefix=socat_v2026_box1deg
    # yearly chunks (multi-year requests time out with HTTP 408 on this server);
    # skip years already covered by an earlier multi-year chunk file.
    periods=""
    for y in $(seq 1990 2026); do
      covered=0
      for f in "$d"/${prefix}_*-*.csv; do
        [ -e "$f" ] || continue
        r=${f##*_}; r=${r%.csv}; a=${r%-*}; b=${r#*-}
        [ "$y" -ge "$a" ] && [ "$y" -lt "$b" ] && covered=1
      done
      [ $covered = 0 ] && periods="$periods $y-01-01/$((y+1))-01-01"
    done
  else
    V='platform_number,cycle_number,direction,time,latitude,longitude,position_qc'
    for v in pres temp psal doxy nitrate chla ph_in_situ_total bbp700; do V="$V,${v}_adjusted,${v}_adjusted_qc"; done
    box=$(python3 -I -c "import math;d=0.9/math.cos(math.radians($lat));print(f'latitude%3E={$lat-0.9:.3f}&latitude%3C={$lat+0.9:.3f}&longitude%3E={$lon-d:.3f}&longitude%3C={$lon+d:.3f}')")
    base="https://erddap.ifremer.fr/erddap/tabledap/ArgoFloats-synthetic-BGC.csv?$V&$box"
    periods=""; for y in $(seq 2002 2026); do periods="$periods $y-01-01/$((y+1))-01-01"; done
    prefix=argo_sprof_box100km
  fi
  echo "$base&time%3E=<t0>&time%3C<t1>" > "$d/query_url.txt"
  for p in $periods; do
    t0=${p%/*}; t1=${p#*/}; out="$d/${prefix}_${t0:0:4}-${t1:0:4}.csv"
    [ -s "$out" ] && continue
    fetch "$base&time%3E=${t0}T00:00:00Z&time%3C${t1}T00:00:00Z" "$out"; rc=$?
    echo "$site $WHAT $t0-$t1 rc=$rc $( [ -f "$out" ] && wc -l < "$out")"
  done
done
