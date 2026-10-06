#!/bin/bash
# Rebuild references/mitgcm-index from upstream MITgcm master + darwin3 + Dustin's branch clones.
# Read-only on all clones except `git fetch origin` in the upstream reference clone (updates remote refs only).
# Usage: ~/.claude/skills/mitgcm-ecco/scripts/refresh_index.sh
set -euo pipefail
SKILL=$(cd "$(dirname "$0")/.." && pwd)
R=$HOME/Documents/research
UP=$R/ECCO/BBL/MITgcm                      # any clone with origin = MITgcm/MITgcm
TMP=$(mktemp -d /tmp/mitgcm_index.XXXXXX)
trap 'rm -rf "$TMP" /tmp/mitgcm_base_*' EXIT

git -C "$UP" fetch -q origin
mkdir -p "$TMP/src"
git -C "$UP" archive origin/master | tar -x -C "$TMP/src"
LABEL="origin/master $(git -C "$UP" log -1 --format='%h %ad' --date=short origin/master) ($(git -C "$UP" describe --tags --abbrev=0 origin/master) +)"

EXTRA=()
add() { [ -d "$2" ] && EXTRA+=(--extra "$1=$2${3:+:$3}") || echo "skip $1 ($2 missing)"; }
add bbl        "$R/ECCO/BBL/MITgcm"
add bbl_c68g   "$R/ECCO/BBL/MITgcm_c68g"
add wad        "$R/ECCO/wetting_drying/MITgcm"
add wadcheckin "$R/ECCO/MITgcm_wad_checkin"
add seaicebc   "$R/ECCO/sea_ice_BCs/MITgcm"
add d3backport "$R/debug/darwin3" darwin,radtrans

# external docs -> docs-ext.md. ECCO doc repos are shallow clones in the cache (pulled here);
# ecco_darwin and darwin3 are the user's clones, read as-is (never pulled).
ED=$HOME/.cache/mitgcm-ecco-sources/ecco-docs
mkdir -p "$ED"
for r in ECCO-GROUP/ECCO-v4-Configurations ECCO-GROUP/ECCO-v4-Python-Tutorial gaelforget/ECCOv4; do
  n=$(basename $r)
  if [ -d "$ED/$n/.git" ]; then git -C "$ED/$n" pull -q --ff-only || echo "pull failed: $n"
  else git clone -q --depth 1 "https://github.com/$r.git" "$ED/$n" || echo "clone failed: $n"; fi
done
DOCS=()
doc() { [ -d "$2" ] && DOCS+=(--docroot "$1=$2${3:+:$3}") || echo "skip doc $1 ($2 missing)"; }
doc ecco_v4_configs "$ED/ECCO-v4-Configurations"
doc ecco_v4_python  "$ED/ECCO-v4-Python-Tutorial"
doc eccov4_gael     "$ED/ECCOv4"
doc ecco_darwin     "$HOME/Documents/GitHub/ecco_darwin"
doc darwin3_manual  "$HOME/Documents/GitHub/darwin3/doc" darwin,radtrans

python3 "$SKILL/scripts/build_mitgcm_index.py" "$TMP/src" "$TMP/idx" \
    --darwin3 "$HOME/Documents/GitHub/darwin3" --label "$LABEL" "${EXTRA[@]}" "${DOCS[@]}"

rm -rf "$SKILL/references/mitgcm-index"
mv "$TMP/idx" "$SKILL/references/mitgcm-index"
echo "index refreshed: $(head -3 "$SKILL/references/mitgcm-index/README.md" | tail -1)"
if git -C "$SKILL" rev-parse --git-dir >/dev/null 2>&1; then
  git -C "$SKILL" add references/mitgcm-index
  git -C "$SKILL" commit -q -m "refresh MITgcm index: $LABEL" && echo "committed index refresh" || echo "index unchanged"
fi
