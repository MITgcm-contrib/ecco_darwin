#!/bin/bash
# Refresh the raw material behind references/troubleshooting/: re-download the mitgcm-support
# archive + MITgcm GitHub issues/PRs, then rebuild the topic digests (only threads since YEAR).
# The distilled files are NOT regenerated automatically: read the new digest items and add/merge
# entries by hand (or with agents) into references/troubleshooting/<topic>.md, then bump the
# "Distilled through" date in references/troubleshooting/README.md.
# Usage: ~/.claude/skills/mitgcm-ecco/scripts/refresh_troubleshooting.sh [YEAR]   (default: this year)
set -euo pipefail
SKILL=$(cd "$(dirname "$0")/.." && pwd)
CACHE=$HOME/.cache/mitgcm-ecco-sources
SINCE=${1:-$(date +%Y)}
"$SKILL/scripts/fetch_support_sources.sh" "$CACHE"
python3 "$SKILL/scripts/digest_support_sources.py" "$CACHE" --since "$SINCE"
echo "digests for $SINCE+ in $CACHE/digest/ (compare against 'Distilled through' in references/troubleshooting/README.md)"
