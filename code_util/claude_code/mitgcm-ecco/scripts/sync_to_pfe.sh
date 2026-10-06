#!/bin/bash
# Push this skill to Pleiades (~/.claude/skills/mitgcm-ecco on pfe). Needs the ssh ControlMaster:
# run `! ssh -fN pfe` in Claude Code first if `ssh -O check pfe` fails.
set -euo pipefail
SKILL=$(cd "$(dirname "$0")/.." && pwd)
ssh -O check pfe >/dev/null 2>&1 || { echo "no pfe master connection: run 'ssh -fN pfe' first"; exit 1; }
ssh -o BatchMode=yes pfe 'mkdir -p ~/.claude/skills'
rsync -az --delete --exclude .git --exclude __pycache__ -e "ssh -o BatchMode=yes" "$SKILL/" pfe:.claude/skills/mitgcm-ecco/
echo "synced to pfe:~/.claude/skills/mitgcm-ecco"
