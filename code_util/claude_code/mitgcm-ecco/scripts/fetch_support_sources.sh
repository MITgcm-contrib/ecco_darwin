#!/bin/bash
# Download the raw sources behind references/troubleshooting.md into a cache outside the skill:
#   - mitgcm-support mailing-list archive (pipermail monthly .txt.gz, 2003-present)
#   - MITgcm/MITgcm GitHub issues + PRs with their comments (via gh)
# Incremental: months already downloaded are skipped except the latest two; GitHub is re-pulled.
# Read-only on everything except the cache dir. Then run scripts/digest_support_sources.py.
# Usage: ~/.claude/skills/mitgcm-ecco/scripts/fetch_support_sources.sh [cache_dir]
set -euo pipefail
CACHE=${1:-$HOME/.cache/mitgcm-ecco-sources}
LIST=http://mailman.mitgcm.org/pipermail/mitgcm-support
mkdir -p "$CACHE/mail" "$CACHE/github"

# mailing list
months=$(curl -sf "$LIST/" | grep -o 'href="[0-9]\{4\}-[A-Za-z]*\.txt\.gz"' | sed 's/href="//;s/"//')
n=0
for m in $months; do
  n=$((n+1))
  f="$CACHE/mail/$m"
  if [ -s "$f" ] && [ $n -gt 2 ]; then continue; fi
  curl -sf -o "$f" "$LIST/$m" || echo "failed $m"
done
echo "mail: $(ls "$CACHE/mail" | wc -l | tr -d ' ') months in $CACHE/mail"

# GitHub issues + PRs (the issues endpoint returns both) and all comments
gh api --paginate 'repos/MITgcm/MITgcm/issues?state=all&per_page=100' \
  --jq '.[] | {number, title, state, created_at, closed_at, user: .user.login, labels: [.labels[].name], is_pr: (.pull_request != null), body}' \
  > "$CACHE/github/issues.jsonl"
gh api --paginate 'repos/MITgcm/MITgcm/issues/comments?per_page=100' \
  --jq '.[] | {issue: (.issue_url | split("/") | last | tonumber), user: .user.login, created_at, body}' \
  > "$CACHE/github/comments.jsonl"
gh api --paginate 'repos/MITgcm/MITgcm/pulls/comments?per_page=100' \
  --jq '.[] | {issue: (.pull_request_url | split("/") | last | tonumber), user: .user.login, created_at, path, body}' \
  > "$CACHE/github/review_comments.jsonl"
echo "github: $(wc -l < "$CACHE/github/issues.jsonl" | tr -d ' ') issues/PRs, $(wc -l < "$CACHE/github/comments.jsonl" | tr -d ' ') comments, $(wc -l < "$CACHE/github/review_comments.jsonl" | tr -d ' ') review comments"
date -u +%Y-%m-%dT%H:%MZ > "$CACHE/fetched_at"
