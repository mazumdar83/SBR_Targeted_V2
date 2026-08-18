#!/usr/bin/env bash
# Initialise git and push to GitHub (or any remote).
#
#   ./scripts/init_repo.sh https://github.com/<user>/bfgm.git
#   ./scripts/init_repo.sh                       # uses gh to create the repo
set -euo pipefail
cd "$(dirname "$0")/.."

REMOTE="${1:-}"

[ -d .git ] || git init -b main
git add -A
git diff --cached --quiet || git commit -m "bfgm: function-to-gene pipeline with Claude Code agents"

if [ -z "$REMOTE" ]; then
  command -v gh >/dev/null 2>&1 || {
    echo "No remote given and gh CLI not found."
    echo "Usage: $0 https://github.com/<user>/bfgm.git"
    exit 1
  }
  gh repo create bfgm --private --source=. --remote=origin --push
  echo "Created and pushed via gh."
  exit 0
fi

if git remote get-url origin >/dev/null 2>&1; then
  git remote set-url origin "$REMOTE"
else
  git remote add origin "$REMOTE"
fi
echo "Remote set to $REMOTE"
git push -u origin main
