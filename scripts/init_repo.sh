#!/usr/bin/env bash
# Initialise git and push to Bitbucket.
#   ./scripts/init_repo.sh git@bitbucket.org:<workspace>/bfgm.git
set -euo pipefail
REMOTE="${1:?usage: init_repo.sh <bitbucket-remote-url>}"
cd "$(dirname "$0")/.."
[ -d .git ] || git init -b main
git add -A
git diff --cached --quiet || git commit -m "bfgm: function-to-gene mapping pipeline with seed agent skill"
git remote get-url origin >/dev/null 2>&1 && git remote set-url origin "$REMOTE" || git remote add origin "$REMOTE"
echo "Remote set to $REMOTE"
echo "Now run:  git push -u origin main"
