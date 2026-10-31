#!/usr/bin/env bash
# Prints one line per tracked or untracked non-ignored file of the current git worktree: MODE SHA256 PATH
# (symlink: `link SHA256(target)`, tracked but gone: `MISSING -`). The gate takes this before its first step
# and diffs it in tree-clean, so any content, mode, symlink, add or delete change is caught, including a
# change to a file that was already dirty (which `git status --porcelain` cannot see).
set -euo pipefail
git ls-files -z --cached --others --exclude-standard | sort -zu | while IFS= read -r -d '' f; do
  if [[ -L $f ]]; then printf 'link %s %s\n' "$(readlink -- "$f" | sha256sum | cut -c1-64)" "$f"
  elif [[ -f $f ]]; then printf '%s %s %s\n' "$(stat -c %a -- "$f")" "$(sha256sum < "$f" | cut -c1-64)" "$f"
  else printf 'MISSING - %s\n' "$f"; fi
done
