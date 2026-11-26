#!/usr/bin/env bash
# Every manifest row's patch (tools/mutants/manifest.txt and manifest.d/*.txt; patch path relative to tools/mutants) must
# `git apply --check` cleanly against REPO's working tree. Prints `STALE <id>: <patch> does not apply` per failure and exits 1;
# a missing patch file is also STALE. Never modifies the tree. usage: tools/mutants_apply_check.sh [REPO]   (default: this checkout)
set -u
export LC_ALL=C
repo=$(cd "${1:-$(dirname "$0")/..}" && pwd) || { echo "mutants-apply-check: FAIL no repo"; exit 1; }
mdir=$repo/tools/mutants; n=0; bad=0
for m in "$mdir/manifest.txt" "$mdir"/manifest.d/*.txt; do
  [ -f "$m" ] || continue
  while IFS=$'\t' read -r _ id patch _; do
    case "$id" in ''|'#'*) continue ;; esac
    case "$patch" in '') continue ;; esac
    n=$((n + 1))
    if [ ! -f "$mdir/$patch" ] || ! git -C "$repo" apply --check "$mdir/$patch" >/dev/null 2>&1; then
      echo "STALE $id: $patch does not apply"; bad=$((bad + 1))
    fi
  done < <(grep -v '^[[:space:]]*#' "$m")
done
[ "$n" -gt 0 ] || { echo "mutants-apply-check: FAIL no manifest rows under $mdir"; exit 1; }
[ "$bad" -eq 0 ] || { echo "mutants-apply-check: FAIL $bad of $n patches do not apply"; exit 1; }
echo "mutants-apply-check: $n patches apply"
