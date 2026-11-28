#!/usr/bin/env bash
# No shell file under tools/ or scripts/ pipes into `grep -q` (Spence, 2026-10-29, after d8d955f): runlog.sh sets pipefail, and
# grep -q exits at its first match, so the writer upstream can die of SIGPIPE (rc 141) and the pipeline reads as "no match"
# (mutate.sh's equivalent-id check falsely failed 070-L8c and 096-E-op-npd in Paul's 1520-full_gate). A `cmd && fail` on such a
# pipeline can falsely PASS. Use a here-string (`grep -q PAT <<< "$x"`) or read a file (`grep -q PAT FILE`).
# Banned on code lines (full-line comments are skipped): a single `|` (not `||`) followed by grep with any -q-bearing short option
# (-q, -qx, -qF, -qxF, -sq, ...) or --quiet / --silent. Shell files: *.sh, or a #!...sh / #!...bash first line.
# usage: tools/grepq_scan.sh [ROOT] (default: the repo root); scans ROOT/tools and ROOT/scripts. Exit 0 clean, 1 on hits.
set -u
export LC_ALL=C
root=${1:-$(cd "$(dirname "$0")/.." && pwd)}
re='(^|[^|])\|[[:space:]]*grep([[:space:]]+[^|;&)[:space:]]+)*[[:space:]]+(-[a-zA-Z0-9]*q[a-zA-Z0-9]*|--quiet|--silent)([[:space:]]|$)'
files=()
while IFS= read -r f; do
  case "$f" in *.sh) files+=("$f"); continue ;; esac
  IFS= read -r first < "$f" || first=""
  [[ $first =~ ^\#!.*(/|[[:space:]])(ba)?sh([[:space:]]|$) ]] && files+=("$f")
done < <(for d in "$root/tools" "$root/scripts"; do [ -d "$d" ] && find "$d" -type f; done | sort)
hits=""
[ "${#files[@]}" -gt 0 ] && hits=$(grep -HnE -- "$re" "${files[@]}" | grep -vE '^[^:]+:[0-9]+:[[:space:]]*#' | sort -t: -k1,1 -k2,2n)
if [ -n "$hits" ]; then
  n=$(grep -c . <<< "$hits")
  sed -E 's/^([^:]+:[0-9]+):[[:space:]]*/\1: /' <<< "$hits"
  echo "grepq_scan: FAIL $n hits: a pipe into grep -q under pipefail can SIGPIPE its writer (rc 141) and read as no match; use the here-string form grep -q PAT <<< \"\$x\" or read a file"
  exit 1
fi
echo "grepq_scan: ok (${#files[@]} shell files, 0 hits)"
