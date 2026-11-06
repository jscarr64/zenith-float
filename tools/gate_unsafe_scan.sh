#!/usr/bin/env bash
# Fail when `unsafe` appears in library sources (comment lines excepted).
# grep exits 0 on a match, 1 on no match, 2+ on an error; anything but 0/1 fails the step.
set -euo pipefail
command -v grep >/dev/null 2>&1 || { echo "error: grep not found" >&2; exit 2; }
dirs=()
for d in zenith-float-num/src zenith-float-macro/src src; do
  [[ -d "$d" ]] && dirs+=("$d")
done
[[ "${#dirs[@]}" -gt 0 ]] || { echo "error: no source dirs to scan" >&2; exit 2; }
status=0
hits="$(grep -rn --include='*.rs' -E '\bunsafe\b' "${dirs[@]}")" || status=$?
if [[ "$status" -ge 2 ]]; then
  echo "error: grep failed with status $status" >&2
  exit 2
fi
[[ "$status" -eq 1 ]] && { echo "unsafe-scan: ok"; exit 0; }
# drop hits whose source text is a comment: file:line:<ws>// or file:line:<ws>*
status=0
code="$(printf '%s\n' "$hits" | grep -vE '^[^:]+:[0-9]+:[[:space:]]*(//|\*)')" || status=$?
if [[ "$status" -ge 2 ]]; then
  echo "error: grep failed with status $status" >&2
  exit 2
fi
if [[ "$status" -eq 0 ]]; then
  printf '%s\n' "$code"
  echo "error: unsafe in library sources" >&2
  exit 1
fi
echo "unsafe-scan: ok"
