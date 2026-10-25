#!/usr/bin/env bash
# golden_sums_check.sh <goldens-dir> <required-file>...  Verify <dir>/SHA256SUMS with `sha256sum -c`, failing closed:
# a missing or empty SHA256SUMS, a required golden not listed in it, or any hash mismatch / missing listed file fails.
# Regenerating a pinned golden means updating SHA256SUMS in the same commit (review sees both).
set -u
[ $# -ge 2 ] || { echo "usage: golden_sums_check.sh <goldens-dir> <required-file>..." >&2; exit 2; }
dir=$1; shift; sums="$dir/SHA256SUMS"
[ -f "$sums" ] || { echo "golden_sums: FAIL $sums missing"; exit 1; }
[ -n "$(grep -v '^[[:space:]]*$' "$sums")" ] || { echo "golden_sums: FAIL $sums is empty"; exit 1; }
for f in "$@"; do
  grep -Eq "^[0-9a-f]{64}  \*?$f\$" "$sums" || { echo "golden_sums: FAIL $f not listed in $sums"; exit 1; }
done
( cd "$dir" && sha256sum --strict -c SHA256SUMS ) || { echo "golden_sums: FAIL hash mismatch in $sums"; exit 1; }
echo "golden_sums: ok ($sums, $# required)"
