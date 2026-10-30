#!/usr/bin/env bash
# Prints the `--test NAME` arguments for the gate's test-i686 step: every integration test in the workspace
# except the oracle tests, which are left out by explicit name below.
# A test with required-features runs too, with exactly those features (`--features PKG/FEAT`); a required
# feature its package does not declare is an error (exit 1), never a skip.
# An excluded name that no longer exists is an error, so the list cannot go stale quietly.
set -euo pipefail
cd "$(dirname "$0")/.."
# Oracle tests: they check against the host MPFR or the mpmath oracle goldens (tests/oracle/goldens, 64-bit
# limbs; tests/oracle.rs does not build for 32-bit words). The host steps oracle-fast, oracle-deep and mpfr run them.
EXCLUDE=(oracle specials_oracle methods_oracle accuracy_table elliptic_phi)
meta=$(cargo metadata --no-deps --format-version 1 --offline)
names=$(jq -r '.packages[] | .name as $p | (.features // {} | keys) as $have | .targets[] | select(.kind == ["test"])
  | "\(.name) \([(.["required-features"] // [])[] as $f | if ($have | index([$f])) then "\($p)/\($f)" else "MISSING:\($p)/\($f)" end] | join(","))"' <<< "$meta" | sort -u)
for e in "${EXCLUDE[@]}"; do
  grep -q "^$e " <<< "$names" || { echo "i686_test_args: excluded oracle test '$e' is not a test target" >&2; exit 1; }
done
feats=()
while read -r n req; do
  [[ " ${EXCLUDE[*]} " == *" $n "* ]] && continue
  [[ $req != *MISSING:* ]] || { echo "i686_test_args: $n requires an undeclared feature ($req)" >&2; exit 1; }
  printf -- '--test %s ' "$n"
  [[ -z $req ]] || feats+=("${req//,/ }")
done <<< "$names"
[[ ${#feats[@]} -eq 0 ]] || printf -- '--features %s ' "${feats[@]}"
