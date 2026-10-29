#!/usr/bin/env bash
# Prints the `--test NAME` arguments for the gate's test-i686 step: every integration test in the workspace
# except the oracle tests, which are left out by explicit name below.
# A test with required-features is left out the same way plain `cargo test` leaves it out (default features).
# An excluded name that no longer exists is an error, so the list cannot go stale quietly.
set -euo pipefail
cd "$(dirname "$0")/.."
# Oracle tests: they check against the host MPFR or the mpmath oracle goldens (tests/oracle/goldens, 64-bit
# limbs; tests/oracle.rs does not build for 32-bit words). The host steps oracle-fast, oracle-deep and mpfr run them.
EXCLUDE=(oracle specials_oracle methods_oracle accuracy_table elliptic_phi)
names=$(cargo metadata --no-deps --format-version 1 --offline | jq -r '.packages[].targets[] | select(.kind == ["test"]) | "\(.name) \(.["required-features"] | length)"' | sort -u)
for e in "${EXCLUDE[@]}"; do
  grep -q "^$e " <<< "$names" || { echo "i686_test_args: excluded oracle test '$e' is not a test target" >&2; exit 1; }
done
while read -r n req; do
  [[ " ${EXCLUDE[*]} " == *" $n "* ]] && continue
  [[ $req -eq 0 ]] || { echo "i686_test_args: $n needs features, skipped as in plain cargo test" >&2; continue; }
  printf -- '--test %s ' "$n"
done <<< "$names"
