#!/usr/bin/env bash
# Coverage check for a host test step's cargo test log (file $2, step label $1): FAILs (exit 1) when a test
# binary ran 0 tests (none, all filtered out, all ignored), a Running line has no summary, the log has no
# binaries, or a binary named after --expect did not run at all (a test skipped for a missing feature).
# Binaries that are empty by construction are listed in ZERO_OK with the source fact checked each run, so a
set -euo pipefail
label=${1:?usage: tools/host_test_counts.sh LABEL CARGO_TEST_LOG [--expect NAME...]}
out=${2:?usage: tools/host_test_counts.sh LABEL CARGO_TEST_LOG [--expect NAME...]}
shift 2
expect=" "
if [[ $# -gt 0 ]]; then
  [[ $1 == --expect ]] || { echo "$label coverage: unknown argument $1" >&2; exit 2; }
  shift; expect+="$* "
fi
root=${TEST_COUNTS_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}
# name|file|fact: `notest` = no #[test] in that file or dir; otherwise a line that must be in the file
ZERO_OK=(
  "zenith_float|src|notest"
  "hex_golds|src/bin/hex_golds.rs|notest"
  "integration_tests|zenith-float-num/tests/integration_tests.rs|#[cfg(feature = \"mpfr-tests\")]"
)
ok=" "
for e in "${ZERO_OK[@]}"; do
  IFS='|' read -r n f fact <<< "$e"
  [[ -e "$root/$f" ]] || { echo "$label coverage: zero-test entry $n: $f does not exist" >&2; exit 2; }
  case $fact in
    notest) if grep -rqF '#[test]' "$root/$f"; then echo "$label coverage: zero-test entry $n: $f now has #[test]" >&2; exit 2; fi;;
    *) grep -qxF "$fact" "$root/$f" || { echo "$label coverage: zero-test entry $n: '$fact' not in $f" >&2; exit 2; };;
  esac
  ok+="$n "
done
awk -v ok="$ok" -v expect="$expect" -v label="$label" '
  /^ +Running / { if (b != "") { bad = bad " " b "(no-summary)" } b = $NF; gsub(/[()]/, "", b); sub(/.*\//, "", b); sub(/-[^-]*$/, "", b); nb++; seen[b] = 1 }
  /^test result:/ { p = f = 0; for (k = 1; k <= NF; k++) { if ($k ~ /^passed;?/) p = $(k-1); if ($k ~ /^failed;?/) f = $(k-1) }
                    if (b != "" && p + f == 0 && index(ok, " " b " ") == 0) bad = bad " " b "(0-run)"; b = "" }
  END { if (b != "") bad = bad " " b "(no-summary)"
        if (nb == 0) bad = bad " no-binaries"
        n = split(expect, want, " "); for (k = 1; k <= n; k++) if (!(want[k] in seen)) bad = bad " " want[k] "(feature-skipped)"
        if (bad != "") { print label " coverage: FAIL" bad; exit 1 }
        print label " coverage: ok " nb " binaries" }' "$out"
