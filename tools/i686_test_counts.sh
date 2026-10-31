#!/usr/bin/env bash
# Coverage check for the gate's test-i686 output (cargo test log, file $1): FAILs (exit 1) when a test binary
# ran no test (0 passed and 0 failed: no tests, all filtered out, or all ignored), when a `Running` binary has no
# `test result:` line (crash, abort), or when no binary ran at all. The only binaries allowed to run 0 tests are
# listed in ZERO_OK with the source fact that makes them empty on i686; each fact is checked here, so a listed
# entry whose reason is gone is an input error (exit 2), not a quiet pass.
set -euo pipefail
out=${1:?usage: tools/i686_test_counts.sh CARGO_TEST_LOG}
root=${I686_COUNTS_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}
# name|file|fact: `notest` = no #[test] in that file or dir; otherwise a line that must be in the file
ZERO_OK=(
  "zenith_float|src|notest"
  "hex_golds|src/bin/hex_golds.rs|notest"
  "integration_tests|zenith-float-num/tests/integration_tests.rs|#[cfg(feature = \"mpfr-tests\")]"
  "r13_div_exact|zenith-float-num/tests/r13_div_exact.rs|#![cfg(not(target_pointer_width = \"32\"))]"
)
ok=" "
for e in "${ZERO_OK[@]}"; do
  IFS='|' read -r n f fact <<< "$e"
  [[ -e "$root/$f" ]] || { echo "i686 coverage: zero-test entry $n: $f does not exist" >&2; exit 2; }
  case $fact in
    notest) if grep -rqF '#[test]' "$root/$f"; then echo "i686 coverage: zero-test entry $n: $f now has #[test]" >&2; exit 2; fi;;
    *) grep -qxF "$fact" "$root/$f" || { echo "i686 coverage: zero-test entry $n: '$fact' not in $f" >&2; exit 2; };;
  esac
  ok+="$n "
done
awk -v ok="$ok" '
  /^ +Running / { if (b != "") { bad = bad " " b "(no-summary)" } b = $NF; gsub(/[()]/, "", b); sub(/.*\//, "", b); sub(/-[^-]*$/, "", b); nb++ }
  /^test result:/ { p = f = 0; for (k = 1; k <= NF; k++) { if ($k ~ /^passed;?/) p = $(k-1); if ($k ~ /^failed;?/) f = $(k-1) }
                    if (b != "" && p + f == 0 && index(ok, " " b " ") == 0) bad = bad " " b "(0-run)"; b = "" }
  END { if (b != "") bad = bad " " b "(no-summary)"
        if (nb == 0) bad = bad " no-binaries"
        if (bad != "") { print "i686 coverage: FAIL" bad; exit 1 }
        print "i686 coverage: ok " nb " binaries" }' "$out"
