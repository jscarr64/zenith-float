#!/usr/bin/env bash
# Fails when any workspace test target (lib/bin unit harness, integration test,
# or Doc-tests) holds 0 tests. A target with nothing to run must be switched off
# in its Cargo.toml (`test = false` / `doctest = false`). Lists tests only; runs none.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"
out=$(mktemp)
trap 'rm -f "$out"' EXIT
cargo test --workspace -- --list > "$out" 2>&1 || { cat "$out" >&2; echo "check_zero_test_targets: cargo test --list failed" >&2; exit 1; }
awk '/^ *(Running|Doc-tests) /{sub(/^ +/,""); h=$0; next}
     /^[0-9]+ tests?, [0-9]+ benchmarks?$/ { if (h!="") { t++; if ($1==0) { z++; print "ZERO-TEST TARGET: " h } } h="" }
     END { printf "check_zero_test_targets: targets=%d zero=%d\n", t, z; exit (t==0 || z>0) }' "$out"
