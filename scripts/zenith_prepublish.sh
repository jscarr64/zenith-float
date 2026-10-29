#!/usr/bin/env bash
# Pre-publish gate. Exits 0 only when all 12 checks pass.
# Usage: CARGO_TARGET_DIR=/tmp/target-<agent>-<worktree> ./scripts/zenith_prepublish.sh
# Every cargo test runs through tools/fresh_test_run.sh (one target dir per worktree; the test binaries are checked
# before any test runs), so CARGO_TARGET_DIR must be set and owned by this worktree.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"
# shellcheck source=use_rust_1_99.sh
source "$root/scripts/use_rust_1_99.sh"

n=0
ok() { echo "ok  [$n/12] $*"; }
fail() { echo "error: [$n/12] $*" >&2; exit 1; }
step() { n=$((n + 1)); echo "== [$n/12] $* =="; }

# --all-features includes mpfr-tests. Check 1 covers std/random/serde (and
# compare-crate features). Check 2 is the MPFR oracle in release, single-threaded.
step "cargo test --all-features (workspace random,serde; mpfr-tests is [2/12])"
"$root/tools/fresh_test_run.sh" cargo test --workspace --features random,serde
ok "workspace tests with random,serde"

step "cargo test --features mpfr-tests"
if [[ "$(uname -s)" == Linux && "$(uname -m)" == x86_64 ]]; then
  "$root/tools/fresh_test_run.sh" cargo test -p zenith-float-num --features mpfr-tests --release -- --test-threads=1
  ok "mpfr-tests"
else
  fail "mpfr-tests require Linux x86_64"
fi

step "cargo build --no-default-features"
cargo build --no-default-features
cargo build -p zenith-float-num --no-default-features
ok "no_std host"

step "cargo build --no-default-features --target thumbv7em-none-eabihf"
rustup target add --toolchain 1.99 thumbv7em-none-eabihf >/dev/null
cargo build --no-default-features --target thumbv7em-none-eabihf
cargo build -p zenith-float-num --no-default-features --target thumbv7em-none-eabihf
ok "thumbv7em-none-eabihf"

step "scripts/compare-bench.sh --quick"
# output to COMPARE_OUT (the gate sets its log dir), else the target dir; never the tracked doc/compare-results.tsv
COMPARE_OUT="${COMPARE_OUT:-$CARGO_TARGET_DIR/compare-results.tsv}" bash "$root/scripts/compare-bench.sh" --quick
ok "compare-bench --quick"

step "kernel purity grep (no f32/f64, no arch intrinsics, no kernel unsafe)"
bash "$root/scripts/ci_no_intrinsics.sh"
ok "purity"

step "scripts/check_precision_docs.sh"
bash "$root/scripts/check_precision_docs.sh"
ok "precision rustdoc"

step "CAPABILITIES version matches Cargo.toml"
python3 "$root/scripts/prepublish_lib.py" version
ok "version"

step "scripts/check_leaves.sh"
bash "$root/scripts/check_leaves.sh"
ok "expr!/cexpr! leaves"

step "no <!-- verify --> in doc/"
python3 "$root/scripts/prepublish_lib.py" verify
ok "no verify markers"

step "proptest PROPTEST_CASES=1000"
python3 - <<'PY'
from pathlib import Path
text = Path("zenith-float-num/src/defs.rs").read_text()
if "pub const PROPTEST_CASES: u32 = 1000;" not in text:
    raise SystemExit("error: PROPTEST_CASES is not 1000 in defs.rs")
print("PROPTEST_CASES = 1000")
PY
"$root/tools/fresh_test_run.sh" cargo test -p zenith-float-num --test proptest_props -- --test-threads=1
ok "proptest"

step "no TODO/FIXME/HACK; unwrap() has a proof comment"
python3 "$root/scripts/prepublish_lib.py" hygiene
ok "source hygiene"

echo "zenith-float prepublish ok (12/12)"
