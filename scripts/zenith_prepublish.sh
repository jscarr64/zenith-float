#!/usr/bin/env bash
# Pre-publish gate. Exits 0 only when all 11 checks pass (thumbv7em step dropped, 64-bit only).
# Usage: CARGO_TARGET_DIR=/tmp/target-<agent>-<worktree> ./scripts/zenith_prepublish.sh
# Every cargo test runs through tools/fresh_test_run.sh (one target dir per worktree; the test binaries are checked
# before any test runs), so CARGO_TARGET_DIR must be set and owned by this worktree.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"
# shellcheck source=use_rust_1_99.sh
source "$root/scripts/use_rust_1_99.sh"

n=0
ok() { echo "ok  [$n/11] $*"; }
fail() { echo "error: [$n/11] $*" >&2; exit 1; }
step() { n=$((n + 1)); echo "== [$n/11] $* =="; }

# --all-features includes mpfr-tests. Check 1 covers std/random/serde (and
# compare-crate features). Check 2 is the MPFR oracle in release, single-threaded.
step "cargo test --all-features (workspace random,serde; mpfr-tests is [2/11])"
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

# No 32-bit target build (thumbv7em dropped): 2.0 is 64-bit only (P6; Jeff 2026-11-10 19:15 CT).

step "scripts/compare-bench.sh --quick"
# The gate never rewrites a tracked file: the comparison goes to scratch under CARGO_TARGET_DIR (outside the worktree),
# not to the tracked doc/compare-results.tsv (refresh that by hand with scripts/compare-bench.sh and commit it).
mkdir -p "${CARGO_TARGET_DIR:?}/prepublish"
case "$(cd "$CARGO_TARGET_DIR/prepublish" && pwd -P)/" in "$(cd "$root" && pwd -P)/"*) fail "CARGO_TARGET_DIR is inside the worktree; the compare output must go outside it" ;; esac
COMPARE_OUT="$CARGO_TARGET_DIR/prepublish/compare-results.tsv" bash "$root/scripts/compare-bench.sh" --quick
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

echo "zenith-float prepublish ok (11/11)"
