#!/usr/bin/env bash
# Pre-publish gate. Exits 0 only when all 12 checks pass.
# Usage: ./scripts/zenith_prepublish.sh
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"

n=0
ok() { echo "ok  [$n/12] $*"; }
fail() { echo "error: [$n/12] $*" >&2; exit 1; }
step() { n=$((n + 1)); echo "== [$n/12] $* =="; }

# --all-features includes mpfr-tests. Check 1 covers std/random/serde (and
# compare-crate features). Check 2 is the MPFR oracle in release, single-threaded.
step "cargo test --all-features (workspace random,serde; mpfr-tests is [2/12])"
cargo test --workspace --features random,serde
ok "workspace tests with random,serde"

step "cargo test --features mpfr-tests"
if [[ "$(uname -s)" == Linux && "$(uname -m)" == x86_64 ]]; then
  cargo test -p zenith-float-num --features mpfr-tests --release -- --test-threads=1
  ok "mpfr-tests"
else
  fail "mpfr-tests require Linux x86_64"
fi

step "cargo build --no-default-features"
cargo build --no-default-features
cargo build -p zenith-float-num --no-default-features
ok "no_std host"

step "cargo build --no-default-features --target thumbv7em-none-eabihf"
cargo build --no-default-features --target thumbv7em-none-eabihf
cargo build -p zenith-float-num --no-default-features --target thumbv7em-none-eabihf
ok "thumbv7em-none-eabihf"

step "scripts/compare-bench.sh --quick"
bash "$root/scripts/compare-bench.sh" --quick
ok "compare-bench --quick"

step "kernel purity grep (no f32/f64 tokens in zenith-float-num/src)"
if rg -n 'f32|f64' --glob '*.rs' zenith-float-num/src; then
  fail "hardware IEEE type tokens in zenith-float-num/src"
fi
ok "purity"

step "scripts/check_precision_docs.sh"
bash "$root/scripts/check_precision_docs.sh"
ok "precision rustdoc"

step "CAPABILITIES version/date match Cargo.toml"
python3 "$root/scripts/prepublish_lib.py" version
ok "version"

step "scripts/check_leaves.sh"
bash "$root/scripts/check_leaves.sh"
ok "expr!/cexpr! leaves"

step "no <!-- verify --> in CAPABILITIES §24"
python3 "$root/scripts/prepublish_lib.py" verify
ok "comparison table verified"

step "proptest PROPTEST_CASES=1000"
python3 - <<'PY'
from pathlib import Path
text = Path("zenith-float-num/src/defs.rs").read_text()
if "pub const PROPTEST_CASES: u32 = 1000;" not in text:
    raise SystemExit("error: PROPTEST_CASES is not 1000 in defs.rs")
print("PROPTEST_CASES = 1000")
PY
cargo test -p zenith-float-num --test proptest_props -- --test-threads=1
ok "proptest"

step "no TODO/FIXME/HACK; unwrap() has a proof comment"
python3 "$root/scripts/prepublish_lib.py" hygiene
ok "source hygiene"

echo "zenith-float prepublish ok (12/12)"
