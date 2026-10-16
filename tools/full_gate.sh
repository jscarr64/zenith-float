#!/usr/bin/env bash
# Canonical full gate. usage: tools/full_gate.sh TAG
# Every step runs under `timeout` (STEP_TIMEOUT, default 20m); the rc of each step is written to $T/TAG.rc and its
# output to $T/TAG-<step>.log (T defaults to /tmp/zf-gate-logs). REPO selects the checkout (default: this one);
# CARGO_TARGET_DIR is honoured.
#
# Clone ratchet "in progress" mode: `rules-scan` scans the added lines since the ratchet started (RATCHET_BASE,
# default c74e834^ = the parent of the first commit with rules_scan and the baseline) and enforces the baseline
# count; `no-clone-progress` prints the remaining no_clone.sh violations (a tracked count, it never fails the gate).
# The zero-clone goal is unchanged: the baseline in tools/clone_count_baseline.txt only goes down.
# `rules-scan` also enforces (via rules_scan.py count checks) the helper-copy baseline (tools/helper_copy_count_baseline.txt),
# the exception allowlist (stale entries fail, tools/exception_count_baseline.txt) and the library-scoped ref ratchets;
# and the neg_via_mul use-count baseline (tools/neg_via_mul_count_baseline.txt, all scanned code incl. tests; it only goes down).
# the clone count includes test code. `no-clone-progress` also prints the helper-copy count. All baselines only go down.
set -u
T=${T:-/tmp/zf-gate-logs}; mkdir -p "$T"
cd "${REPO:-$(cd "$(dirname "$0")/.." && pwd)}" || exit 2
tag=${1:?usage: tools/full_gate.sh TAG}; : > "$T/$tag.rc"
run() { name=$1; shift; echo "start $name $(date +%T)" >> "$T/$tag.rc"; timeout "${STEP_TIMEOUT:-20m}" bash -c "$*" > "$T/$tag-$name.log" 2>&1; echo "$name=$?" >> "$T/$tag.rc"; }
ORACLES="--test oracle --test specials_oracle --test methods_oracle --test accuracy_table --test elliptic_phi"
run build-debug "cargo test --workspace"
run build-release "cargo test --workspace --release"
run build-nostd "cargo build --no-default-features && cargo build -p zenith-float-num --no-default-features"
run build-thumbv7em "cargo build --no-default-features --target thumbv7em-none-eabihf && cargo build -p zenith-float-num --no-default-features --target thumbv7em-none-eabihf"
run features-random-serde "cargo test --workspace --features random,serde"
run oracle-fast "cargo test -p zenith-float-num --release --features mpfr-tests $ORACLES -- --test-threads=1"
run oracle-deep "ORACLE_TIER=deep cargo test -p zenith-float-num --release --features mpfr-tests $ORACLES -- --test-threads=1"
run mpfr "cargo test -p zenith-float-num --release --features mpfr-tests -- --test-threads=1"
# ratchet mode: the clone ratchet started at c74e834 (first commit with rules_scan + baseline); older added lines are tracked by no-clone-progress, not flagged here
run rules-scan "tools/rules_scan_selftest.sh && tools/rules_scan.sh ${RATCHET_BASE:-c74e834^}..HEAD"
# informational: remaining violations are a tracked count (must only go down; hard fail is for the final step)
run no-clone-progress "tools/no_clone_selftest.sh && { tools/no_clone.sh | tail -1; tools/no_clone.sh --helpers | tail -1; true; }"
run fmt "cargo fmt --all -- --check"
run clippy-x86_64 "cargo clippy --workspace --all-targets --all-features -- -D warnings"
run clippy-i686 "cargo clippy -p zenith-float-num --all-targets --target i686-unknown-linux-musl -- -D warnings"
run audit "cargo audit --deny warnings"
run no-intrinsics "bash scripts/ci_no_intrinsics.sh"
run hex-32bit "bash scripts/ci_hex_32bit.sh"
run unsafe-scan "! grep -rn --include='*.rs' -E '\bunsafe\b' zenith-float-num/src zenith-float-macro/src src 2>/dev/null | grep -v '^\s*//' | grep -vE ':\s*(//|\*)'"
run prepublish "bash scripts/zenith_prepublish.sh"
echo done >> "$T/$tag.rc"
