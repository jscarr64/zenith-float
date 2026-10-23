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
gate_dir=$(cd "$(dirname "$(readlink -f "$0")")" && pwd)
T=${T:-/tmp/zf-gate-logs}; mkdir -p "$T"
cd "${REPO:-$(cd "$(dirname "$0")/.." && pwd)}" || exit 2
tag=${1:?usage: tools/full_gate.sh TAG}; : > "$T/$tag.rc"
# shellcheck disable=SC1091
if [[ -f "$HOME/pop-logs/bin/poplog.sh" ]]; then source "$HOME/pop-logs/bin/poplog.sh"; else source "$gate_dir/../scripts/poplog.sh"; fi
poplog_init "${POPLOG_AGENT:-rusty}" "full_gate-$tag" "$PWD" "tools/full_gate.sh $tag"
echo "poplog: $POPLOG_FILE" | tee -a "$T/$tag.rc" >&2
# Pop scheduling (Paul, 2026-10-23): heavy gates share 3 counting-lock slots (/tmp/pop-gate.slot1..3, replaces the
# single flock /tmp/pop-gate.lock), 8 jobs each (24 of 32 cores). The slot is held on an fd until this script exits;
# waiting shows as WAITING in ~/pop-logs/RUNNING.tsv. Skipped when a caller already holds a slot (POPSLOT_HELD).
# CARGO_BUILD_JOBS / RUST_TEST_THREADS default to the slot value 6 (pop_slot_acquire) unless the caller sets them.
if [[ -z "${POPSLOT_HELD:-}" ]]; then pop_slot_acquire "${POPSLOT_SLOTS:-4}"; fi
export CARGO_BUILD_JOBS=${CARGO_BUILD_JOBS:-6} RUST_TEST_THREADS=${RUST_TEST_THREADS:-6}
TT="-- --test-threads=$RUST_TEST_THREADS"
# Preflight: every external tool a step needs must exist; a missing tool fails the gate loudly instead of
# letting a step's check read "command not found" as "no matches".
missing=()
for tool in bash timeout cargo rustup git grep python3 tail; do command -v "$tool" >/dev/null 2>&1 || missing+=("$tool"); done
cargo audit --version >/dev/null 2>&1 || missing+=("cargo-audit")
if [[ "${#missing[@]}" -gt 0 ]]; then
  echo "error: missing tools: ${missing[*]}" | tee -a "$T/$tag.rc" >&2
  echo "preflight=2" >> "$T/$tag.rc"; POPLOG_WARNINGS+=("preflight: missing ${missing[*]}"); poplog_finish 2; exit 2
fi
# pipefail: a failing command inside a step's pipeline fails the step
# Jeff's logging rule (2026-10-23): every step is logged via poplog (~/pop-logs/bin/poplog.sh, else the vendored
# scripts/poplog.sh): header (cmd, repo, SHA, rustc, CARGO_TARGET_DIR), per-step start/end/duration/rc/test counts,
# 5-min heartbeats, RUNNING.tsv row, PASS/FAIL footer. A test step reporting zero tests FAILs. Exit is nonzero on FAIL.
TEST_STEPS=" build-debug build-release features-random-serde oracle-fast oracle-deep mpfr "
run() {
  name=$1; shift; echo "start $name $(date +%T)" >> "$T/$tag.rc"
  log_step_start "$name"; ln -sf "$T/$tag-$name.log" "$(_pl_out)"  # nested: own out file, never the parent's
  timeout "${STEP_TIMEOUT:-20m}" bash -o pipefail -c "$*" > "$T/$tag-$name.log" 2>&1 &
  POPLOG_CHILD_PID=$!; echo "$POPLOG_CHILD_PID" > "$POPLOG_STATE_DIR/child"
  wait "$POPLOG_CHILD_PID"; rc=$?; POPLOG_CHILD_PID=""
  echo "$name=$rc" >> "$T/$tag.rc"
  tests=""; [[ "$TEST_STEPS" == *" $name "* ]] && tests=--tests
  log_step_end "$name" "$rc" "$T/$tag-$name.log" "$tests" || true
  echo "  output: $T/$tag-$name.log" >> "$POPLOG_FILE"
}
ORACLES="--test oracle --test specials_oracle --test methods_oracle --test accuracy_table --test elliptic_phi"
run build-debug "cargo test --workspace $TT"
# early compile of the mpfr-tests integration tests (6ab75ce prepublish broke here only at the last step)
run mpfr-tests-build "cargo test -p zenith-float-num --release --features mpfr-tests --no-run"
run build-release "cargo test --workspace --release $TT"
run build-nostd "cargo build --no-default-features && cargo build -p zenith-float-num --no-default-features"
run build-thumbv7em "cargo build --no-default-features --target thumbv7em-none-eabihf && cargo build -p zenith-float-num --no-default-features --target thumbv7em-none-eabihf"
run features-random-serde "cargo test --workspace --features random,serde $TT"
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
run unsafe-scan "bash tools/gate_unsafe_scan.sh"
run prepublish "bash scripts/zenith_prepublish.sh"
echo "done" >> "$T/$tag.rc"
poplog_finish
exit "$POPLOG_OVERALL"
