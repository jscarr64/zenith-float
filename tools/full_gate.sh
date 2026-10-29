#!/usr/bin/env bash
# Canonical full gate. usage: tools/full_gate.sh TAG
# Every step runs under `timeout` (STEP_TIMEOUT, default 20m; the mpfr step uses MPFR_STEP_TIMEOUT, default 40m); the rc of each step is written to $T/TAG.rc and its
# output to $T/TAG-<step>.log (T, else ZF_GATE_LOG_DIR, else $HOME/rusty-work/logs/zf-gate/<sha>-<time>; never /tmp). REPO selects the checkout (default: this one).
# CARGO_TARGET_DIR is REQUIRED and belongs to one worktree (Paul, 2026-10-23: /tmp/target-<agent>-<worktree>). The first
# step, target-dir, runs the tools/fresh_test_run.sh selftest and --claim: it fails if CARGO_TARGET_DIR is unset or owned
# by another worktree (marker $CARGO_TARGET_DIR/.worktree-owner). Every cargo test step (mpfr-tests-build included) runs
# through tools/fresh_test_run.sh: build with --no-run, then fail on zero test binaries or one outside CARGO_TARGET_DIR,
# before any test runs; only then the tests run.
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
T=${T:-${ZF_GATE_LOG_DIR:-$HOME/rusty-work/logs/zf-gate/$(git -C "$gate_dir" rev-parse --short HEAD)-$(date +%Y%m%d-%H%M%S)}}; mkdir -p "$T"
cd "${REPO:-$(cd "$(dirname "$0")/.." && pwd)}" || exit 2
tag=${1:?usage: tools/full_gate.sh TAG}; : > "$T/$tag.rc"
# shellcheck disable=SC1091
if [[ -f "$HOME/pop-logs/bin/poplog.sh" ]]; then source "$HOME/pop-logs/bin/poplog.sh"; else source "$gate_dir/../scripts/poplog.sh"; fi
poplog_init "${POPLOG_AGENT:-rusty}" "full_gate-$tag" "$PWD" "tools/full_gate.sh $tag"
echo "poplog: $POPLOG_FILE" | tee -a "$T/$tag.rc" >&2
# Pop scheduling (Paul, 2026-10-23): heavy gates share 3 counting-lock slots (/tmp/pop-gate.slot1..3, replaces the
# single flock /tmp/pop-gate.lock), 8 jobs each (24 of 32 cores). The slot is held on an fd until this script exits;
# waiting shows as WAITING in ~/pop-logs/RUNNING.tsv. Skipped when a caller already holds a slot (POPSLOT_HELD).
# CARGO_BUILD_JOBS / RUST_TEST_THREADS default to the slot value POPSLOT_JOBS (4) (pop_slot_acquire) unless the caller sets them.
if [[ -z "${POPSLOT_HELD:-}" ]]; then pop_slot_acquire "${POPSLOT_SLOTS:-8}"; fi
export CARGO_BUILD_JOBS=${CARGO_BUILD_JOBS:-${POPSLOT_JOBS:-4}} RUST_TEST_THREADS=${RUST_TEST_THREADS:-${POPSLOT_JOBS:-4}}
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
TEST_STEPS=" build-debug build-release test-i686-lib features-random-serde oracle-fast oracle-deep mpfr "
# GATE_STEPS (space- or comma-separated): run only the steps with exactly these names; an unknown name fails the gate.
failed=(); ran=()
run() {
  name=$1; shift; if [[ -n ${GATE_STEPS:-} && " ${GATE_STEPS//,/ } " != *" $name "* ]]; then return 0; fi; ran+=("$name")
  echo "start $name $(date +%T)" >> "$T/$tag.rc"
  rm -f "$(_pl_out)"  # the last step left it a symlink to its log; log_step_start truncates it
  log_step_start "$name"; ln -sf "$T/$tag-$name.log" "$(_pl_out)"  # nested: own out file, never the parent's
  to=${STEP_TIMEOUT:-20m}; [[ $name == mpfr ]] && to=${MPFR_STEP_TIMEOUT:-40m}  # mpfr: integration_tests alone is ~1090 s
  timeout "$to" bash -o pipefail -c "$*" > "$T/$tag-$name.log" 2>&1 &
  POPLOG_CHILD_PID=$!; echo "$POPLOG_CHILD_PID" > "$POPLOG_STATE_DIR/child"
  wait "$POPLOG_CHILD_PID"; rc=$?; POPLOG_CHILD_PID=""
  echo "$name=$rc" >> "$T/$tag.rc"
  tests=""; [[ "$TEST_STEPS" == *" $name "* ]] && tests=--tests
  log_step_end "$name" "$rc" "$T/$tag-$name.log" "$tests" || failed+=("$name")
  echo "  output: $T/$tag-$name.log" >> "$POPLOG_FILE"
  # forward the step's test-binary and summary lines to stdout (full output stays in $T) so a wrapping
  # fresh_test_run.sh (CMD... mode) sees which binaries ran and checks they lie inside CARGO_TARGET_DIR
  grep -E '^ *Running .*\(/[^)]*\)$|^test result: ' "$T/$tag-$name.log" || true
}
ORACLES="--test oracle --test specials_oracle --test methods_oracle --test accuracy_table --test elliptic_phi"
git status --porcelain > "$T/$tag.pre-status"
run target-dir "tools/fresh_test_run_selftest.sh && tools/fresh_test_run.sh --claim"
# goldens pinned by SHA256SUMS next to them (fails closed: missing/empty sums, unlisted golden, drift)
run golden-sums "tools/golden_sums_selftest.sh && tools/golden_sums_check.sh zenith-float-num/tests/oracle/goldens borrow_ops.jsonl"
run build-debug "tools/fresh_test_run.sh cargo test --workspace $TT"
# early compile of the mpfr-tests integration tests (6ab75ce prepublish broke here only at the last step)
run mpfr-tests-build "tools/fresh_test_run.sh cargo test -p zenith-float-num --release --features mpfr-tests --no-run"
run build-release "tools/fresh_test_run.sh cargo test --workspace --release $TT"
# i686 lib tests (merge-blocking from 9d66035): the 32-bit target clippy-i686 checks, run natively (static musl), own target subdir
run test-i686-lib "CARGO_TARGET_DIR=$CARGO_TARGET_DIR/i686 tools/fresh_test_run.sh cargo test -p zenith-float-num --lib --target i686-unknown-linux-musl $TT"
run build-nostd "cargo build --no-default-features && cargo build -p zenith-float-num --no-default-features"
run build-thumbv7em "cargo build --no-default-features --target thumbv7em-none-eabihf && cargo build -p zenith-float-num --no-default-features --target thumbv7em-none-eabihf"
run features-random-serde "tools/fresh_test_run.sh cargo test --workspace --features random,serde $TT"
run oracle-fast "tools/fresh_test_run.sh cargo test -p zenith-float-num --release --features mpfr-tests $ORACLES -- --test-threads=1"
run oracle-deep "ORACLE_TIER=deep tools/fresh_test_run.sh cargo test -p zenith-float-num --release --features mpfr-tests $ORACLES -- --test-threads=1"
# MPFR_TIER selects the MPFR tier (tests/mpfr/common.rs): the gate runs `full` (every run count, the
# 100032-bit tests) unless the caller sets MPFR_TIER=fast.
run mpfr "MPFR_TIER=${MPFR_TIER:-full} tools/fresh_test_run.sh cargo test -p zenith-float-num --release --features mpfr-tests -- --test-threads=1"
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
# compare-bench output goes to the run log dir, never the tracked doc/compare-results.tsv
run prepublish "COMPARE_OUT=$T/compare-results.tsv bash scripts/zenith_prepublish.sh"
# the gate must leave the worktree as it found it: compare against the status taken before the first step
run tree-clean "git status --porcelain | diff -u $T/$tag.pre-status - && echo 'tree unchanged'"
for want in ${GATE_STEPS:+${GATE_STEPS//,/ }}; do [[ " ${ran[*]} " == *" $want "* ]] || { failed+=("GATE_STEPS:unknown-step:$want"); POPLOG_OVERALL=1; }; done
# GATE_STEPS set but matching no step (e.g. "," or " "): nothing ran, which is a FAIL, never a PASS
if [[ -n ${GATE_STEPS:-} && ${#ran[@]} -eq 0 ]]; then failed+=("GATE_STEPS='$GATE_STEPS': no step matched"); POPLOG_OVERALL=1; fi
# a selective run is PARTIAL (stdout, .rc, poplog footer) so it can never be cited as a merge gate
part=""; [[ -n ${GATE_STEPS:-} ]] && part="PARTIAL "
if [[ ${#failed[@]} -gt 0 ]]; then res="RESULT: ${part}FAIL failed steps: ${failed[*]}"; else res="RESULT: ${part}PASS steps: ${ran[*]}"; fi
echo "$res" | tee -a "$T/$tag.rc"; echo "gate $res" >> "$POPLOG_FILE"
[[ -n $part ]] && POPLOG_WARNINGS+=("$res (selective GATE_STEPS run, not a merge gate)")
echo "done" >> "$T/$tag.rc"
poplog_finish
exit "$POPLOG_OVERALL"
