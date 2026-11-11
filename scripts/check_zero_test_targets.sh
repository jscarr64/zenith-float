#!/usr/bin/env bash
# Fails when any workspace test target (lib/bin unit harness, integration test,
# or Doc-tests) holds 0 tests. A target with nothing to run must be switched off
# in its Cargo.toml (`test = false` / `doctest = false`). Lists tests only; runs none.
#
# Fail-closed parsing (Spence review of PR #18): the old awk passed silently on
# every one of these, so each is now a hard failure:
#   - a header with no following count line (dangling header, including a header
#     followed by an error line and then nothing);
#   - a count line with no preceding header (orphan count);
#   - a target that vanishes entirely, caught by asserting the expected target
#     count for each pass (update the count in the same PR that adds or removes
#     a target; a count change without a target change is a failure, not a pass).
# While a header is pending, only test-name lines (": test") and blank lines may
# appear before its count; anything else (an error line, a second header, build
# chatter) fails, so an error cannot be swallowed between header and count.
#
# Expected counts on this branch (default features, no mpfr-tests): the
# workspace list has 66 targets. The mpfr-tests pass below is only for
# zenith-float-num, which has 51 targets with mpfr-tests on x86_64-linux,
# including integration_tests.
#
# The counts moved 52 -> 53 (46 -> 47 for mpfr-tests) when main brought in
# the crlf_inputs.rs test target from zf #19 (merged 2026-11-03), then
# 53 -> 54 (47 -> 48) with mp16_modf.rs, then 54 -> 55 (48 -> 49) with
# mp18_int_getters.rs on this stacked land tip, then 55 -> 56 (49 -> 50)
# with fp16_bf16.rs (zenith-float 2.0 Slice 1, Ieee16 / Bfloat16), then 56 -> 65 (mpfr-tests
# unchanged at 50) with the bound-layer fold (fold PR 1): zenith-float-bound's lib, its doc-tests and
# its 6 integration test targets (array, ball, enclosure, harness, mag, types), plus the facade's
# tests/fold_pr1_facade.rs, then 65 -> 66 (50 -> 51) with zenith-float-num/tests/top_precision_edges.rs
# (8da300e, a valid precision near the top is never InvalidArgument); the count moves in the
# same PR that adds or removes a target; a count change without a target
# change is a failure, not a pass).
# integration_tests minor (Spence review of PR #18): that target's tests live
# behind cfg(all(target_arch = "x86_64", target_os = "linux")) in
# zenith-float-num/tests/integration_tests.rs:7, so with mpfr-tests it is 0
# tests off x86_64-linux, and the default pass never sees it at all (it only
# builds under mpfr-tests). The mpfr pass below therefore runs only on
# x86_64-linux, where integration_tests actually has tests and is checked;
# off that platform it prints an explicit skip naming the cfg gate, instead of
# silently never checking it. The mpfr pass is a different feature set from the
# default pass, so its list is not a duplicate rebuild of the default list.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"

EXPECTED_TARGETS=66
EXPECTED_MPFR_TARGETS=51

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

check_list() {
  local label="$1" expected="$2"
  shift 2
  local out="$tmpdir/list-$label.txt"
  # A list file override exists only so the parser's fail-closed cases can be
  # demonstrated without a real workspace: gates never set it.
  if [[ -n "${ZERO_TEST_LIST_FILE:-}" ]]; then
    cp "$ZERO_TEST_LIST_FILE" "$out"
  elif ! cargo test "$@" -- --list > "$out" 2>&1; then
    cat "$out" >&2
    echo "check_zero_test_targets: cargo test --list failed ($label)" >&2
    return 1
  fi
  awk -v expected="$expected" -v label="$label" '
    /^ *(Running|Doc-tests) / {
      if (pending != "") {
        print "FAIL (" label "): header without a count: " pending > "/dev/stderr"
        fail = 1
      }
      pending = $0; sub(/^ +/, "", pending); next
    }
    /^[0-9]+ tests?, [0-9]+ benchmarks?$/ {
      if (pending == "") {
        print "FAIL (" label "): count with no header: " $0 > "/dev/stderr"
        fail = 1; next
      }
      total++
      if ($1 == 0) { zero++; print "ZERO-TEST TARGET (" label "): " pending }
      pending = ""; next
    }
    {
      if (pending != "") {
        if ($0 ~ /: test$/ || $0 == "") next
        print "FAIL (" label "): unexpected line after header " pending ": " $0 > "/dev/stderr"
        fail = 1; pending = ""
      }
      next
    }
    END {
      if (pending != "") {
        print "FAIL (" label "): header without a count: " pending > "/dev/stderr"
        fail = 1
      }
      printf "check_zero_test_targets (%s): targets=%d zero=%d\n", label, total, zero
      if (fail || total == 0 || zero > 0) exit 1
      if (total != expected) {
        print "FAIL (" label "): expected " expected " targets, got " total > "/dev/stderr"
        exit 1
      }
    }' "$out"
}

check_list "workspace" "$EXPECTED_TARGETS" --workspace

if [[ "$(uname -s)" == Linux && "$(uname -m)" == "x86_64" ]]; then
  check_list "mpfr-tests" "$EXPECTED_MPFR_TARGETS" -p zenith-float-num --features mpfr-tests
else
  echo "check_zero_test_targets: mpfr-tests pass skipped off x86_64-linux (integration_tests is 0 there by cfg at zenith-float-num/tests/integration_tests.rs:7) - not checked on this platform"
fi
