#!/usr/bin/env bash
# Case counter for the tools/*_selftest.sh scripts (Paul, 2026-10-23: a selftest step ran 0 cases and exited 0).
# Source it, then: st_init DIR (a fresh temp dir of the selftest), `ck;` at the start of every case line (one call per case executed,
# subshells included, because the count is a file), and st_done NAME WANT last: it FAILS unless exactly WANT cases ran.
# WANT is the number of cases the script runs; change it in the same commit as the cases.
st_init() { ST_CASES="$1/.selftest-cases"; : > "$ST_CASES" || { echo "selftest_cases: cannot create $ST_CASES"; exit 2; }; }
ck() { printf 'x\n' >> "${ST_CASES:?st_init was not called}"; }
st_done() { # st_done NAME WANT
  local n; n=$(grep -c x "${ST_CASES:?st_init was not called}" 2>/dev/null) || n=0
  [[ $2 =~ ^[0-9]+$ ]] && [ "$2" -gt 0 ] || { echo "$1 selftest FAIL: expected case count '$2' is not a positive integer"; exit 1; }
  [ "$n" = "$2" ] || { echo "$1 selftest FAIL: ran $n cases, expected $2 (a case was skipped, or added without updating the count)"; exit 1; }
}
