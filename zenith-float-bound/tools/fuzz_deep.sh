#!/usr/bin/env bash
# S4 deep fuzz run (TEST_PLAN: the gate runs 1e6 seeded cases; this runs many more, OUTSIDE the gate).
# Runs every target in tools/fuzz_deep_targets.txt (FUZZ_DEEP_TARGETS overrides the file; default fuzz_mag::fuzz_mag_gate) in release, each
# split into independent shards (one process each, one distinct base seed each), with ZB_FUZZ_SEED / ZB_FUZZ_CASES set per shard. A shard
# passes only if its log says '1 passed' (a renamed or missing test runs 0 tests under --exact and must fail). Every failure prints its
# replay line; rerun one case with
#   tools/fuzz_deep.sh --replay <seed hex>:<case> [--target MODULE::TEST]     (default target: the first in the list)
# Deep runs take `flock /tmp/pop-deep.lock` (FUZZ_DEEP_LOCK overrides the path, FUZZ_DEEP_LOCK_WAIT=<s> bounds the wait, default 3600),
# never the gate lock. Timeouts: FUZZ_DEEP_BUILD_TIMEOUT (default 30m) on the release build, FUZZ_DEEP_SHARD_TIMEOUT (default 6h) per shard.
# The case count comes from ZB_FUZZ_CASES or --cases (the flag wins): a plain integer or <digits>e<digits> (1e8); garbage is refused.
#   tools/fuzz_deep.sh [--cases N|NeK] [--shards K] [--seed HEX] [--out DIR]     defaults: 100000000 cases in total, 7 shards, seed 0x5eed
#   tools/fuzz_deep.sh --replay <seed hex>:<case>                           one case, the replay line's seed=... and case=...
#   tools/fuzz_deep.sh --selftest                                             proves argument checks, shard seeds, a small run and replay
# Deep runs are env-driven: CARGO_TARGET_DIR (default $(zb_scratch)/target-fuzz-deep), CARGO_BUILD_JOBS (default 7). No repo file is written.
# shellcheck source=tools/runlog.sh
. "${BASH_SOURCE[0]%/*}/runlog.sh"; runlog_wrap fuzz_deep "$@"   # run log (tools/runlog.sh)
set -u
export LC_ALL=C
root=$(cd "$(dirname "$0")/.." && pwd)
export CARGO_TARGET_DIR=${CARGO_TARGET_DIR:-$(zb_scratch)/target-fuzz-deep}
export CARGO_BUILD_JOBS=${CARGO_BUILD_JOBS:-7}
fail() { echo "fuzz_deep FAIL: $*" >&2; exit 2; }
# count_of <text> <label>: a case count, plain decimal digits or <digits>e<digits> (1e8, 25e5); at most 10^18; anything else fails naming <label>
count_of() {
  local t=$1 m e n
  case "$t" in
    ''|*[!0-9e]*|e*|*e|*e*e*) fail "$2 must be a plain integer or <digits>e<digits> such as 1e8, got '$t'" ;;
    *e*) m=${t%%e*}; e=${t#*e} ;;
    *) m=$t; e=0 ;;
  esac
  m=$(( 10#$m )); e=$(( 10#$e ))
  if ! { [ "${#m}" -le 18 ] && [ $(( ${#m} + e )) -le 19 ]; }; then fail "$2 is too large (at most 10^18), got '$t'"; fi
  n=$m; while [ "$e" -gt 0 ]; do n=$(( n * 10 )); e=$(( e - 1 )); done
  [ "$n" -le 1000000000000000000 ] || fail "$2 is too large (at most 10^18), got '$t'"
  echo "$n"
}
lock=${FUZZ_DEEP_LOCK:-/tmp/pop-deep.lock}; lock_wait=${FUZZ_DEEP_LOCK_WAIT:-3600}
build_to=${FUZZ_DEEP_BUILD_TIMEOUT:-30m}; shard_to=${FUZZ_DEEP_SHARD_TIMEOUT:-6h}
tfile=${FUZZ_DEEP_TARGETS:-$root/tools/fuzz_deep_targets.txt}
[ -f "$tfile" ] || fail "no target list at $tfile"
mapfile -t targets < <(grep -v -E '^[[:space:]]*(#|$)' "$tfile")
[ "${#targets[@]}" -ge 1 ] || fail "target list $tfile names no target"
cases=100000000; cases_text=; shards=7; seed=5eed; out=; replay=; selftest=0; rtarget=
if [ -n "${ZB_FUZZ_CASES+x}" ]; then cases=$(count_of "$ZB_FUZZ_CASES" ZB_FUZZ_CASES) || exit 2; fi
while [ $# -gt 0 ]; do
  case "$1" in
    --cases) cases_text=${2:-}; [ -n "$cases_text" ] || fail "--cases needs a value"; shift 2 ;;
    --shards) shards=${2:-}; shift 2 ;;
    --seed) seed=${2:-}; shift 2 ;;
    --out) out=${2:-}; shift 2 ;;
    --replay) replay=${2:-}; shift 2 ;;
    --target) rtarget=${2:-}; shift 2 ;;
    --selftest) selftest=1; shift ;;
    *) fail "unknown argument $1" ;;
  esac
done
seed=${seed#0x}
if [ -n "$cases_text" ]; then cases=$(count_of "$cases_text" --cases) || exit 2; fi
case "$shards" in ''|*[!0-9]*) fail "--shards must be a decimal number" ;; esac
case "$seed" in ''|*[!0-9a-fA-F]*) fail "--seed must be hex" ;; esac
[ "${#seed}" -le 16 ] || fail "--seed has more than 16 hex digits"
[ "$cases" -ge 1 ] || fail "--cases must be at least 1"
if ! { [ "$shards" -ge 1 ] && [ "$shards" -le 64 ]; }; then fail "--shards must be 1..64"; fi
[ "$cases" -ge "$shards" ] || fail "--cases must be at least --shards"

shard_seed() { # shard_seed <base hex> <shard index>: base + index * 0x9E3779B97F4A7C15 mod 2^64 (bash arithmetic is 64-bit modular; the Rust side only sees the hex)
  printf '%016x' $(( 0x$1 + $2 * 0x9E3779B97F4A7C15 ))
}
run_test() { (cd "$root" && timeout "$shard_to" cargo test --release --test harness "$1" -- --exact --test-threads=1 --nocapture 2>&1); }   # run_test MODULE::TEST
one_passed() { grep -q 'test result: ok\. 1 passed;' "$1"; }

if [ "$selftest" = 1 ]; then
  export RUNLOG_PLAIN=1 RUNLOG_SELFTEST=1   # selftest-only marker (runlog refuses RUNLOG_PLAIN alone); fixture runs below (many expected to fail) stay out of the run log; this selftest's own log captures their output
  me=$(cd "$(dirname "$0")" && pwd)/$(basename "$0")
  bad() { "$me" "$@" >/dev/null 2>&1 && fail "selftest: expected failure for: $*"; return 0; }
  bad --cases 0; bad --shards 0; bad --shards 65; bad --cases 3 --shards 7; bad --seed zz; bad --cases x; bad --bogus; bad --replay nocolon
  a=$(shard_seed 5eed 0); b=$(shard_seed 5eed 1); c=$(shard_seed 5eed 2)
  [ "$a" = 0000000000005eed ] || fail "selftest: shard 0 seed is $a"
  [ "$b" = 9e3779b97f4adb02 ] || fail "selftest: shard 1 seed is $b"
  if ! { [ "$a" != "$b" ] && [ "$b" != "$c" ]; }; then fail "selftest: shard seeds collide"; fi
  t=$(mktemp -d "$(zb_scratch)"/fuzz-deep-selftest.XXXXXX); runlog_on_exit "rm -rf \"$t\""
  "$me" --cases 4000 --shards 2 --out "$t" >"$t/run.txt" 2>&1 || { cat "$t/run.txt"; fail "selftest: small run failed"; }
  grep -q 'fuzz_deep: 2 shards ok, 4000 cases' "$t/run.txt" || { cat "$t/run.txt"; fail "selftest: no summary line"; }
  if ! { [ -s "$t/fuzz_mag.fuzz_mag_gate/shard-0.log" ] && [ -s "$t/fuzz_mag.fuzz_mag_gate/shard-1.log" ]; }; then fail "selftest: shard logs missing"; fi
  "$me" --replay 5eed:17 >"$t/replay.txt" 2>&1 || { cat "$t/replay.txt"; fail "selftest: replay failed"; }
  grep -q 'fuzz_deep: replay 5eed:17 ok' "$t/replay.txt" || fail "selftest: no replay summary"
  # case counts: ZB_FUZZ_CASES and --cases take plain integers or <digits>e<digits>; garbage is refused with a clear error
  for g in garbage 1e 1e99 e8 -5 1.5e3 0x10 ''; do
    msg=$(ZB_FUZZ_CASES="$g" "$me" --shards 2 2>&1) && fail "selftest: ZB_FUZZ_CASES=$g was accepted"
    case "$msg" in *"ZB_FUZZ_CASES"*) ;; *) fail "selftest: ZB_FUZZ_CASES=$g error does not name the variable: $msg" ;; esac
  done
  msg=$("$me" --cases 1e 2>&1) && fail "selftest: --cases 1e was accepted"
  case "$msg" in *"--cases"*) ;; *) fail "selftest: --cases 1e error does not name the flag: $msg" ;; esac
  [ "$(count_of 1e8 x)" = 100000000 ] || fail "selftest: 1e8 is not 100000000"
  [ "$(count_of 25e2 x)" = 2500 ] || fail "selftest: 25e2 is not 2500"
  [ "$(count_of 4000 x)" = 4000 ] || fail "selftest: 4000 is not 4000"
  ZB_FUZZ_CASES=2e3 "$me" --shards 2 --out "$t/env" >"$t/env.txt" 2>&1 || { cat "$t/env.txt"; fail "selftest: ZB_FUZZ_CASES=2e3 run failed"; }
  grep -q 'fuzz_deep: 2 shards ok, 2000 cases' "$t/env.txt" || { cat "$t/env.txt"; fail "selftest: ZB_FUZZ_CASES=2e3 was not honoured"; }
  ZB_FUZZ_CASES=2e3 "$me" --cases 3000 --shards 2 --out "$t/flag" >"$t/flag.txt" 2>&1 || fail "selftest: --cases over env run failed"
  grep -q '3000 cases' "$t/flag.txt" || fail "selftest: --cases does not override ZB_FUZZ_CASES"
  # deep runs take the deep lock (default /tmp/pop-deep.lock, never the gate lock); a busy lock blocks the run
  grep -q 'lock /tmp/pop-deep.lock' "$t/run.txt" || { cat "$t/run.txt"; fail "selftest: the small run did not report the deep lock"; }
  grep -q 'pop-gate.lock' "$t/run.txt" && fail "selftest: the run mentions the gate lock"
  lk="$t/deep.lock"
  ( exec 8>"$lk"; flock 8; sleep 30 8>&- ) & holder=$!   # only the subshell holds the lock fd, so killing it frees the lock
  sleep 0.3
  FUZZ_DEEP_LOCK="$lk" FUZZ_DEEP_LOCK_WAIT=1 "$me" --cases 1000 --shards 1 --out "$t/busy" >"$t/busy.txt" 2>&1 && { kill "$holder"; fail "selftest: the run did not wait for a held deep lock"; }
  grep -q 'deep lock .* busy' "$t/busy.txt" || { kill "$holder"; cat "$t/busy.txt"; fail "selftest: no busy-lock message"; }
  pkill -P "$holder" sleep 2>/dev/null; kill "$holder" 2>/dev/null; wait "$holder" 2>/dev/null
  FUZZ_DEEP_LOCK="$lk" FUZZ_DEEP_LOCK_WAIT=5 "$me" --cases 1000 --shards 1 --out "$t/free" >"$t/free.txt" 2>&1 || { cat "$t/free.txt"; fail "selftest: the run failed with a free lock"; }
  # Spence review: a renamed or missing fuzz test runs 0 tests under --exact and exits 0; every shard log must say '1 passed'
  printf '# fixture\nfuzz_mag::fuzz_mag_gate_renamed\n' > "$t/tg.txt"
  FUZZ_DEEP_TARGETS="$t/tg.txt" FUZZ_DEEP_LOCK="$t/l1" "$me" --cases 100 --shards 1 --out "$t/ren" >"$t/ren.txt" 2>&1 && { cat "$t/ren.txt"; fail "selftest: a renamed fuzz test passed"; }
  grep -q "1 passed" "$t/ren.txt" || { cat "$t/ren.txt"; fail "selftest: renamed-test failure does not say '1 passed' is missing"; }
  # timeouts on the build and on each shard
  FUZZ_DEEP_BUILD_TIMEOUT=0.001 FUZZ_DEEP_LOCK="$t/l2" "$me" --cases 100 --shards 1 --out "$t/bto" >"$t/bto.txt" 2>&1 && fail "selftest: build timeout not enforced"
  grep -q 'build timed out' "$t/bto.txt" || { cat "$t/bto.txt"; fail "selftest: no build-timeout message"; }
  FUZZ_DEEP_SHARD_TIMEOUT=1 FUZZ_DEEP_LOCK="$t/l3" "$me" --cases 1000000000 --shards 1 --out "$t/sto" >"$t/sto.txt" 2>&1 && fail "selftest: shard timeout not enforced"
  grep -q 'timed out' "$t/sto.txt" || { cat "$t/sto.txt"; fail "selftest: no shard-timeout message"; }
  # the deep lock has a default wait (no unbounded flock)
  grep -q 'waiting up to 3600s for deep lock' "$t/run.txt" || { cat "$t/run.txt"; fail "selftest: no default deep-lock wait"; }
  echo "fuzz_deep selftest ok"; exit 0
fi

if [ -n "$replay" ]; then
  case "$replay" in *:*) ;; *) fail "--replay must be <seed hex>:<case>" ;; esac
  rt=${rtarget:-${targets[0]}}; rl=$(mktemp "$(zb_scratch)"/fuzz-deep-replay.XXXXXX); runlog_on_exit "rm -f \"$rl\""
  ZB_REPLAY="$replay" run_test "$rt" > "$rl"; rc=$?
  grep -E 'replay|why:|test result' "$rl" | head -20
  [ "$rc" = 0 ] || { echo "fuzz_deep: replay $replay FAILED (rc $rc)"; exit 1; }
  one_passed "$rl" || { echo "fuzz_deep: replay $replay of $rt did not report '1 passed' (renamed or missing test?)"; exit 1; }
  echo "fuzz_deep: replay $replay ok"; exit 0
fi

# the deep lock (never the gate lock): one deep run at a time on this machine
exec 9>"$lock" || fail "cannot open the deep lock $lock"
echo "fuzz_deep: waiting up to ${lock_wait}s for deep lock $lock"
flock -w "$lock_wait" 9 || fail "deep lock $lock is busy (waited ${lock_wait}s)"
echo "fuzz_deep: holding deep lock $lock"
[ -n "$out" ] || out=$(mktemp -d "$(zb_scratch)"/fuzz-deep.XXXXXX)
mkdir -p "$out" || fail "cannot create $out"
(cd "$root" && timeout "$build_to" cargo test --release --test harness --no-run >"$out/build.log" 2>&1); brc=$?
[ "$brc" != 124 ] || fail "release build timed out after $build_to"
[ "$brc" = 0 ] || { cat "$out/build.log"; fail "release build failed"; }
per=$(( (cases + shards - 1) / shards ))
start=$(date +%s)
pids=(); names=(); logs=(); seeds=()
for tg in "${targets[@]}"; do
  d="$out/${tg//::/.}"; mkdir -p "$d" || fail "cannot create $d"
  i=0
  while [ "$i" -lt "$shards" ]; do
    s=$(shard_seed "$seed" "$i")
    ( ZB_FUZZ_SEED="$s" ZB_FUZZ_CASES="$per" run_test "$tg" >"$d/shard-$i.log" ) &
    pids+=($!); names+=("$tg shard $i"); logs+=("$d/shard-$i.log"); seeds+=("$s")
    i=$(( i + 1 ))
  done
done
bad=0; k=0
for p in "${pids[@]}"; do
  wait "$p"; rc=$?
  if [ "$rc" = 124 ]; then bad=$(( bad + 1 )); echo "fuzz_deep: ${names[$k]} timed out after $shard_to (seed ${seeds[$k]}); log ${logs[$k]}"
  elif [ "$rc" != 0 ]; then bad=$(( bad + 1 )); echo "fuzz_deep: ${names[$k]} FAILED (seed ${seeds[$k]}); log ${logs[$k]}"; grep -E 'suite=|why:' "${logs[$k]}" | head -6
  elif ! one_passed "${logs[$k]}"; then bad=$(( bad + 1 )); echo "fuzz_deep: ${names[$k]} did not report '1 passed' (renamed or missing test?); log ${logs[$k]}"; fi
  k=$(( k + 1 ))
done
secs=$(( $(date +%s) - start ))
nsh=$(( shards * ${#targets[@]} ))
total=$(( per * nsh ))
[ "$bad" = 0 ] || { echo "fuzz_deep: $bad of $nsh shards failed after ${secs}s; replay with: tools/fuzz_deep.sh --replay <seed>:<case> --target <target> (the failing line's seed= and case=)"; exit 1; }
echo "fuzz_deep: $nsh shards ok, $total cases in ${secs}s, targets: ${targets[*]} (logs in $out)"
