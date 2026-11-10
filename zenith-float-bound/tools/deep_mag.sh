#!/usr/bin/env bash
# Deep run of the Mag property fuzz (tests/harness/fuzz_mag.rs): ZB_FUZZ_CASES cases per property (default 100000000) on the DEEP seed stream.
# Never part of tools/full_gate.sh. It takes its OWN lock, /tmp/pop-deep.lock, never the gate lock /tmp/pop-gate.lock, and refuses to start
# while another deep run holds it (exit 75), so a deep run and the gate can overlap; use a separate CARGO_TARGET_DIR (default below).
#
# usage: tools/deep_mag.sh [TAG]            env: ZB_FUZZ_CASES=N  ZB_FUZZ_SEED=S  DEEP_TARGET=DIR  T=LOGDIR  REPO=DIR  DEEP_JOBS=N
# output: $T/deep-mag-TAG.log (per-property case counts and outcome counts, test result, wall time), $T/deep-mag-TAG.rc (cargo exit code)
# A failing property prints a replay line (`ZB_REPLAY=<seed>:<case>`) with the failing input: report it, do not patch src/.
# Replay a deep failure with ZB_FUZZ_DEEP=1 ZB_REPLAY=<seed>:<case> cargo test --test harness fuzz_mag:: (the seed already encodes the stream).
set -u
. "${BASH_SOURCE[0]%/*}/scratch.sh"
export PYTHONNOUSERSITE=1
tag=${1:-$(date +%Y%m%d-%H%M%S)}
T=${T:-$(zb_scratch)/zb-gate-logs}; mkdir -p "$T"
cd "${REPO:-$(cd "$(dirname "$0")/.." && pwd)}" || exit 2
exec 9> /tmp/pop-deep.lock
flock -n 9 || { echo "deep_mag: /tmp/pop-deep.lock is held by another deep run" >&2; exit 75; }
export CARGO_TARGET_DIR=${DEEP_TARGET:-$(zb_scratch)/target-deep-mag}
export CARGO_BUILD_JOBS=${DEEP_JOBS:-7}
export ZB_FUZZ_DEEP=1
export ZB_FUZZ_CASES=${ZB_FUZZ_CASES:-100000000}
log="$T/deep-mag-$tag.log"; rcf="$T/deep-mag-$tag.rc"
rm -f "$rcf"
start=$(date +%s)
{
  echo "deep_mag start $(date '+%F %T %Z') cases=$ZB_FUZZ_CASES head=$(git rev-parse --short HEAD) target=$CARGO_TARGET_DIR"
  cargo test --release --test harness "fuzz_mag::" -- --nocapture --test-threads="${DEEP_JOBS:-7}"
  echo "cargo_rc=$?" 
  echo "deep_mag wall=$(( $(date +%s) - start ))s"
} > "$log" 2>&1
rc=$(grep -o 'cargo_rc=[0-9]*' "$log" | tail -1 | cut -d= -f2)
echo "${rc:-1}" > "$rcf"
exit "${rc:-1}"
