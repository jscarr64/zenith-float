#!/usr/bin/env bash
# Regenerates zenith-float-num/tests/mpfr/golden (r6 item 6, Paul call 2026-11-24): every MPFR comparison test runs with
# `--features mpfr-golden-gen`, sends each value query to gen.py (python-flint, pinned team env oracle-2026-10-30) and
# writes <module>__<test>.gold (binary zfgold2, parts <name>.gold.<k> above 90 MB; the readable text form goes to
# golden-full/, not committed); then gen.py --finish writes SHA256SUMS and PROVENANCE. Run under popslot, with the pinned
# env's bin first on PATH (gen.py fails closed otherwise, through goldlib):
#   PATH="$HOME/team-env/oracle-2026-10-30/bin:$PATH" popslot --agent <you> -- tools/mpfr_golden/gen.sh [TEST...]
# With TEST names (exact libtest names) only those streams are regenerated; SHA256SUMS/PROVENANCE are rewritten for all.
set -u
repo="$(cd "$(dirname "$0")/../.." && pwd)"; cd "$repo" || exit 1
python3 tools/mpfr_golden/gen.py --check-env || exit 2
dir=zenith-float-num/tests/mpfr/golden
threads=${MPFR_GOLDEN_THREADS:-$(nproc)}
if [ $# -eq 0 ]; then
  # every stream from scratch
  rm -f "$dir"/*.gold "$dir"/*.gold.[0-9]* "$dir/SHA256SUMS" "$dir/PROVENANCE"; rm -rf "$dir-full"
  cargo test -p zenith-float-num --release --features mpfr-golden-gen --test integration_tests -- --test-threads="$threads" || exit 1
else
  # only the named tests (exact libtest names, e.g. mpfr::compare_ops_test::mpfr_compare_ops); their streams are rewritten
  rm -f "$dir/SHA256SUMS" "$dir/PROVENANCE"
  cargo test -p zenith-float-num --release --features mpfr-golden-gen --test integration_tests -- --test-threads="$threads" --exact "$@" || exit 1
fi
python3 tools/mpfr_golden/gen.py --finish "$dir"
