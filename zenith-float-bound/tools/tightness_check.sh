#!/usr/bin/env bash
# Gate step (S12): the tightness ratchet. Reads the results file the tightness sets write (default
# ${CARGO_TARGET_DIR:-target}/tmp/tightness_results.txt, i.e. CARGO_TARGET_TMPDIR of the harness target; override with TIGHTNESS_RESULTS),
# prints the report and checks it against tools/tightness_baseline.txt. Reads only, never regenerates. See tools/tightness_report.py.
# shellcheck source=tools/runlog.sh
. "${BASH_SOURCE[0]%/*}/runlog.sh"; runlog_wrap tightness_check "$@"   # run log (tools/runlog.sh)
set -u
export LC_ALL=C PYTHONNOUSERSITE=1
root=${TIGHTNESS_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}
results=${TIGHTNESS_RESULTS:-${CARGO_TARGET_DIR:-$root/target}/tmp/tightness_results.txt}
tool=/usr/bin/python3
"$tool" "$root/tools/tightness_report.py" report "$results" || exit $?
"$tool" "$root/tools/tightness_report.py" check "$results" "$root/tools/tightness_baseline.txt"; exit $?
