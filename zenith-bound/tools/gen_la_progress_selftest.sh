#!/usr/bin/env bash
# Self-test of the heavy-run progress contract of scripts/gen_la_goldens.py (Jeff, 2026-10-24, after la4-heavy2 sat 90 min at
# 100% CPU with zero output): whenever heavy rows are computed, every row prints one stderr line
#   progress: op=<op> case=<case_id> tier=<gate|heavy> row_s=<s> elapsed_s=<s>
# after the row is on disk (heavy rows: --heavy DIR/<slice>/la_<op>.rows<k>.partial, flushed + fsynced per row), and a per-case
# wall budget (--case-budget SECONDS, default 1800) fails loudly with rc 5 naming the op and the case_id being computed.
# Needs the oracle venv (sympy/mpmath/gmpy2 pins) and a popslot scope (popguard); run via poprun/popslot.
# shellcheck source=tools/runlog.sh
. "${BASH_SOURCE[0]%/*}/runlog.sh"; runlog_wrap gen_la_progress_selftest "$@"   # run log (tools/runlog.sh)
runlog_detach_children
set -u
here=$(cd "$(dirname "$0")" && pwd); gen="$here/../scripts/gen_la_goldens.py"
fail() { echo "gen_la_progress selftest FAIL: $*"; exit 1; }
. "$here/scratch.sh"; t=$(mktemp -d "$(zb_scratch)/gen-la-progress-selftest.XXXXXX") || fail "mktemp"; runlog_on_exit "rm -rf \"$t\""; . "${BASH_SOURCE[0]%/*}/selftest_cases.sh"; st_init "$t"
py=${ZB_ORACLE_PY:-$HOME/.local/oracle-venv/bin/python}; [[ -x $py ]] || fail "oracle python missing: $py"
run() { "$@" > "$t/out" 2> "$t/err"; echo $? > "$t/rc"; }
rc() { cat "$t/rc"; }
# 1. budget: a 1 ms per-case budget must trip on some case, rc 5, naming op and case_id
ck; run "$py" "$gen" --heavy "$t/h1" --case-budget 0.001; [[ $(rc) == 5 ]] || fail "1 ms budget exited $(rc), want 5: $(tail -c 400 "$t/err")"
ck; grep -Eq '^gen_la_goldens: CASE BUDGET exceeded \(0\.001 s\): op=[a-z_]+ case=[a-z_]+\.[^ ]+' "$t/err" || fail "budget failure does not name op and case: $(tail -c 400 "$t/err")"
# 2. a non-positive budget is a usage error (argparse rc 2), never "no budget"
ck; run "$py" "$gen" --heavy "$t/h0" --case-budget 0; [[ $(rc) == 2 ]] || fail "--case-budget 0 exited $(rc), want 2"
# 3. progress + fsync: run --heavy until the first heavy progress line, SIGKILL it, then that row must already be in a .partial file
"$py" "$gen" --heavy "$t/h2" > "$t/out2" 2> "$t/err2" & pid=$!
for _ in $(seq 1 6000); do grep -q 'tier=heavy' "$t/err2" && break; kill -0 "$pid" 2>/dev/null || break; sleep 0.05; done
kill -9 "$pid" 2>/dev/null; wait "$pid" 2>/dev/null
ck; line=$(grep -m1 'tier=heavy' "$t/err2") || fail "no heavy progress line within 300 s: $(tail -c 400 "$t/err2")"
ck; [[ $line =~ ^progress:\ op=[a-z_]+\ case=([a-z_]+\.[^ ]+)\ tier=heavy\ row_s=[0-9]+\.[0-9]+\ elapsed_s=[0-9]+\.[0-9]+$ ]] || fail "bad progress line: $line"
cid=${BASH_REMATCH[1]}
ck; grep -Eq '^progress: op=[a-z_]+ case=[^ ]+ tier=gate row_s=' "$t/err2" || fail "no gate progress lines"
ck; grep -qs "^$cid " "$t/h2"/*/*.partial "$t/h2"/*/*_heavy.txt || fail "row $cid was reported but is on disk in neither a .partial nor a finished heavy file after SIGKILL"
st_done gen_la_progress 7   # WANT: cases this script runs (tools/selftest_cases.sh)
echo "gen_la_progress selftest: ok"
