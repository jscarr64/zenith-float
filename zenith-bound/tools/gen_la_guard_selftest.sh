#!/usr/bin/env bash
# Self-test of the popguard gate at the top of scripts/gen_la_goldens.py (Harper/Paul, after the 11:48 58 GB OOM): the generator runs
# ~/pop-logs/bin/popguard.sh gen_la_goldens before any heavy import and exits with its nonzero rc, or 3 if the guard is missing or
# cannot run. Outside a popslot scope (no POPSLOT_SCOPE, or a copied/forged one) it must exit 3 without doing any work.
# shellcheck source=tools/runlog.sh
. "${BASH_SOURCE[0]%/*}/runlog.sh"; runlog_wrap gen_la_guard_selftest "$@"   # run log (tools/runlog.sh)
runlog_detach_children
set -u
here=$(cd "$(dirname "$0")" && pwd); gen="$here/../scripts/gen_la_goldens.py"
fail() { echo "gen_la_guard selftest FAIL: $*"; exit 1; }
. "$here/scratch.sh"; t=$(mktemp -d "$(zb_scratch)/gen-la-guard-selftest.XXXXXX") || fail "mktemp"; runlog_on_exit "rm -rf \"$t\""; . "${BASH_SOURCE[0]%/*}/selftest_cases.sh"; st_init "$t"
py=$(command -v python3) || fail "python3 missing"
run() { "$@" > "$t/out" 2>&1; echo $? > "$t/rc"; }
rc() { cat "$t/rc"; }
# fake HOME whose ~/pop-logs/bin/popguard.sh is a stub (exit code $1), or absent ($1 = none)
fakehome() { local h="$t/home-$1"; mkdir -p "$h/pop-logs/bin"; [[ $1 == none ]] || { printf "#!/bin/sh\nexit %s\n" "$1" > "$h/pop-logs/bin/popguard.sh"; chmod +x "$h/pop-logs/bin/popguard.sh"; }; echo "$h"; }
[[ -x $HOME/pop-logs/bin/popguard.sh ]] || fail "real popguard.sh missing at $HOME/pop-logs/bin"
ck; run env -u POPSLOT_SCOPE "$py" "$gen" --check --skip-heavy; [[ $(rc) == 3 ]] || fail "bare run (no POPSLOT_SCOPE) exited $(rc), want 3: $(head -c 400 "$t/out")"
ck; grep -q "popguard" "$t/out" || fail "bare run did not name popguard: $(head -c 400 "$t/out")"
ck; grep -q "^check:" "$t/out" && fail "bare run did work before the guard"
ck; run env POPSLOT_SCOPE=popslot-forged-1.scope "$py" "$gen" --check --skip-heavy; [[ $(rc) == 3 ]] || fail "forged POPSLOT_SCOPE exited $(rc), want 3"
h=$(fakehome none); ck; run env HOME="$h" "$py" "$gen" --check --skip-heavy; [[ $(rc) == 3 ]] || fail "missing popguard.sh exited $(rc), want 3"
ck; grep -q "missing\|cannot run" "$t/out" || fail "missing guard not named: $(head -c 400 "$t/out")"
h=$(fakehome 7); ck; run env HOME="$h" "$py" "$gen" --check --skip-heavy; [[ $(rc) == 7 ]] || fail "guard rc 7 propagated as $(rc), want 7"
h=$(fakehome 0); ck; run env HOME="$h" "$py" "$gen" --no-such-flag; [[ $(rc) != 3 ]] || fail "guard rc 0 still refused (rc 3)"
ck; grep -q "gen_la_goldens: popguard" "$t/out" && fail "guard rc 0 still printed a refusal"
st_done gen_la_guard 9   # WANT: cases this script runs (tools/selftest_cases.sh)
echo "gen_la_guard selftest: ok"
