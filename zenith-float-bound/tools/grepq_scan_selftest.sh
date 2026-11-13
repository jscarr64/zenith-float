#!/usr/bin/env bash
# Self-test for tools/grepq_scan.sh: a fixture tree with planted violations (every pipe-into-grep -q form the scan bans) must be
# flagged at exactly their file:line, a clean fixture (here-string, file read, comments, grep -c / grep -v pipes) must pass, and a
# shebang shell file without .sh is scanned too. usage: tools/grepq_scan_selftest.sh
# shellcheck source=tools/runlog.sh
. "${BASH_SOURCE[0]%/*}/runlog.sh"; runlog_wrap grepq_scan_selftest "$@"   # run log (tools/runlog.sh)
set -u
export LC_ALL=C
here="$(cd "$(dirname "$0")" && pwd)"; . "$here/scratch.sh"
fail() { echo "grepq-scan selftest FAIL: $*"; exit 1; }
t=$(mktemp -d "$(zb_scratch)/grepq-selftest.XXXXXX") || fail "mktemp failed"; runlog_on_exit "rm -rf \"$t\""; . "$here/selftest_cases.sh"; st_init "$t"
p='|'   # the pipe is spliced in at run time, so this file never contains the banned form itself
mkdir -p "$t/ok/tools" "$t/ok/scripts"
{ printf '#!/usr/bin/env bash\n# a comment naming printf x %s grep -q y is fine\n' "$p"
  printf 'grep -q y <<< "$x" && echo hit\ngrep -qxF -- y f.txt\nfalse %s%s grep -q y f.txt\nn=$(printf x %s grep -c y)\nprintf x %s grep -v y > f\n' "$p" "$p" "$p" "$p"; } > "$t/ok/tools/good.sh"
ck; out=$("$here/grepq_scan.sh" "$t/ok" 2>&1) || fail "clean fixture flagged: $out"
mkdir -p "$t/bad/tools" "$t/bad/scripts"
{ printf '#!/usr/bin/env bash\nok=1\n'
  printf 'printf x %s grep -q y && echo a\n' "$p"
  printf 'echo "$x" %sgrep -qxF y || fail b\n' "$p"
  printf 'echo "$x" %s grep -A1 k %s grep -qF y\n' "$p" "$p"
  printf 'echo "$x" %s grep -E -q y\n' "$p"
  printf '  echo "$x" %s grep --quiet y\n' "$p"; } > "$t/bad/tools/bad.sh"
{ printf '#!/bin/sh\n'; printf 'cat f %s grep -sq y\n' "$p"; } > "$t/bad/scripts/noext"
printf 'just text %s grep -q y\n' "$p" > "$t/bad/scripts/notes.txt"   # not a shell file: not scanned
ck; out=$("$here/grepq_scan.sh" "$t/bad" 2>&1) && fail "planted violations passed: $out"
want="$t/bad/scripts/noext:2
$t/bad/tools/bad.sh:3
$t/bad/tools/bad.sh:4
$t/bad/tools/bad.sh:5
$t/bad/tools/bad.sh:6
$t/bad/tools/bad.sh:7"
ck; got=$(sed -n 's/^\([^ ]*:[0-9]*\): .*/\1/p' <<< "$out")
[ "$got" = "$want" ] || fail "wrong hits: want
$want
got
$got
(full output: $out)"
ck; grep -q '^grepq_scan: FAIL 6 hits' <<< "$out" || fail "no FAIL line with the hit count: $out"
ck; grep -q 'here-string' <<< "$out" || fail "the FAIL message does not point to the here-string form: $out"
echo "grepq-scan selftest: ok"
st_done grepq_scan 5   # WANT: cases this script runs (tools/selftest_cases.sh)
