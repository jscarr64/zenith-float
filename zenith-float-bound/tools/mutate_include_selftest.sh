#!/usr/bin/env bash
# Self-test for the per-slice include files of tools/mutate.sh (S0, 2026-10-24): tools/mutants/manifest.d/<slice>.txt (rows),
# tools/mutants/equivalent.d/<slice>.txt (documented equivalents) and tools/mutants/floor.d/<slice>.txt (floor shares).
# Runs only `mutate.sh --list` and `mutate.sh --floor` (no cargo, no worktree) on a fixture under $(zb_scratch).
# Proves: include rows are listed after the central manifest; an id outside its file's `# range: LO-HI` FAILS; a file without a
# range line FAILS; two files with overlapping ranges FAIL; an id duplicated across files FAILS; floor shares add to floor.txt.
# Hardening (Paul, 2026-10-24): a central id inside a slice's range FAILS; an include row whose patch is not under <slice>/ FAILS;
# suffixed ids sharing a number (205-a, 205-b) FAIL the duplicate check; a row without 4 tab-separated columns FAILS; equivalent.d
# and floor.d entries are validated (known id in the slice's range with a reason; a matching manifest.d file; each key once).
# shellcheck source=tools/runlog.sh
. "${BASH_SOURCE[0]%/*}/runlog.sh"; runlog_wrap mutate_include_selftest "$@"   # run log (tools/runlog.sh)
set -u
here="$(cd "$(dirname "$0")" && pwd)"; t="$(mktemp -d "$(zb_scratch)/mutate-include-selftest.XXXXXX")"; runlog_on_exit "rm -rf \"$t\""
. "${BASH_SOURCE[0]%/*}/selftest_cases.sh"; st_init "$t"
fail() { echo "SELFTEST FAIL: $*"; exit 1; }
m="$t/mutants"; mkdir -p "$m/manifest.d" "$m/floor.d" "$m/equivalent.d"
tab=$'\t'
printf '# central\nsrc/a.rs:f%s001-A%s001-A.patch%ssuite\n' "$tab" "$tab" "$tab" > "$m/manifest.txt"
printf 'killed\t1\nkilled+equivalent\t1\n' > "$m/floor.txt"
printf '# range: 200-299\nsrc/b.rs:g%s201-B%ss1/201-B.patch%ssuite\n' "$tab" "$tab" "$tab" > "$m/manifest.d/s1.txt"
printf '# range: 300-399\n' > "$m/manifest.d/s2.txt"
printf 'killed\t2\nkilled+equivalent\t3\n' > "$m/floor.d/s1.txt"
printf 'killed\t0\nkilled+equivalent\t0\n' > "$m/floor.d/s2.txt"
mkdir -p "$m/s1"; touch "$m/001-A.patch" "$m/s1/201-B.patch"
mu() { MUTATE_MANIFEST="$m/manifest.txt" RUNLOG_PLAIN=1 RUNLOG_SELFTEST=1 "$here/mutate.sh" "$@" 2>&1; }
reset() { printf '# central\nsrc/a.rs:f%s001-A%s001-A.patch%ssuite\n' "$tab" "$tab" "$tab" > "$m/manifest.txt"
  printf '# range: 200-299\nsrc/b.rs:g%s201-B%ss1/201-B.patch%ssuite\n' "$tab" "$tab" "$tab" > "$m/manifest.d/s1.txt"
  rm -f "$m"/equivalent.d/*.txt; printf 'killed\t2\nkilled+equivalent\t3\n' > "$m/floor.d/s1.txt"; rm -f "$m/floor.d/s9.txt"; }
reset; ck; out=$(mu --list) || fail "reset fixture fails: $out"
printf 'src/a.rs:f%s250-Z%s001-A.patch%ssuite\n' "$tab" "$tab" "$tab" >> "$m/manifest.txt"
ck; out=$(mu --list) && fail "central id 250 inside s1's range passed: $out"
ck; grep -q 'owned by slice' <<< "$out" || fail "central-in-range message missing: $out"
reset; printf 'src/b.rs:g%s202-E%s001-A.patch%ssuite\n' "$tab" "$tab" "$tab" >> "$m/manifest.d/s1.txt"
ck; out=$(mu --list) && fail "s1 row whose patch is outside s1/ passed: $out"
ck; grep -q 'not under s1/' <<< "$out" || fail "patch-ownership message missing: $out"
reset; printf 'src/b.rs:g%s205-a%ss1/201-B.patch%ssuite\nsrc/b.rs:g%s205-b%ss1/201-B.patch%ssuite\n' "$tab" "$tab" "$tab" "$tab" "$tab" "$tab" >> "$m/manifest.d/s1.txt"
ck; out=$(mu --list) && fail "suffixed ids 205-a and 205-b passed the duplicate check: $out"
ck; grep -q 'duplicate manifest id number 205' <<< "$out" || fail "suffix-dup message missing: $out"
reset; printf 'src/a.rs:f 003-N 001-A.patch suite\n' >> "$m/manifest.txt"
ck; out=$(mu --list) && fail "central row without tabs passed: $out"
ck; grep -q '4 tab-separated columns' <<< "$out" || fail "no-tab message missing: $out"
reset; printf 'src/b.rs:g 203-N s1/201-B.patch suite\n' >> "$m/manifest.d/s1.txt"
ck; out=$(mu --list) && fail "include row without tabs passed: $out"
reset; printf '201-B\tunobservable by construction\n' > "$m/equivalent.d/s1.txt"
ck; out=$(mu --list) || fail "a valid equivalent.d entry failed: $out"
printf '299-X\tno such mutant\n' > "$m/equivalent.d/s1.txt"
ck; out=$(mu --list) && fail "equivalent.d id not in the manifest passed: $out"
printf '201-B\n' > "$m/equivalent.d/s1.txt"
ck; out=$(mu --list) && fail "equivalent.d entry without a reason passed: $out"
reset; printf '201-B\treason\n' > "$m/equivalent.d/s2.txt"
ck; out=$(mu --list) && fail "equivalent.d/s2 claiming s1's id passed: $out"
reset; printf '201-B\treason\n' > "$m/equivalent.d/s9.txt"
ck; out=$(mu --list) && fail "equivalent.d file without a manifest.d slice passed: $out"
reset; printf 'killed\t0\nkilled+equivalent\t0\n' > "$m/floor.d/s9.txt"
ck; out=$(mu --floor) && fail "floor.d file without a manifest.d slice passed: $out"
reset; printf 'killed\t1\nkilled\t2\nkilled+equivalent\t3\n' > "$m/floor.d/s1.txt"
ck; out=$(mu --floor) && fail "floor.d with a duplicate key passed: $out"
reset; printf 'killed\t1\nkilled+equivalent\t3\nkiled\t5\n' > "$m/floor.d/s1.txt"
ck; out=$(mu --floor) && fail "floor.d with an unknown key passed: $out"
reset
ck; out=$(mu --list) || fail "--list with includes failed: $out"
ck; [ "$(printf '%s\n' "$out" | tr '\n' ' ')" = "001-A 201-B " ] || fail "--list gave '$out', want 001-A then 201-B"
ck; out=$(mu --floor) || fail "--floor failed: $out"
ck; [ "$(printf '%s\n' "$out" | tr '\n' ' ')" = "killed 3 killed+equivalent 4 " ] || fail "--floor gave '$out', want killed 3, killed+equivalent 4"
printf 'src/c.rs:h%s305-C%s201-B.patch%ssuite\n' "$tab" "$tab" "$tab" >> "$m/manifest.d/s1.txt"
ck; out=$(mu --list) && fail "id 305 outside s1's range 200-299 passed: $out"
ck; grep -q 'outside' <<< "$out" || fail "out-of-range message missing: $out"
printf '# range: 200-299\nsrc/b.rs:g%s201-B%ss1/201-B.patch%ssuite\n' "$tab" "$tab" "$tab" > "$m/manifest.d/s1.txt"
printf 'src/b.rs:g%s301-D%s201-B.patch%ssuite\n' "$tab" "$tab" "$tab" > "$m/manifest.d/s3.txt"
ck; out=$(mu --list) && fail "include file without a range line passed: $out"
printf '# range: 250-320\n' > "$m/manifest.d/s3.txt"
ck; out=$(mu --list) && fail "overlapping ranges passed: $out"
printf '# range: 400-499\nsrc/b.rs:g%s201-B%ss1/201-B.patch%ssuite\n' "$tab" "$tab" "$tab" > "$m/manifest.d/s3.txt"
ck; out=$(mu --list) && fail "id duplicated across include files passed: $out"
rm "$m/manifest.d/s3.txt"
# Spence (2026-10-24, 1422-zb-s2-const-steps: 368 `ERR rc=1 line=92 owner=$(inrange ...)`): an id outside every slice range is the
# normal case and must not trip the poplog ERR trap; a real inrange failure (malformed range row) must still be loud and stop the run.
pl() { local n=$1; shift   # mutate.sh with its own real run log (ERR trap on); prints output, rc, then the log path
  ( unset POPLOG_FILE POPLOG_PID POPLOG_STATE_DIR POPLOG_PROGRESS_FILE POPLOG_RUNNING RUNLOG_INNER RUNLOG_ACTIVE RUNLOG_PLAIN
    MUTATE_MANIFEST="$m/manifest.txt" POPLOG_ROOT="$t/pl-$n" RUNLOG_AGENT=selftest "$here/mutate.sh" "$@" 2>&1; echo "rc=$?" )
  find "$t/pl-$n" -name '*-mutate.log' | head -1; }
reset; printf 'src/a.rs:f%s002-A%s001-A.patch%ssuite\nsrc/a.rs:f%s950-A%s001-A.patch%ssuite\n' "$tab" "$tab" "$tab" "$tab" "$tab" "$tab" >> "$m/manifest.txt"
printf '201-B\tunobservable by construction\n' > "$m/equivalent.d/s1.txt"
out=$(pl quiet --list); lg=$(tail -1 <<< "$out")
ck; grep -q '^rc=0$' <<< "$out" && [ -f "$lg" ] || fail "err-quiet: --list with out-of-range central ids failed or wrote no run log: $out"
ck; grep -q ' ERR ' "$lg" && fail "err-quiet: out-of-range ids (no owning slice) logged ERR trap lines: $(grep -c ' ERR ' "$lg") e.g. $(grep -m1 ' ERR ' "$lg")"
# inrange contract, on the function itself (ranges are validated before it runs, so a malformed row is an internal error)
ir=$(awk '/^inrange\(\) \{/{f=1;s=NR} f{print} f&&(/^\}$/||(NR==s&&/; \}$/)){exit}' "$here/mutate.sh")
ck; [ -n "$ir" ] || fail "inrange() not found in mutate.sh"
r=$(bash -c 'eval "$1"; ranges="200 299 s1.txt"$'"'"'\n'"'"'; o=$(inrange 5); echo "rc=$? o=[$o]"; o=$(inrange 250); echo "rc=$? o=[$o]"' _ "$ir" 2>&1)
ck; [ "$r" = "rc=0 o=[]"$'\n'"rc=0 o=[s1.txt]" ] || fail "inrange: a miss must be rc 0 with empty output, a hit rc 0 + file; got: $r"
# inrange forces base 10 itself (Spence nit on b5ac93f): a leading-zero id is decimal, never octal (0213 is 213, not 139; 0219 is no error)
r=$(bash -c 'eval "$1"; ranges="200 299 s1.txt"$'"'"'\n'"'"'; o=$(inrange 0213); echo "rc=$? o=[$o]"; o=$(inrange 0219); echo "rc=$? o=[$o]"; o=$(inrange 0139); echo "rc=$? o=[$o]"' _ "$ir" 2>&1)
ck; [ "$r" = "rc=0 o=[s1.txt]"$'\n'"rc=0 o=[s1.txt]"$'\n'"rc=0 o=[]" ] || fail "inrange: leading-zero ids must be base 10 (0213, 0219 in 200-299; 0139 not); got: $r"
r=$(bash -c 'eval "$1"; ranges="x 299 s1.txt"$'"'"'\n'"'"'; o=$(inrange 5); echo "rc=$?"' _ "$ir" 2>&1)
ck; grep -q '^rc=0$' <<< "$r" && fail "inrange: a malformed range row returned rc 0 (a real failure went quiet): $r"
ck; grep -q 'inrange' <<< "$r" || fail "inrange: a malformed range row gave no message: $r"
ck; grep -nE 'inrange .*\|\| *true' "$here/mutate.sh" && fail "a call site masks inrange's rc with || true (a real failure would be silent)"
reset
printf 'killed\tx\n' > "$m/floor.d/s2.txt"
ck; out=$(mu --floor) && fail "non-integer floor share passed: $out"
echo "mutate_include selftest: ok"
st_done mutate_include 36   # WANT: cases this script runs (tools/selftest_cases.sh)
