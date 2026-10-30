#!/usr/bin/env bash
# Self-test for tools/mutate.sh: a tiny fixture crate (a toy function and a toy suite) is generated under $(zb_scratch) at run time as its own git repo.
# Proves the runner KILLS a real mutant, reports a surviving mutant as SURVIVED (exit nonzero), accepts it only when documented as
# equivalent, flags a non-building mutant and a non-applying patch as INVALID, counts a timeout as a kill, FAILS with 0 mutants,
# and ALWAYS removes its worktree. The repo under test is never touched.
# shellcheck source=tools/runlog.sh
. "${BASH_SOURCE[0]%/*}/runlog.sh"; runlog_wrap mutate_selftest "$@"   # run log (tools/runlog.sh)
runlog_detach_children
set -u
here="$(cd "$(dirname "$0")" && pwd)"; t="$(mktemp -d "$(zb_scratch)"/mutate-selftest.XXXXXX)"; runlog_on_exit "rm -rf \"$t\""
fail() { echo "SELFTEST FAIL: $*"; exit 1; }
r="$t/fix"; mkdir -p "$r/src" "$r/tests" "$r/tools/mutants"; cd "$r" || exit 1
git init -q . && git config user.email a@b && git config user.name t
cat > Cargo.toml <<'TOML'
[package]
name = "toyfix"
version = "0.0.0"
edition = "2021"
[dependencies]
TOML
cat > src/lib.rs <<'RS'
pub fn add(a: u32, b: u32) -> u32 {
    a + b
}
pub fn is_small(x: u32) -> bool {
    x < 10
}
RS
cat > tests/toy.rs <<'RS'
use toyfix::{add, is_small};
#[test]
fn add_works() {
    assert_eq!(add(2, 3), 5);
}
#[test]
fn small_works() {
    assert!(is_small(9));
    assert!(!is_small(10));
}
RS
printf '/Cargo.lock\n/target\n' > .gitignore   # mutate.sh runs `cargo fetch` in the repo, which writes Cargo.lock
git add -A; git commit -qm fixture
mk() { # mk <name> <sed-expr>: a patch made by editing src/lib.rs, then reverting
  sed -i "$2" src/lib.rs; git diff > "tools/mutants/$1.patch"; git checkout -q src/lib.rs
}
mk 001-add-mul 's/a + b/a * b/'                 # killed: add(2,3) becomes 6
mk 002-add-commute 's/a + b/b + a/'             # equivalent: addition commutes, the suite cannot tell
mk 003-nobuild 's/a + b/a + undefined_name/'    # does not build
mk 004-hang 's/a + b/loop {}/'                  # never returns: timeout counts as a kill
printf 'garbage that is not a patch\n' > tools/mutants/005-badpatch.patch
printf 'src/lib.rs:add\t001-add-mul\t001-add-mul.patch\ttoy\n' > m.killed
printf 'src/lib.rs:add\t002-add-commute\t002-add-commute.patch\ttoy\n' > m.surv
printf 'src/lib.rs:add\t003-nobuild\t003-nobuild.patch\ttoy\n' > m.nobuild
printf 'src/lib.rs:add\t004-hang\t004-hang.patch\ttoy\n' > m.hang
printf 'src/lib.rs:add\t005-badpatch\t005-badpatch.patch\ttoy\n' > m.bad
printf '# header only\n' > tools/mutants/EQUIVALENT.txt
git add -A; git commit -qm mutants
export MUTATE_REPO="$r" MUTATE_TMP="$t"
# Every run() is one case; the count is asserted at the end, so a selftest that silently skips its cases cannot read ok (Paul, 2026-10-23:
# the zb-green-mutate-dirty2 step ran 0 tests and exited 0).
printf '# fixture floor\nkilled\t0\nkilled+equivalent\t1\n' > "$t/mfloor.txt"
run() { printf 'x\n' >> "$t/cases"; MUTATE_FLOOR_ZERO_OK="${MZERO-1}" MUTATE_FLOOR="${MFLOOR:-$t/mfloor.txt}" MUTATE_MANIFEST="$r/tools/mutants/manifest.txt" "$here/mutate.sh" "$@" 2>&1; }   # mutate.sh FAIL lines go to stderr (9e8612d); the cases grep both streams
leftovers() { git -C "$r" worktree list | grep -c 'zb-mutate\.' ; find "$t" -maxdepth 1 \( -name "zb-mutate.*" -o -name "zb-mutate-target.*" \) | wc -l; }
cm() { git -C "$r" add -A; git -C "$r" commit -qm fixture-edit >/dev/null 2>&1 || :; }   # commit a fixture edit: mutate.sh refuses a dirty tree
setm() { cp "$r/$1" "$r/tools/mutants/manifest.txt"; cm; }   # committed: mutate.sh refuses a dirty tree

# manifest integrity (checked before any mutant runs): an id listed twice, or a patch file that does not exist, FAILS the whole run
printf 'src/lib.rs:add\t001-add-mul\t001-add-mul.patch\ttoy\nsrc/lib.rs:add\t001-add-mul\t002-add-commute.patch\ttoy\n' > "$r/tools/mutants/manifest.txt"; cm
out="$(run)" && fail "duplicate manifest id passed: $out"; echo "$out" | grep -c "^mutate: FAIL duplicate manifest id 001-add-mul" > /dev/null || fail "duplicate id: wrong message: $out"
printf 'src/lib.rs:add\t001-add-mul\t001-add-mul.patch\ttoy\nsrc/lib.rs:add\t006-gone\t006-gone.patch\ttoy\n' > "$r/tools/mutants/manifest.txt"; cm
out="$(run)" && fail "missing patch file passed: $out"; echo "$out" | grep -c "^mutate: FAIL manifest id 006-gone names missing patch 006-gone.patch" > /dev/null || fail "missing patch: wrong message: $out"

# 0 mutants: a header-only manifest FAILS (Spence review: the gate's mutants step must run at least one mutant)
printf '# empty\n' > "$r/tools/mutants/manifest.txt"; cm
out="$(run)" && fail "empty manifest passed (0 mutants proves nothing): $out"; echo "$out" | grep -c "0 mutants ran" > /dev/null || fail "empty manifest: wrong message: $out"

# a real mutant must be KILLED
setm m.killed; out="$(run)" || fail "killable mutant not accepted: $out"
echo "$out" | grep -c "^KILLED 001-add-mul" > /dev/null || fail "killed: wrong report: $out"

# a surviving mutant makes the runner exit nonzero and is reported as SURVIVED
setm m.surv; out="$(run)" && fail "surviving mutant passed (runner exited 0)"
echo "$out" | grep -c "^SURVIVED 002-add-commute" > /dev/null || fail "survivor: wrong report: $out"
# ... unless it is documented as equivalent (reported differently, exit 0)
printf '002-add-commute\taddition commutes\n' > "$r/tools/mutants/EQUIVALENT.txt"; cm
out="$(run)" || fail "documented-equivalent survivor failed the run: $out"
echo "$out" | grep -c "^SURVIVED-EQUIVALENT 002-add-commute" > /dev/null || fail "equivalent: wrong report: $out"
printf '# header only\n' > "$r/tools/mutants/EQUIVALENT.txt"; cm

# a mutant that does not build and a patch that does not apply are INVALID (nonzero)
setm m.nobuild; out="$(run)" && fail "non-building mutant passed"; echo "$out" | grep -c "^INVALID 003-nobuild" > /dev/null || fail "nobuild: wrong report: $out"
setm m.bad; out="$(run)" && fail "non-applying patch passed"; echo "$out" | grep -c "^INVALID 005-badpatch" > /dev/null || fail "badpatch: wrong report: $out"

# a hang is killed by the timeout
setm m.hang; out="$(MUTATE_TIMEOUT=5s run)" || fail "hanging mutant not accepted as killed: $out"
echo "$out" | grep -c "^KILLED 004-hang.*timed out" > /dev/null || fail "hang: wrong report: $out"

# --only selects one row; an unknown id is an error
cat m.killed m.surv > "$r/tools/mutants/manifest.txt"; cm
out="$(run --only 001-add-mul)" || fail "--only killed row failed: $out"; echo "$out" | grep -c "1 mutants" > /dev/null || fail "--only: wrong count: $out"
run --only nope >/dev/null 2>&1 && fail "unknown --only id passed"

# MUTATE_TARGET: a persistent target dir is reused and kept (same verdicts), the fixed worktree is still removed; bad values are refused
pt="$t/persist-target"; setm m.killed
out="$(MUTATE_TARGET="$pt" run)" || fail "persistent target: killable mutant not accepted: $out"
echo "$out" | grep -c "^KILLED 001-add-mul" > /dev/null || fail "persistent target: wrong report: $out"
[ -d "$pt/debug" ] || fail "persistent target: $pt/debug was not kept"
touch "$pt/marker"
setm m.surv; out="$(MUTATE_TARGET="$pt" run)" && fail "persistent target: surviving mutant passed"
echo "$out" | grep -c "^SURVIVED 002-add-commute" > /dev/null || fail "persistent target: survivor wrong report: $out"
[ -f "$pt/marker" ] || fail "persistent target: the dir was not reused (marker gone)"
[ -e "$t/zb-mutate.persist" ] && fail "persistent target: fixed worktree left behind"
MUTATE_TARGET="relative/dir" run >/dev/null 2>&1 && fail "relative MUTATE_TARGET accepted"
MUTATE_TARGET="$r/target" run >/dev/null 2>&1 && fail "MUTATE_TARGET inside the repo accepted"
# the fixed persistent worktree lives beside the persistent target (or under MUTATE_PERSIST), never under MUTATE_TMP, whose default was
# /tmp (wiped on reboot): an unwritable MUTATE_TMP must not matter in persistent mode
ro="$t/ro"; mkdir -p "$ro"; chmod 555 "$ro"; setm m.killed
out="$(MUTATE_TMP="$ro" MUTATE_PERSIST="$t/pp" MUTATE_TARGET="$pt" run 2>&1)" || fail "MUTATE_PERSIST: worktree not placed under it: $out"
[ -d "$t/pp" ] && [ ! -e "$t/pp/zb-mutate.persist" ] || fail "MUTATE_PERSIST: dir missing or worktree left behind"
out="$(MUTATE_TMP="$ro" MUTATE_TARGET="$pt" run 2>&1)" || fail "persistent worktree not placed beside MUTATE_TARGET: $out"
chmod 755 "$ro"
rm -rf "$pt"

# mutate.sh tests the COMMITTED HEAD, so a dirty tree would test code other than what is on disk (2102-zb-mut180: SURVIVED because the
# uncommitted mag_sort tests were not in HEAD 7191a1f). It must refuse a modified or untracked file, and print the SHA it tested.
setm m.killed; head_sha=$(git -C "$r" rev-parse HEAD)
out="$(run)" || fail "clean tree: killable mutant not accepted: $out"
echo "$out" | grep -qx "mutate: testing committed HEAD $head_sha" || fail "clean tree: tested SHA not printed: $out"
echo '// edit' >> "$r/tests/toy.rs"
out="$(run)" && fail "modified tracked file passed: $out"; echo "$out" | grep -q "^mutate: FAIL dirty tree" || fail "modified: wrong message: $out"
git -C "$r" checkout -q tests/toy.rs; printf '#[test]\nfn new_one() {}\n' > "$r/tests/new_untracked.rs"
out="$(run)" && fail "untracked file passed: $out"; echo "$out" | grep -q "^mutate: FAIL dirty tree" || fail "untracked: wrong message: $out"
rm -f "$r/tests/new_untracked.rs"
# Paul (2026-10-23): a full run FAILS when killed or killed+equivalent is below the checked-in floor (tools/mutants/floor.txt; MUTATE_FLOOR
# overrides it). --only runs skip the floor but still need one mutant run. A missing, garbage or CRLF floor FAILS (fail closed).
setm m.killed
printf 'killed\t2\nkilled+equivalent\t2\n' > "$t/f.high"; out="$(MFLOOR="$t/f.high" run)" && fail "below the mutant floor passed: $out"
echo "$out" | grep -q "^mutate: FAIL BELOW floor: 1 killed < 2" || fail "below floor: wrong message: $out"
printf 'killed\t1\nkilled+equivalent\t1\n' > "$t/f.at"; out="$(MFLOOR="$t/f.at" run)" || fail "at the mutant floor failed: $out"
echo "$out" | grep -q "^mutate: floor ok" || fail "at floor: no floor line: $out"
out="$(MFLOOR="$t/f.high" run --only 001-add-mul)" || fail "--only must skip the floor: $out"
out="$(MFLOOR="$t/f.none" run)" && fail "a missing mutant floor passed: $out"; echo "$out" | grep -q "^mutate: FAIL NOFLOOR" || fail "missing floor: wrong message: $out"
printf 'killed\tlots\nkilled+equivalent\t1\n' > "$t/f.bad"; out="$(MFLOOR="$t/f.bad" run)" && fail "a garbage mutant floor passed: $out"
echo "$out" | grep -q "^mutate: FAIL BADCOUNT" || fail "garbage floor: wrong message: $out"
printf 'killed\t1\r\nkilled+equivalent\t1\r\n' > "$t/f.crlf"; out="$(MFLOOR="$t/f.crlf" run)" && fail "a CRLF mutant floor passed: $out"
echo "$out" | grep -q "^mutate: FAIL BADCOUNT.*CRLF" || fail "CRLF floor: wrong message: $out"
# Spence (2026-10-23): a real floor of 0 FAILS; only this fixture (MUTATE_FLOOR_ZERO_OK=1, set by run()) may use 0.
out="$(MZERO= run)" && fail "a 0 mutant floor without MUTATE_FLOOR_ZERO_OK passed: $out"
echo "$out" | grep -q "^mutate: FAIL BADCOUNT.*must be 1 or more" || fail "0 floor: wrong message: $out"
printf '# empty\n' > "$r/tools/mutants/manifest.txt"; cm
out="$(MFLOOR="$t/f.at" run)" && fail "0 mutants with a floor passed: $out"; echo "$out" | grep -c "0 mutants ran" > /dev/null || fail "0 mutants: wrong message: $out"
# 0934-zb-slice-a-full-gate3: Cargo.lock is gitignored, so the mutant worktree had none; every cargo call re-resolved the git dependency
# over the network (DNS retries per mutant: the 20m budget ran out) and a resolve error in the test run (rc 101) was counted as KILLED.
# mutate.sh must copy the repo's Cargo.lock into the worktree, run cargo --offline --frozen, and call only a real test failure a kill.
lockok() { [ "$(grep -c '|fetch$' "$1")" = 1 ] && ! grep -vE '^samelock\|(fetch$|.*--offline.*--frozen)' "$1" | grep -q .; }
realcargo=$(command -v cargo); sh="$t/shim"; mkdir -p "$sh"
cat > "$sh/cargo" <<SHIM
#!/usr/bin/env bash
l=nolock; [ -f Cargo.lock ] && cmp -s Cargo.lock "$r/Cargo.lock" && l=samelock
printf '%s|%s\n' "\$l" "\$*" >> "$t/shim.log"
case " \$* " in *" --no-run "*) ;; *" test "*) if [ "\${SHIM_NETFAIL:-}" = 1 ]; then echo "error: failed to get \\\`x\\\` as a dependency: Could not resolve host: index.crates.io"; exit 101; fi ;; esac
"$realcargo" "\$@"
SHIM
chmod +x "$sh/cargo"; setm m.killed; : > "$t/shim.log"
out="$(PATH="$sh:$PATH" run)" || fail "lock/offline: killable mutant not accepted: $out"
grep -q . "$t/shim.log" || fail "lock/offline: cargo shim never ran"
lockok "$t/shim.log" || fail "lock/offline: a cargo call without the repo Cargo.lock or without --offline --frozen: $(cat "$t/shim.log")"
# Spence (820f061 block): the allow-list let any number of online fetches through. A second fetch line must fail the check.
printf 'x\n' >> "$t/cases"; { cat "$t/shim.log"; echo 'samelock|fetch'; } > "$t/shim2.log"
lockok "$t/shim2.log" && fail "lock/offline: a log with two fetch lines was accepted"
out="$(SHIM_NETFAIL=1 PATH="$sh:$PATH" run)" && fail "a test run that failed before any test ran was accepted: $out"
echo "$out" | grep -q "^INVALID 001-add-mul: suite run failed without a test failure" || fail "netfail: wrong report: $out"
echo "$out" | grep -q "^KILLED" && fail "netfail: counted as KILLED: $out"
# Paul (2026-10-24): a dependency fetch failure must never count as KILLED. mutate.sh fetches once up front (FAIL FETCH, exit 2, no
# mutant runs) and builds/tests every mutant with --offline --frozen. The fixture's unfetchable git dependency is a missing local repo.
setm m.killed; cp "$r/Cargo.toml" "$t/Cargo.toml.keep"
printf 'nofetch = { git = "file://%s/no-such-repo" }\n' "$t" >> "$r/Cargo.toml"; cm
out="$(run)" && fail "an unfetchable dependency passed: $out"
echo "$out" | grep -q "^mutate: FAIL FETCH" || fail "fetch failure: no FAIL FETCH line: $out"
echo "$out" | grep -q "^KILLED" && fail "fetch failure was counted as KILLED: $out"
cp "$t/Cargo.toml.keep" "$r/Cargo.toml"; cm
# Spence (2026-10-24, 1207-zb-popguard-full-gate1: 184 `ERR rc=1 line=92 owner=$(inrange ...)` and one `ERR rc=101 line=214` per
# killed mutant): expected failures must not trip the poplog ERR trap; an unexpected one still must. Each runs with a real run log
# (fresh POPLOG_ROOT, not nested in this selftest's log).
logged() { # logged NAME [ENV...]: mutate.sh with its own run log under $t/pl-NAME; prints the output, then the log path on the last line
  local n=$1; shift
  ( unset POPLOG_FILE POPLOG_PID POPLOG_STATE_DIR POPLOG_PROGRESS_FILE POPLOG_RUNNING RUNLOG_INNER RUNLOG_ACTIVE
    env "$@" POPLOG_ROOT="$t/pl-$n" RUNLOG_AGENT=selftest MUTATE_FLOOR_ZERO_OK=1 MUTATE_FLOOR="$t/mfloor.txt" \
      MUTATE_MANIFEST="$r/tools/mutants/manifest.txt" "$here/mutate.sh" 2>&1; echo "rc=$?" )
  find "$t/pl-$n/selftest" -name '*-mutate.log' | head -1
}
setm m.killed
printf 'x\n' >> "$t/cases"; out="$(logged quiet)"; lg=$(tail -1 <<< "$out")
grep -q '^rc=0$' <<< "$out" && grep -q '^KILLED 001-add-mul' <<< "$out" || fail "err-quiet: killable mutant not accepted: $out"
[ -f "$lg" ] || fail "err-quiet: no run log: $out"
grep -q ' ERR ' "$lg" && fail "err-quiet: expected failures (inrange miss, killed rc 101) logged ERR trap lines: $(grep ' ERR ' "$lg")"
printf 'x\n' >> "$t/cases"; out="$(logged loud SHIM_NETFAIL=1 PATH="$sh:$PATH")"; lg=$(tail -1 <<< "$out")
grep -q '^rc=0$' <<< "$out" && fail "err-loud: an unexpected suite-run failure passed: $out"
grep -q '^INVALID 001-add-mul: suite run failed without a test failure' <<< "$out" || fail "err-loud: wrong report: $out"
[ -f "$lg" ] && grep -q ' ERR rc=101 ' "$lg" || fail "err-loud: an unexpected failure (rc 101, no test failed) did not log an ERR trap line: $(cat "$lg" 2>/dev/null | grep ERR)"
# every worktree and target dir was removed, in every case above
n="$(leftovers)"; [ "$(printf '%s' "$n" | tr -d '\n 0')" = "" ] || fail "leftover worktrees or dirs: $n"
want_cases=34; ncases=$(grep -c x "$t/cases" 2>/dev/null || echo 0)
[ "$ncases" = "$want_cases" ] || fail "ran $ncases cases, expected $want_cases (a case was skipped or added without updating want_cases)"
echo "mutate selftest: ok ($ncases cases)"
