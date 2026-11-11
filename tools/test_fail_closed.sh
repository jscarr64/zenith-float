#!/usr/bin/env bash
# ZF-T04 self-test: the gate's search steps and the bench baseline refresh fail closed.
# - The scans that scripts/zenith_prepublish.sh and scripts/ci.sh run (scripts/ci_no_intrinsics.sh,
#   scripts/ci_no_name.sh) see a stub git/rg/grep on PATH: rc 0 = hit (must FAIL), rc 1 = clean (must PASS),
#   rc 2 or 127 = tool error or missing tool (must FAIL with an "error:" line).
# - scripts/bench-compare.sh --update-from: a log that parses to zero rows (empty, or no Criterion lines) is
#   refused and the baseline is left byte-for-byte unchanged; a one-row log still refreshes it.
# - scripts/bench-compare.sh --check and compare mode: the same zero-row inputs are refused and the baseline
#   is left unchanged; a one-row log that matches the baseline still passes (compare mode uses a stub cargo,
#   no real bench run).
# usage: tools/test_fail_closed.sh   (no cargo build; scratch in a mktemp dir; compare mode uses a stub cargo)
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"
t="$(mktemp -d)"
trap 'rm -rf "$t"' EXIT
bad=0
fail() { echo "FAIL: $*" >&2; bad=1; }
ok() { echo "ok   $*"; }

# The callers must still run the scans this test covers (a scan moved back inline would be untested).
for pair in "scripts/zenith_prepublish.sh|scripts/ci_no_intrinsics.sh" "scripts/ci.sh|scripts/ci_no_intrinsics.sh" "scripts/ci.sh|scripts/ci_no_name.sh"; do
  IFS='|' read -r caller scan <<< "$pair"
  rc=0; grep -qF "$scan" "$caller" || rc=$?
  case $rc in
    0) ok "$caller runs $scan" ;;
    1) fail "$caller no longer runs $scan" ;;
    *) fail "grep on $caller exited $rc" ;;
  esac
done

# stub <dir> <rc>: fake git, rg and grep that print nothing and exit rc
stub() { mkdir -p "$1"; for tool in git rg grep; do printf '#!/bin/sh\nexit %s\n' "$2" > "$1/$tool"; chmod +x "$1/$tool"; done; }
for rc in 0 1 2 127; do stub "$t/bin$rc" "$rc"; done

for scan in scripts/ci_no_intrinsics.sh scripts/ci_no_name.sh; do
  for rc in 0 1 2 127; do
    got=0; PATH="$t/bin$rc:$PATH" bash "$scan" > "$t/out" 2>&1 || got=$?
    case $rc in
      1) if [[ $got -eq 0 ]]; then ok "$scan tool rc 1 (clean) passes"; else fail "$scan tool rc 1 (clean): rc=$got, want 0: $(tail -1 "$t/out")"; fi ;;
      0) if [[ $got -ne 0 ]]; then ok "$scan tool rc 0 (hit) fails rc=$got"; else fail "$scan tool rc 0 (hit): rc=0, want fail"; fi ;;
      2|127)
        if [[ $got -eq 0 ]]; then fail "$scan tool rc $rc (error/missing): rc=0, want fail"
        else
          mrc=0; grep -q 'error:' "$t/out" || mrc=$?
          case $mrc in
            0) ok "$scan tool rc $rc fails rc=$got with an error line" ;;
            1) fail "$scan tool rc $rc fails rc=$got but prints no error: line" ;;
            *) fail "grep on $scan output exited $mrc" ;;
          esac
        fi ;;
      *) fail "unexpected stub rc $rc" ;;
    esac
  done
done

# bench-compare.sh --update-from
base="$t/base.tsv"
printf 'benchmark_id\ttime_ns\tthroughput_per_s\nkeep/1\t5\t-\n' > "$base"
before="$(sha256sum < "$base")"
: > "$t/empty.log"
printf 'Compiling zenith-float-num\nno criterion lines here\n' > "$t/junk.log"
for log in empty junk; do
  got=0; REGRESSION_FILE="$base" bash scripts/bench-compare.sh --update-from "$t/$log.log" > "$t/out" 2>&1 || got=$?
  if [[ $got -eq 0 ]]; then fail "bench-compare.sh --update-from $log log (zero rows): rc=0, want fail"; else ok "bench-compare.sh --update-from $log log refused rc=$got"; fi
  if [[ "$(sha256sum < "$base")" == "$before" ]]; then ok "baseline untouched after $log log"; else fail "bench-compare.sh --update-from $log log overwrote the baseline"; fi
done
# bench-compare.sh --check: zero-row logs must be refused with the parse-time error; baseline untouched.
# A later awk "no overlapping benchmarks" exit is not enough: that is a different failure.
for log in empty junk; do
  got=0; REGRESSION_FILE="$base" bash scripts/bench-compare.sh --check "$t/$log.log" > "$t/out" 2>&1 || got=$?
  if [[ $got -eq 0 ]]; then fail "bench-compare.sh --check $log log (zero rows): rc=0, want fail"
  else
    mrc=0; grep -q 'error: no Criterion benchmark rows parsed' "$t/out" || mrc=$?
    case $mrc in
      0) ok "bench-compare.sh --check $log log refused rc=$got with zero-row error" ;;
      1) fail "bench-compare.sh --check $log log fails rc=$got but not with zero-row error: $(tail -1 "$t/out")" ;;
      *) fail "grep on --check $log output exited $mrc" ;;
    esac
  fi
  if [[ "$(sha256sum < "$base")" == "$before" ]]; then ok "baseline untouched after --check $log log"; else fail "bench-compare.sh --check $log log overwrote the baseline"; fi
done

# bench-compare.sh compare mode: stub cargo emits a zero-row log (no real bench run)
mkdir -p "$t/cargoempty" "$t/cargojunk"
printf '#!/bin/sh\nexit 0\n' > "$t/cargoempty/cargo"
chmod +x "$t/cargoempty/cargo"
printf '#!/bin/sh\necho "Compiling zenith-float-num"\necho "no criterion lines here"\nexit 0\n' > "$t/cargojunk/cargo"
chmod +x "$t/cargojunk/cargo"
for kind in empty junk; do
  got=0; PATH="$t/cargo$kind:$PATH" REGRESSION_FILE="$base" bash scripts/bench-compare.sh > "$t/out" 2>&1 || got=$?
  if [[ $got -eq 0 ]]; then fail "bench-compare.sh compare $kind cargo log (zero rows): rc=0, want fail"
  else
    mrc=0; grep -q 'error: no Criterion benchmark rows parsed' "$t/out" || mrc=$?
    case $mrc in
      0) ok "bench-compare.sh compare $kind cargo log refused rc=$got with zero-row error" ;;
      1) fail "bench-compare.sh compare $kind cargo log fails rc=$got but not with zero-row error: $(tail -1 "$t/out")" ;;
      *) fail "grep on compare $kind output exited $mrc" ;;
    esac
  fi
  if [[ "$(sha256sum < "$base")" == "$before" ]]; then ok "baseline untouched after compare $kind cargo log"; else fail "bench-compare.sh compare $kind cargo log overwrote the baseline"; fi
done

# one-row matching log still passes in --check and compare; baseline untouched
printf 'keep/1                  time:   [4.0000 ns 5.0000 ns 6.0000 ns]\n' > "$t/keep.log"
got=0; REGRESSION_FILE="$base" bash scripts/bench-compare.sh --check "$t/keep.log" > "$t/out" 2>&1 || got=$?
if [[ $got -eq 0 && "$(sha256sum < "$base")" == "$before" ]]; then ok "bench-compare.sh --check one-row matching log passes, baseline untouched"
else fail "bench-compare.sh --check one-row matching log: rc=$got, baseline unchanged: $([[ "$(sha256sum < "$base")" == "$before" ]] && echo yes || echo no): $(tail -1 "$t/out")"; fi
mkdir -p "$t/cargokeep"
printf '#!/bin/sh\necho "keep/1                  time:   [4.0000 ns 5.0000 ns 6.0000 ns]"\nexit 0\n' > "$t/cargokeep/cargo"
chmod +x "$t/cargokeep/cargo"
got=0; PATH="$t/cargokeep:$PATH" REGRESSION_FILE="$base" bash scripts/bench-compare.sh > "$t/out" 2>&1 || got=$?
if [[ $got -eq 0 && "$(sha256sum < "$base")" == "$before" ]]; then ok "bench-compare.sh compare one-row matching cargo log passes, baseline untouched"
else fail "bench-compare.sh compare one-row matching cargo log: rc=$got, baseline unchanged: $([[ "$(sha256sum < "$base")" == "$before" ]] && echo yes || echo no): $(tail -1 "$t/out")"; fi

printf 'add/64                  time:   [1.0000 ns 2.0000 ns 3.0000 ns]\n' > "$t/one.log"
got=0; REGRESSION_FILE="$base" bash scripts/bench-compare.sh --update-from "$t/one.log" > "$t/out" 2>&1 || got=$?
if [[ $got -eq 0 && "$(sha256sum < "$base")" != "$before" ]]; then ok "bench-compare.sh --update-from one-row log refreshes the baseline"
else fail "bench-compare.sh --update-from one-row log: rc=$got, baseline changed: $([[ "$(sha256sum < "$base")" != "$before" ]] && echo yes || echo no): $(tail -1 "$t/out")"; fi

if [[ "$bad" -ne 0 ]]; then echo "test_fail_closed: FAIL" >&2; exit 1; fi
echo "test_fail_closed: ok"
