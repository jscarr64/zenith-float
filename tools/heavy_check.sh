#!/usr/bin/env bash
# Heavy-tier check (ZB-P03, Paul 2026-10-30). libtest has no runtime skip: an early-return SKIP reports "passed", so the gate reads the output.
#   --heavy LOG          PASS only with exactly one `HEAVY s6_heavy_rows_replay_through_zb: rows run=30 passed=30 ` line and no SKIP line.
#   --logs DIR TAG [SELF] FAIL (NOT-RUN) when any DIR/TAG-*.log (except TAG-SELF.log) has a SKIP word; FAIL when there are no logs at all.
#   --pin PINFILE FILE   PASS only when sha256(FILE) equals PINFILE's `  matrix_rows_heavy.txt` line.
# Exit 0 pass, 1 fail, 2 usage or grep error.
set -u
export LC_ALL=C
skipgrep() { local rc=0; grep -nw SKIP "$@" || rc=$?; [ "$rc" -le 1 ] || { echo "heavy_check: grep failed (rc $rc)"; return 2; }; return "$rc"; }
case ${1:-} in
  --heavy)
    [ $# -eq 2 ] || { echo "usage: $0 --heavy LOG"; exit 2; }
    [ -f "$2" ] || { echo "heavy_check: NOT-RUN: no heavy log $2"; exit 1; }
    s=0; skipgrep "$2" || s=$?
    case $s in 0) echo "heavy_check: NOT-RUN: SKIP line in $2"; exit 1 ;; 1) ;; *) exit 2 ;; esac
    n=0; n=$(grep -c '^HEAVY s6_heavy_rows_replay_through_zb: rows run=30 passed=30 ' "$2") || true
    [ "$n" = 1 ] || { echo "heavy_check: FAIL: want exactly 1 'rows run=30 passed=30' line in $2, got ${n:-0}"; exit 1; }
    echo "heavy_check: PASS rows run=30 passed=30 ($2)" ;;
  --logs)
    [ $# -ge 3 ] && [ $# -le 4 ] || { echo "usage: $0 --logs DIR TAG [SELF]"; exit 2; }
    files=()
    for f in "$2/$3"-*.log "$2/$3"-*.out; do [ -f "$f" ] || continue; [ "$f" = "$2/$3-${4:-}.log" ] && continue; files+=("$f"); done
    [ ${#files[@]} -gt 0 ] || { echo "heavy_check: FAIL: no gate logs $2/$3-*.log (nothing scanned)"; exit 1; }
    s=0; skipgrep "${files[@]}" /dev/null || s=$?
    case $s in 0) echo "heavy_check: NOT-RUN: SKIP line(s) above in gate logs"; exit 1 ;; 1) ;; *) exit 2 ;; esac
    echo "heavy_check: PASS no SKIP in ${#files[@]} gate logs" ;;
  --pin)
    [ $# -eq 3 ] || { echo "usage: $0 --pin PINFILE FILE"; exit 2; }
    [ -f "$3" ] || { echo "heavy_check: FAIL: no heavy file $3"; exit 1; }
    want=$(awk '$2 == "matrix_rows_heavy.txt" { print $1 }' "$2") || exit 2
    got=$(sha256sum "$3" | awk '{ print $1 }') || exit 2
    [ -n "$want" ] && [ "$want" = "$got" ] || { echo "heavy_check: FAIL: $3 sha256 $got != pin '${want}' ($2)"; exit 1; }
    echo "heavy_check: PASS pin $got" ;;
  *) echo "usage: $0 --heavy LOG | --logs DIR TAG [SELF] | --pin PINFILE FILE"; exit 2 ;;
esac
