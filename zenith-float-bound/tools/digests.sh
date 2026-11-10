#!/usr/bin/env bash
# S13 cross-target grid-digest job (TEST_PLAN P10, section 13). The harness test grid_digests_write_the_digest_file_for_the_script writes one
# digest row per suite; this script runs it for every target in tools/digest_targets.txt, in debug and release, with 1 and 7 test threads,
# and requires every run to produce exactly the rows of tools/digest_expected.txt (so also identical to each other).
#   tools/digests.sh run [MANIFEST] OUTDIR       one file OUTDIR/<target>-<profile>-t<threads>.txt per run
#   tools/digests.sh compare OUTDIR [EXPECTED]   exit 1 and print the first differing line of every run that differs from EXPECTED
#   tools/digests.sh all                         run into a temp dir, compare, remove it (the gate step)
# The '#' lines of a run file (arch, pointer width) are not compared. Env: DIGESTS_CARGO (default cargo), DIGESTS_TMP (where the test writes;
# default <target dir>/tmp, or <target dir>/<triple>/tmp with --target). A manifest line is a target triple or `host`; a target needs its std installed and a way to
# execute its binaries (a native one, or CARGO_TARGET_<TRIPLE>_RUNNER); a run that fails makes the job fail, it is never skipped.
# shellcheck source=tools/runlog.sh
. "${BASH_SOURCE[0]%/*}/runlog.sh"; runlog_wrap digests "$@"   # run log (tools/runlog.sh)
set -u
export LC_ALL=C
root="$(cd "$(dirname "$0")/.." && pwd)"
cargo_cmd=${DIGESTS_CARGO:-cargo}
tbase=${CARGO_TARGET_DIR:-$root/target}   # the harness writes to $CARGO_TARGET_TMPDIR: <target dir>/tmp for the host, <target dir>/<triple>/tmp for --target
threads_list="1 7"
test_filter=grid_digests_write

rows() { grep -v '^#' "$1"; }

do_run() { # do_run MANIFEST OUTDIR
  local manifest="$1" out="$2" n=0 target profile th label args
  mkdir -p "$out" || return 1
  while IFS= read -r target; do
    case "$target" in ''|'#'*) continue ;; esac
    n=$((n + 1))
    for profile in debug release; do
      for th in $threads_list; do
        label="$target-$profile-t$th"
        args=(test)
        [ "$profile" = release ] && args+=(--release)
        [ "$target" != host ] && args+=(--target "$target")
        tmpdir=${DIGESTS_TMP:-$tbase/tmp}
        [ "$target" != host ] && tmpdir=${DIGESTS_TMP:-$tbase/$target/tmp}
        args+=(--test harness "$test_filter" -- "--test-threads=$th")
        rm -f "$tmpdir/grid_digests.txt"
        echo "digests: $label"
        (cd "$root" && "$cargo_cmd" "${args[@]}") >"$out/$label.log" 2>&1 || { echo "digests: RUN FAILED: $label (see $out/$label.log)"; tail -5 "$out/$label.log"; return 1; }
        [ -s "$tmpdir/grid_digests.txt" ] || { echo "digests: RUN WROTE NO DIGEST FILE: $label"; return 1; }
        cp "$tmpdir/grid_digests.txt" "$out/$label.txt" || return 1
      done
    done
  done < "$manifest"
  [ "$n" -gt 0 ] || { echo "digests: the manifest $manifest lists no target"; return 1; }
}

do_compare() { # do_compare OUTDIR EXPECTED
  local dir="$1" want="$2" f bad=0 n=0 first
  [ -f "$want" ] || { echo "digests: no expected file $want"; return 1; }
  [ -n "$(rows "$want")" ] || { echo "digests: the expected file has no rows"; return 1; }
  for f in "$dir"/*.txt; do
    [ -e "$f" ] || break
    n=$((n + 1))
    if [ -z "$(rows "$f")" ]; then echo "digests: DIFF $(basename "$f"): no digest rows"; bad=$((bad + 1)); continue; fi
    if ! diff <(rows "$want") <(rows "$f") >/dev/null; then
      first=$(diff <(rows "$want") <(rows "$f") | grep '^[<>]' | head -2 | tr '\n' ' ')
      echo "digests: DIFF $(basename "$f"): first differing lines (< expected, > got): $first"
      bad=$((bad + 1))
    fi
  done
  [ "$n" -gt 0 ] || { echo "digests: no run files in $dir"; return 1; }
  if [ "$bad" -gt 0 ]; then echo "digests: FAILED ($bad of $n runs differ)"; return 1; fi
  echo "digests: ok ($n runs, all equal to $(basename "$want"))"
}

case "${1:-}" in
  run) if [ "$#" = 3 ]; then do_run "$2" "$3"; else do_run "$root/tools/digest_targets.txt" "${2:?usage: digests.sh run [MANIFEST] OUTDIR}"; fi ;;
  compare) do_compare "${2:?usage: digests.sh compare OUTDIR [EXPECTED]}" "${3:-$root/tools/digest_expected.txt}" ;;
  all) w=$(mktemp -d "$(zb_scratch)"/zb-digests.XXXXXX); runlog_on_exit "rm -rf \"$w\""; do_run "$root/tools/digest_targets.txt" "$w" && do_compare "$w" "$root/tools/digest_expected.txt" ;;
  *) echo "usage: tools/digests.sh run [MANIFEST] OUTDIR | compare OUTDIR [EXPECTED] | all"; exit 2 ;;
esac
