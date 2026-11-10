#!/usr/bin/env bash
# History: fba8ec9 (slice C) first added this step because `mutate.sh --list`/`--floor` alone never apply a patch, so a stale patch passed
# the gate (slice C merge 8597039, 2026-10-24); Spence chose the design below, which folds in fba8ec9's ID_REGEX (merge resolution).
# Per-slice mutant run for the gate (Spence-approved design, 2026-10-24). It does NOT trust tools/mutate.sh's exit code: it checks the
# result lines itself and FAILS on any of
#   TIMEOUT     the slice ran past its budget: names every id that never ran
#   INCOMPLETE  a manifest id produced no result line (e.g. 257 results for 268 ids): names the missing ids
#   INVALID     an INVALID mutant
#   SURVIVED    a survivor with no equivalent entry (mutate.sh prints SURVIVED-EQUIVALENT for documented ones)
#   BELOW floor the slice's KILLED / KILLED+SURVIVED-EQUIVALENT counts are under its floor (mutate.sh --floor --slice)
# and on a duplicate result id, a result for an id other than the run's --only id, an unreadable id list or floor, or a nonzero runner exit.
# Before a full slice run it clears that slice's logdir of stale *.out files once (clear_outs below; mutate.sh --only keeps other ids' files).
# A result is an UNINDENTED line `KILLED|SURVIVED-EQUIVALENT|SURVIVED|INVALID <id>` (mutate.sh echoes suite output indented; that never counts).
#   TIMEOUT <id> one mutant ran past MUTANT_TIMEOUT_S (each id runs as its own `mutate.sh --slice S --only ID` under that limit)
#   KILLED-OOM <id> the run was SIGKILLed (exit 137) before its limit: not timeout's own -k kill, so the OOM killer or another SIGKILL
# usage: tools/mutate_slice.sh SLICE     one slice (central, or a manifest.d/<slice>.txt name)
#        tools/mutate_slice.sh SLICE ID_REGEX   PARTIAL: only the slice ids matching the ERE; floor and completeness are NOT checked (says so,
#                                        never prints a full pass); a regex matching no id fails. --all takes no regex (the gate never runs partial).
#        tools/mutate_slice.sh --all     every slice from `mutate.sh --slices`, then the whole-manifest completeness check
#        tools/mutate_slice.sh --step-timeout   prints the seconds the gate's mutants step may take (sum of slice timeouts + 600)
# Per-slice timeout = ids x MUTANT_BUDGET_S + MUTATE_SETUP_S (fetch, plus the first build of a slice's suite into the persistent target).
# BUDGET (measured, see CHANGELOG 2026-10-24): BUDGET_S and SETUP_S below come from a real tiny-slice run; the values and their source are
# in the CHANGELOG entry for this file. env: MUTANT_BUDGET_S, MUTATE_SETUP_S (positive integers; override for selftests only), MUTATE_SH;
# tools/full_gate.sh unsets all of them (and MUTANT_TIMEOUT_S, MUTATE_FLOOR_ZERO_OK) for the gate.
# shellcheck source=tools/runlog.sh
. "${BASH_SOURCE[0]%/*}/runlog.sh"; runlog_wrap mutate_slice "$@"   # run log (tools/runlog.sh)
set -u
export LC_ALL=C
M=${MUTATE_SH:-${BASH_SOURCE[0]%/*}/mutate.sh}
# Measured 2026-10-24 (poprun log harper/2026-10-24/1644-zb-mslice-measure, f854a9d, tiny 4-mutant slice from a fresh MUTATE_TARGET,
# one per suite mag/array/ball/harness): per-mutant max 14s (178-U-underflow-zero), whole run 21s (so setup ~1s). Also 1520-zb-mutq-merge-full-gate:
# 184 mutants in a 4m19s step (~1.4s mean). BUDGET_S = max 14s x ~4 margin = 60s; SETUP_S = 300s (fetch + first suite build under load).
BUDGET_S=60; SETUP_S=300
# Per-mutant limit: the same measurement (max 14s, 178-U-underflow-zero, 1644-zb-mslice-measure) x ~4 = 60s. A mutant past it is
# `TIMEOUT <id>` and fails the slice (mutate.sh's own 20m test timeout would count a hang as KILLED). env MUTANT_TIMEOUT_S (selftests only).
MUTANT_TIMEOUT_DEFAULT_S=60
mto=${MUTANT_TIMEOUT_S:-$MUTANT_TIMEOUT_DEFAULT_S}
[[ $mto =~ ^[0-9]{1,6}$ ]] && [ "$mto" -gt 0 ] || { echo "mutate_slice: FAIL MUTANT_TIMEOUT_S '$mto' is not a positive integer"; exit 2; }
budget=${MUTANT_BUDGET_S:-$BUDGET_S}; setup=${MUTATE_SETUP_S:-$SETUP_S}
[[ $budget =~ ^[0-9]{1,6}$ ]] && [ "$budget" -gt 0 ] || { echo "mutate_slice: FAIL MUTANT_BUDGET_S '$budget' is not a positive integer"; exit 2; }
[[ $setup =~ ^[0-9]{1,6}$ ]] || { echo "mutate_slice: FAIL MUTATE_SETUP_S '$setup' is not an integer"; exit 2; }
w=$(mktemp -d "$(zb_scratch)/zb-mutate-slice.XXXXXX") || { echo "mutate_slice: FAIL cannot create scratch dir"; exit 2; }
runlog_on_exit "rm -rf \"$w\""
words() { printf '%s\n' "$1" | awk 'NF' | tr '\n' ' ' | sed 's/ $//'; }
slices() { local sl; sl=$("$M" --slices) || { echo "mutate_slice: FAIL mutate.sh --slices failed" >&2; return 1; }
  [ -n "$sl" ] || { echo "mutate_slice: FAIL mutate.sh --slices listed nothing" >&2; return 1; }; printf '%s\n' "$sl"; }
# slice_ids SLICE: the slice's manifest ids (fails closed)
slice_ids() { "$M" --list --slice "$1" || { echo "mutate_slice: FAIL mutate.sh --list --slice $1 failed" >&2; return 1; }; }
# clear_outs SLICE (Spence, 2026-10-29): each id runs as `mutate.sh --only ID`, which removes only its own <id>.out, so a full slice run
# first clears the slice logdir ONCE, by mutate.sh's full-run rule (regular *.out files directly in <logdir>/<slice>; other files, subdirs
# and other slices are left alone). An id that has left the manifest cannot leave a stale .out behind. The logdir is resolved as mutate.sh
# does: MUTATE_LOGDIR, else MUTATE_PERSIST (default the parent of MUTATE_TARGET) with MUTATE_TARGET, else MUTATE_TMP (default zb_scratch),
# then /mutate-out. A dir that does not exist yet is fine (mutate.sh creates it). A PARTIAL (ID_REGEX) run does not clear.
clear_outs() { local lb d
  if [ -n "${MUTATE_TARGET:-}" ]; then lb=${MUTATE_PERSIST:-$(dirname "$MUTATE_TARGET")}; else lb=${MUTATE_TMP:-$(zb_scratch)}; fi
  d=${MUTATE_LOGDIR:-$lb/mutate-out}/$1; [ -d "$d" ] || return 0
  find "$d" -mindepth 1 -maxdepth 1 -type f -name '*.out' -delete || { echo "mutate_slice: FAIL cannot clear stale .out files in $d"; return 1; }
}
# one SLICE: run and check one slice; its result ids go to $w/SLICE.got
one() {
  local s=$1 re=${2-} partial=0 ids n fl id left lim t0 t1 el irc tmo="" oom="" cut="" nov="" why  fk fke to rc=0 out=$w/$1.out res got nres miss dup oth inv sur nk ne bad=0 tmax
  ids=$(slice_ids "$s") || return 1; ids=$(printf '%s\n' "$ids" | awk 'NF')
  if [ $# -ge 2 ]; then
    partial=1
    irc=0; ids=$(printf '%s\n' "$ids" | grep -E -- "$re") || irc=$?
    [ "$irc" -le 1 ] || { echo "mutate_slice: FAIL ID_REGEX '$re' is not a valid ERE"; return 1; }
    [ -n "$ids" ] || { echo "mutate_slice: FAIL ID_REGEX '$re' matches no id in slice $s"; return 1; }
  fi
  n=$(printf '%s\n' "$ids" | awk 'NF' | wc -l)
  fl=$("$M" --floor --slice "$s") || { echo "mutate_slice: FAIL slice $s: mutate.sh --floor --slice failed"; return 1; }
  fk=$(printf '%s\n' "$fl" | awk '$1 == "killed" { print $2 }'); fke=$(printf '%s\n' "$fl" | awk '$1 == "killed+equivalent" { print $2 }')
  [[ $fk =~ ^[0-9]+$ && $fke =~ ^[0-9]+$ ]] || { echo "mutate_slice: FAIL slice $s: unreadable floor: $(words "$fl")"; return 1; }
  to=$((n * budget + setup))
  echo "mutate_slice: slice $s: $n ids, timeout ${to}s ($n x ${budget}s + ${setup}s setup), floor killed $fk killed+equivalent $fke"
  : > "$w/$s.got"; : > "$w/$s.res"; : > "$w/$s.oth"
  [ "$partial" = 1 ] || clear_outs "$s" || return 1
  if [ "$n" -gt 0 ]; then
    : > "$out"; t0=$SECONDS
    for id in $ids; do
      left=$((to - (SECONDS - t0)))
      if [ "$left" -le 0 ]; then cut+=" $id"; continue; fi
      lim=$mto; [ "$left" -lt "$lim" ] && lim=$left
      t1=$SECONDS; irc=0; timeout -k 30 "$lim" "$M" --slice "$s" --only "$id" > "$out.1" 2>&1 < /dev/null || irc=$?; el=$((SECONDS - t1))
      cat "$out.1"; cat "$out.1" >> "$out"
      # results of this run: unindented lines only, and only for the --only id (anything else is named and fails the slice)
      awk -v want="$id" -v res="$w/$s.res" -v oth="$w/$s.oth" '/^(KILLED|SURVIVED-EQUIVALENT|SURVIVED|INVALID) / { r = $2; sub(/:$/, "", r)
        if (r == want) print $1, r >> res; else print want, r >> oth }' "$out.1"
      # timeout exits 124 at the limit, or 137 when its -k 30 KILL follows (so not before the limit); a 137 before it is another SIGKILL
      if [ "$irc" = 124 ] || { [ "$irc" = 137 ] && [ "$el" -ge "$lim" ]; }; then
        if [ "$lim" = "$mto" ]; then echo "TIMEOUT $id: ran past the ${mto}s per-mutant limit"; tmo+=" $id"; else cut+=" $id"; fi
      elif [ "$irc" = 137 ]; then echo "KILLED-OOM $id: the runner was SIGKILLed (exit 137) after ${el}s, under its ${lim}s limit (OOM killer?)"; oom+=" $id"
      # NO-VERDICT (Harper, 2026-10-29): the --only run ended without a result line for its id (e.g. a manifest check failed first);
      # name the id, its exit code and its own FAIL line here, so INCOMPLETE is never the only trace
      elif ! awk -v want="$id" '$2 == want { f = 1 } END { exit !f }' "$w/$s.res"; then
        why=$(grep -m1 '^mutate: FAIL' "$out.1" || awk 'NF && !/^runlog: / { l = $0 } END { print l }' "$out.1")
        echo "NO-VERDICT $id: mutate.sh --only $id exited $irc before a verdict: ${why:-no output}"; nov+=" $id"
        if [ "$irc" != 0 ] && [ "$rc" = 0 ]; then rc=$irc; fi
      elif [ "$irc" != 0 ] && [ "$rc" = 0 ]; then rc=$irc; fi
    done
  else : > "$out"; fi
  res=$(cat "$w/$s.res")
  got=$(printf '%s\n' "$res" | awk 'NF { print $2 }'); printf '%s\n' "$got" | awk 'NF' > "$w/$s.got"
  nres=$(awk 'NF' "$w/$s.got" | wc -l)
  miss=$(printf '%s\n' "$ids" | awk -v f="$w/$s.got" 'BEGIN { while ((getline l < f) > 0) g[l] = 1 } NF && !($1 in g)')
  dup=$(sort "$w/$s.got" | uniq -d)
  oth=$(awk '{ printf "the --only %s run printed a result for %s; ", $1, $2 }' "$w/$s.oth")
  [ -z "$tmo" ] || { echo "mutate_slice: FAIL TIMEOUT (per-mutant ${mto}s): $(words "$(printf '%s\n' $tmo)") (slice $s)"; bad=1
    miss=$(printf '%s\n' "$miss" | grep -vxF -f <(printf '%s\n' $tmo) || true); }
  # a NO-VERDICT id stays in the INCOMPLETE missing list below: it has no result, and full_gate_selftest (mcase missing, dcase forged) pins that line
  [ -z "$nov" ] || { echo "mutate_slice: FAIL NO-VERDICT: $(words "$(printf '%s\n' $nov)") (slice $s)"; bad=1; }
  [ -z "$oom" ] || { echo "mutate_slice: FAIL KILLED-OOM (SIGKILL before the limit): $(words "$(printf '%s\n' $oom)") (slice $s)"; bad=1
    miss=$(printf '%s\n' "$miss" | grep -vxF -f <(printf '%s\n' $oom) || true); }
  if [ "$partial" = 1 ]; then
    [ -z "$cut" ] || { echo "mutate_slice: FAIL TIMEOUT: slice $s ran past ${to}s; never ran: $(words "$(printf '%s\n' $cut)")"; bad=1; }
  elif [ -n "$cut" ]; then
    echo "mutate_slice: FAIL TIMEOUT: slice $s ran past ${to}s ($n x ${budget}s + ${setup}s): $nres results for $n ids; never ran: $(words "$miss")"; bad=1
  elif [ -n "$miss" ]; then
    echo "mutate_slice: FAIL INCOMPLETE: $nres results for $n ids in slice $s; missing: $(words "$miss")"; bad=1
  fi
  [ -z "$dup" ] || { echo "mutate_slice: FAIL slice $s: duplicate result for $(words "$dup")"; bad=1; }
  [ -z "$oth" ] || { echo "mutate_slice: FAIL slice $s: ${oth%; }"; bad=1; }
  inv=$(printf '%s\n' "$res" | awk '$1 == "INVALID" { print $2 }'); sur=$(printf '%s\n' "$res" | awk '$1 == "SURVIVED" { print $2 }')
  [ -z "$inv" ] || { echo "mutate_slice: FAIL INVALID: $(words "$inv") (slice $s)"; bad=1; }
  [ -z "$sur" ] || { echo "mutate_slice: FAIL SURVIVED (no equivalent entry): $(words "$sur") (slice $s)"; bad=1; }
  nk=$(printf '%s\n' "$res" | awk '$1 == "KILLED"' | wc -l); ne=$(printf '%s\n' "$res" | awk '$1 == "SURVIVED-EQUIVALENT"' | wc -l)
  if [ "$partial" = 0 ] && { [ "$nk" -lt "$fk" ] || [ "$((nk + ne))" -lt "$fke" ]; }; then
    echo "mutate_slice: FAIL BELOW floor: slice $s: $nk killed (floor $fk), $((nk + ne)) killed+equivalent (floor $fke)"; bad=1
  fi
  [ "$rc" = 0 ] || { echo "mutate_slice: FAIL slice $s: mutate.sh exited $rc"; bad=1; }
  tmax=$(awk '/^TIME / { v = $3; sub(/s$/, "", v); if (v + 0 > m) { m = v + 0; id = $2 } } END { if (id != "") print m "s (" id ")"; else print "-" }' "$out")
  if [ "$partial" = 1 ]; then
    echo "mutate_slice: PARTIAL run (ID_REGEX=$re): floor/completeness not checked"
    echo "mutate_slice: slice $s PARTIAL: $nres results for $n matching ids, $nk killed, $ne survived-equivalent, max per-mutant time $tmax (limit ${mto}s) -> $([ "$bad" = 0 ] && echo 'PARTIAL, not a gate pass' || echo FAIL)"
    return "$bad"
  fi
  echo "mutate_slice: slice $s: $nres/$n results, $nk killed, $ne survived-equivalent, max per-mutant time $tmax (budget ${budget}s) -> $([ "$bad" = 0 ] && echo ok || echo FAIL)"
  return "$bad"
}
case ${1:-} in
  --step-timeout)
    sl=$(slices) || exit 1; tot=600
    for s in $sl; do ids=$(slice_ids "$s") || exit 1; n=$(printf '%s\n' "$ids" | awk 'NF' | wc -l); [ "$n" = 0 ] || tot=$((tot + n * budget + setup)); done
    echo "$tot"; exit 0 ;;
  --all)
    [ $# = 1 ] || { echo "mutate_slice: FAIL --all takes no ID_REGEX (the gate never runs a partial slice)"; exit 2; }
    sl=$(slices) || exit 1; frc=0; failed=""
    for s in $sl; do one "$s" || { frc=1; failed+=" $s"; }; done
    # whole-manifest completeness: every id of `mutate.sh --list` belongs to a slice and produced exactly one result
    all=$("$M" --list) || { echo "mutate_slice: FAIL mutate.sh --list failed"; exit 1; }
    cat "$w"/*.got 2>/dev/null | awk 'NF' > "$w/all.got"
    na=$(printf '%s\n' "$all" | awk 'NF' | wc -l); ng=$(awk 'NF' "$w/all.got" | wc -l)
    amiss=$(printf '%s\n' "$all" | awk -v f="$w/all.got" 'BEGIN { while ((getline l < f) > 0) g[l] = 1 } NF && !($1 in g)')
    if [ -n "$amiss" ] || [ "$na" = 0 ]; then echo "mutate_slice: FAIL INCOMPLETE: $ng results for $na manifest ids; missing: $(words "$amiss")"; frc=1; fi
    if [ "$frc" = 0 ]; then echo "mutate_slice: all slices ok ($ng results for $na manifest ids)"; else echo "mutate_slice: FAIL slices:$failed ($ng results for $na manifest ids)"; fi
    exit "$frc" ;;
  ''|-*) echo "usage: tools/mutate_slice.sh SLICE [ID_REGEX] | --all | --step-timeout"; exit 2 ;;
  *) [ $# -le 2 ] || { echo "usage: tools/mutate_slice.sh SLICE [ID_REGEX] | --all | --step-timeout"; exit 2; }; one "$@"; exit $? ;;
esac
