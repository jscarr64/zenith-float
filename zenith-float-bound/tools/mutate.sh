#!/usr/bin/env bash
# Hand-written mutant runner (TEST_PLAN section 11, default plan D3). No `cfg` switch is added to src/; a mutant is a patch.
#
# usage: tools/mutate.sh [--only MUTANT_ID] [--slice SLICE] [--list] [--floor] [--slices]
#   --slices prints every slice name (central = manifest.txt, then each manifest.d/<slice>.txt) and exits.
#   --slice SLICE limits --list, --floor and a run to that slice's ids; --floor/the run's floor check then use only that slice's floor
#   (floor.txt for central, floor.d/SLICE.txt for a slice, 0 if it has none). An empty slice runs 0 mutants and passes a 0 floor.
#   Every mutant prints `TIME <id> <seconds>s` after its result (tools/mutate_slice.sh reports the max; the gate budget comes from it).
#   --list prints every mutant id (central manifest, then the per-slice includes) and exits; --floor prints the summed floor and exits.
# PER-SLICE INCLUDES (S0, 2026-10-24), so parallel slices never edit a shared registry:
#   manifest.d/<slice>.txt    more manifest rows (same 4 columns; the patch path is relative to the manifest's directory, e.g.
#                             arith/201-x.patch). Its FIRST non-blank line must be `# range: LO-HI`: every id's leading number must lie in
#                             LO..=HI, and no two include files may have overlapping ranges. Ids stay unique over all files.
#                             Each row's patch must be under <slice>/ (the slice owns it); a central-manifest id may not sit in any slice's
#                             range; ids sharing a NUMBER (205-a, 205-b) are duplicates; a row without 4 tab-separated columns FAILS.
#   equivalent.d/<slice>.txt  more EQUIVALENT.txt entries (same format: ID<TAB>REASON; the id must be in the manifest and in the slice's range).
#   floor.d/<slice>.txt       (needs manifest.d/<slice>.txt; keys killed and killed+equivalent, each exactly once) the slice's share of the floor (same KEY<TAB>MIN rows as floor.txt; 0 allowed here). The floor checked is
#                             floor.txt plus every share. A slice raises its own share, never floor.txt.
#   env: MUTATE_REPO=DIR            repo to mutate (default: the parent of tools/). Mutants apply to its COMMITTED HEAD.
#        MUTATE_MANIFEST=FILE       default $MUTATE_REPO/tools/mutants/manifest.txt (patches are read from that file's directory)
#        MUTATE_EQUIVALENT=FILE     default: EQUIVALENT.txt next to the manifest
#        MUTATE_TIMEOUT=DURATION    per `cargo test` run of a suite (default 20m); a timeout counts as a kill
#        MUTATE_BUILD_TIMEOUT=DUR   per build (default 20m)
#        MUTATE_TMP=DIR             parent of the throwaway worktree and target dir (default "$(zb_scratch)", tools/scratch.sh)
#        MUTATE_PERSIST=DIR         parent of that fixed worktree (default: the parent of MUTATE_TARGET, so not /tmp, which a reboot wipes)
#        MUTATE_TARGET=DIR          optional PERSISTENT cargo target dir, kept across mutants and runs (outside the repo). The worktree then has a
#                                   fixed path ($MUTATE_PERSIST/zb-mutate.persist, still created fresh from committed HEAD and removed after each
#                                   mutant), so cargo rebuilds only what the patch touched. Unset: a fresh target dir per mutant, removed after it.
#
# For each manifest row it creates a throwaway `git worktree` under /tmp (never inside the repo), applies ONE patch, builds the suite,
# runs `cargo test --test <suite>` under `timeout`, and REQUIRES the suite to fail. The worktree is always removed (EXIT/INT/TERM).
#   KILLED               a test failed ("error: test failed, to rerun pass"), or the run timed out; any other run error is INVALID
#   SURVIVED             the suite passed on the mutant: exit 1 (unless the id is in EQUIVALENT.txt: SURVIVED-EQUIVALENT, not a failure)
#   INVALID              the patch does not apply or the mutant does not build (proves nothing): exit 1
# Exit 0 only if at least one mutant ran and every mutant was KILLED or SURVIVED-EQUIVALENT (0 mutants is a failure), and, for a full
# run (no --only), killed and killed+equivalent reach tools/mutants/floor.txt (MUTATE_FLOOR overrides).
# shellcheck source=tools/runlog.sh
. "${BASH_SOURCE[0]%/*}/runlog.sh"; runlog_wrap mutate "$@"   # run log (tools/runlog.sh)
set -u
export LC_ALL=C
repo=${MUTATE_REPO:-$(cd "$(dirname "$0")/.." && pwd)}
manifest=${MUTATE_MANIFEST:-$repo/tools/mutants/manifest.txt}
equiv=${MUTATE_EQUIVALENT:-$(dirname "$manifest")/EQUIVALENT.txt}
mfloor=${MUTATE_FLOOR:-$(dirname "$manifest")/floor.txt}
ttl=${MUTATE_TIMEOUT:-20m}; bttl=${MUTATE_BUILD_TIMEOUT:-20m}
only=""; list=0; slice=
while [ $# -gt 0 ]; do
  case "$1" in
    --only) only=${2:-}; shift 2 ;;
    --list) list=1; shift ;;
    --floor) list=2; shift ;;
    --slices) list=3; shift ;;
    --slice) slice=${2:-}; [ -n "$slice" ] || { echo "mutate: --slice needs a name"; exit 2; }; shift 2 ;;
    *) echo "usage: tools/mutate.sh [--only MUTANT_ID] [--slice SLICE] [--list] [--floor] [--slices]"; exit 2 ;;
  esac
done
[ -f "$manifest" ] || { echo "mutate: no manifest at $manifest"; exit 2; }
rows=$(grep -v -E '^[[:space:]]*(#|$)' "$manifest" || true)
mdir=$(dirname "$manifest")
# badrows FILE ROWS: every row needs exactly 4 non-empty tab-separated columns (a row with no tab used to slip through --list)
badrows() { local bad; bad=$(printf '%s\n' "$2" | awk -F'\t' 'NF && (NF != 4 || $1 == "" || $2 == "" || $3 == "" || $4 == "") && n++ < 3')   # n<3 in awk, not | head -3: SIGPIPE rc 141 under pipefail
  [ -z "$bad" ] || { printf 'mutate: FAIL %s: row needs 4 tab-separated columns: %s\n' "$1" "$bad"; exit 1; }; }
badrows "$manifest" "$rows"
slmap=$(printf '%s\n' "$rows" | awk -F'\t' 'NF { print $2 "\tcentral" }')   # id<TAB>slice, for --slice
# per-slice includes (header above): append each file's rows after checking its `# range: LO-HI` line and its ids against it
ranges=""
for inc in "$mdir"/manifest.d/*.txt; do
  [ -f "$inc" ] || continue
  rl=$(grep -m1 -v -E '^[[:space:]]*$' "$inc" || true)   # -m1, not | head -1: SIGPIPE rc 141 under pipefail (Paul 1146-full_gate)
  if [[ ! $rl =~ ^#\ range:\ ([0-9]{1,6})-([0-9]{1,6})$ ]]; then echo "mutate: FAIL $inc: first line must be '# range: LO-HI'"; exit 1; fi
  lo=$((10#${BASH_REMATCH[1]})); hi=$((10#${BASH_REMATCH[2]}))
  [ "$lo" -le "$hi" ] || { echo "mutate: FAIL $inc: range $lo-$hi is empty"; exit 1; }
  while read -r olo ohi oname; do
    [ -n "$olo" ] || continue
    if [ "$lo" -le "$ohi" ] && [ "$olo" -le "$hi" ]; then echo "mutate: FAIL $inc: range $lo-$hi overlaps $oname ($olo-$ohi)"; exit 1; fi
  done <<< "$ranges"
  ranges+="$lo $hi $(basename "$inc")"$'\n'
  irows=$(grep -v -E '^[[:space:]]*(#|$)' "$inc" || true)
  badrows "$inc" "$irows"
  islice=$(basename "$inc" .txt)
  while IFS=$'\t' read -r _ iid ipatch _; do
    [ -n "$iid" ] || continue
    num=${iid%%-*}
    if [[ ! $num =~ ^[0-9]{1,6}$ ]] || [ "$((10#$num))" -lt "$lo" ] || [ "$((10#$num))" -gt "$hi" ]; then
      echo "mutate: FAIL $inc: id $iid is outside the file's range $lo-$hi"; exit 1
    fi
    # owned by the slice: its patch lives under the slice's own patch dir
    case "$ipatch" in "$islice"/*) ;; *) echo "mutate: FAIL $inc: id $iid: patch $ipatch is not under $islice/ (not owned by slice $islice)"; exit 1 ;; esac
  done <<< "$irows"
  [ -n "$irows" ] && rows+=$'\n'"$irows" && slmap+=$'\n'$(printf '%s\n' "$irows" | awk -F'\t' -v s="$islice" 'NF { print $2 "\t" s }')
done
rows=$(printf '%s\n' "$rows" | grep -v -E '^[[:space:]]*$' || true)
# inrange NUM: prints the include file whose range holds NUM (empty if none)
# Contract: rc 0 on a hit AND on a miss (a miss is the normal case for central ids, so it must not trip the poplog
# ERR trap); rc 2 + message only on a real failure (a malformed range row), which call sites leave unmasked so it is loud.
inrange() {
  local olo ohi oname
  while read -r olo ohi oname; do
    [ -n "$olo" ] || continue
    [[ $olo =~ ^[0-9]+$ && $ohi =~ ^[0-9]+$ && $1 =~ ^[0-9]+$ ]] || { echo "mutate: internal: inrange $1: bad range row '$olo $ohi $oname'" >&2; return 2; }
    # base 10 forced here: a leading-zero id (0213) is never octal
    if (( 10#$1 >= 10#$olo && 10#$1 <= 10#$ohi )); then echo "$oname"; return 0; fi
  done <<< "$ranges"
  return 0
}
# a central-manifest id must not sit inside a slice's range: that number is owned by the slice's include file
while IFS=$'\t' read -r _ cid _ _; do
  [ -n "$cid" ] || continue
  cnum=${cid%%-*}; [[ $cnum =~ ^[0-9]{1,6}$ ]] || { echo "mutate: FAIL $manifest: id $cid does not start with a number"; exit 1; }
  owner=$(inrange "$((10#$cnum))"); irc=$?   # a miss is rc 0 (empty); nonzero is a real failure: left unmasked so the ERR trap logs it
  [ "$irc" -eq 0 ] || { echo "mutate: FAIL inrange $cid: rc $irc"; exit 1; }
  [ -z "$owner" ] || { echo "mutate: FAIL $manifest: id $cid is inside the range of manifest.d/$owner (owned by slice ${owner%.txt}); move it there"; exit 1; }
done <<< "$(grep -v -E '^[[:space:]]*(#|$)' "$manifest" || true)"
# manifest integrity, before any mutant runs: every id exactly once, every referenced patch file present
dups=$(printf '%s\n' "$rows" | cut -f2 | grep -v '^$' | sort | uniq -d)
if [ -n "$dups" ]; then printf 'mutate: FAIL duplicate manifest id %s\n' $dups; exit 1; fi
# ... and every id NUMBER exactly once: suffixed ids (205-a, 205-b) must not slip past the full-id check
ndups=$(printf '%s\n' "$rows" | cut -f2 | grep -v '^$' | sed -E 's/-.*//; s/^0+([0-9])/\1/' | sort | uniq -d)
if [ -n "$ndups" ]; then printf 'mutate: FAIL duplicate manifest id number %s (ids sharing a number, e.g. a suffix, count as duplicates)\n' $ndups; exit 1; fi
missing=0
while IFS=$'\t' read -r _ mid mpatch _; do
  if [ -n "$mpatch" ] && [ ! -f "$(dirname "$manifest")/$mpatch" ]; then echo "mutate: FAIL manifest id $mid names missing patch $mpatch"; missing=1; fi
done <<< "$rows"
[ "$missing" = 0 ] || exit 1
ids=$(printf '%s\n' "$rows" | cut -f2)
equivs=()
[ -f "$equiv" ] && equivs+=("$equiv")
for e in "$mdir"/equivalent.d/*.txt; do [ -f "$e" ] && equivs+=("$e"); done
# equivalent entries: ID<TAB>REASON, a known id, a non-empty reason; an equivalent.d/<slice>.txt needs manifest.d/<slice>.txt and only
# lists ids in that slice's range
for e in "${equivs[@]}"; do
  eslice=""; erange=""
  if [ "$e" != "$equiv" ]; then
    eslice=$(basename "$e" .txt)
    [ -f "$mdir/manifest.d/$eslice.txt" ] || { echo "mutate: FAIL $e: no manifest.d/$eslice.txt for this slice"; exit 1; }
  fi
  while IFS= read -r eline; do
    case "$eline" in ''|'#'*|[[:space:]]*) continue ;; esac
    eid=${eline%%$'\t'*}; ereason=${eline#*$'\t'}
    { [ "$eid" != "$eline" ] && [ -n "${ereason//[[:space:]]/}" ]; } || { echo "mutate: FAIL $e: entry '$eline' needs ID<TAB>REASON"; exit 1; }
    # here-string, not `printf | grep -q`: under pipefail grep -q exiting at its match SIGPIPEs printf (rc 141), a false miss (Harper, 2026-10-29)
    grep -qxF -- "$eid" <<< "$ids" || { echo "mutate: FAIL $e: equivalent id $eid is not in the manifest"; exit 1; }
    if [ -n "$eslice" ]; then
      enum=${eid%%-*}; erange=$(inrange "$((10#$enum))"); irc=$?
      [ "$irc" -eq 0 ] || { echo "mutate: FAIL inrange $eid: rc $irc"; exit 1; }
      [ "$erange" = "$eslice.txt" ] || { echo "mutate: FAIL $e: id $eid is not in slice $eslice's range"; exit 1; }
    fi
  done < "$e"
done
# floor.d shares: a matching manifest.d slice, only the keys killed and killed+equivalent, each exactly once
for f in "$mdir"/floor.d/*.txt; do
  [ -f "$f" ] || continue
  [ -f "$mdir/manifest.d/$(basename "$f")" ] || { echo "mutate: FAIL $f: no manifest.d/$(basename "$f") for this slice"; exit 1; }
  fbad=$(awk -F'\t' '/^[[:space:]]*(#|$)/ { next } $1 != "killed" && $1 != "killed+equivalent" { print "unknown key " $1; next } { n[$1]++ }
    END { for (k in n) if (n[k] > 1) print "duplicate key " k; if (!n["killed"]) print "no killed row"; if (!n["killed+equivalent"]) print "no killed+equivalent row" }' "$f")
  [ -z "$fbad" ] || { echo "mutate: FAIL $f: $(printf '%s' "$fbad" | tr '\n' ';')"; exit 1; }
done
if [ "$list" = 3 ]; then echo central; for inc in "$mdir"/manifest.d/*.txt; do [ -f "$inc" ] && basename "$inc" .txt; done; exit 0; fi
if [ -n "$slice" ]; then
  [ "$slice" = central ] || [ -f "$mdir/manifest.d/$slice.txt" ] || { echo "mutate: FAIL no slice $slice (central or manifest.d/$slice.txt)"; exit 2; }
  sids=$(printf '%s\n' "$slmap" | awk -F'\t' -v s="$slice" '$2 == s { print $1 }')
  rows=$(printf '%s\n' "$rows" | awk -F'\t' -v keep="$sids" 'BEGIN { n = split(keep, a, "\n"); for (i = 1; i <= n; i++) k[a[i]] = 1 } NF && ($2 in k)')
  ids=$sids
fi
if [ "$list" = 1 ]; then [ -z "$ids" ] || printf '%s\n' "$ids"; exit 0; fi
# floor_total KEY: floor.txt's row plus every floor.d share (each a plain integer, a share may be 0); fails closed
floor_total() {
  local f v tot=0 first=1
  for f in "$mfloor" "$mdir"/floor.d/*.txt; do
    [ -f "$f" ] || { [ "$first" = 1 ] && { echo "mutate: FAIL NOFLOOR: mutant floor $f is missing or unreadable" >&2; return 1; }; continue; }
    v=$(awk -F'\t' -v k="$1" '$1 == k { print $2; found = 1; exit } END { exit !found }' "$f") || { echo "mutate: FAIL NOFLOOR: no '$1' row in $f" >&2; return 1; }
    case "$v" in *$'\r') echo "mutate: FAIL BADCOUNT: floor row ends in CR (CRLF line ending): '${v%$'\r'}' in $f; convert it to LF" >&2; return 1 ;; esac
    [[ $v =~ ^[0-9]{1,18}$ ]] || { echo "mutate: FAIL BADCOUNT: floor value '$v' in $f is not a plain integer" >&2; return 1; }
    if [ "$first" = 1 ] && [[ $v =~ ^0+$ ]] && [ "${MUTATE_FLOOR_ZERO_OK:-}" != 1 ]; then echo "mutate: FAIL BADCOUNT: floor value '$v' in $f must be 1 or more" >&2; return 1; fi
    tot=$((tot + 10#$v)); first=0
  done
  printf '%s\n' "$tot"
}
# floor_slice KEY: the --slice share (floor.txt for central, with its 1-or-more rule; floor.d/SLICE.txt for a slice, 0 if absent)
floor_slice() {
  local f v
  if [ "$slice" = central ]; then f=$mfloor; else f="$mdir/floor.d/$slice.txt"; [ -f "$f" ] || { echo 0; return 0; }; fi
  [ -f "$f" ] || { echo "mutate: FAIL NOFLOOR: mutant floor $f is missing or unreadable" >&2; return 1; }
  v=$(awk -F'\t' -v k="$1" '$1 == k { print $2; found = 1; exit } END { exit !found }' "$f") || { echo "mutate: FAIL NOFLOOR: no '$1' row in $f" >&2; return 1; }
  [[ $v =~ ^[0-9]{1,18}$ ]] || { echo "mutate: FAIL BADCOUNT: floor value '$v' in $f is not a plain integer" >&2; return 1; }
  if [ "$slice" = central ] && [[ $v =~ ^0+$ ]] && [ "${MUTATE_FLOOR_ZERO_OK:-}" != 1 ]; then echo "mutate: FAIL BADCOUNT: floor value '$v' in $f must be 1 or more" >&2; return 1; fi
  printf '%s\n' "$((10#$v))"
}
if [ -n "$slice" ]; then floor_total() { floor_slice "$@"; }; fi
if [ "$list" = 2 ]; then
  fk=$(floor_total killed) || exit 1; fke=$(floor_total killed+equivalent) || exit 1
  printf 'killed %s\nkilled+equivalent %s\n' "$fk" "$fke"; exit 0
fi

base=${MUTATE_TMP:-$(zb_scratch)}; case "$base" in /home/jsc64/MathWorks*|"$repo"*) echo "mutate: MUTATE_TMP must be outside the repo"; exit 2 ;; esac
persist=${MUTATE_TARGET:-}
if [ -n "$persist" ]; then
  case "$persist" in /*) ;; *) echo "mutate: MUTATE_TARGET must be an absolute path"; exit 2 ;; esac
  case "$persist" in /home/jsc64/MathWorks*|"$repo"|"$repo"/*) echo "mutate: MUTATE_TARGET must be outside the repo"; exit 2 ;; esac
  mkdir -p "$persist" || { echo "mutate: cannot create MUTATE_TARGET $persist"; exit 2; }
  pbase=${MUTATE_PERSIST:-$(dirname "$persist")}
  case "$pbase" in /*) ;; *) echo "mutate: MUTATE_PERSIST must be an absolute path"; exit 2 ;; esac
  case "$pbase" in /home/jsc64/MathWorks*|"$repo"|"$repo"/*) echo "mutate: MUTATE_PERSIST must be outside the repo"; exit 2 ;; esac
  mkdir -p "$pbase" || { echo "mutate: cannot create MUTATE_PERSIST $pbase"; exit 2; }
fi
newtgt() { if [ -n "$persist" ]; then tgt=$persist; else tgt=$(mktemp -d "$base/zb-mutate-target.XXXXXX"); fi; }
newwork() { if [ -n "$persist" ]; then work="$pbase/zb-mutate.persist"; git -C "$repo" worktree remove --force "$work" >/dev/null 2>&1; rm -rf "$work"; git -C "$repo" worktree prune >/dev/null 2>&1; else work=$(mktemp -d "$base/zb-mutate.XXXXXX"); rmdir "$work"; fi; }
work=""; tgt=""
cleanup() {
  if [ -n "$work" ]; then
    git -C "$repo" worktree remove --force "$work" >/dev/null 2>&1
    rm -rf "$work"; git -C "$repo" worktree prune >/dev/null 2>&1
  fi
  [ -n "$tgt" ] && [ -z "$persist" ] && rm -rf "$tgt"
  work=""; tgt=""
}
runlog_on_exit cleanup; trap 'cleanup; exit 130' INT TERM

head_sha=$(git -C "$repo" rev-parse HEAD 2>/dev/null) || { echo "mutate: FAIL cannot read HEAD of $repo"; exit 2; }
# Mutants run on a worktree of the COMMITTED HEAD, so a dirty tree (edits or untracked files, e.g. new tests not yet committed) would
# test other code than the tree on disk. Refuse it (Spence, 2026-10-23: 2102-zb-mut180 SURVIVED on HEAD 7191a1f, before mag_sort was committed).
dirty=$(git -C "$repo" status --porcelain 2>&1) || { echo "mutate: FAIL git status failed in $repo: $dirty"; exit 2; }
if [ -n "$dirty" ]; then echo "mutate: FAIL dirty tree in $repo (mutants test the committed HEAD $head_sha; commit or stash first):"; sed -n '1,20s/^/    /p' <<< "$dirty"; exit 2; fi
# Cargo.lock is gitignored, so a worktree of HEAD has none: without it every cargo call re-resolves the git dependency over the network
# (0934-zb-slice-a-full-gate3: DNS retries on every mutant ran out the 20m budget, and resolve errors were miscounted). Use the repo's lock,
# offline and locked; no lock FAILS (fail closed).
echo "mutate: testing committed HEAD $head_sha"
# Fetch every dependency ONCE, up front (Paul, 2026-10-24): a network or registry failure inside a mutant's `cargo test` used to exit
# nonzero and read as KILLED. Here it FAILS the run before any mutant; each mutant then builds and runs --offline --frozen against the
# fetched Cargo.lock (copied into its worktree), so it cannot fetch at all. Cargo.lock is .gitignored, so the fetch keeps the tree clean.
if ! fout=$(cd "$repo" && timeout "$bttl" cargo fetch 2>&1 </dev/null); then
  echo "mutate: FAIL FETCH: cargo fetch failed in $repo; no mutant ran (a fetch failure is never a kill):"; printf '%s\n' "$fout" | tail -5 | sed 's/^/    /'; exit 2
fi
lock="$repo/Cargo.lock"; [ -f "$lock" ] || { echo "mutate: FAIL FETCH: cargo fetch left no Cargo.lock in $repo"; exit 2; }
# Per-ID test output (Harper, 2026-10-29): each mutant's `cargo test` output is kept as $logdir/<slice>/<id>.out (it used to be one
# run.log the next mutant overwrote), and a KILLED result is followed by ONE indented line naming the failing tests (cargo's
# `---- NAME stdout ----` blocks, capped at 5 plus `+N more`; `(timed out)` for a timeout kill, so EVERY KILLED line has one). The result line and the summary line are unchanged; mutate_slice.sh counts
# only unindented result lines, so the extra line is ignored by every caller. env MUTATE_LOGDIR (default MUTATE_TMP/mutate-out, or MUTATE_PERSIST/mutate-out with MUTATE_TARGET).
lbase=$base; [ -z "$persist" ] || lbase=$pbase   # persistent mode never writes under MUTATE_TMP (it may be unwritable)
logdir=${MUTATE_LOGDIR:-$lbase/mutate-out}/${slice:-all}; mkdir -p "$logdir" || { echo "mutate: cannot create log dir $logdir"; exit 2; }
# stale outputs (Spence review, 2026-10-29): a full run first removes every *.out in THIS slice's logdir (only regular files directly in it;
# other files, subdirs and other slices are left alone), so an .out from an earlier HEAD cannot pass for a current one. An --only run
# removes just its own <id>.out, since mutate_slice.sh runs a slice as one --only run per id and must keep the earlier ids' files.
if [ -n "$only" ]; then [[ $only =~ ^[A-Za-z0-9._-]+$ && $only != *..* ]] || { echo "mutate: bad --only id $only"; exit 2; }; rm -f -- "$logdir/$only.out" || { echo "mutate: cannot clear $logdir/$only.out"; exit 2; }
else find "$logdir" -mindepth 1 -maxdepth 1 -type f -name '*.out' -delete || { echo "mutate: cannot clear stale .out files in $logdir"; exit 2; }; fi
# failnames FILE: the failing test names in a cargo test log, at most 5, then `+N more`
failnames() { awk '/^---- .* stdout ----$/ { n++; if (n <= 5) s = s (s == "" ? "" : " ") $2 } END { if (n > 5) s = s " +" (n - 5) " more"; print (n ? s : "(none parsed)") }' "$1"; }
nk=0; ns=0; ne=0; ni=0; total=0
newtgt
if [ -n "$rows" ]; then
while IFS="$(printf '\t')" read -r fnpath id patch suite; do
  [ -n "$only" ] && [ "$only" != "$id" ] && continue
  total=$((total + 1)); t0=$SECONDS
  if [ -z "$fnpath" ] || [ -z "$id" ] || [ -z "$patch" ] || [ -z "$suite" ]; then echo "INVALID $id: manifest row needs 4 tab-separated columns"; ni=$((ni + 1)); continue; fi
  pfile=$(dirname "$manifest")/$patch
  [ -f "$pfile" ] || { echo "INVALID $id: patch file $pfile not found"; ni=$((ni + 1)); continue; }
  newwork
  if ! git -C "$repo" worktree add --detach -q "$work" HEAD >/dev/null 2>&1; then echo "INVALID $id: cannot create worktree"; ni=$((ni + 1)); work=""; continue; fi
  src=${fnpath%%:*}
  if [ ! -f "$work/$src" ]; then echo "INVALID $id: $src does not exist at HEAD"; ni=$((ni + 1)); cleanup; newtgt; continue; fi
  if ! git -C "$work" apply "$pfile" >/dev/null 2>&1; then echo "INVALID $id: patch does not apply"; ni=$((ni + 1)); cleanup; newtgt; continue; fi
  cp "$lock" "$work/Cargo.lock" || { echo "INVALID $id: cannot copy Cargo.lock into the worktree"; ni=$((ni + 1)); cleanup; newtgt; continue; }
  if ! (cd "$work" && CARGO_TARGET_DIR="$tgt" timeout "$bttl" cargo test --offline --frozen --test "$suite" --no-run >"$tgt/build.log" 2>&1 </dev/null); then
    echo "INVALID $id: mutant does not build suite $suite (see: first lines below)"; head -5 "$tgt/build.log" | sed 's/^/    /'; ni=$((ni + 1))
  else
    # rc captured explicitly: a failing suite is the EXPECTED outcome for a mutant, so it must not trip the poplog ERR trap
    rc=0; (cd "$work" && CARGO_TARGET_DIR="$tgt" timeout "$ttl" cargo test --offline --frozen --test "$suite" >"$tgt/run.log" 2>&1 </dev/null) || rc=$?
    cp "$tgt/run.log" "$logdir/$id.out" 2>/dev/null || :
    if [ "$rc" = 0 ]; then
      if [ "${#equivs[@]}" -gt 0 ] && awk -F'\t' -v id="$id" '!/^[[:space:]]*(#|$)/ && $1 == id { f = 1 } END { exit !f }' "${equivs[@]}"; then
        echo "SURVIVED-EQUIVALENT $id ($fnpath, suite $suite): documented in EQUIVALENT.txt or equivalent.d"; ne=$((ne + 1))
      else
        echo "SURVIVED $id ($fnpath, suite $suite): the suite passed on the mutant"; ns=$((ns + 1))
      fi
    elif [ "$rc" = 124 ]; then echo "KILLED $id ($fnpath, suite $suite): timed out after $ttl"; nk=$((nk + 1))
      echo "    failing tests: (timed out) (output: $logdir/$id.out)"   # every KILLED line is followed by one failing-tests line
    # a kill is a test binary that ran and failed: cargo then prints "error: test failed, to rerun pass". Anything else (resolve,
    # network, lock, I/O error) proves nothing about the mutant.
    elif grep -q '^error: test failed, to rerun pass' "$tgt/run.log"; then echo "KILLED $id ($fnpath, suite $suite): suite failed (rc $rc)"; nk=$((nk + 1))
      echo "    failing tests: $(failnames "$tgt/run.log") (output: $logdir/$id.out)"
    else (exit "$rc")   # unexpected: re-raise rc so the ERR trap still records it
      echo "INVALID $id: suite run failed without a test failure (rc $rc; see: first lines below)"; head -5 "$tgt/run.log" | sed 's/^/    /'; ni=$((ni + 1)); fi
  fi
  cleanup; newtgt
  echo "TIME $id $((SECONDS - t0))s"
done <<< "$rows"
fi
if [ -n "$only" ] && [ "$total" = 0 ]; then echo "mutate: no mutant with id $only"; exit 2; fi
echo "mutate: $total mutants: $nk killed, $ns survived, $ne survived-equivalent, $ni invalid"
[ "$total" -gt 0 ] || [ -n "$slice" ] || { echo "mutate: FAIL 0 mutants ran (an empty manifest proves nothing)"; exit 1; }
[ "$ns" = 0 ] && [ "$ni" = 0 ] || exit 1
# Mutant floor (Paul, 2026-10-23): a full run must kill at least the checked-in number, so a manifest that loses rows, or a suite that
# stops killing, cannot pass. --only runs one row by design and skips it (it still needs its one mutant run, above). Fails closed.
[ -n "$only" ] && exit 0
# floor.txt plus the per-slice shares in floor.d (floor_total, above; a real floor.txt value of 0 still fails)
fk=$(floor_total killed) || exit 1
fke=$(floor_total killed+equivalent) || exit 1
nke=$((nk + ne)); rc=0
if [ "$nk" -ge "$fk" ]; then :; else echo "mutate: FAIL BELOW floor: $nk killed < $fk ($mfloor)"; rc=1; fi
if [ "$nke" -ge "$fke" ]; then :; else echo "mutate: FAIL BELOW floor: $nke killed+equivalent < $fke ($mfloor)"; rc=1; fi
[ "$rc" = 0 ] && echo "mutate: floor ok ($nk killed >= $fk, $nke killed+equivalent >= $fke)"
exit "$rc"
