#!/usr/bin/env bash
# Selftest for tools/pr_done_check.py (Jeff 2026-11-25 2:53 PM CT: "Make sure PR's are closed once merged").
# A PR is DONE only when it is merged into its base, the PR is closed, and its head branch is deleted.
# FAIL cases: an OPEN PR whose head is already on the base (compare status identical or behind), and a landed PR
# (merged, or closed with its head on the base) whose head branch still exists. Read-only: only gh GET calls.
# gh is mocked: a stub on PATH answers from fixture files, so this never touches GitHub.
# PR_DONE_CHECK can point at another copy of pr_done_check.py (the mutation runner uses this).
# Each section below landed as its own red commit before the code that makes it pass.
set -u
here="$(cd "$(dirname "$0")" && pwd)"
pc="${PR_DONE_CHECK:-$here/pr_done_check.py}"
root="${ZF_SCRATCH:-$HOME/.cache/zenith-float}/pr-done-selftest"
mkdir -p "$root" || exit 2
work=$(mktemp -d "$root/run.XXXXXX") || exit 2
trap 'rm -rf "$work"' EXIT
pass=0; fail=0
ok()  { pass=$((pass+1)); echo "ok   $1"; }
bad() { fail=$((fail+1)); echo "FAIL $1"; }

# --- mock gh -------------------------------------------------------------------------------------------------------
mkdir -p "$work/bin"
cat > "$work/bin/gh" <<'STUB'
#!/usr/bin/env bash
# fixture layout: $MOCK/prs.json ; $MOCK/compare/<head sha> (status word) ; $MOCK/branches/<name with / as __>
echo "$*" >> "$MOCK/calls.log"
case "$1 $2" in
  "pr list")
    [ -f "$MOCK/fail-list" ] && { echo "gh: HTTP 502" >&2; exit 1; }
    repo=""; prev=""; for a in "$@"; do [ "$prev" = "-R" ] && repo="$a"; prev="$a"; done
    f="$MOCK/prs-${repo//\//__}.json"; [ -f "$f" ] || f="$MOCK/prs.json"
    [ -f "$MOCK/empty-list" ] && exit 0
    [ -f "$MOCK/garbage-list" ] && { echo "not json"; exit 0; }
    [ -f "$MOCK/object-list" ] && { echo '{"number":1}'; exit 0; }
    [ -f "$MOCK/scalar-list" ] && { echo '[1]'; exit 0; }
    st=all; lim=30; prev=""; for a in "$@"; do [ "$prev" = "--state" ] && st="$a"; [ "$prev" = "--limit" ] && lim="$a"; prev="$a"; done
    # like gh: --state filters, --limit caps the rows (gh's own default is 30)
    python3 -c 'import json,sys; st=sys.argv[2].upper(); d=json.load(open(sys.argv[1])); print(json.dumps([p for p in d if st=="ALL" or p["state"]==st][:int(sys.argv[3])]))' "$f" "$st" "$lim"; exit 0;;
  "api "*)
    path="$2"
    case "$path" in
      repos/*/compare/*)
        sha="${path##*...}"
        [ -f "$MOCK/compare/$sha" ] || { echo "gh: Not Found (HTTP 404)" >&2; exit 1; }
        st=$(cat "$MOCK/compare/$sha"); [ "$st" = nostatus ] && { echo '{"ahead_by":0}'; exit 0; }
        printf '{"status":"%s","ahead_by":0,"behind_by":0}\n' "$st"; exit 0;;
      repos/*/git/ref/heads/*)
        name="${path#*/git/ref/heads/}"; f="$MOCK/branches/${name//\//__}"
        [ -f "$MOCK/fail-branch" ] && { echo "gh: HTTP 500" >&2; exit 1; }
        [ -f "$f" ] || { echo "gh: Not Found (HTTP 404)" >&2; exit 1; }
        printf '{"ref":"refs/heads/%s","object":{"sha":"%s"}}\n' "$name" "$(cat "$f")"; exit 0;;
    esac;;
esac
echo "mock gh: unexpected call: $*" >&2; exit 3
STUB
chmod +x "$work/bin/gh"
export PATH="$work/bin:$PATH"

# pr JSON helper: pr NUM STATE HEADREF HEADSHA [MERGED_AT]
pr() { printf '{"number":%s,"state":"%s","headRefName":"%s","headRefOid":"%s","baseRefName":"main","mergedAt":%s,"isCrossRepository":false,"url":"https://example/pull/%s"}' \
  "$1" "$2" "$3" "$4" "${5:-null}" "$1"; }
# scenario NAME WANT_RC 'json array' [compare:sha=status ...] [branch:name=sha ...] [flag:file] -- [expected output regex]
scenario() {
  local name=$1 want=$2 json=$3; shift 3
  local m="$work/m-$name"; mkdir -p "$m/compare" "$m/branches"; printf '%s\n' "$json" > "$m/prs.json"
  local re=""
  while [ $# -gt 0 ]; do
    case "$1" in
      compare:*) kv=${1#compare:}; echo "${kv#*=}" > "$m/compare/${kv%%=*}";;
      branch:*)  kv=${1#branch:}; n=${kv%%=*}; echo "${kv#*=}" > "$m/branches/${n//\//__}";;
      flag:*)    touch "$m/${1#flag:}";;
      re:*)      re=${1#re:};;
    esac; shift
  done
  local out rc
  out=$(MOCK="$m" python3 -B "$pc" --repo o/r ${EXTRA:-} 2>&1); rc=$?
  if [ "$rc" = "$want" ] && { [ -z "$re" ] || grep -Eq "$re" <<<"$out"; }; then ok "$name (rc=$rc)"
  else bad "$name (want rc=$want${re:+ and /$re/}, got rc=$rc)"; printf '%s\n' "$out" | sed 's/^/     | /'; fi
  if grep -Eqv '^(pr list|api repos/)' "$m/calls.log" 2>/dev/null || grep -Eq -- '(-X|--method|POST|PATCH|DELETE|PUT)' "$m/calls.log" 2>/dev/null; then
    bad "$name: non-read gh call"; sed 's/^/     | /' "$m/calls.log"; fi
}

# ---- merged-open: an OPEN PR whose head is already on its base; every gh error fails closed (rc 2) ----
scenario "no PRs"                         0 '[]'
scenario "open, head ahead of main"       0 "[$(pr 1 OPEN f/a aaaa)]" compare:aaaa=ahead
grep -q -- '--limit 1000' "$work/m-open, head ahead of main/calls.log" && ok "gh pr list asks for --limit 1000 by default" || { bad "gh pr list --limit"; sed 's/^/     | /' "$work/m-open, head ahead of main/calls.log"; }
scenario "open, diverged"                 0 "[$(pr 2 OPEN f/b bbbb)]" compare:bbbb=diverged
scenario "open, head identical to main (merged-open)" 1 "[$(pr 3 OPEN f/c cccc)]" compare:cccc=identical re:'#3.*merged-open'
scenario "open, head behind main (merged-open)"       1 "[$(pr 4 OPEN f/d dddd)]" compare:dddd=behind re:'#4.*merged-open'
scenario "mix: one bad among good"        1 "[$(pr 11 OPEN f/l llll),$(pr 13 OPEN f/n nnnn)]" compare:llll=ahead compare:nnnn=identical re:'#13.*merged-open'
scenario "gh pr list error fails closed"  2 '[]' flag:fail-list
scenario "open PR, compare 404 fails closed" 2 "[$(pr 15 OPEN f/p pppp)]"
# Cassie B10: every gh error class gives rc 2
scenario "gh pr list empty output fails closed"  2 '[]' flag:empty-list re:'empty output'
scenario "gh pr list JSON that is not a list fails closed" 2 '[]' flag:object-list
scenario "gh pr list output not JSON fails closed" 2 '[]' flag:garbage-list
scenario "compare JSON without status fails closed" 2 "[$(pr 20 OPEN f/u uuuu)]" compare:uuuu=nostatus
m="$work/m-nogh"; mkdir -p "$m"; out=$(GH="$work/no-such-gh" MOCK="$m" python3 -B "$pc" --repo o/r 2>&1); rc=$?
[ $rc = 2 ] && ok "gh binary missing fails closed (rc=2)" || { bad "gh binary missing (rc=$rc)"; printf '%s\n' "$out" | sed 's/^/     | /'; }
# Cassie B8: a list that reaches --limit may be truncated -> rc 2
many() { python3 -c 'import json,sys; print(json.dumps([{"number":i,"state":"OPEN","headRefName":"f/%d"%i,"headRefOid":"x%d"%i,"baseRefName":"main","mergedAt":None,"isCrossRepository":False,"url":"u"} for i in range(int(sys.argv[1]))]))' "$1"; }
EXTRA="--limit 3" scenario "3 PRs with --limit 3: possibly truncated" 2 "$(many 3)" compare:x0=ahead compare:x1=ahead compare:x2=ahead re:'truncated'
EXTRA="--limit 4" scenario "3 PRs with --limit 4: ok" 0 "$(many 3)" compare:x0=ahead compare:x1=ahead compare:x2=ahead
m="$work/m-1000"; mkdir -p "$m/compare" "$m/branches"; many 1000 > "$m/prs.json"
out=$(MOCK="$m" python3 -B "$pc" --repo o/r 2>&1); rc=$?
[ $rc = 2 ] && grep -q truncated <<<"$out" && ok "1000 PRs at the default --limit: rc 2 (truncated)" || { bad "1000 PRs at the default --limit (rc=$rc)"; printf '%s\n' "$out" | tail -3 | sed 's/^/     | /'; }
EXTRA="--limit 0" scenario "--limit 0 always reads as truncated" 2 '[]'
# two repos: a bad first repo is not hidden by a good second one
m="$work/m-two"; mkdir -p "$m/compare" "$m/branches"; echo "[$(pr 19 OPEN f/t tttt)]" > "$m/prs-o__bad.json"; echo '[]' > "$m/prs-o__good.json"; echo identical > "$m/compare/tttt"
MOCK="$m" python3 -B "$pc" --repo o/bad --repo o/good >/dev/null 2>&1; rc=$?; [ $rc = 1 ] && ok "bad repo then good repo (rc=1)" || bad "bad repo then good repo (rc=$rc)"
# usage
python3 -B "$pc" >/dev/null 2>&1; rc=$?; [ $rc = 2 ] && ok "no --repo is a usage error (rc=2)" || bad "no --repo (rc=$rc)"
# ---- branch-not-deleted: a landed PR (merged, or closed with its head on the base) whose branch still exists ----
scenario "merged and closed, branch deleted"          0 "[$(pr 5 MERGED f/e eeee '"2026-11-07T16:23:25Z"')]"
scenario "merged and closed, branch kept"             1 "[$(pr 6 MERGED f/g gggg '"2026-11-07T16:23:25Z"')]" branch:f/g=gggg re:'#6.*branch-not-deleted'
scenario "merged, branch name reused at other sha"    1 "[$(pr 7 MERGED f/h hhhh '"2026-11-07T16:23:25Z"')]" branch:f/h=zzzz re:'#7.*branch-not-deleted'
scenario "closed unmerged, head not on main, branch kept" 0 "[$(pr 8 CLOSED f/i iiii)]" compare:iiii=diverged branch:f/i=iiii
scenario "closed unmerged, head on main, branch kept (landed outside GitHub)" 1 "[$(pr 9 CLOSED f/j jjjj)]" compare:jjjj=behind branch:f/j=jjjj re:'#9.*branch-not-deleted'
scenario "closed unmerged, head on main, branch deleted" 0 "[$(pr 10 CLOSED f/k kkkk)]" compare:kkkk=identical
scenario "mix: merged with branch kept among good"   1 "[$(pr 21 OPEN f/v vvvv),$(pr 12 MERGED f/m mmmm '"2026-11-01T00:00:00Z"')]" compare:vvvv=ahead branch:f/m=mmmm re:'#12.*branch-not-deleted'
scenario "branch lookup error fails closed" 2 "[$(pr 14 MERGED f/o oooo '"2026-11-07T16:23:25Z"')]" flag:fail-branch
scenario "cross-repo PR is skipped, not failed" 0 "$(printf '[{"number":16,"state":"MERGED","headRefName":"f/q","headRefOid":"qqqq","baseRefName":"main","mergedAt":"2026-11-07T16:23:25Z","isCrossRepository":true,"url":"u"},{"number":22,"state":"OPEN","headRefName":"f/w","headRefOid":"wwww","baseRefName":"main","mergedAt":null,"isCrossRepository":true,"url":"u"}]')" branch:f/q=qqqq compare:wwww=identical
# --since: landed PRs merged before the date (UTC) are not checked for branch deletion
m="$work/m-since"; mkdir -p "$m/compare" "$m/branches"; echo "[$(pr 17 MERGED f/r rrrr '"2026-10-01T00:00:00Z"')]" > "$m/prs.json"; echo rrrr > "$m/branches/f__r"
MOCK="$m" python3 -B "$pc" --repo o/r --since 2026-11-01 >/dev/null 2>&1; rc=$?; [ $rc = 0 ] && ok "--since skips older merges (rc=0)" || bad "--since skips older merges (rc=$rc)"
MOCK="$m" python3 -B "$pc" --repo o/r --since 2026-09-01 >/dev/null 2>&1; rc=$?; [ $rc = 1 ] && ok "--since includes newer merges (rc=1)" || bad "--since includes newer merges (rc=$rc)"
m="$work/m-since-eq"; mkdir -p "$m/compare" "$m/branches"; echo "[$(pr 18 MERGED f/s ssss '"2026-11-01T09:00:00Z"')]" > "$m/prs.json"; echo ssss > "$m/branches/f__s"
MOCK="$m" python3 -B "$pc" --repo o/r --since 2026-11-01 >/dev/null 2>&1; rc=$?; [ $rc = 1 ] && ok "--since includes a merge on that date (rc=1)" || bad "--since boundary (rc=$rc)"
# the UTC date decides: merged 2026-11-24 9 PM CST = 2026-11-25T03:00Z counts as 11-25
m="$work/m-utc"; mkdir -p "$m/compare" "$m/branches"; echo "[$(pr 23 MERGED f/x xxxx '"2026-11-25T03:00:00Z"')]" > "$m/prs.json"; echo xxxx > "$m/branches/f__x"
MOCK="$m" python3 -B "$pc" --repo o/r --since 2026-11-25 >/dev/null 2>&1; rc=$?; [ $rc = 1 ] && ok "--since uses the UTC date of mergedAt (rc=1)" || bad "--since UTC date (rc=$rc)"
# Cassie B9: a --since that is not a real YYYY-MM-DD date is a usage error, not a silent skip
for bd in yesterday 2026-1-5 2026-02-30 20261105 ""; do
  MOCK="$m" python3 -B "$pc" --repo o/r --since "$bd" >/dev/null 2>&1; rc=$?; [ $rc = 2 ] && ok "--since '$bd' is a usage error (rc=2)" || bad "--since '$bd' (want rc=2, got rc=$rc)"
done
# ---- round 2 (Cassie R2-B3): malformed entries and dates fail closed; a landed PR with no date is still checked ----
scenario "list entry that is not an object fails closed (R2-B3a)" 2 '[]' flag:scalar-list
scenario "mergedAt that is a number fails closed (TypeError, Cassie Y07)" 2 "[$(pr 30 MERGED f/w wwww 5)]" branch:f/w=wwww
EXTRA="--since 2026-11-01" scenario "landed PR with no merge/close date under --since: branch still checked (Cassie Y08)" 1 "[$(pr 31 MERGED f/y yyyy)]" branch:f/y=yyyy re:'#31.*branch-not-deleted'
scenario "JSON object instead of a list: says so (G03)" 2 '[]' flag:object-list re:'not a JSON list'
echo "pr_done_check_selftest: $([ $fail = 0 ] && echo PASS || echo FAIL) ($pass ok, $fail failed)"
[ $fail = 0 ]
