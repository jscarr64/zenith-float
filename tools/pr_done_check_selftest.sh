#!/usr/bin/env bash
# Selftest for tools/pr_done_check.py (Jeff 2026-11-25 2:53 PM CT: "Make sure PR's are closed once merged").
# A PR is DONE only when it is merged into its base, the PR is closed, and its head branch is deleted.
# FAIL cases: an OPEN PR whose head is already on the base (compare status identical or behind), and a landed PR
# (merged, or closed with its head on the base) whose head branch still exists. Read-only: only gh GET calls.
# gh is mocked: a stub on PATH answers from fixture files, so this never touches GitHub.
# PR_DONE_CHECK can point at another copy of pr_done_check.py (the mutation runner uses this).
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
[ -f "$pc" ] || { echo "pr_done_check_selftest: $pc missing"; echo "pr_done_check_selftest: FAIL (0 ok, 1 failed)"; exit 1; }

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
    st=all; prev=""; for a in "$@"; do [ "$prev" = "--state" ] && st="$a"; prev="$a"; done
    python3 -c 'import json,sys; st=sys.argv[2].upper(); d=json.load(open(sys.argv[1])); print(json.dumps([p for p in d if st=="ALL" or p["state"]==st]))' "$f" "$st"; exit 0;;
  "api "*)
    path="$2"
    case "$path" in
      repos/*/compare/*)
        sha="${path##*...}"
        [ -f "$MOCK/compare/$sha" ] || { echo "gh: Not Found (HTTP 404)" >&2; exit 1; }
        printf '{"status":"%s","ahead_by":0,"behind_by":0}\n' "$(cat "$MOCK/compare/$sha")"; exit 0;;
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
  out=$(MOCK="$m" python3 -B "$pc" --repo o/r 2>&1); rc=$?
  if [ "$rc" = "$want" ] && { [ -z "$re" ] || grep -Eq "$re" <<<"$out"; }; then ok "$name (rc=$rc)"
  else bad "$name (want rc=$want${re:+ and /$re/}, got rc=$rc)"; printf '%s\n' "$out" | sed 's/^/     | /'; fi
  if grep -Eqv '^(pr list|api repos/)' "$m/calls.log" 2>/dev/null || grep -Eq -- '(-X|--method|POST|PATCH|DELETE|PUT)' "$m/calls.log" 2>/dev/null; then
    bad "$name: non-read gh call"; sed 's/^/     | /' "$m/calls.log"; fi
}

scenario "no PRs"                         0 '[]'
scenario "open, head ahead of main"       0 "[$(pr 1 OPEN f/a aaaa)]" compare:aaaa=ahead
scenario "open, diverged"                 0 "[$(pr 2 OPEN f/b bbbb)]" compare:bbbb=diverged
scenario "open, head identical to main (merged-open)" 1 "[$(pr 3 OPEN f/c cccc)]" compare:cccc=identical re:'#3.*merged-open'
scenario "open, head behind main (merged-open)"       1 "[$(pr 4 OPEN f/d dddd)]" compare:dddd=behind re:'#4.*merged-open'
scenario "merged and closed, branch deleted"          0 "[$(pr 5 MERGED f/e eeee '"2026-11-07T16:23:25Z"')]"
scenario "merged and closed, branch kept"             1 "[$(pr 6 MERGED f/g gggg '"2026-11-07T16:23:25Z"')]" branch:f/g=gggg re:'#6.*branch-not-deleted'
scenario "merged, branch name reused at other sha"    1 "[$(pr 7 MERGED f/h hhhh '"2026-11-07T16:23:25Z"')]" branch:f/h=zzzz re:'#7.*branch-not-deleted'
scenario "closed unmerged, head not on main, branch kept" 0 "[$(pr 8 CLOSED f/i iiii)]" compare:iiii=diverged branch:f/i=iiii
scenario "closed unmerged, head on main, branch kept (landed outside GitHub)" 1 "[$(pr 9 CLOSED f/j jjjj)]" compare:jjjj=behind branch:f/j=jjjj re:'#9.*branch-not-deleted'
scenario "closed unmerged, head on main, branch deleted" 0 "[$(pr 10 CLOSED f/k kkkk)]" compare:kkkk=identical
scenario "mix: one bad among good"        1 "[$(pr 11 OPEN f/l llll),$(pr 12 MERGED f/m mmmm '"2026-11-01T00:00:00Z"'),$(pr 13 OPEN f/n nnnn)]" compare:llll=ahead compare:nnnn=identical re:'#13.*merged-open'
scenario "gh pr list error fails closed"  2 '[]' flag:fail-list
scenario "branch lookup error fails closed" 2 "[$(pr 14 MERGED f/o oooo '"2026-11-07T16:23:25Z"')]" flag:fail-branch
scenario "open PR, compare 404 fails closed" 2 "[$(pr 15 OPEN f/p pppp)]"
scenario "cross-repo PR is skipped, not failed" 0 "$(printf '[{"number":16,"state":"MERGED","headRefName":"f/q","headRefOid":"qqqq","baseRefName":"main","mergedAt":"2026-11-07T16:23:25Z","isCrossRepository":true,"url":"u"}]')" branch:f/q=qqqq
# --since: landed PRs merged before the date are not checked for branch deletion
m="$work/m-since"; mkdir -p "$m/compare" "$m/branches"; echo "[$(pr 17 MERGED f/r rrrr '"2026-10-01T00:00:00Z"')]" > "$m/prs.json"; echo rrrr > "$m/branches/f__r"
MOCK="$m" python3 -B "$pc" --repo o/r --since 2026-11-01 >/dev/null 2>&1; rc=$?; [ $rc = 0 ] && ok "--since skips older merges (rc=0)" || bad "--since skips older merges (rc=$rc)"
MOCK="$m" python3 -B "$pc" --repo o/r --since 2026-09-01 >/dev/null 2>&1; rc=$?; [ $rc = 1 ] && ok "--since includes newer merges (rc=1)" || bad "--since includes newer merges (rc=$rc)"
# --since boundary: merged ON the date is checked
m="$work/m-since-eq"; mkdir -p "$m/compare" "$m/branches"; echo "[$(pr 18 MERGED f/s ssss '"2026-11-01T09:00:00Z"')]" > "$m/prs.json"; echo ssss > "$m/branches/f__s"
MOCK="$m" python3 -B "$pc" --repo o/r --since 2026-11-01 >/dev/null 2>&1; rc=$?; [ $rc = 1 ] && ok "--since includes a merge on that date (rc=1)" || bad "--since boundary (rc=$rc)"
# two repos: a bad first repo is not hidden by a good second one
m="$work/m-two"; mkdir -p "$m/compare" "$m/branches"; echo "[$(pr 19 OPEN f/t tttt)]" > "$m/prs-o__bad.json"; echo '[]' > "$m/prs-o__good.json"; echo identical > "$m/compare/tttt"
MOCK="$m" python3 -B "$pc" --repo o/bad --repo o/good >/dev/null 2>&1; rc=$?; [ $rc = 1 ] && ok "bad repo then good repo (rc=1)" || bad "bad repo then good repo (rc=$rc)"
# usage
python3 -B "$pc" >/dev/null 2>&1; rc=$?; [ $rc = 2 ] && ok "no --repo is a usage error (rc=2)" || bad "no --repo (rc=$rc)"
echo "pr_done_check_selftest: $([ $fail = 0 ] && echo PASS || echo FAIL) ($pass ok, $fail failed)"
[ $fail = 0 ]
