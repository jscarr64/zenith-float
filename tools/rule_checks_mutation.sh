#!/usr/bin/env bash
# Hand mutation test for the shell/Python rule checks (no cargo-mutants for these; plan v2 step 1).
# Each row flips one comparison, threshold, guard or exit code in a scratch copy of a check and runs that
# check's selftest against the mutant (via DRIFT_CHECK / PR_DONE_CHECK / PRE_COMMIT_HOOK / PRE_PUSH_HOOK).
# KILLED = the selftest fails on the mutant. Any SURVIVED row fails this script (rc 1).
# Usage: tools/rule_checks_mutation.sh [-j N]      Scratch: ${ZF_SCRATCH:-$HOME/.cache/zenith-float} (never /tmp)
set -u
cd "$(dirname "$0")/.." || exit 2
jobs=4; [ "${1:-}" = "-j" ] && jobs=$2
work=$(mktemp -d "${ZF_SCRATCH:-$HOME/.cache/zenith-float}/mutation.XXXXXX") || exit 2
trap 'rm -rf "$work"' EXIT
# id | file | occurrence | old ==> new
rows() { cat <<'T'
C01|.githooks/pre-commit|1| --pre-commit --repo ==>  --repo
C02|.githooks/pre-commit|1|exec python3 -B "$top/tools/drift_check.py" --pre-commit --repo "$top" ==> python3 -B "$top/tools/drift_check.py" --pre-commit --repo "$top"; exit 0
D01|tools/drift_check.py|1|if n > max_ahead: ==> if n >= max_ahead:
D02|tools/drift_check.py|1|if n > max_ahead: ==> if n > max_ahead + 1:
D03|tools/drift_check.py|2|if n > max_ahead: ==> if n >= max_ahead:
D04|tools/drift_check.py|2|if n > max_ahead: ==> if n > max_ahead + 1:
D05|tools/drift_check.py|1|elif n and now - oldest > max_age: ==> elif n and now - oldest >= max_age:
D06|tools/drift_check.py|1|elif n and now - oldest > max_age: ==> elif n and now - oldest > max_age + 60:
D07|tools/drift_check.py|2|elif n and now - oldest > max_age: ==> elif n and now - oldest >= max_age:
D08|tools/drift_check.py|2|elif n and now - oldest > max_age: ==> elif n and now - oldest > max_age + 60:
D09|tools/drift_check.py|1|if d is not None and now - d > max_age: ==> if d is not None and now - d >= max_age:
D10|tools/drift_check.py|1|if d is not None and now - d > max_age: ==> if d is not None and now - d > max_age + 60:
D11|tools/drift_check.py|1|        if n >= max_ahead: ==>         if n > max_ahead:
D12|tools/drift_check.py|1|"--max-ahead", type=int, default=1) ==> "--max-ahead", type=int, default=2)
D13|tools/drift_check.py|1|"--max-age-min", type=float, default=30) ==> "--max-age-min", type=float, default=31)
D14|tools/drift_check.py|1|max_age = a.max_age_min * 60 ==> max_age = a.max_age_min * 61
D15|tools/drift_check.py|1|"--untracked-files=all" ==> "--untracked-files=no"
D16|tools/drift_check.py|1|spec = [ref, "--not", "--remotes"] ==> spec = [ref, "--not", "--all"]
D17|tools/drift_check.py|1|e.pop(k, None) ==> pass
D18|tools/drift_check.py|1|if xy[0] in "RC": ==> if False:
D19|tools/drift_check.py|1|m = index_mtime ==> m = 0
D20|tools/drift_check.py|1|    return 1 if bad else 0 ==>     return 0
D21|tools/drift_check.py|1|bad |= check_repo ==> bad = check_repo
D23|tools/drift_check.py|1|if (common, name) in seen: ==> if False:
D24|tools/drift_check.py|1|seen.add((common, name)) ==> pass
D25|tools/drift_check.py|1|lambda s: print("drift_check: " + s), seen) ==> lambda s: print("drift_check: " + s))
D22|tools/drift_check.py|1|if detached and head not in tips: ==> if False:
H01|.githooks/pre-push|1|unset GIT_DIR GIT_WORK_TREE GIT_INDEX_FILE GIT_COMMON_DIR GIT_OBJECT_DIRECTORY GIT_ALTERNATE_OBJECT_DIRECTORIES ==> :
H02|.githooks/pre-push|1|(cd "$t" && bash tools/rules_scan.sh ==> (cd "$src" && bash tools/rules_scan.sh
H03|.githooks/pre-push|1|no_suppress_scan.py --root "$t" ==> no_suppress_scan.py --root "$src"
H04|.githooks/pre-push|1|if [ "$lsha" = "$zero" ] ==> if false
H05|.githooks/pre-push|1|else echo "pre-push: FAIL rules-scan"; rc=1; fi ==> else echo "pre-push: FAIL rules-scan"; fi
H06|.githooks/pre-push|1|else echo "pre-push: FAIL no-suppress"; rc=1; fi ==> else echo "pre-push: FAIL no-suppress"; fi
H07|.githooks/pre-push|1|has no tools/rules_scan.sh"; rc=1 ==> has no tools/rules_scan.sh"
H08|.githooks/pre-push|1|checkout -q --detach "$lsha" ==> checkout -q --detach HEAD
H09|.githooks/pre-push|1|scratch="${ZF_SCRATCH:-$HOME/.cache/zenith-float}/pre-push" ==> scratch="$HOME/.cache/zenith-float/pre-push"
H10|.githooks/pre-push|1|exit $rc ==> exit 0
H12|.githooks/pre-push|1|rc=1; continue ==> continue
H13|.githooks/pre-push|1|n=$((n + 1)); t="$tmp/$n" ==> t="$tmp/$lsha"
P01|tools/pr_done_check.py|1|LANDED = ("identical", "behind") ==> LANDED = ("identical",)
P02|tools/pr_done_check.py|1|LANDED = ("identical", "behind") ==> LANDED = ("behind",)
P03|tools/pr_done_check.py|1|LANDED = ("identical", "behind") ==> LANDED = ("identical", "behind", "ahead")
P04|tools/pr_done_check.py|1|if state == "OPEN": ==> if state == "MERGED":
P05|tools/pr_done_check.py|1|landed = state == "MERGED" or ( ==> landed = state == "MERGED" and (
P06|tools/pr_done_check.py|1|if landed and branch_exists ==> if branch_exists
P07|tools/pr_done_check.py|1|if pr.get("isCrossRepository"): ==> if False:
P08|tools/pr_done_check.py|1|when[:10] < since ==> when[:10] <= since
P09|tools/pr_done_check.py|1|if allow_404 and "HTTP 404" in err: ==> if allow_404:
P10|tools/pr_done_check.py|1|    return 1 if bad else 0 ==>     return 0
P11|tools/pr_done_check.py|1|bad |= b ==> bad = b
P12|tools/pr_done_check.py|1|"--state", "all" ==> "--state", "open"
T
}
one() {  # id file occ old new
  local id=$1 f=$2 occ=$3 old=$4 new=$5 d="$work/$1"
  mkdir -p "$d/scratch"; local m="$d/$(basename "$f")"
  python3 -B - "$f" "$m" "$occ" "$old" "$new" <<'PY' || { echo "$id ERROR pattern not found: $old"; return; }
import sys; src, dst, occ, old, new = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5]
s = open(src).read(); i = -1
for _ in range(occ):
    i = s.find(old, i + 1)
    if i < 0: sys.exit(1)
open(dst, "w").write(s[:i] + new + s[i + len(old):])
PY
  chmod +x "$m"
  local st env
  case "$f" in
    tools/drift_check.py) st=tools/drift_check_selftest.sh; env="DRIFT_CHECK=$m";;
    .githooks/pre-commit) st=tools/drift_check_selftest.sh; env="PRE_COMMIT_HOOK=$m";;
    tools/pr_done_check.py) st=tools/pr_done_check_selftest.sh; env="PR_DONE_CHECK=$m";;
    .githooks/pre-push) st=tools/pre_push_selftest.sh; env="PRE_PUSH_HOOK=$m";;
  esac
  if env "$env" ZF_SCRATCH="$d/scratch" bash "$st" >"$d/log" 2>&1; then
    echo "$id SURVIVED $f: '$old' -> '$new'"
  else
    echo "$id KILLED   $f: '$old' -> '$new' ($(grep -m1 '^FAIL' "$d/log" | cut -c1-90))"
  fi
}
export -f one; export work
out=$(rows | while IFS='|' read -r id f occ rest; do old=${rest%% ==> *}; new=${rest#* ==> }; printf '%s\0%s\0%s\0%s\0%s\0' "$id" "$f" "$occ" "$old" "$new"; done \
  | xargs -0 -n5 -P "$jobs" bash -c 'one "$@"' _ | sort)
echo "$out"
total=$(echo "$out" | grep -c .); killed=$(echo "$out" | grep -c ' KILLED '); surv=$(echo "$out" | grep -c ' SURVIVED '); err=$(echo "$out" | grep -c ' ERROR ')
echo "rule_checks_mutation: $killed/$total killed, $surv survived, $err errors"
[ "$surv" = 0 ] && [ "$err" = 0 ]
