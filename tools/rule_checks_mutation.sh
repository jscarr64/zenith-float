#!/usr/bin/env bash
# Hand mutation test for the shell/Python rule checks (no cargo-mutants for these; plan v2 step 1).
# Each row flips one comparison, threshold, guard or exit code in a scratch copy of a check and runs that
# check's selftest against the mutant (via DRIFT_CHECK / PR_DONE_CHECK / PRE_COMMIT_HOOK / PRE_PUSH_HOOK).
# KILLED = the selftest fails on the mutant. Any SURVIVED row fails this script (rc 1).
# Usage: tools/rule_checks_mutation.sh [-j N]      Scratch: ${ZF_SCRATCH:-$HOME/.cache/zenith-float} (never /tmp)
set -u
cd "$(dirname "$0")/.." || exit 2
jobs=4; [ "${1:-}" = "-j" ] && jobs=$2
base="${ZF_SCRATCH:-$HOME/.cache/zenith-float}"; mkdir -p "$base" || exit 2   # a fresh ZF_SCRATCH may not exist yet (Cassie N2)
work=$(mktemp -d "$base/mutation.XXXXXX") || exit 2
trap 'rm -rf "$work"' EXIT
# id | file | occurrence | old ==> new
rows() { cat <<'T'
A01|tools/drift_check.py|1|if now - last <= window: ==> if now - last < window:
A02|tools/drift_check.py|1|if now - last <= window: ==> if now - last <= window + 1:
A03|tools/drift_check.py|1|"--active-window-min", type=float, default=1440) ==> "--active-window-min", type=float, default=1441)
A04|tools/drift_check.py|1|if newest is not None and now - newest < max_age: ==> if newest is not None and now - newest <= max_age:
A05|tools/drift_check.py|1|if newest is not None and now - newest < max_age: ==> if False:
A06|tools/drift_check.py|1|times.append(newest) ==> pass
A07|tools/drift_check.py|1|times.append(os.stat(p).st_mtime) ==> pass
A08|tools/drift_check.py|1|for f in ("index", os.path.join("logs", "HEAD")): ==> for f in ("index",):
A09|tools/drift_check.py|1|for f in ("index", os.path.join("logs", "HEAD")): ==> for f in (os.path.join("logs", "HEAD"),):
A10|tools/drift_check.py|1|return True, "activity unknown ==> return False, "activity unknown
A11|tools/drift_check.py|1|if how == "scan" and not a.all: ==> if how == "scan":
A12|tools/drift_check.py|1|if how == "scan" and not a.all: ==> if not a.all:
A13|tools/drift_check.py|1|explicit = list(a.repo) + list(a.include) ==> explicit = list(a.repo)
A14|tools/drift_check.py|1|if l.strip() and not l.strip().startswith("#")] ==> if l.strip()]
A15|tools/drift_check.py|1|                win = window ==>                 win = None
A16|tools/drift_check.py|1|and now - when > window: ==> and now - when >= window:
A17|tools/drift_check.py|1|and now - when > window: ==> and now - when > window + 1:
A18|tools/drift_check.py|1|if window is not None and name != head_ref and ==> if window is not None and
A19|tools/drift_check.py|1|and when is not None and now ==> and now
A20|tools/drift_check.py|1|    if window < 0: ==>     if window < -60:
A21|tools/drift_check.py|1|            if key in done: ==>             if False:
A22|tools/drift_check.py|1|        explicit += [l.strip() for l in open(f) ==>         explicit = [l.strip() for l in open(f)
A23|tools/drift_check.py|1|times = [int(git(repo, "log", "-1", "--format=%ct", "HEAD").strip())] ==> times = [0]
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
D22|tools/drift_check.py|1|if detached and head not in tips: ==> if False:
D23|tools/drift_check.py|1|if (common, name) in seen: ==> if False:
D24|tools/drift_check.py|1|seen.add((common, name)) ==> pass
D25|tools/drift_check.py|1|a.pre_commit, seen, win) ==> a.pre_commit, None, win)
G01|tools/pr_done_check.py|1|except OSError as e: ==> except ImportError as e:
G02|tools/pr_done_check.py|1|    if not out: ==>     if False:
G03|tools/pr_done_check.py|1|if not isinstance(prs, list): ==> if False:
G04|tools/pr_done_check.py|1|if len(prs) >= limit: ==> if len(prs) > limit:
G05|tools/pr_done_check.py|1|"--limit", str(limit) ==> "--limit", str(limit + 1)
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
N01|tools/drift_check.py|1|if not os.path.lexists(os.path.join(c, ".git")): ==> if not os.path.isdir(os.path.join(c, ".git")):
N02|tools/drift_check.py|1|unreadable.append(c) ==> pass
N03|tools/drift_check.py|1|    if unreadable: ==>     if False:
P01|tools/pr_done_check.py|1|LANDED = ("identical", "behind") ==> LANDED = ("identical",)
P02|tools/pr_done_check.py|1|LANDED = ("identical", "behind") ==> LANDED = ("behind",)
P03|tools/pr_done_check.py|1|LANDED = ("identical", "behind") ==> LANDED = ("identical", "behind", "ahead")
P04|tools/pr_done_check.py|1|if state == "OPEN": ==> if state == "MERGED":
P05|tools/pr_done_check.py|1|landed = state == "MERGED" or ( ==> landed = state == "MERGED" and (
P06|tools/pr_done_check.py|1|if landed and branch_exists ==> if branch_exists
P07|tools/pr_done_check.py|1|if pr.get("isCrossRepository"): ==> if False:
P08|tools/pr_done_check.py|1|when < since ==> when <= since
P09|tools/pr_done_check.py|1|if allow_404 and "HTTP 404" in err: ==> if allow_404:
P10|tools/pr_done_check.py|1|    return 1 if bad else 0 ==>     return 0
P11|tools/pr_done_check.py|1|bad |= b ==> bad = b
P12|tools/pr_done_check.py|1|"--state", "all" ==> "--state", "open"
V01|tools/pr_done_check.py|1|if a.since is not None and not valid_date(a.since): ==> if False:
V02|tools/pr_done_check.py|1|r"\d{4}-\d{2}-\d{2}" ==> r".*"
V03|tools/pr_done_check.py|1|        datetime.date.fromisoformat(s) ==>         pass
X01|tools/drift_check.py|1|return len(out), min(int(t) for t in out) ==> return len(out), max(int(t) for t in out)
X02|tools/drift_check.py|1|oldest = m if oldest is None else min(oldest, m) ==> oldest = m if oldest is None else max(oldest, m)
X03|tools/drift_check.py|1|os.path.realpath(git(path, "rev-parse", "--show-toplevel").strip()) == os.path.realpath(path) ==> True
X04|tools/drift_check.py|1|n, _ = unpushed(repo, ref, up) if head else (0, None) ==> n, _ = unpushed(repo, ref, up)
X05|tools/drift_check.py|1|if not os.path.isdir(d): ==> if False:
X06|tools/drift_check.py|1|except GitError as e: ==> except ImportError as e:
X07|tools/drift_check.py|1|if upstream.startswith("refs/remotes/") and _ref_exists(repo, upstream): ==> if upstream.startswith("refs/remotes/"):
X07b|tools/drift_check.py|1|if upstream.startswith("refs/remotes/") and _ref_exists(repo, upstream): ==> if _ref_exists(repo, upstream):
X08|tools/drift_check.py|1|if detached and head not in tips: ==> if detached:
X09|tools/drift_check.py|1|"--ignore-submodules=none" ==> "--ignore-submodules=all"
X11|tools/pr_done_check.py|1|except (GhError, KeyError, ValueError, TypeError) as e: ==> except (GhError,) as e:
X12|tools/pr_done_check.py|1|ap.add_argument("--limit", type=int, default=1000 ==> ap.add_argument("--limit", type=int, default=1
X13|.githooks/pre-push|1|[ -z "${lsha:-}" ] && continue ==> :
X14|.githooks/pre-push|1|base="${RATCHET_BASE:-c74e834^}" ==> base="${RATCHET_BASE:-HEAD}"
X15|.githooks/pre-push|1|trap 'rm -rf "$tmp"' EXIT ==> :
X16|.githooks/pre-commit|1||| { echo "pre-commit: not in a git work tree" >&2; exit 2; } ==> || true
X18|tools/pr_done_check.py|1|    if p.returncode != 0: ==>     if p.returncode > 1:
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
