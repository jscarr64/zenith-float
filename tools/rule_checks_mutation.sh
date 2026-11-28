#!/usr/bin/env bash
# Hand mutation test for the shell/Python rule checks (no cargo-mutants for these; plan v2 step 1).
# Each row flips one comparison, threshold, guard or exit code in a scratch copy of a check and runs that
# check's selftest against the mutant (via DRIFT_CHECK / PR_DONE_CHECK / PRE_COMMIT_HOOK / PRE_PUSH_HOOK /
# RULE_CHECKS_DAILY / SCOPED_GATE).
# KILLED = the selftest fails on the mutant. Any SURVIVED row fails this script (rc 1).
# Usage: tools/rule_checks_mutation.sh [-j N]      Scratch: ${ZF_SCRATCH:-$HOME/.cache/zenith-float} (never /tmp)
set -u
cd "$(dirname "$0")/.." || exit 2
jobs=4; [ "${1:-}" = "-j" ] && jobs=$2
base="${ZF_SCRATCH:-$HOME/.cache/zenith-float}"; mkdir -p "$base" || exit 2   # a fresh ZF_SCRATCH may not exist yet (Cassie N2)
work=$(mktemp -d "$base/mutation.XXXXXX") || exit 2
trap 'rm -rf "$work"' EXIT
# id | file | occurrence | old ==> new        (a literal \n in old or new is a newline)
# Equivalent, so not a row: Cassie r8 M11, OWNER_NAME += '$'. fullmatch after strip already anchors the end.
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
N01|tools/drift_check.py|1|except (FileNotFoundError, NotADirectoryError): ==> except OSError:
N02|tools/drift_check.py|1|unreadable.append(c) ==> pass
N03|tools/drift_check.py|1|    if unreadable or badorigin or nomatch: ==>     if badorigin or nomatch:
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
Y01|tools/drift_check.py|1|except (GitError, OSError, ValueError) as e: ==> except GitError as e:
Y02|tools/drift_check.py|1|if os.path.exists(p): ==> if False:
Y03|tools/drift_check.py|1|last = max(times) ==> last = min(times)
Y04|tools/drift_check.py|1|if head and not detached else None ==> if False else None
Y05|tools/drift_check.py|1|except OSError as e: ==> except ImportError as e:
Y05b|tools/drift_check.py|2|except OSError as e: ==> except ImportError as e:
Y06|tools/drift_check.py|1|if not act: ==> if act:
Y07|tools/pr_done_check.py|1|except (GhError, KeyError, ValueError, TypeError) as e: ==> except (GhError, KeyError, ValueError) as e:
Y08|tools/pr_done_check.py|1|if since and when and when < since: ==> if since and when < since:
Y10|tools/drift_check.py|1|if window < 0: ==> if window <= 0:
Y12|tools/drift_check.py|1|int(when) if when.isdigit() else None ==> None
Y16|tools/drift_check.py|1|for c in [d] + sorted( ==> for c in sorted(
Y19|tools/pr_done_check.py|1|or pr.get("closedAt") or "")[:10] ==> or pr.get("closedAt") or "")[:9]
Z01|tools/drift_check.py|1|unreadable.append(c)  # EACCES ==> pass  # EACCES
Z02|tools/drift_check.py|1|if not math.isfinite(a.max_age_min) or a.max_age_min < 0: ==> if a.max_age_min < 0:
Z03|tools/drift_check.py|1|if not math.isfinite(a.max_age_min) or a.max_age_min < 0: ==> if not math.isfinite(a.max_age_min):
Z04|tools/drift_check.py|1|or a.max_age_min < 0: ==> or a.max_age_min <= 0:
Z05|tools/drift_check.py|1|if not math.isfinite(a.active_window_min): ==> if False:
Z06|tools/drift_check.py|1|if oldest is not None and now - oldest > max_age: ==> if oldest is not None and now - oldest >= max_age:
Z07|tools/drift_check.py|1|if oldest is not None and now - oldest > max_age: ==> if False:
Z08|tools/pr_done_check.py|1|if not isinstance(pr, dict): ==> if False:
Z11|tools/rule_checks_daily.sh|1|--all "${args[@]}" ==> "${args[@]}"
Z15|tools/rule_checks_daily.sh|1|[ -d "$m" ] && { nd=$((nd+1)); roots+=("$m"); } ==> { nd=$((nd+1)); roots+=("$m"); }
Z16|tools/rule_checks_daily.sh|1|[ -d "$p" ] || { ==> true || {
Z17|tools/rule_checks_daily.sh|1|[ ${#roots[@]} -gt 0 ] || { ==> true || {
Z19|tools/rule_checks_daily.sh|1|nf=$(grep -Ev "$sumre" <<<"$out" | grep -c '^drift_check: FAIL ') ==> nf=$(grep -c '^drift_check: FAIL ' <<<"$out")
Z20|tools/rule_checks_daily.sh|1|1) v=FAIL;; ==> 1) v=PASS;;
Z21|tools/rule_checks_daily.sh|1|exit $rc ==> exit 0
Z22|tools/rule_checks_daily.sh|1|*[\*\?\[]*) ==> *[\?\[]*)
Z30|tools/scoped_gate.sh|1|[ "$k" = 1 ] || ==> [ "$k" -ge 1 ] ||
Z31|tools/scoped_gate.sh|1|failed+=("${names[$i]}") ==> :
Z32|tools/scoped_gate.sh|1|exit 1; fi ==> exit 0; fi
Z33|tools/scoped_gate.sh|1|[ -x "$poprun" ] || ==> true ||
Z34|tools/scoped_gate.sh|1|"rules-scan|@full_gate" ==> "rules-scan|echo skipped"
Z35|tools/scoped_gate.sh|1|cd "$repo" || ==> cd "$here" ||
Z36|tools/scoped_gate.sh|1|--step "${names[$i]}" ==> --step main
Z37|tools/scoped_gate.sh|1|[ -n "$c" ] || ==> true ||
W01|tools/drift_check.py|1|os.lstat(os.path.join(c, ".git")) ==> os.listdir(os.path.join(c, ".git"))
W02|tools/drift_check.py|1|except (FileNotFoundError, NotADirectoryError): ==> except FileNotFoundError:
W03|tools/drift_check.py|1|except (FileNotFoundError, NotADirectoryError): ==> except NotADirectoryError:
W04|tools/drift_check.py|1|if oldest is not None and now - oldest > max_age: ==> if oldest is not None and now - oldest > max_age + 60:
W05|tools/drift_check.py|1|oldest, newest = dirty_span(repo) ==> newest, oldest = dirty_span(repo)
W06|tools/drift_check.py|1|if not math.isfinite(a.active_window_min): ==> if math.isinf(a.active_window_min):
W07|tools/pr_done_check.py|1|if not isinstance(pr, dict): ==> if isinstance(pr, (int, float, str)):
W11|tools/rule_checks_daily.sh|1|unknown option $1" >&2; exit 2;; ==> unknown option $1" >&2; shift;;
W12|tools/rule_checks_daily.sh|1|[ -f "$dc" ] || ==> true ||
W13|tools/rule_checks_daily.sh|1|*) v=ERROR; rc=2;; ==> *) v=ERROR; rc=0;;
W14|tools/rule_checks_daily.sh|1|compgen -G "$p" | sort) ==> compgen -G "$p" | sort | head -1)
W15|tools/rule_checks_daily.sh|1|exit 2; }; roots+=("$p");; ==> exit 2; };;
W17|tools/rule_checks_daily.sh|1|case $rc in 0) v=PASS;; ==> case $rc in 0) v=PASS; rc=1;;
W20|tools/scoped_gate.sh|1|[ $# = 1 ] || usage ==> true || usage
W21|tools/scoped_gate.sh|1|unknown scope '$scope'" >&2; usage;; ==> unknown scope '$scope'" >&2; exit 0;;
W22|tools/scoped_gate.sh|1|bash -o pipefail -c "${cmds[$i]}" ==> bash -c "${cmds[$i]}"
W23|tools/scoped_gate.sh|1|if [ $rc = 0 ]; then echo ==> if true; then echo
W24|tools/scoped_gate.sh|1| "gate-tool-checks|tools/test_gate_tool_checks.sh" ==>  
W25|tools/scoped_gate.sh|1|agent="${POPLOG_AGENT:-spence}" ==> agent=spence
W26|tools/scoped_gate.sh|1|repo="${REPO:-$here/..}" ==> repo="$here/.."
W27|tools/scoped_gate.sh|1|[ "$k" = 1 ] || ==> [ "$k" != 0 ] ||
W28|tools/scoped_gate.sh|1|"full-gate-poplog|tools/test_full_gate_poplog.sh" ==> "full-gate-poplog|true"
W29|tools/scoped_gate.sh|1|"full-gate-syntax|bash -n tools/full_gate.sh" ==> "full-gate-syntax|true"
V15|tools/rule_checks_daily.sh|1|[ $nd -gt 0 ] || ==> true ||
V16|tools/rule_checks_daily.sh|1|[ -d "$m" ] && { nd=$((nd+1)); roots+=("$m"); } ==> nd=$((nd+1)); [ -d "$m" ] && roots+=("$m")
R4X10|tools/rule_checks_daily.sh|1|[ -d "$m" ] && { nd= ==> [ -e "$m" ] && { nd=
R4X11|tools/rule_checks_daily.sh|1|matches no directory" >&2; exit 2; } ==> matches no directory" >&2; exit 1; }
R4X12|tools/rule_checks_daily.sh|1|[ $nd -gt 0 ] || ==> [ $nd -ge 0 ] ||
R4X15|tools/rule_checks_daily.sh|1|                nd=0\n ==>                 nd=1\n
R4X17|tools/rule_checks_daily.sh|1|matches no directory" >&2; exit 2; } ==> matches no directory"; exit 2; }
R4X19|tools/rule_checks_daily.sh|1|*[\*\?\[]*) mapfile ==> *[\*\?]*) mapfile
R4X20|tools/rule_checks_daily.sh|1|roots+=("$m"); }; done ==> roots+=("$m"); }; done; nd=${#ms[@]}
R4X21|tools/rule_checks_daily.sh|1|matches no directory" >&2; exit 2; };; ==> matches no directory" >&2; };;
R4D17|tools/rule_checks_daily.sh|1|  k=$(readlink -f -- "$r"); ==>   k=$r;
R4D18|tools/rule_checks_daily.sh|1|[ -n "${seen[$k]+x}" ] && continue ==> false && continue
R4D19|tools/rule_checks_daily.sh|1|  seen[$k]=1; uniq+=("$r") ==>   uniq+=("$r")
R4D20|tools/rule_checks_daily.sh|1|roots=("${uniq[@]+"${uniq[@]}"}") ==> :
R5Q21|tools/rule_checks_daily.sh|1|  k=$(readlink -f -- "$r"); ==>   k=$(realpath -s -- "$r");
R5Q25|tools/rule_checks_daily.sh|1|  k=$(readlink -f -- "$r"); ==>   k=${r%/};
CR01|tools/drift_check.py|1|    if host not in HOSTS: ==>     if not any(h in u.lower() for h in HOSTS):
CR02|tools/drift_check.py|1|p = p[:-4] if p.lower().endswith(".git") else p ==> p = p.replace(".git", "")
CR03|tools/drift_check.py|1|    u = url.strip() ==>     u = url
CR04|tools/drift_check.py|1|        return ""  # a local path ==>         return u.lower()  # a local path
CR05|tools/drift_check.py|1|    return p.lower()\n ==>     return "/".join(p.split("/")[:2]).lower()\n
CR06|tools/drift_check.py|1|            if not is_work_tree(c):  # before the origin check: a broken .git is this ERROR (Cassie r7 CR06)\n                unreadable.append(c)\n                continue\n            trees += 1\n            if only is not None:\n                try:\n                    r = origin_repo(c)\n                except GitError as e:\n                    badorigin.append((c, str(e)))\n                    continue\n ==>             trees += 1\n            if only is not None:\n                try:\n                    r = origin_repo(c)\n                except GitError as e:\n                    badorigin.append((c, str(e)))\n                    continue\n            if not is_work_tree(c):\n                unreadable.append(c)\n                continue\n
CR07|tools/drift_check.py|1|                if r != only: ==>                 if r and r != only:
CR08|tools/drift_check.py|1|report("skip %s (origin is not GitHub %s)" % (c, only)) ==> pass
CR09|tools/drift_check.py|1|p = p[:-4] if p.lower().endswith(".git") else p ==> p = p[:-4] if p.lower().endswith("git") else p
CR10|tools/drift_check.py|1|for c in [d] + sorted( ==> for c in sorted(
CR11|tools/drift_check.py|1|        raise GitError("origin url cannot be read") ==>         return "jscarr64/zenith-float"
CR12|tools/drift_check.py|1|HOSTS = ("github.com", "www.github.com") ==> HOSTS = ("github.com",)
R701|tools/drift_check.py|1|only = a.only_origin.strip().lower() ==> only = a.only_origin
R702|tools/drift_check.py|1|only = a.only_origin.strip().lower() ==> only = a.only_origin.lower()
R703|tools/drift_check.py|1|only = a.only_origin.strip().lower() ==> only = a.only_origin.strip()
R704|tools/drift_check.py|1|if not OWNER_NAME.fullmatch(only) or only.endswith(".git"): ==> if not OWNER_NAME.fullmatch(only):
R705|tools/drift_check.py|1|if not OWNER_NAME.fullmatch(only) or only.endswith(".git"): ==> if only.endswith(".git"):
R706|tools/drift_check.py|1|OWNER_NAME = re.compile(r"[a-z0-9-]+/[a-z0-9._-]+") ==> OWNER_NAME = re.compile(r"[a-z0-9-]*/[a-z0-9._-]*")
R707|tools/drift_check.py|1|OWNER_NAME = re.compile(r"[a-z0-9-]+/[a-z0-9._-]+") ==> OWNER_NAME = re.compile(r"[a-z0-9-]+/[a-z0-9._ -]+")
R708|tools/drift_check.py|1|    if a.only_origin is not None: ==>     if a.only_origin:
R709|tools/drift_check.py|1|nomatch = only is not None and trees and not matched ==> nomatch = False
R710|tools/drift_check.py|1|nomatch = only is not None and trees and not matched ==> nomatch = only is not None and not matched
R711|tools/drift_check.py|1|            matched += 1 ==>             pass
R712|tools/drift_check.py|1|            trees += 1 ==>             pass
R713|tools/drift_check.py|1|    if unreadable or badorigin or nomatch: ==>     if unreadable or nomatch:
R714|tools/drift_check.py|1|    if unreadable or badorigin or nomatch: ==>     if unreadable or badorigin:
R715|tools/drift_check.py|1|    if n != 1: ==>     if n > 2:
R716|tools/drift_check.py|1|    if n != 1: ==>     if n < 1:
R717|tools/drift_check.py|1|u = git(repo, "ls-remote", "--get-url", "origin") ==> u = git(repo, "config", "--get", "remote.origin.url")
R718|tools/drift_check.py|1|if "origin" not in git(repo, "remote").split(): ==> if False:
R719|tools/drift_check.py|1|        if scheme == "file": ==>         if False:
R720|tools/drift_check.py|1|if scheme not in SCHEMES or not host: ==> if not host:
R721|tools/drift_check.py|1|if scheme not in SCHEMES or not host: ==> if scheme not in SCHEMES:
R722|tools/drift_check.py|1|    elif ":" not in u.split("/", 1)[0]: ==>     elif ":" not in u:
R723|tools/drift_check.py|1|        if not m:\n            return None ==>         if not m:\n            return ""
R724|tools/drift_check.py|1|    if host not in HOSTS: ==>     if False:
R725|tools/drift_check.py|1|p = path[1:] if path.startswith("/") else path ==> p = path
R726|tools/drift_check.py|1|p = p[:-1] if p.endswith("/") else p ==> p = p
R727|tools/drift_check.py|1|p = p[:-4] if p.lower().endswith(".git") else p ==> p = p[:-4] if p.endswith(".git") else p
R728|tools/drift_check.py|1|URL = re.compile(r"([A-Za-z][A-Za-z0-9+.-]*)://(?:[^@/]*@)? ==> URL = re.compile(r"([A-Za-z][A-Za-z0-9+.-]*)://
R729|tools/drift_check.py|1|(?::[0-9]*)?(/.*)?") ==> (/.*)?")
R730|tools/drift_check.py|1|SCP = re.compile(r"(?:[^@/:]*@)?([^/:@]+):(.*)") ==> SCP = re.compile(r"([^/:@]+):(.*)")
R731|tools/drift_check.py|1|        raise GitError("origin url cannot be parsed") ==>         raise GitError("origin url cannot be parsed: " + u)
R732|tools/drift_check.py|1|        raise GitError("origin has %d url values" % n) ==>         return ""
R733|tools/drift_check.py|1|                    report("skip %s (origin is not GitHub %s)" % (c, only))\n                    continue ==>                     report("skip %s (origin is not GitHub %s)" % (c, only))
R734|tools/drift_check.py|1|        print("drift_check: ERROR %s: %s" % (c, why), file=sys.stderr) ==>         pass
M01|tools/drift_check.py|1|if not OWNER_NAME.fullmatch(only) or only.endswith(".git"): ==> if not OWNER_NAME.match(only) or only.endswith(".git"):
M02|tools/drift_check.py|1|if not OWNER_NAME.fullmatch(only) or only.endswith(".git"): ==> if not OWNER_NAME.search(only) or only.endswith(".git"):
M03|tools/drift_check.py|1|m = URL.fullmatch(u) ==> m = URL.match(u)
M04|tools/drift_check.py|1|m = URL.fullmatch(u) ==> m = URL.search(u)
M05|tools/drift_check.py|1|        m = SCP.fullmatch(u) ==>         m = SCP.search(u)
M06|tools/drift_check.py|1|if "origin" not in git(repo, "remote").split(): ==> if "origin" not in git(repo, "remote"):
M07|tools/drift_check.py|1|    if host not in HOSTS: ==>     if host.rstrip(".") not in HOSTS:
M08|tools/drift_check.py|1|    p = path[1:] if path.startswith("/") else path ==>     p = (path.split("?")[0].split("#")[0])[1:] if path.startswith("/") else path.split("?")[0].split("#")[0]
M09|tools/drift_check.py|1|scheme, host, path = m.group(1).lower(), m.group(2).lower(), m.group(3) or "" ==> scheme, host, path = m.group(1), m.group(2).lower(), m.group(3) or ""
M10|tools/drift_check.py|1|print("drift_check: ERROR --only-origin %s matched none of %d work tree(s) under --scan" % (only, trees), ==> pass  #
M12|tools/drift_check.py|1|u = git(repo, "ls-remote", "--get-url", "origin") ==> u = git(repo, "ls-remote", "--get-url", "origin").split("#")[0]
T
}
one() {  # id file occ old new
  local id=$1 f=$2 occ=$3 old=$4 new=$5 d="$work/$1"
  mkdir -p "$d/scratch"; local m="$d/$(basename "$f")"
  python3 -B - "$f" "$m" "$occ" "$old" "$new" <<'PY' || { echo "$id ERROR pattern not found: $old"; return; }
import sys; src, dst, occ, old, new = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5]
old, new = old.replace("\\n", "\n"), new.replace("\\n", "\n")  # a row's \n is a newline (Cassie r4 harness)
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
    tools/rule_checks_daily.sh) st=tools/rule_checks_daily_selftest.sh; env="RULE_CHECKS_DAILY=$m";;
    tools/scoped_gate.sh) st=tools/scoped_gate_selftest.sh; env="SCOPED_GATE=$m";;
  esac
  # a mutant that hangs (a walk that never ends) counts as killed after 15 min
  if timeout 900 env "$env" ZF_SCRATCH="$d/scratch" bash "$st" >"$d/log" 2>&1; then
    echo "$id SURVIVED $f: '$old' -> '$new'"
  else
    echo "$id KILLED   $f: '$old' -> '$new' ($(grep -m1 '^FAIL' "$d/log" | cut -c1-90))"
  fi
}
export -f one; export work
out=$(rows | while IFS='|' read -r id f occ rest; do old=${rest%% ==> *}; new=${rest#* ==> }; printf '%s\0%s\0%s\0%s\0%s\0' "$id" "$f" "$occ" "$old" "$new"; done \
  | xargs -0 -n5 -P "$jobs" bash -c 'one "$@"' _ | sort)
echo "$out"
total=$(echo "$out" | grep -c .); killed=$(echo "$out" | grep -c '^[^ ]* KILLED '); surv=$(echo "$out" | grep -c '^[^ ]* SURVIVED '); err=$(echo "$out" | grep -c '^[^ ]* ERROR ')
echo "rule_checks_mutation: $killed/$total killed, $surv survived, $err errors"
[ "$surv" = 0 ] && [ "$err" = 0 ]
