# Rule checks: drift, unpushed cap, pre-push, PR done

Plan v2 step 1 (Jeff 2026-11-25: "push and merge as we go", "make sure PRs are closed once merged",
"limit the check to active work trees"). Piled-up local work and forgotten PRs fail loudly instead of growing.

| Check | File | FAILS when | rc |
|---|---|---|---|
| Drift | `tools/drift_check.py` | a local branch or detached HEAD is **more than 1** commit ahead of what is pushed; the oldest unpushed commit is **older than 30 min**; the oldest dirty path (modified, deleted, untracked, not ignored, submodules included) is **older than 30 min**. Exactly 30 min passes. | 0 ok, 1 drift, 2 usage or git error |
| Unpushed cap + library rules | `.githooks/pre-commit` -> `drift_check.py --pre-commit`, then #40's `rules_scan.sh --working`, `rules_scan.sh --count-only --no-rise HEAD`, `cargo fmt --check`, `cargo clippy -D warnings` | the current branch (or detached HEAD) already has **1** unpushed commit: push it first (an unborn branch has 0; a refusal or error stops the hook before the rules run); or any of #40's four rule steps fails (all four always run and report). | 0 / 1 / 2 (not in a work tree, git error) |
| Pre-push | `.githooks/pre-push` | for each pushed ref, in a clean `--shared` clone checked out at the **pushed sha** (tags peeled): `tools/rules_scan.sh "${RATCHET_BASE:-c74e834^}..<sha>"` fails or is missing; `tools/no_suppress_scan.py` fails or is missing (required since #40 landed; it was reported and skipped while #40 was open); the clone cannot be made. No tests run, so a red-first commit pushes. | 0 / 1 |
| PR done | `tools/pr_done_check.py` | an **OPEN** PR whose head is already on its base (merged-open); a **landed** PR (MERGED, or CLOSED with head on base) whose head branch still exists (branch-not-deleted). Cross-repo PRs skipped. | 0 ok, 1 not done, 2 any gh error, or a list that reaches `--limit` |

**"Ahead of what is pushed."** Commits on the branch that are not on its upstream, when the upstream is an existing
remote-tracking ref (`refs/remotes/*`). With no upstream, an upstream that is gone, or an upstream that is another
**local** branch, it is the commits on no remote-tracking ref.

**Active mode (`--scan`).** Only active work trees are checked. A tree is active when:
- its newest change is at most **24 h** old (`--active-window-min 1440`; exactly 24 h is active). The changes
  counted are the HEAD commit time, the index mtime, the `logs/HEAD` mtime and the newest dirty file; or
- it has dirty work **under 30 min** old, whatever the window; or
- it has dirty work **over 30 min** old (that is the violation itself, so it is never skipped as idle), whatever the
  window (Cassie R2-B2d); or
- it is listed with `--include DIR` or in `--include-file FILE` (one path per line, `#` comments); or
- its activity cannot be read (fail closed: the tree is checked, and a git error then gives rc 2).

In an active scanned tree the HEAD branch is checked, plus other local branches whose tip is at most 24 h old.
Idle trees and branches print `skip`. `--all` is the full sweep; `--repo` and `--include` trees are always checked
in full. `--scan DIR` checks DIR itself and its direct children. A child that has a `.git` but cannot be read
(missing gitdir, corrupt HEAD), or whose `.git` cannot even be looked up (no permission, e.g. mode 000), is named and
gives rc 2; a child without `.git` (a plain dir, a bare repo) is skipped. An unreadable DIR gives rc 2. Linked worktrees share branches, so each branch
of a repository is reported once.

Options: `drift_check.py --repo DIR | --scan DIR [--all] [--include DIR] [--include-file FILE]
[--active-window-min 1440] [--max-ahead 1] [--max-age-min 30] [--only-origin OWNER/NAME] [--pre-commit] [--now EPOCH]`.
`--only-origin OWNER/NAME` is the allowlist (Jeff 2026-11-27 4:42 PM CT: check only our own zenith-float folders).
It fails closed (Cassie r7 B1-B3):
- The value is trimmed and lower-cased and must be exactly OWNER/NAME (`[a-z0-9-]+/[a-z0-9._-]+`, no `.git`), else rc 2.
- A `--scan` entry (DIR or a direct child) is checked only when origin is GitHub OWNER/NAME. The URL is the one git
  fetches from (`git ls-remote --get-url origin`, so `insteadOf` counts): `ssh://`, `git://`, `http://`, `https://`
  (any user, token or port) or scp form `[user@]host:path`; host `github.com` or `www.github.com`; any case; an
  optional trailing `/` and `.git`. The URL is never printed (it may hold a token).
- Skipped, as `skip PATH (origin is not GitHub OWNER/NAME)`: no origin remote, a local path (git's rule: no `:`
  before the first `/`) or `file://`, another host (look-alikes included), or a GitHub path that is not exactly
  OWNER/NAME (extra or empty parts, near names, zenith-bound).
- ERROR, rc 2, never a skip: an origin that cannot be read or parsed (an unknown scheme, an empty host, an origin
  remote with no url) or that has two url values; and work trees found under `--scan` with none matching
  (`ERROR --only-origin X matched none of N work tree(s) under --scan`). An empty agent folder is not a work tree.
- The work-tree check runs first, so a broken `.git` is `ERROR PATH has a .git but is not a readable work tree`.
  That message also covers a child that is a symlink loop (rc 2); a `--scan` DIR or daily root that is a loop is
  `not a directory` (rc 2).
- `--repo` and `--include` are never filtered (the pre-commit hook uses `--repo`).
`pr_done_check.py --repo OWNER/NAME [--since YYYY-MM-DD] [--limit 1000]` (GH=<gh binary>; GET calls only).
`--max-age-min` and `--active-window-min` must be finite and >= 0 (`nan`, `inf` or a negative value gives rc 2).
`--include-file` entries are paths as written, resolved against the current directory: use absolute paths.
`pr_done_check.py` fails closed (rc 2) on any malformed gh output, including a list entry that is not an object.
`--since` must be a real `YYYY-MM-DD` date. It is compared with the **UTC** date of mergedAt / closedAt, so a PR
merged 2026-11-24 9 PM CST (03:00Z on 11-25) counts as 2026-11-25.

Known limits (Cassie N4, N5):
- git gives a pre-commit hook no way to tell an amend apart, so `git commit --amend` of the one unpushed commit is
  refused too. Push it first; then the amend is a new commit.
- Dirty-file age is the file mtime, so a `chmod` or `touch` of an old file counts as a change now. This errs
  toward failing.

All git calls use `--no-optional-locks` and drop `GIT_DIR`/`GIT_WORK_TREE`/`GIT_INDEX_FILE`/`GIT_COMMON_DIR`,
so the drift check only reads other agents' worktrees. Scratch goes to `${ZF_SCRATCH:-$HOME/.cache/zenith-float}`,
never `/tmp`. Enable the hooks once per clone with `git config core.hooksPath .githooks`. Fix a refusal; never
bypass it with `--no-verify`.

## Why the pre-push hook was rewritten (2:11 PM CT refusal of d41d938)

Reproduced in `tools/pre_push_selftest.sh` against #40's hook:
1. git exports `GIT_DIR=<gitdir>` to hooks. From a linked worktree it is an absolute path, so #40's
   `git -C <scratch clone> checkout` ran against the pusher's repo. The scratch clone stayed empty, and the hook
   reported "could not check out ... or it has no tools/no_suppress_scan.py". The new hook unsets the git env first.
2. #40's hook ran `rules_scan.sh` in the pusher's worktree. Its count checks read working-tree files, so dirty
   files decided the verdict for a clean commit. The new hook runs every scanner in the clean clone at the pushed sha.

#40 also added `.githooks/pre-commit` (rules_scan, fmt, clippy). The #40 merge (Jeff 2026-11-27 10:42 PM CT) keeps
both: the unpushed cap runs first, then #40's four steps. `drift_check_selftest.sh` stubs `tools/rules_scan.sh` and
`cargo` in its fixture (it tests the cap); `tools/test_rules_enforcement.sh` part 4 tests the rules part and now pushes
each accepted fixture commit (the cap allows one unpushed commit). The mutation runner runs both selftests for a
pre-commit mutant. `tools/pre_push_selftest.sh` has two parts: part 1 (#42, wiring with stub scanners) and part 2
(#40, the real scanners on fixture commits of a clone of this repo's HEAD).

## Tests

`tools/drift_check_selftest.sh`, `tools/pr_done_check_selftest.sh` (mock gh that honours `--state` and `--limit`),
`tools/pre_push_selftest.sh` (part 1: real pushes to a local bare origin, stub scanners; part 2: the real scanners), `tools/rule_checks_daily_selftest.sh`,
`tools/scoped_gate_selftest.sh` (stub poprun), and `tools/rule_checks_mutation.sh` (243 mutants, including the #40-merge rows C02 (adapted), P01-P07 and H14, Cassie's review rows, her r7 CR01-CR12, the R701-R734 allowlist rows and her r8 M01-M10 and M12; every one must be killed; a row's `\n` is a newline, and a
mutant that hangs is killed after 15 min). Every selftest run is capped at 60 s, so a hang fails its case fast.
The rows for the retired `--exclude` mechanism (the step-by-step symlink walk, the git status pathspec and the
tracked-file skip) were removed with it; D25, N03, W05 and X20 were adapted to the code that is left. The r7 rewrite
of the origin check replaced L01-L14 with R701-R734, and Cassie's CR rows are adapted to it with the same intent
(CR01 host substring, CR02/CR09 the .git strip, CR04 the local-path return, CR05 extra parts, CR06 the work-tree
check order, CR07 the skip test, CR11 a read error, CR12 www.). Cassie's r8 M11 (`OWNER_NAME` += `$`) is equivalent,
not a row: `fullmatch` after `strip()` already anchors the end. Each behavior
landed as one red commit, then one fix commit. Test-only commits add cases for code that was already right; their
red proof is the mutant each one kills, all in the sweep. full_gate.sh runs all six as step `rule-checks`.

## Follow-ups (Cassie r8, non-blocking; left alone)

- `ssh://github.com:org/repo` and `@` inside userinfo are unparseable ERROR, not skip (Cassie r8 N1).
- A daily sweep of a Team whose every `*-work` folder has no work trees is rc 0 with 0 repos checked (Cassie r8,
  empty-Team). A typo'd Team with no `*-work` dirs is already rc 2 at the glob.

## Scoped gate (in the repo)

`tools/scoped_gate.sh rules` is the scoped gate for a PR that changes only these tools (Jeff Q6: a scoped gate per PR,
the full gate at tier boundaries). Steps: `full-gate-syntax` (`bash -n tools/full_gate.sh`), `rules-scan` and
`rule-checks` (their commands are read verbatim from `tools/full_gate.sh`, so the two cannot drift),
`gate-tool-checks` and `full-gate-poplog`. Each step is one `poprun --step`; all steps run; rc 1 names the failed
ones. `--list` / `--dry-run` print the plan. Run it from a clean checkout of the PR head:

    PYTHONDONTWRITEBYTECODE=1 ~/pop-logs/bin/poprun --agent spence --job scoped-gate-rules -- bash tools/scoped_gate.sh rules

The poprun footer must say `RESULT: PASS` and `gate-commit: <PR head> (clean at start and end)`. It is not
full_gate evidence: it covers 2 of full_gate's 25 run steps (no Rust builds or tests).

## Daily full sweep (Jeff 2026-11-25 5:20 PM CT)

"Yes, add a daily full sweep over every worktree that fails on any drift." Active mode hides idle drift (Cassie N1),
so once a day `tools/rule_checks_daily.sh` runs `drift_check.py --all` (every tree and every branch, nothing skipped)
over every root and FAILs (rc 1) on any drift; rc 2 on no roots, a bad root or a git error. Arguments after `--` go
to drift_check (`run-daily.sh` passes `--only-origin jscarr64/zenith-float`).
- Roots are dirs or quoted globs of dirs. Non-dir matches are ignored, but a glob that matches no directory at all is
  an error (rc 2), so a mistyped root, or a Team with no agent folders, cannot silently pass.
- The same dir given twice (a literal and a glob, `./`, a symlink) is swept once, first spelling kept. A root nested
  inside another root is still listed as its own root; drift_check checks each repo once, so findings are not
  doubled.
- The summary line counts the FAIL (finding) and ERROR lines and the roots.

## 15-minute scheduler on Pop

A systemd **user** timer, the same pattern as `poprotate.timer`. No cron, no system units. Both timers are disabled
until Jeff re-enables them.
- Scope (Jeff 2026-11-27 4:42 PM CT): only our own zenith-float clones and worktrees. The roots are
  `/home/jsc64/MathWorks/Team/<agent>-work` (the symlink resolved with `readlink -f` at run time) and their direct
  children; a work tree counts only when its `remote.origin.url` is GitHub `jscarr64/zenith-float`
  (`--only-origin`). Other repos, zenith-bound clones, local-path origins and anything deeper are not checked.
- `~/.config/systemd/user/zf-rule-checks.timer`: `OnCalendar=*:0/15`, `Persistent=true`.
- `~/.config/systemd/user/zf-rule-checks.service`: oneshot,
  `ExecStart=/home/jsc64/MathWorks/Team/spence-work/rule-checks-bin/run.sh`.
- `run.sh` runs copies of `drift_check.py` and `pr_done_check.py` pinned to one commit, kept in
  `/home/jsc64/MathWorks/Team/spence-work/rule-checks-bin/<sha>/` with `SHA256SUMS`. It refuses to run if a hash
  does not match. It runs `drift_check.py --scan` (active mode) over every `Team/*-work` with
  `--only-origin jscarr64/zenith-float` and `--include-file rule-checks-bin/include.txt`, then `pr_done_check.py` on
  `jscarr64/zenith-float`. Output is appended to `/home/jsc64/pop-logs/spence/rule-checks/YYYY-MM-DD.log`.
- Daily: `zf-rule-checks-daily.timer` (`OnCalendar=*-*-* 06:00`, `Persistent=true`) starts
  `zf-rule-checks-daily.service` (oneshot, `ExecStart=/home/jsc64/MathWorks/Team/spence-work/rule-checks-bin/run-daily.sh`).
  It checks the same pinned `SHA256SUMS` and runs the pinned `rule_checks_daily.sh` over `Team/*-work` with the same
  `--only-origin`. Output: `/home/jsc64/pop-logs/spence/rule-checks/YYYY-MM-DD-daily.log`.
- To repin after a merge, copy `drift_check.py`, `pr_done_check.py` and `rule_checks_daily.sh` at the new sha into
  `rule-checks-bin/<sha>/`, write `SHA256SUMS`, and change `PIN` in `run.sh` and `run-daily.sh`. To stop them:
  `systemctl --user disable --now zf-rule-checks.timer zf-rule-checks-daily.timer`.
- Linger is on for jsc64 (Jeff 2026-11-25), so both timers run without a session on Pop.
