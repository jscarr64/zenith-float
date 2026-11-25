# Rule checks: drift, unpushed cap, pre-push, PR done

Plan v2 step 1 (Jeff 2026-11-25: "push and merge as we go", "make sure PRs are closed once merged").
These checks make piled-up local work and forgotten PRs fail loudly instead of growing silently.

| Check | File | FAILS when | rc |
|---|---|---|---|
| Drift | `tools/drift_check.py` | a local branch is **more than 1** commit ahead of its upstream (no upstream: commits on no remote ref); a detached HEAD has more than 1 commit on no remote ref; the oldest unpushed commit is **older than 30 min**; the oldest dirty file (modified, deleted, untracked, not ignored) is **older than 30 min**. Exactly 30 min passes. | 0 ok, 1 drift, 2 usage |
| Unpushed cap | `.githooks/pre-commit` → `drift_check.py --pre-commit` | the current branch (or detached HEAD) already has **1** unpushed commit: push it before the next commit | 0 / 1 |
| Pre-push | `.githooks/pre-push` | for each pushed ref, in a clean `--shared` clone checked out at the **pushed sha**: `tools/rules_scan.sh` fails or is missing; `tools/no_suppress_scan.py` (when the pushed tree has it) fails; the clone cannot be made. It runs no tests, so a red-first commit pushes. | 0 / 1 |
| PR done | `tools/pr_done_check.py` | an **OPEN** PR whose head is already on its base (`compare` = identical or behind): merged-open; a **landed** PR (MERGED, or CLOSED with head on base) whose head branch still exists: branch-not-deleted. Cross-repo PRs skipped. gh errors fail closed. | 0 ok, 1 not done, 2 error |

Options: `drift_check.py --repo DIR | --scan DIR [--max-ahead 1] [--max-age-min 30] [--pre-commit] [--now EPOCH]`;
`pr_done_check.py --repo OWNER/NAME [--since YYYY-MM-DD]` (GH=<gh binary>; GET calls only).
All git calls use `--no-optional-locks` and drop `GIT_DIR`/`GIT_WORK_TREE`/`GIT_INDEX_FILE`, so the drift check is
read-only in other agents' worktrees. Scratch goes to `${ZF_SCRATCH:-$HOME/.cache/zenith-float}`, never `/tmp`.

Enable the hooks once per clone: `git config core.hooksPath .githooks`. A refusal is fixed, never bypassed with `--no-verify`.

## Why the pre-push hook was rewritten (2:11 PM CT refusal of d41d938)

Reproduced in `tools/pre_push_selftest.sh` against #40's hook:
1. git exports `GIT_DIR=<gitdir>` to hooks. From a linked worktree it is an absolute path, so #40's
   `git -C <scratch clone> checkout` ran against the pusher's repo, the scratch clone stayed empty, and the hook
   reported "could not check out ... or it has no tools/no_suppress_scan.py". The new hook unsets the git env first.
2. #40's hook ran `rules_scan.sh` in the pusher's worktree; its count checks read working-tree files, so dirty
   files decided the verdict for a clean commit. The new hook runs every scanner in the clean clone at the pushed sha.

## Tests

`tools/drift_check_selftest.sh`, `tools/pr_done_check_selftest.sh` (mock gh), `tools/pre_push_selftest.sh` (real pushes
to a local bare origin, stub scanners), and `tools/rule_checks_mutation.sh` (each comparison, threshold, guard and exit
code flipped in a scratch copy; every mutant must be killed). full_gate.sh runs all four as step `rule-checks`.

## Follow-up for Jeff (not installed)

A 15-minute scheduler on Pop (cron or a systemd user timer; neither was touched) running, for example:

    cd /mnt/hub/spence-work/zf-step1 && python3 -B tools/drift_check.py --scan /mnt/hub/paul-work \
      --scan /mnt/hub/rusty-work --scan /mnt/hub/spence-work --scan /mnt/hub/vera-work
    python3 -B tools/pr_done_check.py --repo jscarr64/zenith-float --repo jscarr64/zenith-bound

with a FAIL reported to Jeff. Installing it needs Jeff's go (it changes cron/systemd).
