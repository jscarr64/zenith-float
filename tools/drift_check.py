#!/usr/bin/env python3
"""drift_check: no work piles up unpushed (Jeff 2026-11-25: "We should be able to push and merge as we go").

Plan v2 C-DRIFT / C-AHEAD. FAIL (rc 1) when, in a checked repo,
  * a local branch, or a detached HEAD, is more than MAX_AHEAD (default 1) commits ahead of what is pushed, or
  * an unpushed commit is older than MAX_AGE (default 30 min; exactly 30 min passes, 30 min + 1 s fails), or
  * a dirty file (modified, deleted, untracked and not ignored) is older than MAX_AGE (file mtime; a deleted
    file uses the index mtime).
"Ahead" = commits on the branch that are not on its upstream; with no upstream, or an upstream ref that no longer
exists, commits not on any remote-tracking ref (`--not --remotes`). A detached HEAD is checked the same way unless it
sits on a local branch tip.

--pre-commit: refuse (rc 1) a new commit on the current branch or detached HEAD when MAX_AHEAD commits are already
unpushed (the new one would exceed the cap); ages are not checked there. Used by .githooks/pre-commit.
--scan DIR: check every git work tree that is DIR itself or a direct child of DIR (bare repos are skipped).
Read-only: every git call uses --no-optional-locks; nothing is written. rc 0 ok, 1 drift, 2 usage or git error.
"""
import argparse
import os
import subprocess
import sys
import time


class GitError(Exception):
    pass


def git(repo, *args):
    p = subprocess.run(["git", "--no-optional-locks", "-C", repo] + list(args),
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=_env())
    if p.returncode != 0:
        raise GitError("git %s: %s" % (" ".join(args), p.stderr.decode(errors="replace").strip()))
    return p.stdout.decode(errors="replace")


def _env():
    e = dict(os.environ)
    for k in ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE", "GIT_COMMON_DIR"):
        e.pop(k, None)  # hooks export these; -C <repo> must win
    return e


def is_work_tree(path):
    try:
        return git(path, "rev-parse", "--is-inside-work-tree").strip() == "true" and \
            os.path.realpath(git(path, "rev-parse", "--show-toplevel").strip()) == os.path.realpath(path)
    except GitError:
        return False


def unpushed(repo, ref, upstream):
    """(count, oldest committer time or None) of commits on ref not pushed."""
    if upstream and _ref_exists(repo, upstream):
        spec = ["%s..%s" % (upstream, ref)]
    else:
        spec = [ref, "--not", "--remotes"]
    out = git(repo, "log", "--format=%ct", *spec).split()
    if not out:
        return 0, None
    return len(out), min(int(t) for t in out)


def _ref_exists(repo, ref):
    try:
        git(repo, "rev-parse", "--verify", "--quiet", ref + "^{commit}")
        return True
    except GitError:
        return False


def branches(repo):
    rows = []
    for line in git(repo, "for-each-ref", "--format=%(refname)%00%(upstream)%00%(objectname)", "refs/heads").splitlines():
        if line:
            name, up, sha = line.split("\0")
            rows.append((name, up, sha))
    return rows


def head_state(repo):
    """(detached?, head sha or None)."""
    try:
        sha = git(repo, "rev-parse", "--verify", "--quiet", "HEAD").strip()
    except GitError:
        return False, None  # unborn branch
    try:
        git(repo, "symbolic-ref", "-q", "HEAD")
        return False, sha
    except GitError:
        return True, sha


def dirty_oldest(repo):
    """oldest mtime among dirty paths, or None when clean."""
    out = git(repo, "status", "--porcelain=v1", "-z", "--untracked-files=all", "--ignore-submodules=none")
    entries = out.split("\0")
    oldest = None
    i = 0
    index_mtime = None
    while i < len(entries):
        e = entries[i]
        i += 1
        if not e:
            continue
        xy, path = e[:2], e[3:]
        if xy[0] in "RC":
            i += 1  # the original path of a rename/copy follows
        full = os.path.join(repo, path)
        try:
            m = os.lstat(full).st_mtime
        except FileNotFoundError:
            if index_mtime is None:
                ix = git(repo, "rev-parse", "--git-path", "index").strip()
                index_mtime = os.stat(ix if os.path.isabs(ix) else os.path.join(repo, ix)).st_mtime
            m = index_mtime
        oldest = m if oldest is None else min(oldest, m)
    return oldest


def check_repo(repo, now, max_ahead, max_age, pre_commit, report, seen=None):
    """seen: set of (common git dir, branch) already checked; linked worktrees share branches, report each once."""
    bad = 0
    seen = set() if seen is None else seen
    common = os.path.realpath(os.path.join(repo, git(repo, "rev-parse", "--git-common-dir").strip()))
    detached, head = head_state(repo)
    if pre_commit:
        if detached:
            n, _ = unpushed(repo, "HEAD", None)
            what = "detached HEAD"
        else:
            ref = git(repo, "symbolic-ref", "-q", "HEAD").strip()
            up = git(repo, "for-each-ref", "--format=%(upstream)", ref).strip() if head else ""
            n, _ = unpushed(repo, ref, up) if head else (0, None)
            what = ref[len("refs/heads/"):]
        if n >= max_ahead:
            report("FAIL %s %s: %d unpushed commit(s) already; push before committing again (cap %d)"
                   % (repo, what, n, max_ahead))
            return 1
        return 0
    tips = set()
    for name, up, sha in branches(repo):
        tips.add(sha)
        if (common, name) in seen:
            continue
        seen.add((common, name))
        n, oldest = unpushed(repo, name, up)
        short = name[len("refs/heads/"):]
        if n > max_ahead:
            report("FAIL %s %s: %d commits ahead of what is pushed (max %d)" % (repo, short, n, max_ahead)); bad = 1
        elif n and now - oldest > max_age:
            report("FAIL %s %s: unpushed commit %d s old (max %d s)" % (repo, short, now - oldest, max_age)); bad = 1
    if detached and head not in tips:
        n, oldest = unpushed(repo, "HEAD", None)
        if n > max_ahead:
            report("FAIL %s detached HEAD: %d commits not on any remote ref (max %d)" % (repo, n, max_ahead)); bad = 1
        elif n and now - oldest > max_age:
            report("FAIL %s detached HEAD: unpushed commit %d s old (max %d s)" % (repo, now - oldest, max_age)); bad = 1
    d = dirty_oldest(repo)
    if d is not None and now - d > max_age:
        report("FAIL %s: dirty work %d s old (max %d s)" % (repo, int(now - d), max_age)); bad = 1
    if not bad:
        report("ok   %s" % repo)
    return bad


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--repo", action="append", default=[], help="a work tree to check (repeatable)")
    ap.add_argument("--scan", action="append", default=[], help="check DIR and its direct children (repeatable)")
    ap.add_argument("--max-ahead", type=int, default=1)
    ap.add_argument("--max-age-min", type=float, default=30)
    ap.add_argument("--now", type=float, default=None, help="epoch seconds (tests)")
    ap.add_argument("--pre-commit", action="store_true")
    a = ap.parse_args(argv)
    now = time.time() if a.now is None else a.now
    max_age = a.max_age_min * 60
    repos = []
    for r in a.repo:
        if not is_work_tree(r):
            print("drift_check: %s is not the top of a git work tree" % r, file=sys.stderr)
            return 2
        repos.append(r)
    for d in a.scan:
        if not os.path.isdir(d):
            print("drift_check: %s is not a directory" % d, file=sys.stderr)
            return 2
        cands = [d] + sorted(os.path.join(d, c) for c in os.listdir(d))
        repos += [c for c in cands if os.path.isdir(c) and is_work_tree(c)]
    if not a.repo and not a.scan:
        print("drift_check: give --repo or --scan", file=sys.stderr)
        return 2
    bad = 0
    seen = set()
    try:
        for r in repos:
            bad |= check_repo(r, now, a.max_ahead, max_age, a.pre_commit, lambda s: print("drift_check: " + s), seen)
    except GitError as e:
        print("drift_check: %s" % e, file=sys.stderr)
        return 2
    print("drift_check: %s (%d repo(s), max %d ahead, max age %d min)"
          % ("FAIL" if bad else "ok", len(repos), a.max_ahead, a.max_age_min))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
