#!/usr/bin/env python3
"""drift_check: no work piles up unpushed (Jeff 2026-11-25: "We should be able to push and merge as we go").

Plan v2 C-DRIFT / C-AHEAD. FAIL (rc 1) when, in a checked work tree,
  * a local branch, or a detached HEAD, is more than MAX_AHEAD (default 1) commits ahead of what is pushed, or
  * its OLDEST unpushed commit (committer time) is older than MAX_AGE (default 30 min; exactly 30 min passes), or
  * its OLDEST dirty path (modified, deleted, untracked and not ignored, submodules included) is older than MAX_AGE
    (file mtime, so a chmod or touch counts as a change; a deleted file uses the index mtime).
"Ahead of what is pushed" = commits on the branch that are not on its upstream, when the upstream is a
remote-tracking ref (refs/remotes/*) that exists; otherwise (no upstream, an upstream that is gone, or an upstream
that is another LOCAL branch) commits that are on no remote-tracking ref (`--not --remotes`). A detached HEAD is
checked the same way unless it sits on a local branch tip (then that branch covers it).
Read-only: every git call uses --no-optional-locks and drops GIT_DIR/GIT_WORK_TREE/GIT_INDEX_FILE/GIT_COMMON_DIR.
rc 0 ok, 1 drift, 2 usage or git error (fail closed).
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


def _ref_exists(repo, ref):
    try:
        git(repo, "rev-parse", "--verify", "--quiet", ref + "^{commit}")
        return True
    except GitError:
        return False


def unpushed(repo, ref, upstream):
    """(count, oldest committer time or None) of commits on ref that are not pushed."""
    if upstream.startswith("refs/remotes/") and _ref_exists(repo, upstream):
        spec = ["%s..%s" % (upstream, ref)]
    else:
        spec = [ref, "--not", "--remotes"]
    out = git(repo, "log", "--format=%ct", *spec).split()
    if not out:
        return 0, None
    return len(out), min(int(t) for t in out)


def branches(repo):
    rows = []
    for line in git(repo, "for-each-ref", "--format=%(refname)%00%(upstream)%00%(objectname)", "refs/heads").splitlines():
        if line:
            name, up, sha = line.split("\0")
            rows.append((name, up, sha))
    return rows


def dirty_span(repo):
    """(oldest, newest) mtime among dirty paths, (None, None) when clean."""
    out = git(repo, "status", "--porcelain=v1", "-z", "--untracked-files=all", "--ignore-submodules=none")
    entries = out.split("\0")
    oldest = newest = None
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
        newest = m if newest is None else max(newest, m)
    return oldest, newest


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


def check_repo(repo, now, max_ahead, max_age, report):
    bad = 0
    detached, head = head_state(repo)
    tips = set()
    for name, up, sha in branches(repo):
        tips.add(sha)
        n, oldest = unpushed(repo, name, up)
        short = name[len("refs/heads/"):]
        if n > max_ahead:
            report("FAIL %s %s: %d commits ahead of what is pushed (max %d)" % (repo, short, n, max_ahead)); bad = 1
        elif n and now - oldest > max_age:
            report("FAIL %s %s: unpushed commit %d s old (max %d s)" % (repo, short, now - oldest, max_age)); bad = 1
    if detached and head not in tips:
        n, oldest = unpushed(repo, "HEAD", "")
        if n > max_ahead:
            report("FAIL %s detached HEAD: %d commits not on any remote ref (max %d)" % (repo, n, max_ahead)); bad = 1
        elif n and now - oldest > max_age:
            report("FAIL %s detached HEAD: unpushed commit %d s old (max %d s)" % (repo, now - oldest, max_age)); bad = 1
    d, _ = dirty_span(repo)
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
    ap.add_argument("--all", action="store_true", help="check every scanned work tree and branch (full sweep)")
    ap.add_argument("--include", action="append", default=[], help="a work tree always checked (repeatable)")
    ap.add_argument("--include-file", action="append", default=[], help="file of work trees always checked")
    ap.add_argument("--active-window-min", type=float, default=1440)
    a = ap.parse_args(argv)
    now = time.time() if a.now is None else a.now
    max_age = a.max_age_min * 60
    report = lambda s: print("drift_check: " + s)
    repos = []
    for r in a.repo:
        if not is_work_tree(r):
            print("drift_check: %s is not the top of a git work tree" % r, file=sys.stderr)
            return 2
        repos.append(r)
    if not repos:
        print("drift_check: give --repo", file=sys.stderr)
        return 2
    bad = 0
    try:
        for r in repos:
            bad |= check_repo(r, now, a.max_ahead, max_age, report)
    except GitError as e:
        print("drift_check: %s" % e, file=sys.stderr)
        return 2
    print("drift_check: %s (%d repo(s), max %d ahead, max age %g min)"
          % ("FAIL" if bad else "ok", len(repos), a.max_ahead, a.max_age_min))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
