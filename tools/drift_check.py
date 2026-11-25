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
--pre-commit: refuse (rc 1) a new commit on the current branch or detached HEAD when MAX_AHEAD commits are already
unpushed (the new one would exceed the cap); ages are not checked there. An unborn branch has 0. Used by
.githooks/pre-commit. git gives a pre-commit hook no way to tell an amend apart, so `git commit --amend` of the one
unpushed commit is refused too: push it, then amend is a new commit.
--scan DIR: DIR itself and each direct child of DIR that has a .git (file or directory). One whose .git cannot be
read as the top of a work tree (missing gitdir, corrupt HEAD, ...) is named and gives rc 2 (fail closed); a
directory without .git (a plain dir, a bare repo) is skipped; a child whose .git cannot even be looked up (no
permission) is unreadable too, and an unreadable DIR gives rc 2. Linked worktrees share branches: each branch of a
repository is reported once per run.
Active mode (Jeff 2026-11-25 3:55 PM CT: "limit the check to active work trees"): a --scan tree is checked only when
  * its newest change is at most WINDOW old (--active-window-min, default 1440 = 24 h; exactly 24 h is active):
    HEAD commit time, index mtime, logs/HEAD mtime, newest dirty file mtime; or
  * it has dirty work less than MAX_AGE (30 min) old, whatever the window; or
  * it is listed with --include PATH or in --include-file FILE (one path per line, # comments); or
  * its activity cannot be read (fail closed: it is checked, and a git error then gives rc 2).
  In an active scanned tree the HEAD branch is checked, plus other local branches whose tip commit is at most WINDOW
  old (a tip with no commit date is checked). Idle trees and branches are named as "skip".
  --all checks every scanned tree and every branch (the full sweep). --repo and --include trees are always checked
  in full.
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
    fmt = "--format=%(refname)%00%(upstream)%00%(objectname)%00%(committerdate:unix)"
    for line in git(repo, "for-each-ref", fmt, "refs/heads").splitlines():
        if line:
            name, up, sha, when = line.split("\0")
            rows.append((name, up, sha, int(when) if when.isdigit() else None))
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


def activity(repo, now, window, max_age):
    """(active?, reason) for a --scan tree. Anything unreadable counts as active (fail closed)."""
    try:
        gd = git(repo, "rev-parse", "--absolute-git-dir").strip()
        times = [int(git(repo, "log", "-1", "--format=%ct", "HEAD").strip())]
        for f in ("index", os.path.join("logs", "HEAD")):
            p = os.path.join(gd, f)
            if os.path.exists(p):
                times.append(os.stat(p).st_mtime)
        _, newest = dirty_span(repo)
    except (GitError, OSError, ValueError) as e:
        return True, "activity unknown, checked (%s)" % str(e).splitlines()[0][:120]
    if newest is not None and now - newest < max_age:
        return True, "dirty work under %d s old" % max_age
    if newest is not None:
        times.append(newest)
    last = max(times)
    if now - last <= window:
        return True, "last change %d s ago" % (now - last)
    return False, "idle: last change %d s ago (window %d s)" % (now - last, window)


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


def check_repo(repo, now, max_ahead, max_age, report, pre_commit=False, seen=None, window=None):
    """seen: set of (common git dir, branch) already checked; linked worktrees share branches, report each once.
    window: None checks every branch; else only the HEAD branch and branches whose tip is at most window old."""
    bad = 0
    seen = set() if seen is None else seen
    common = os.path.realpath(os.path.join(repo, git(repo, "rev-parse", "--git-common-dir").strip()))
    detached, head = head_state(repo)
    if pre_commit:
        if detached:
            n, _ = unpushed(repo, "HEAD", "")
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
    head_ref = git(repo, "symbolic-ref", "-q", "HEAD").strip() if head and not detached else None
    for name, up, sha, when in branches(repo):
        tips.add(sha)
        if window is not None and name != head_ref and when is not None and now - when > window:
            report("skip %s %s: branch tip %d s old (window %d s)" % (repo, name[len("refs/heads/"):], now - when, window))
            continue
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
    window = a.active_window_min * 60
    if window < 0:
        print("drift_check: --active-window-min must be >= 0", file=sys.stderr)
        return 2
    report = lambda s: print("drift_check: " + s)
    explicit = list(a.repo) + list(a.include)
    for f in a.include_file:
        try:
            explicit += [l.strip() for l in open(f) if l.strip() and not l.strip().startswith("#")]
        except OSError as e:
            print("drift_check: --include-file %s: %s" % (f, e), file=sys.stderr)
            return 2
    repos = []
    for r in explicit:
        if not is_work_tree(r):
            print("drift_check: %s is not the top of a git work tree" % r, file=sys.stderr)
            return 2
        repos.append((r, None))
    unreadable = []
    for d in a.scan:
        if not os.path.isdir(d):
            print("drift_check: %s is not a directory" % d, file=sys.stderr)
            return 2
        try:
            names = os.listdir(d)
        except OSError as e:
            print("drift_check: --scan %s cannot be read: %s" % (d, e), file=sys.stderr)
            return 2
        for c in [d] + sorted(os.path.join(d, c) for c in names):
            try:
                os.lstat(os.path.join(c, ".git"))
            except (FileNotFoundError, NotADirectoryError):
                continue  # no .git: a plain dir, a file or a bare repo
            except OSError:
                unreadable.append(c)  # EACCES and the like: cannot tell, fail closed
                continue
            if is_work_tree(c):
                repos.append((c, "scan"))
            else:
                unreadable.append(c)
    if not explicit and not a.scan:
        print("drift_check: give --repo or --scan", file=sys.stderr)
        return 2
    bad = 0
    seen = set()
    for c in unreadable:
        print("drift_check: ERROR %s has a .git but is not a readable work tree" % c, file=sys.stderr)
    done = set()
    checked = skipped = 0
    try:
        for r, how in repos:
            key = os.path.realpath(r)
            if key in done:
                continue
            done.add(key)
            win = None
            if how == "scan" and not a.all:
                act, why = activity(r, now, window, max_age)
                if not act:
                    report("skip %s (%s)" % (r, why)); skipped += 1
                    continue
                win = window
            checked += 1
            bad |= check_repo(r, now, a.max_ahead, max_age, report, a.pre_commit, seen, win)
    except GitError as e:
        print("drift_check: %s" % e, file=sys.stderr)
        return 2
    print("drift_check: %s (%d repo(s) checked, %d idle skipped, max %d ahead, max age %g min, %s)"
          % ("FAIL" if bad else "ok", checked, skipped, a.max_ahead, a.max_age_min,
             "all work trees" if a.all else "active window %g min" % a.active_window_min))
    if unreadable:
        return 2
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
