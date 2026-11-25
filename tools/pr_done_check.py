#!/usr/bin/env python3
"""pr_done_check: a PR is done only when merged, closed and its branch deleted (Jeff 2026-11-25 2:53 PM CT:
"Make sure PR's are closed once merged"; plan v2 C-PRDONE).

FAIL (rc 1):
  merged-open         an OPEN PR whose head is already on its base branch (compare base...head is identical or
                      behind), i.e. it landed but the PR was left open;
  branch-not-deleted  a landed PR (MERGED, or CLOSED with its head on the base) whose head branch still exists.
Cross-repository PRs are skipped (their branch is not ours to delete). --since YYYY-MM-DD limits the branch check
to PRs merged or closed on or after that date. Read-only: only `gh pr list` and `gh api` GET calls.
rc 0 ok, 1 not done, 2 usage or gh error (fails closed). GH=<path> selects the gh binary (default: gh on PATH).
"""
import argparse
import json
import os
import subprocess
import sys

LANDED = ("identical", "behind")


class GhError(Exception):
    pass


def gh(*args, allow_404=False):
    p = subprocess.run([os.environ.get("GH", "gh")] + list(args), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if p.returncode != 0:
        err = p.stderr.decode(errors="replace")
        if allow_404 and "HTTP 404" in err:
            return None
        raise GhError("gh %s: %s" % (" ".join(args), err.strip()))
    return json.loads(p.stdout.decode() or "null")


def head_on_base(repo, base, sha):
    r = gh("api", "repos/%s/compare/%s...%s" % (repo, base, sha))
    return r["status"] in LANDED


def branch_exists(repo, name):
    return gh("api", "repos/%s/git/ref/heads/%s" % (repo, name), allow_404=True) is not None


def check(repo, since, report):
    prs = gh("pr", "list", "-R", repo, "--state", "all", "--limit", "1000", "--json",
             "number,state,headRefName,headRefOid,baseRefName,mergedAt,closedAt,isCrossRepository,url")
    bad = 0
    for pr in prs:
        if pr.get("isCrossRepository"):
            continue
        n, state, base = pr["number"], pr["state"], pr["baseRefName"]
        if state == "OPEN":
            if head_on_base(repo, base, pr["headRefOid"]):
                report("FAIL %s#%d merged-open: head %s is on %s but the PR is still open (%s)"
                       % (repo, n, pr["headRefOid"][:9], base, pr["url"])); bad = 1
            continue
        when = (pr.get("mergedAt") or pr.get("closedAt") or "")[:10]
        if since and when and when[:10] < since:
            continue
        landed = state == "MERGED" or (state == "CLOSED" and head_on_base(repo, base, pr["headRefOid"]))
        if landed and branch_exists(repo, pr["headRefName"]):
            report("FAIL %s#%d branch-not-deleted: %s landed (%s) but branch %s still exists (%s)"
                   % (repo, n, state.lower(), when or "?", pr["headRefName"], pr["url"])); bad = 1
    return bad, len(prs)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--repo", action="append", default=[], help="OWNER/NAME (repeatable)")
    ap.add_argument("--since", default=None, help="YYYY-MM-DD: skip branch checks for PRs landed before this date")
    a = ap.parse_args(argv)
    if not a.repo:
        print("pr_done_check: give --repo OWNER/NAME", file=sys.stderr)
        return 2
    bad = 0
    try:
        for r in a.repo:
            b, n = check(r, a.since, lambda s: print("pr_done_check: " + s))
            print("pr_done_check: %s %s (%d PRs)" % ("FAIL" if b else "ok  ", r, n))
            bad |= b
    except (GhError, KeyError, ValueError) as e:
        print("pr_done_check: error, failing closed: %s" % e, file=sys.stderr)
        return 2
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
