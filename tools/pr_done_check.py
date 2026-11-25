#!/usr/bin/env python3
"""pr_done_check: a PR is done only when merged, closed and its branch deleted (Jeff 2026-11-25 2:53 PM CT:
"Make sure PR's are closed once merged"; plan v2 C-PRDONE).

FAIL (rc 1):
  merged-open         an OPEN PR whose head is already on its base branch (compare base...head is identical or
                      behind), i.e. it landed but the PR was left open;
  (branch-not-deleted lands next, with its red test first; --since is accepted and not used yet.)
--limit N (default 1000) is the most PRs read per repo; a list that reaches N may be truncated, so it fails (rc 2).
Read-only: only `gh pr list` and `gh api` GET calls.
rc 0 ok, 1 not done, 2 usage or gh error of any kind (fails closed). GH=<path> selects the gh binary.
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
    try:
        p = subprocess.run([os.environ.get("GH", "gh")] + list(args), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except OSError as e:
        raise GhError("cannot run gh: %s" % e)
    if p.returncode != 0:
        err = p.stderr.decode(errors="replace")
        if allow_404 and "HTTP 404" in err:
            return None
        raise GhError("gh %s: %s" % (" ".join(args), err.strip()))
    out = p.stdout.decode(errors="replace").strip()
    if not out:
        raise GhError("gh %s: empty output" % " ".join(args))
    return json.loads(out)


def head_on_base(repo, base, sha):
    r = gh("api", "repos/%s/compare/%s...%s" % (repo, base, sha))
    return r["status"] in LANDED


def check(repo, since, limit, report):
    prs = gh("pr", "list", "-R", repo, "--state", "all", "--limit", str(limit), "--json",
             "number,state,headRefName,headRefOid,baseRefName,mergedAt,closedAt,isCrossRepository,url")
    if not isinstance(prs, list):
        raise GhError("gh pr list -R %s: not a JSON list" % repo)
    if len(prs) >= limit:
        raise GhError("gh pr list -R %s returned %d PRs = --limit %d: possibly truncated; raise --limit"
                      % (repo, len(prs), limit))
    bad = 0
    for pr in prs:
        n, state, base = pr["number"], pr["state"], pr["baseRefName"]
        if state == "OPEN":
            if head_on_base(repo, base, pr["headRefOid"]):
                report("FAIL %s#%d merged-open: head %s is on %s but the PR is still open (%s)"
                       % (repo, n, pr["headRefOid"][:9], base, pr["url"])); bad = 1
            continue
    return bad, len(prs)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--repo", action="append", default=[], help="OWNER/NAME (repeatable)")
    ap.add_argument("--since", default=None, help="YYYY-MM-DD (UTC): skip branch checks for PRs landed before it")
    ap.add_argument("--limit", type=int, default=1000, help="most PRs read per repo; reaching it fails (rc 2)")
    a = ap.parse_args(argv)
    if not a.repo:
        print("pr_done_check: give --repo OWNER/NAME", file=sys.stderr)
        return 2
    bad = 0
    try:
        for r in a.repo:
            b, n = check(r, a.since, a.limit, lambda s: print("pr_done_check: " + s))
            print("pr_done_check: %s %s (%d PRs)" % ("FAIL" if b else "ok  ", r, n))
            bad |= b
    except (GhError, KeyError, ValueError, TypeError) as e:
        print("pr_done_check: error, failing closed: %s" % e, file=sys.stderr)
        return 2
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
