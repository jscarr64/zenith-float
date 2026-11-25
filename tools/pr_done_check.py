#!/usr/bin/env python3
"""pr_done_check: a PR is done only when merged, closed and its branch deleted (Jeff 2026-11-25 2:53 PM CT:
"Make sure PR's are closed once merged"; plan v2 C-PRDONE).

FAIL (rc 1):
  merged-open         an OPEN PR whose head is already on its base branch (compare base...head is identical or
                      behind), i.e. it landed but the PR was left open;
  branch-not-deleted  a landed PR (MERGED, or CLOSED with its head on the base) whose head branch still exists.
Cross-repository PRs are skipped (their branch is not ours to delete). --since YYYY-MM-DD limits the branch check
to PRs merged or closed on or after that date. The date compared is the UTC date of mergedAt / closedAt (GitHub
gives them in UTC), so a PR merged 2026-11-24 9 PM CST (03:00Z on 11-25) counts as 2026-11-25.
--limit N (default 1000) is the most PRs read per repo; a list that reaches N may be truncated, so it fails (rc 2).
Read-only: only `gh pr list` and `gh api` GET calls.
rc 0 ok, 1 not done, 2 usage or gh error of any kind, a list entry that is not an object included (fails closed). GH=<path> selects the gh binary.
"""
import argparse
import datetime
import json
import os
import re
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


def branch_exists(repo, name):
    return gh("api", "repos/%s/git/ref/heads/%s" % (repo, name), allow_404=True) is not None


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
        if not isinstance(pr, dict):
            raise GhError("gh pr list -R %s: entry %r is not a JSON object" % (repo, pr))
        if pr.get("isCrossRepository"):
            continue
        n, state, base = pr["number"], pr["state"], pr["baseRefName"]
        if state == "OPEN":
            if head_on_base(repo, base, pr["headRefOid"]):
                report("FAIL %s#%d merged-open: head %s is on %s but the PR is still open (%s)"
                       % (repo, n, pr["headRefOid"][:9], base, pr["url"])); bad = 1
            continue
        when = (pr.get("mergedAt") or pr.get("closedAt") or "")[:10]
        if since and when and when < since:
            continue
        landed = state == "MERGED" or (state == "CLOSED" and head_on_base(repo, base, pr["headRefOid"]))
        if landed and branch_exists(repo, pr["headRefName"]):
            report("FAIL %s#%d branch-not-deleted: %s landed (%s) but branch %s still exists (%s)"
                   % (repo, n, state.lower(), when or "?", pr["headRefName"], pr["url"])); bad = 1
    return bad, len(prs)


def valid_date(s):
    if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", s):
        return False
    try:
        datetime.date.fromisoformat(s)
    except ValueError:
        return False
    return True


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--repo", action="append", default=[], help="OWNER/NAME (repeatable)")
    ap.add_argument("--since", default=None, help="YYYY-MM-DD (UTC): skip branch checks for PRs landed before it")
    ap.add_argument("--limit", type=int, default=1000, help="most PRs read per repo; reaching it fails (rc 2)")
    a = ap.parse_args(argv)
    if not a.repo:
        print("pr_done_check: give --repo OWNER/NAME", file=sys.stderr)
        return 2
    if a.since is not None and not valid_date(a.since):
        print("pr_done_check: --since must be a real date YYYY-MM-DD (UTC), got %r" % a.since, file=sys.stderr)
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
