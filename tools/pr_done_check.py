#!/usr/bin/env python3
"""pr_done_check: a PR is done only when merged, closed and its branch deleted (Jeff 2026-11-25 2:53 PM CT:
"Make sure PR's are closed once merged"; plan v2 C-PRDONE).

Scaffold: the command line only (every option is declared, none does anything yet). Each check lands as a red
selftest commit, then the commit that makes it pass.
"""
import argparse
import sys


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--repo", action="append", default=[], help="OWNER/NAME (repeatable)")
    ap.add_argument("--since", default=None)
    ap.add_argument("--limit", type=int, default=1000)
    ap.parse_args(argv)
    return 0


if __name__ == "__main__":
    sys.exit(main())
