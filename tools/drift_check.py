#!/usr/bin/env python3
"""drift_check: no work piles up unpushed (Jeff 2026-11-25: "We should be able to push and merge as we go").

Scaffold: the command line only (every option is declared, none does anything yet). Each check lands as a red selftest commit, then the commit that makes it pass.
"""
import argparse
import sys


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--repo", action="append", default=[], help="a work tree to check (repeatable)")
    ap.add_argument("--scan", action="append", default=[], help="check DIR and its direct children (repeatable)")
    ap.add_argument("--max-ahead", type=int, default=1)
    ap.add_argument("--max-age-min", type=float, default=30)
    ap.add_argument("--now", type=float, default=None, help="epoch seconds (tests)")
    ap.add_argument("--pre-commit", action="store_true")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--include", action="append", default=[])
    ap.add_argument("--include-file", action="append", default=[])
    ap.add_argument("--active-window-min", type=float, default=1440)
    ap.parse_args(argv)
    return 0


if __name__ == "__main__":
    sys.exit(main())
