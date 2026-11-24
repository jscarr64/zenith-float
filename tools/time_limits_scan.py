#!/usr/bin/env python3
"""List every wall-clock budget literal in the tests: file, kind, value, count (sorted, tab-separated).

Pinned in tools/time_limits.tsv and checked by tools/test_full_gate_debug_suite_timeouts.sh. Jeff
2026-11-23 23:04 CT: timing limits are enforced in the release step only; that moves where a limit
applies, never its value, so this list must not change with it.
"""
import collections
import os
import re
import sys

ROOT = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# The helper itself (its unit tests use literal durations as test data, not budgets).
SKIP = {"zenith-float-num/src/test_time.rs"}
DIRS = ("zenith-float-num/src", "zenith-float-num/tests", "tests")
PATTERNS = [
    ("secs", re.compile(r"Duration::from_secs\(\s*([\d_]+)\s*\)")),
    ("millis", re.compile(r"Duration::from_millis\(\s*([\d_]+)\s*\)")),
    ("within_ms", re.compile(r"within_ms\(\s*[^,()]+,\s*([\d_]+)\s*\)")),
    ("ratio", re.compile(r"saturating_mul\(\s*([\d_]+)\s*\)\s*<\s*(?:dense|td)\b")),
    ("ode_max_ms", re.compile(r"cos_closed_form_bits\(\s*[^,]+,\s*[^,]+,\s*[^,]+,\s*([\d_]+)\s*,")),
    ("debug_budget_x", re.compile(r"if cfg!\(debug_assertions\) \{ (\d+) \} else \{ 1 \}")),
]


def main() -> int:
    counts = collections.Counter()
    for d in DIRS:
        for dirpath, _, files in os.walk(os.path.join(ROOT, d)):
            for f in files:
                if not f.endswith(".rs"):
                    continue
                path = os.path.join(dirpath, f)
                rel = os.path.relpath(path, ROOT)
                if rel in SKIP:
                    continue
                text = open(path, encoding="utf-8").read()
                for kind, rx in PATTERNS:
                    for m in rx.finditer(text):
                        counts[(rel, kind, m.group(1).replace("_", ""))] += 1
    for (rel, kind, value), n in sorted(counts.items()):
        print(f"{rel}\t{kind}\t{value}\t{n}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
