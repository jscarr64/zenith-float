#!/usr/bin/env python3
"""List every place a test branches on debug_assertions to run less (file:line: code), one per line.

Jeff 2026-11-23 23:04 CT: "timing limits are enforced in the release step only; debug still runs every
test for correctness". So no test may skip rows, return early, shrink a range, or be compiled out in a
debug build. Any `debug_assertions` in code (cfg!, #[cfg], #[cfg_attr]) is listed, except:
  - zenith-float-num/src/test_time.rs, the one shared timing switch (LIMITS_ENFORCED) and its tests;
  - hang-watchdog budgets, whose debug values are pinned in tools/time_limits.tsv:
      `const BUDGET: std::time::Duration = if cfg!(debug_assertions) {`
      `.saturating_mul(if cfg!(debug_assertions) { N } else { 1 });`
Also listed: a raw per-row wait (`.recv_timeout(` / `.wait_timeout(`) outside test_time.rs. A row's wall-clock
limit fails it only in the release step; in a debug build it is a hang watchdog only (Spence R1 on 82e4afb). Waits go
through test_time::wait_row (release: the limit; debug: limit x DEBUG_WATCHDOG_FACTOR). Exception: s9, whose row and
table budgets are chosen through test_time::LIMITS_ENFORCED before its wait (WAIT_ALLOWED).
Comments are ignored. Exit 1 if anything is listed. usage: debug_skips_scan.py [ROOT]
"""
import os
import re
import subprocess
import sys

ROOT = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SWITCH = "zenith-float-num/src/test_time.rs"
WATCHDOG = [
    re.compile(r"^\s*const BUDGET: std::time::Duration = if cfg!\(debug_assertions\) \{$"),
    re.compile(r"\.saturating_mul\(if cfg!\(debug_assertions\) \{ \d+ \} else \{ 1 \}\);$"),
]


RAW_WAIT = re.compile(r"\.(recv_timeout|wait_timeout|wait_timeout_while)\(")
WAIT_ALLOWED = {
    ("zenith-float-num/tests/s9_specials_golds.rs", "match rx.recv_timeout(row_budget.min(left)) {"),
}


def rust_files():
    out = subprocess.run(
        ["git", "-C", ROOT, "ls-files", "-co", "--exclude-standard", "--", "*.rs"],
        check=True, capture_output=True, text=True,
    ).stdout
    return sorted(p for p in out.splitlines() if p != SWITCH and os.path.isfile(os.path.join(ROOT, p)))


def main() -> int:
    found = []
    for rel in rust_files():
        with open(os.path.join(ROOT, rel), encoding="utf-8") as fh:
            for n, line in enumerate(fh, 1):
                code = line.split("//", 1)[0].rstrip()
                if RAW_WAIT.search(code) and (rel, code.strip()) not in WAIT_ALLOWED:
                    found.append(f"{rel}:{n}: raw per-row wait (use test_time::wait_row): {code.strip()}")
                if "debug_assertions" not in code:
                    continue
                if any(rx.search(code) for rx in WATCHDOG):
                    continue
                found.append(f"{rel}:{n}: {code.strip()}")
    for f in found:
        print(f)
    return 1 if found else 0


if __name__ == "__main__":
    sys.exit(main())
