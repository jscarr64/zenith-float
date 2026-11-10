#!/usr/bin/env python3
"""Tightness report and ratchet (TEST_PLAN section 9 'Metric', checklist S12). Integers only, no floats; run with /usr/bin/python3.

  tools/tightness_report.py report RESULTS                 print the per-op, per-precision table (min, median, p99, max of rad / R_ref)
  tools/tightness_report.py check  RESULTS [BASELINE]      ratchet: exit 1 if any (op, prec) maximum is WORSE than the baseline, or a cell is
                                                           UNRECORDED (not in the baseline) or MISSING (in the baseline, not in the results)
  tools/tightness_report.py update RESULTS [BASELINE]      lower baseline rows and add new cells; refuses (exit 1) to raise a row or drop one

RESULTS (written by the tightness sets through tests/common/tightness.rs): '#' lines, then one row per cell, sorted by (OP, PREC):
    OP PREC COUNT MIN MEDIAN P99 MAX        ratio = rad / R_ref in units of 1/65536 (q16, rounded up); MIN <= MEDIAN <= P99 <= MAX
BASELINE (tools/tightness_baseline.txt): '#' lines, then rows  OP PREC MAX_Q16  (sorted). Exit codes: 0 ok, 1 ratchet failure, 2 malformed input.

An absent RESULTS file is the empty result set: it passes against an empty baseline and fails (MISSING) against a recorded one, so a
baseline row cannot be silenced by deleting the sets that produce it. The ratchet only fails when a maximum GROWS (an improvement is a note).
TAU below holds the PLACEHOLDER tight factors of TEST_PLAN section 9 (decision D2 is Rusty's): the report marks a maximum above tau as OVER;
it is informational, the assertion lives in the tightness sets.
"""
import os
import sys

Q16 = 65536
# PROVISIONAL (D2): tau per op from TEST_PLAN section 9; parametric ones (powi, sum, dot, matrices) and the L2 specials have no entry here.
TAU = {"add": 2, "sub": 2, "mul": 4, "div": 4, "sqrt": 4, "sqr": 4, "pow": 32}
for _op in ("exp", "exp2", "expm1", "log", "log2", "log10", "log1p"):
    TAU[_op] = 8
for _op in ("sin", "cos", "tan", "sinh", "cosh", "tanh", "asin", "acos", "atan", "asinh", "acosh", "atanh"):
    TAU[_op] = 16

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_BASELINE = os.path.join(HERE, "tightness_baseline.txt")


class Malformed(Exception):
    pass


def num(tok, what, line):
    if not tok.isdigit() or str(int(tok)) != tok:
        raise Malformed("line %d: %s %r is not a canonical non-negative integer" % (line, what, tok))
    return int(tok)


def data_lines(path):
    with open(path, "r") as f:
        text = f.read()
    if "\r" in text or (text and not text.endswith("\n")):
        raise Malformed("%s: needs LF line endings and a final newline" % path)
    return [(i + 1, l) for i, l in enumerate(text.split("\n")[:-1]) if not l.startswith("#")]


def read_results(path):
    if not os.path.exists(path):
        return {}
    cells, prev = {}, None
    for n, l in data_lines(path):
        f = l.split(" ")
        if len(f) != 7 or any(x == "" for x in f):
            raise Malformed("line %d: want 'OP PREC COUNT MIN MEDIAN P99 MAX', got %r" % (n, l))
        op = f[0]
        prec, count, mn, med, p99, mx = (num(f[i], "field %d" % (i + 1), n) for i in range(1, 7))
        if prec == 0 or count == 0:
            raise Malformed("line %d: prec and count must be positive" % n)
        if not (mn <= med <= p99 <= mx):
            raise Malformed("line %d: need MIN <= MEDIAN <= P99 <= MAX" % n)
        key = (op, prec)
        if key in cells:
            raise Malformed("line %d: duplicate cell %s %d" % (n, op, prec))
        if prev is not None and key < prev:
            raise Malformed("line %d: rows are not sorted by (OP, PREC)" % n)
        prev = key
        cells[key] = (count, mn, med, p99, mx)
    return cells


def read_baseline(path):
    rows, prev = {}, None
    for n, l in data_lines(path):
        f = l.split(" ")
        if len(f) != 3 or any(x == "" for x in f):
            raise Malformed("%s line %d: want 'OP PREC MAX_Q16', got %r" % (path, n, l))
        key = (f[0], num(f[1], "PREC", n))
        if key in rows:
            raise Malformed("%s line %d: duplicate cell" % (path, n))
        if prev is not None and key < prev:
            raise Malformed("%s line %d: rows are not sorted" % (path, n))
        prev = key
        rows[key] = num(f[2], "MAX_Q16", n)
    return rows


def dec(q):
    return "%d.%04d" % divmod((q * 10000 + Q16 - 1) // Q16, 10000)


def report(cells):
    print("%-8s %5s %7s %9s %9s %9s %9s %5s %s" % ("op", "prec", "n", "min", "median", "p99", "max", "tau", "status"))
    for (op, prec), (count, mn, med, p99, mx) in sorted(cells.items()):
        tau = TAU.get(op)
        status = "ok" if tau is None or mx <= tau * Q16 else "OVER"
        print("%-8s %5d %7d %9s %9s %9s %9s %5s %s" % (op, prec, count, dec(mn), dec(med), dec(p99), dec(mx), "-" if tau is None else tau, status))
    print("tightness: %d cells (tau = provisional placeholders, decision D2)" % len(cells))


def check(cells, base):
    bad = 0
    for key, (count, mn, med, p99, mx) in sorted(cells.items()):
        if key not in base:
            print("tightness: UNRECORDED cell %s %d (max %d); record it with 'update'" % (key[0], key[1], mx))
            bad += 1
        elif mx > base[key]:
            print("tightness: WORSE %s %d: max %d (%s) > baseline %d (%s)" % (key[0], key[1], mx, dec(mx), base[key], dec(base[key])))
            bad += 1
        elif mx < base[key]:
            print("tightness: improved %s %d: max %d < baseline %d (lower it with 'update')" % (key[0], key[1], mx, base[key]))
    for key in sorted(base):
        if key not in cells:
            print("tightness: MISSING cell %s %d (in the baseline, not in the results)" % key)
            bad += 1
    print("tightness: %d cells, baseline %d rows, %s" % (len(cells), len(base), "FAILED" if bad else "ok"))
    return 1 if bad else 0


def update(cells, base, path):
    new = dict(base)
    for key, (count, mn, med, p99, mx) in cells.items():
        if key in base and mx > base[key]:
            print("tightness: refusing to raise %s %d from %d to %d" % (key[0], key[1], base[key], mx))
            return 1
        new[key] = mx if key not in base else min(mx, base[key])
    for key in base:
        if key not in cells:
            print("tightness: refusing to drop %s %d (it is in the baseline, not in the results)" % key)
            return 1
    with open(path, "r") as f:
        header = [l for l in f.read().split("\n") if l.startswith("#")]
    with open(path, "w") as f:
        for l in header:
            f.write(l + "\n")
        for (op, prec), mx in sorted(new.items()):
            f.write("%s %d %d\n" % (op, prec, mx))
    print("tightness: baseline updated (%d rows)" % len(new))
    return 0


def main(argv):
    if len(argv) < 3 or argv[1] not in ("report", "check", "update") or len(argv) > 4:
        print(__doc__)
        return 2
    try:
        cells = read_results(argv[2])
        if argv[1] == "report":
            report(cells)
            return 0
        bpath = argv[3] if len(argv) == 4 else DEFAULT_BASELINE
        base = read_baseline(bpath)
        return check(cells, base) if argv[1] == "check" else update(cells, base, bpath)
    except (Malformed, OSError) as e:
        print("tightness: malformed or unreadable input: %s" % e)
        return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv))
