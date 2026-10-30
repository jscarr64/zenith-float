#!/home/jsc64/.local/oracle-venv/bin/python
"""
Exact golds for slice S6 (src/matrix.rs, CONTRACT 16.S6.6 to 16.S6.10): dot, matvec, matmul, norm_inf, norm_1, norm_frob.
Truth is exact: Python ints and fractions.Fraction only (no floats, no mpmath values); the generator never calls zenith-bound.
    scripts/gen_matrix_goldens.py                  write the gate file tests/goldens/slices/matrix/matrix_rows.txt, SHA256SUMS and
                                                   HEAVY.SHA256 (sha256 of the heavy file, computed while streaming; not committed)
    scripts/gen_matrix_goldens.py --heavy DIR      as above and also write DIR/matrix_rows_heavy.txt
    scripts/gen_matrix_goldens.py --check DIR      regenerate the gate file into DIR (never tests/goldens), byte-compare it with the
                                                   committed one, and compare the sha256 of the regenerated heavy rows with HEAVY.SHA256
    scripts/gen_matrix_goldens.py --selftest       the generator's own checks only
--skip-heavy may follow any mode: the heavy tier is not computed (write mode leaves HEAVY.SHA256 as is). tools/goldens_regen.sh passes it in
gate mode. Rows stream to disk one at a time (one row is at most a few tens of MB at 32x32x32, 16384 bits).

ROWS (one case per line, fields separated by single spaces; a value is 'SIGN EXP MANT': (-1)^SIGN * MANT * 2^(EXP - bitlen(MANT)),
MANT lowercase hex and odd, zero is '0 0 0'; radii use the same form and are nonnegative with at most 8 mantissa bits):
  products   ID OP PREC M K N  A[M*K]  B[K*N]  OUT[M*N]
             every A/B entry is 'MID RAD' (6 tokens), row-major; dot has M = N = 1 (A is the row vector x, B the column vector y),
             matvec has N = 1 (B is the vector); every OUT entry is 'P0 P1 P2 TREF' (12 tokens):
             P0 = sum_k m(a_ik) m(b_kj) exactly; P1 = the exact product of (A low corner) and (B high corner); P2 = (A high corner)(B low corner)
             (low corner = mid - rad entrywise, high = mid + rad); TREF = R1 + 2^(1-PREC) * S rounded up to 64 bits, with
             R1 = sum_k |m(a)| r(b) + r(a)(|m(b)| + r(b)) and S = sum_k |m(a) m(b)| (CONTRACT 16.S6.9). No midpoint is pinned (Spence ruling 2,
             2026-10-25: the implementation accumulates sequentially at PREC + guard bits). A consumer checks containment of P0, P1, P2 and
             rad <= 2 * TREF.
  norms      ID OP PREC M 0 N  A[M*N]  S HI
             norm_inf / norm_1: S = HI = the exact max row / column sum of |m| + r. norm_frob: S = sum (|m| + r)^2 exactly and HI = sqrt(S)
             rounded up to 64 bits. vnorm_inf / vnorm_1 / vnorm_2 (BallVec norms, M = 1, the vector is the row): max / sum of |m| + r, and
             S = sum of squares, HI = sqrt(S) up. A consumer checks result >= S (norm_frob, vnorm_2: result^2 >= S) and
             result <= HI * (1 + 2^-56). PREC is the draw precision of the mids (norms take no prec).
  spread     the row dot.spread (gate tier, CONTRACT 16.S6.10): x = (2^E, 1, -2^E), y = (1, 1, 1), all exact, E = 2^41, so the terms are
             2^(2^41) apart. P0 = P1 = P2 = 1 analytically and TREF = 2^(1-PREC)(2^(E+1) + 1) rounded up is written by exponent bookkeeping
             (spread_out); no Fraction near 2^E is ever built. The selftest checks the same text against the generic Fraction path at E = 200.
Draws: goldlib.SplitMix64 from goldlib.stream_seed(17, 1) (gate) and stream_seed(17, 2) (heavy), in case-table order (suite 17 = s6_matrix,
registered in tests/common/seeds.rs).
"""
# popguard gate (as gen_la_goldens.py): before any heavy import or work. Outside a popslot/poprun memory scope ~/pop-logs/bin/popguard.sh exits
# nonzero and so do we, with its rc; 3 if the guard is missing or cannot run.
import os as _os, subprocess as _subprocess, sys as _sys
def _popguard():
    guard = _os.path.expanduser("~/pop-logs/bin/popguard.sh")
    try:
        rc = _subprocess.run([guard, "gen_matrix_goldens"], check=False).returncode
    except OSError as e:
        print(f"gen_matrix_goldens: popguard missing or cannot run ({guard}: {e}); refusing to start. "
              "Run it via ~/pop-logs/bin/poprun --agent NAME --job JOB -- ~/pop-logs/bin/popslot -- CMD.", file=_sys.stderr)
        _sys.exit(3)
    if rc != 0:
        print(f"gen_matrix_goldens: popguard refused (rc {rc}); not running outside a popslot memory scope.", file=_sys.stderr)
        _sys.exit(rc if rc > 0 else 3)
_popguard()

import hashlib
import os
import sys
from fractions import Fraction
from math import isqrt

import goldlib
from gen_ball_goldens import vfrac, from_fraction, venc

sys.set_int_max_str_digits(0)

GENERATOR_VERSION = 2
GOLDEN_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tests", "goldens", "slices", "matrix")
GATE_NAME = "matrix_rows.txt"
HEAVY_NAME = "matrix_rows_heavy.txt"
SUITE = 17
RAD_BITS = 8
EXP_BAND = 8
SLACK = Fraction(1, 1 << 56)
SPREAD_E = 1 << 41

# Case tables: (op, prec, M, K, N, kind). kind: "rand" draws mids and radii; "zero", "ident", "cancel", "tie", "exact", "spread" are fixed shapes.
GATE_PRECS = (1, 2, 53, 64, 113, 256)
PRODUCT_SHAPES = {"dot": ((1, 1, 1), (1, 3, 1), (1, 8, 1)),
                  "matvec": ((1, 1, 1), (2, 3, 1), (4, 4, 1)),
                  "matmul": ((1, 1, 1), (2, 3, 2), (3, 3, 3), (2, 4, 3))}
NORM_SHAPES = ((1, 0, 1), (2, 0, 3), (4, 0, 4))
MAT_NORMS = ("norm_inf", "norm_1", "norm_frob")
VEC_NORMS = ("vnorm_inf", "vnorm_1", "vnorm_2")
NORMS = MAT_NORMS + VEC_NORMS


def gate_cases():
    cases = []
    for prec in GATE_PRECS:
        for op in ("dot", "matvec", "matmul"):
            for (m, k, n) in PRODUCT_SHAPES[op]:
                cases.append((op, prec, m, k, n, "rand"))
        for op in MAT_NORMS:
            for (m, k, n) in NORM_SHAPES:
                cases.append((op, prec, m, k, n, "rand"))
    # fixed rows at prec 53 (and 2 for the tie): empty inner dimension, zero, identity, exact cancellation, ToEven tie, exact (no radius)
    cases += [("matmul", 53, 2, 0, 2, "rand"), ("dot", 53, 1, 0, 1, "rand"), ("matvec", 53, 3, 0, 1, "rand"),
              ("matmul", 53, 3, 3, 3, "zero"), ("matmul", 53, 3, 3, 2, "ident"), ("dot", 53, 1, 4, 1, "cancel"),
              ("dot", 2, 1, 2, 1, "tie"), ("matmul", 64, 3, 5, 2, "exact"),
              ("matmul", 8, 8, 8, 8, "rand"), ("norm_inf", 53, 8, 0, 8, "rand"), ("norm_1", 53, 8, 0, 8, "rand"),
              ("norm_frob", 53, 8, 0, 8, "rand"), ("norm_inf", 53, 2, 0, 0, "rand"), ("norm_frob", 53, 0, 0, 3, "rand")]
    # version 2 (Spence rulings 2 and 5): the huge-spread row, then BallVec norms (length 4 per prec, and empty vectors)
    cases.append(("dot", 53, 1, 3, 1, "spread"))
    for prec in GATE_PRECS:
        cases += [(op, prec, 1, 0, 4, "rand") for op in VEC_NORMS]
    cases += [(op, 53, 1, 0, 0, "rand") for op in VEC_NORMS]
    return cases


def heavy_cases():
    cases = []
    for prec in (1024, 4096, 16384):
        cases += [("matmul", prec, 16, 16, 16, "rand"), ("matmul", prec, 32, 32, 32, "rand"),
                  ("matvec", prec, 32, 32, 1, "rand"), ("dot", prec, 1, 32, 1, "rand")]
        cases += [(op, prec, 32, 0, 32, "rand") for op in MAT_NORMS]
    for prec in (1024, 4096, 16384):
        cases += [(op, prec, 1, 0, 32, "rand") for op in VEC_NORMS]
    return cases


EXPECTED_ROWS = {"gate": 150, "heavy": 30}


# ---- dyadic helpers -------------------------------------------------------------------------------------------------------
def pow2(e):
    return Fraction(1 << e) if e >= 0 else Fraction(1, 1 << -e)


def binade(v):
    """e with 2^(e-1) <= |v| < 2^e, v != 0."""
    a = abs(v)
    e = a.numerator.bit_length() - a.denominator.bit_length()
    while pow2(e) <= a:
        e += 1
    while pow2(e - 1) > a:
        e -= 1
    return e


def round_even(x, prec):
    return goldlib.round_fraction(x, prec, "ToEven")[0]


def round_up64(x):
    return goldlib.round_fraction(x, 64, "Up")[0] if x != 0 else Fraction(0)


def ulp(v, prec):
    return Fraction(0) if v == 0 else pow2(binade(v) - prec)


def sqrt_up64(s):
    """sqrt(s) rounded up to 64 bits, s a nonnegative dyadic Fraction (exact integer square roots only)."""
    if s == 0:
        return Fraction(0)
    t = 0
    while True:
        scaled = s * pow2(2 * t)
        if scaled >= pow2(140):
            break
        t += 8
    n = scaled.numerator // scaled.denominator
    r = isqrt(n)
    if r * r < scaled:
        r += 1
    up = Fraction(r) / pow2(t)
    if not (up * up >= s):
        raise AssertionError("sqrt_up64 not an upper bound")
    return round_up64(up)


def enc(x):
    return venc(from_fraction(x))


# ---- draws ------------------------------------------------------------------------------------------------------------------
def draw_mid(rng, prec):
    if rng.next() % 16 == 0:
        return Fraction(0)
    bits = 1 + rng.next() % prec
    mant = 1 << (bits - 1)
    got = 1
    while got < bits:
        take = min(64, bits - got)
        mant |= (rng.next() & ((1 << take) - 1)) << (bits - got - take)
        got += take
    e = (rng.next() % (2 * EXP_BAND + 1)) - EXP_BAND
    v = Fraction(mant) * pow2(e - bits)
    return -v if rng.next() % 2 else v


def draw_rad(rng, prec, mid):
    if rng.next() % 4 == 0:
        return Fraction(0)
    mant = 1 + rng.next() % ((1 << RAD_BITS) - 1)
    base = binade(mid) if mid != 0 else 0
    e = base - prec - (rng.next() % 6) + 2
    return Fraction(mant) * pow2(e - RAD_BITS)


def draw_matrix(rng, prec, rows, cols, kind, which):
    out = []
    for i in range(rows):
        row = []
        for j in range(cols):
            if kind == "rand":
                m = draw_mid(rng, prec)
                row.append((m, draw_rad(rng, prec, m)))
            elif kind == "zero":
                row.append((Fraction(0), Fraction(0)))
            elif kind == "ident":
                # A = I (exact), B = (i*cols + j + 1)/64 +- 2^-60: I*B is B itself
                row.append((Fraction(int(i == j)), Fraction(0)) if which == "A" else (Fraction(i * cols + j + 1, 64), pow2(-60)))
            elif kind == "exact":
                row.append((Fraction(draw_mid(rng, prec)), Fraction(0)))
            elif kind == "cancel":
                # x = (1, 2^-60, 3, -1), y = (1, 1, 1, 3): sum = 1 + 2^-60 + 3 - 3 = exact 1 + 2^-60 needs 61 bits; radii on one entry
                vals = (Fraction(1), pow2(-60), Fraction(3), Fraction(-1)) if which == "A" else (Fraction(1), Fraction(1), Fraction(1), Fraction(3))
                row.append((vals[i * cols + j], pow2(-70) if (which == "B" and j == 0 and i == 0) else Fraction(0)))
            elif kind == "tie":
                # x = (5/8, 0), y = (1, 1): sum 5/8 = 0.101b, a tie at 2 bits between 1/2 (mantissa 10b, even) and 3/4; ToEven gives 1/2
                vals = (Fraction(5, 8), Fraction(0)) if which == "A" else (Fraction(1), Fraction(1))
                row.append((vals[i * cols + j], Fraction(0)))
            elif kind == "spread_small":
                # selftest only: the spread row at a small exponent (SPREAD_SMALL_E), through the generic Fraction path
                vals = (pow2(SPREAD_SMALL_E), Fraction(1), -pow2(SPREAD_SMALL_E)) if which == "A" else (Fraction(1),) * 3
                row.append((vals[i * cols + j], Fraction(0)))
            else:
                raise ValueError(kind)
        out.append(row)
    return out


SPREAD_SMALL_E = 200


# ---- truth ------------------------------------------------------------------------------------------------------------------
def product_entry(arow, bcol, prec):
    p0 = sum((a[0] * b[0] for a, b in zip(arow, bcol)), Fraction(0))
    p1 = sum(((a[0] - a[1]) * (b[0] + b[1]) for a, b in zip(arow, bcol)), Fraction(0))
    p2 = sum(((a[0] + a[1]) * (b[0] - b[1]) for a, b in zip(arow, bcol)), Fraction(0))
    first = sum((abs(a[0]) * b[1] + a[1] * (abs(b[0]) + b[1]) for a, b in zip(arow, bcol)), Fraction(0))
    sabs = sum((abs(a[0] * b[0]) for a, b in zip(arow, bcol)), Fraction(0))
    tref = round_up64(first + pow2(1 - prec) * sabs)
    return p0, p1, p2, tref, first, sabs


def spread_mids(e):
    """The spread row's A entries (x = (2^e, 1, -2^e), exact) and B entries (y = (1, 1, 1)) as text, by exponent bookkeeping."""
    one = "0 1 1"
    zero = "0 0 0"
    a = ["0 %d 1 %s" % (e + 1, zero), "%s %s" % (one, zero), "1 %d 1 %s" % (e + 1, zero)]
    b = ["%s %s" % (one, zero)] * 3
    return a, b


def spread_out(e, prec):
    """P0 P1 P2 TREF of the spread row: P = 1 exactly; TREF = round_up64(2^(1-prec) (2^(e+1) + 1)) = (2^63 + 1) 2^(e+2-prec-63)
    for e >= 64 (2^(e+1) + 1 has e + 2 > 64 bits, rounding up at 64 bits gives 2^(e+1) (1 + 2^-63)); 'SIGN EXP MANT' has EXP = e + 3 - prec."""
    if e < 64:
        raise ValueError("spread_out needs e >= 64")
    one = "0 1 1"
    return "%s %s %s 0 %d %x" % (one, one, one, e + 3 - prec, (1 << 63) + 1)


def norm_values(op, a):
    absm = [[abs(m) + r for (m, r) in row] for row in a]
    if op in ("norm_inf", "vnorm_1"):
        v = max((sum(row, Fraction(0)) for row in absm), default=Fraction(0))
        return v, v
    if op in ("norm_1", "vnorm_inf"):
        cols = len(absm[0]) if absm else 0
        v = max((sum((absm[i][j] for i in range(len(absm))), Fraction(0)) for j in range(cols)), default=Fraction(0))
        return v, v
    s = sum((x * x for row in absm for x in row), Fraction(0))
    return s, sqrt_up64(s)


def entries(mat):
    return [enc(m) + " " + enc(r) for row in mat for (m, r) in row]


def make_row(idx, case, rng):
    op, prec, m, k, n, kind = case
    cid = "%s.%s.p%d.%dx%dx%d.%d" % (op, kind, prec, m, k, n, idx)
    if op in NORMS:
        a = draw_matrix(rng, prec, m, n, kind, "A")
        s, hi = norm_values(op, a)
        return " ".join([cid, op, str(prec), str(m), "0", str(n)] + entries(a) + [enc(s), enc(hi)])
    if kind == "spread":
        sa, sb = spread_mids(SPREAD_E)
        return " ".join([cid, op, str(prec), str(m), str(k), str(n)] + sa + sb + [spread_out(SPREAD_E, prec)])
    a = draw_matrix(rng, prec, m, k, kind, "A")
    b = draw_matrix(rng, prec, k, n, kind, "B")
    outs = []
    for i in range(m):
        for j in range(n):
            p0, p1, p2, tref, _, _ = product_entry(a[i], [b[t][j] for t in range(k)], prec)
            outs.append(" ".join(enc(v) for v in (p0, p1, p2, tref)))
    return " ".join([cid, op, str(prec), str(m), str(k), str(n)] + entries(a) + entries(b) + outs)


def row_width(op, m, k, n):
    if op in NORMS:
        return 6 + 6 * m * n + 6
    return 6 + 6 * m * k + 6 * k * n + 12 * m * n


def header(tier, nrows):
    return ["# zenith-bound slice S6 gold rows (%s tier), generated by scripts/gen_matrix_goldens.py (generator version %d)" % (tier, GENERATOR_VERSION),
            "# format: see the generator docstring; products: ID OP PREC M K N A[M*K](MID RAD) B[K*N](MID RAD) OUT[M*N](P0 P1 P2 TREF)",
            "# norms (norm_*, vnorm_*): ID OP PREC M 0 N A[M*N](MID RAD) S HI; values 'SIGN EXP MANT' (MANT hex, odd; zero '0 0 0')",
            "# rows: %d" % nrows]

def stream_tier(tier, sinks):
    """Write the header and every row of a tier to each sink (callables taking one line); returns the row count."""
    cases = gate_cases() if tier == "gate" else heavy_cases()
    want = EXPECTED_ROWS[tier]
    if len(cases) != want:
        raise SystemExit("gen_matrix_goldens: %s tier has %d cases, EXPECTED_ROWS says %d" % (tier, len(cases), want))
    for line in header(tier, want):
        for s in sinks:
            s(line)
    rng = goldlib.SplitMix64(goldlib.stream_seed(SUITE, 1 if tier == "gate" else 2))
    count = 0
    for idx, case in enumerate(cases):
        line = make_row(idx, case, rng)
        op, _, m, k, n, _ = case
        if len(line.split(" ")) != row_width(op, m, k, n):
            raise AssertionError("row %d width" % idx)
        for s in sinks:
            s(line)
        count += 1
    if count != want:
        raise AssertionError("row count")
    return count


def write_tier(tier, path, hasher):
    with open(path, "wb") if path else open(os.devnull, "wb") as f:
        def sink(line):
            data = (line + "\n").encode("ascii")
            f.write(data)
            hasher.update(data)
        return stream_tier(tier, [sink])


def sha_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def pinned_heavy():
    path = os.path.join(GOLDEN_DIR, "HEAVY.SHA256")
    for line in open(path):
        parts = line.split()
        if len(parts) == 2 and parts[1] == HEAVY_NAME:
            return parts[0]
    raise SystemExit("gen_matrix_goldens: no %s line in HEAVY.SHA256" % HEAVY_NAME)


def write_sums(gate_sha, heavy_sha):
    goldlib.write_lines(os.path.join(GOLDEN_DIR, "SHA256SUMS"), [
        "# Slice S6 golden checksums (sha256sum format), checked by tools/goldens_regen.sh. Every *.txt in this directory must be listed.",
        "%s  %s" % (gate_sha, GATE_NAME)])
    if heavy_sha is not None:
        goldlib.write_lines(os.path.join(GOLDEN_DIR, "HEAVY.SHA256"), [
            "# sha256 of the heavy-tier S6 golds (not committed; rebuilt by scripts/gen_matrix_goldens.py --heavy DIR,",
            "# verified by scripts/gen_matrix_goldens.py --check DIR without --skip-heavy). sha256sum format.",
            "%s  %s" % (heavy_sha, HEAVY_NAME)])


# ---- selftest ---------------------------------------------------------------------------------------------------------------
def _check(cond, msg):
    if not cond:
        raise AssertionError("selftest: " + msg)


def selftest():
    rng = goldlib.SplitMix64(goldlib.stream_seed(SUITE, 99))
    for prec in (1, 2, 53, 300):
        for _ in range(40):
            m = draw_mid(rng, prec)
            r = draw_rad(rng, prec, m)
            _check(vfrac(from_fraction(m)) == m and vfrac(from_fraction(r)) == r, "encode round trip")
            if m != 0:
                _check(abs(m).numerator.bit_length() <= prec + 1 or abs(m).denominator & (abs(m).denominator - 1) == 0, "dyadic mid")
                _check(round_even(m, prec) == m, "mid fits prec %d" % prec)
            _check(r >= 0 and from_fraction(r)[2].bit_length() <= RAD_BITS, "radius bits")
    for prec in (2, 53, 113):
        a = draw_matrix(rng, prec, 3, 4, "rand", "A")
        b = draw_matrix(rng, prec, 4, 2, "rand", "B")
        for i in range(3):
            for j in range(2):
                col = [b[t][j] for t in range(4)]
                p0, p1, p2, tref, first, sabs = product_entry(a[i], col, prec)
                true_rad = sum((abs(x[0]) * y[1] + x[1] * abs(y[0]) + x[1] * y[1] for x, y in zip(a[i], col)), Fraction(0))
                _check(true_rad == first, "first-order formula is the exact product radius")
                _check(abs(p1 - p0) <= first and abs(p2 - p0) <= first, "corners inside the exact product box")
                _check(abs(p0) <= sabs, "|P0| <= S")
                bound = first + pow2(1 - prec) * sabs
                _check(tref >= bound and tref <= bound * (1 + pow2(-63)), "TREF rounding")
    # mutant 902 (radius accounting, CONTRACT 16.S6.11): on cancel, tie and exact rows the rounding is unavoidable, and R1 alone
    # (no rounding term) cannot reach P0 from any value representable at prec: dropping the rounding terms breaks containment.
    def r1_misses(kind, prec, m, k, n):
        ra = draw_matrix(rng, prec, m, k, kind, "A")
        rb = draw_matrix(rng, prec, k, n, kind, "B")
        hits = 0
        for i in range(m):
            for j in range(n):
                p0, _, _, _, first, _ = product_entry(ra[i], [rb[t][j] for t in range(k)], prec)
                lo = goldlib.round_fraction(p0, prec, "Down")[0]
                hi = goldlib.round_fraction(p0, prec, "Up")[0]
                if lo != p0 and first < min(p0 - lo, hi - p0):
                    hits += 1
        return hits
    _check(r1_misses("cancel", 53, 1, 4, 1) == 1, "cancel row: R1 alone misses P0")
    _check(r1_misses("tie", 2, 1, 2, 1) == 1, "tie row: R1 alone misses P0")
    _check(r1_misses("exact", 64, 3, 5, 2) >= 1, "exact row: R1 alone misses P0 in some entry")
    ca = draw_matrix(rng, 53, 1, 4, "cancel", "A")
    cb = draw_matrix(rng, 53, 4, 1, "cancel", "B")
    p0 = product_entry(ca[0], [cb[t][0] for t in range(4)], 53)[0]
    _check(p0 == 1 + pow2(-60) and round_even(p0, 53) == 1, "cancel row: 1 + 2^-60 is not representable at 53 bits")
    ta = draw_matrix(rng, 2, 1, 2, "tie", "A")
    tb = draw_matrix(rng, 2, 2, 1, "tie", "B")
    p0 = product_entry(ta[0], [tb[t][0] for t in range(2)], 2)[0]
    _check(p0 == Fraction(5, 8) and round_even(p0, 2) == Fraction(1, 2), "tie row: 5/8 is a tie at 2 bits")
    # spread row: the bookkeeping text equals the generic Fraction path at E = 200
    sa = draw_matrix(rng, 53, 1, 3, "spread_small", "A")
    sb = draw_matrix(rng, 53, 3, 1, "spread_small", "B")
    ga, gb = spread_mids(SPREAD_SMALL_E)
    _check(entries(sa) == ga and entries(sb) == gb, "spread mids text")
    for prec in (1, 53, 256):
        out = product_entry(sa[0], [sb[t][0] for t in range(3)], prec)
        _check(" ".join(enc(v) for v in out[:4]) == spread_out(SPREAD_SMALL_E, prec), "spread out text at prec %d" % prec)
    _check(spread_out(SPREAD_E, 53).split(" ")[10] == str(SPREAD_E + 3 - 53), "spread TREF exponent")
    # CONTRACT 16.S6.9 accounting (Spence option 1): with U = 2^(1-prec) S and w = prec + g, the K products and K partial sums
    # cost < c U / 4 together, the final rounding (one full ulp of |s_K|) costs U (1 + c 2^(-1-prec)), and the Mag evaluation
    # costs (1 + 2^-63)^(4K+8) <= 1 + 2^-43. The ratio bound against TREF = R1 + U is (1.25 + 2^-prec)(1 + 2^-43) < 2 at
    # every prec >= 1 (the earlier "< 1.6 TREF" claim failed at prec 1 to 3). Exact Fractions, U normalised to 1.
    def accounting_ratio(prec, k):
        g = k.bit_length() + 2                              # ceil(log2(k + 1)) + 2
        _check(pow2(g) >= 4 * (k + 1), "2^g >= 4(K+1) at K=%d" % k)
        w = prec + g
        _check(2 * k * pow2(-w) <= Fraction(1, 2), "c bound premise at prec %d K %d" % (prec, k))
        c = 1 + 4 * k * pow2(-w)                            # >= (1 + 2^-w)^(2K) since 2K 2^-w <= 1/2
        steps = 2 * k * pow2(-w) * c / pow2(1 - prec)       # (a) in units of U
        _check(steps < c / 4, "(a) steps < c/4 at prec %d K %d" % (prec, k))
        final = 1 + steps * pow2(1 - prec)                  # (b) 2^(1-prec)(S + steps U) / U
        total = steps + final                               # (c)
        _check(total <= Fraction(5, 4) + pow2(-prec), "(c) rounding total <= 1.25 + 2^-prec at prec %d K %d" % (prec, k))
        n = 4 * k + 8
        _check(n * pow2(-63) <= Fraction(1, 2) and 2 * n * pow2(-63) <= pow2(-43), "(d) Mag loss <= 2^-43 at K %d" % k)
        return max(Fraction(1), total) * (1 + pow2(-43))
    for prec in (1, 2, 3, 4):
        cap = (Fraction(5, 4) + pow2(-prec)) * (1 + pow2(-43))
        for k in (1, 2, 3, 4, 7, 8, 31, 32, 1000, 2 ** 16):
            ratio = accounting_ratio(prec, k)
            _check(ratio <= cap and ratio < 2, "16.S6.9 ratio %s at prec %d K %d" % (ratio, prec, k))
    _check((Fraction(5, 4) + pow2(-1)) * (1 + pow2(-43)) < 2, "16.S6.9 worst case prec 1 below 2")
    # norms: norm_1 of A is norm_inf of A transposed; frob HI is a tight 64-bit upper bound of sqrt(S); vector norms of a row
    a = draw_matrix(rng, 64, 3, 5, "rand", "A")
    at = [[a[i][j] for i in range(3)] for j in range(5)]
    _check(norm_values("norm_1", a)[0] == norm_values("norm_inf", at)[0], "norm_1 = norm_inf of the transpose")
    s, hi = norm_values("norm_frob", a)
    _check(hi * hi >= s and (hi * (1 - pow2(-62))) ** 2 < s, "frob HI tight")
    _check(norm_values("norm_inf", [])[0] == 0 and norm_values("norm_frob", [[], []])[0] == 0, "empty norms are 0")
    v = draw_matrix(rng, 64, 1, 6, "rand", "A")
    vals = [abs(m) + r for (m, r) in v[0]]
    _check(norm_values("vnorm_inf", v)[0] == max(vals) and norm_values("vnorm_1", v)[0] == sum(vals, Fraction(0)), "vector norms")
    _check(norm_values("vnorm_2", v)[0] == sum((x * x for x in vals), Fraction(0)), "vector 2-norm square")
    _check(all(norm_values(op, [[]])[0] == 0 for op in VEC_NORMS), "empty vector norms are 0")
    # case tables and widths
    _check(len(gate_cases()) == EXPECTED_ROWS["gate"] and len(heavy_cases()) == EXPECTED_ROWS["heavy"], "case counts")
    _check(all(max(c[2], c[3], c[4]) <= 8 and c[1] <= 256 for c in gate_cases()), "gate dims <= 8, prec <= 256")
    line = make_row(0, ("matmul", 53, 2, 3, 2, "rand"), goldlib.SplitMix64(1))
    _check(len(line.split(" ")) == row_width("matmul", 2, 3, 2), "row width")
    line = make_row(0, ("dot", 53, 1, 3, 1, "spread"), goldlib.SplitMix64(1))
    _check(len(line.split(" ")) == row_width("dot", 1, 3, 1) and len(line) < 400, "spread row width and size")
    print("gen_matrix_goldens selftest: OK")


def main(argv):
    goldlib.require_mpmath()
    skip_heavy = "--skip-heavy" in argv
    argv = [a for a in argv if a != "--skip-heavy"]
    selftest()
    if argv == ["--selftest"]:
        return 0
    if len(argv) == 2 and argv[0] == "--check":
        os.makedirs(argv[1], exist_ok=True)
        n = write_tier("gate", os.path.join(argv[1], GATE_NAME), hashlib.sha256())
        print("gate rows: %d" % n)
        bad = goldlib.check_against(argv[1], GOLDEN_DIR, [GATE_NAME])
        if not skip_heavy:
            h = hashlib.sha256()
            n = write_tier("heavy", None, h)
            ok = h.hexdigest() == pinned_heavy()
            print("check: %s %s (heavy rows: %d, sha256 %s)" % ("OK" if ok else "DIFF", HEAVY_NAME, n, h.hexdigest()))
            bad += 0 if ok else 1
        return 1 if bad else 0
    heavy_dir = None
    if len(argv) == 2 and argv[0] == "--heavy":
        heavy_dir = argv[1]
        os.makedirs(heavy_dir, exist_ok=True)
    elif argv != []:
        raise SystemExit("usage: scripts/gen_matrix_goldens.py [--check DIR | --heavy DIR | --selftest] [--skip-heavy]")
    gh = hashlib.sha256()
    n = write_tier("gate", os.path.join(GOLDEN_DIR, GATE_NAME), gh)
    print("gate rows: %d" % n)
    heavy_sha = None
    if not skip_heavy:
        hh = hashlib.sha256()
        n = write_tier("heavy", os.path.join(heavy_dir, HEAVY_NAME) if heavy_dir else None, hh)
        heavy_sha = hh.hexdigest()
        print("heavy rows: %d sha256 %s" % (n, heavy_sha))
    write_sums(gh.hexdigest(), heavy_sha)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
