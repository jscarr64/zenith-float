#!/home/jsc64/team-env/oracle-2026-10-30/bin/python
"""Exact linear-algebra reference goldens for zenith-bound (zenith-float 2.0).

Pinned oracle: sympy 1.14.0, mpmath 1.3.0. Run it with ~/team-env/oracle-2026-10-30/bin/python
(the shebang is the team env; the version check below rejects other pins).
Usage:
    gen_la_goldens.py                 write the gate files tests/goldens/slices/<slice>/la_<op>.txt, each slice's
                                      SHA256SUMS (gate files) and HEAVY.SHA256 (sha256 of the uncommitted heavy file)
    gen_la_goldens.py --heavy DIR     also write the heavy files la_<op>_heavy.txt into DIR (never committed)
    gen_la_goldens.py --check [DIR]   regenerate; byte-compare every gate file with the committed one (copies go to DIR
                                      if given, per the goldlib --check DIR convention) and compare the sha256 of every
                                      regenerated heavy file with HEAVY.SHA256. --skip-heavy limits it to the gate files.

Memory (Paul, 2026-10-24, after a 58 GB OOM): rows stream to unlinked spool files under target/la-goldens-spool/ and are hashed
incrementally; run every invocation under poprun/popslot, /usr/bin/time -v and a cap, e.g.
    prlimit --as=4294967296 -- /usr/bin/time -v ~/team-env/oracle-2026-10-30/bin/python scripts/gen_la_goldens.py --check --skip-heavy
--skip-heavy never computes the heavy cosmo_* rows (gate only; HEAVY.SHA256 is left as is in write mode).

Slices (S0, tests/goldens/slices/<slice>/): LA1 cholesky, LA2 lu (+ solve, logdet), LA3 qr, LA4 gls, LA5a eigen (symmetric), LA5b eigen (general, stub: pending complex Ball).
Gate = every case with all dims <= 12 and every input integer < 2^(10^4+64) in size; heavy = the rest
(n = 16..32, 2^+-10^6 scales). case_ids are fixed strings; the random streams are drawn in the same order for
gate and heavy, so splitting never changes a row. The generator is a pure exact oracle and never calls zenith-bound.

All arithmetic is exact (gmpy2.mpq / sympy QQ / Python ints). Irrational
outputs are emitted as rigorous dyadic enclosures computed by exact integer square
roots, exact rational root refinement (sympy, Sturm-based), or mpmath.iv interval
arithmetic. No floats are used anywhere. Every identity and every enclosure is
re-checked here before a row is emitted.
"""
# popguard gate (Harper/Paul, after the 11:48 58 GB OOM): runs before any heavy import or work. Outside a popslot/poprun
# memory scope ~/pop-logs/bin/popguard.sh exits nonzero (3) and so do we, with its rc; 3 if the guard is missing or cannot run.
# Selftest: tools/gen_la_guard_selftest.sh.
import os as _os, subprocess as _subprocess, sys as _sys
def _popguard():
    guard = _os.path.expanduser("~/pop-logs/bin/popguard.sh")
    try:
        rc = _subprocess.run([guard, "gen_la_goldens"], check=False).returncode
    except OSError as e:
        print(f"gen_la_goldens: popguard missing or cannot run ({guard}: {e}); refusing to start. "
              "Run it via ~/pop-logs/bin/poprun --agent NAME --job JOB -- ~/pop-logs/bin/popslot -- CMD.", file=_sys.stderr)
        _sys.exit(3)
    if rc != 0:
        print(f"gen_la_goldens: popguard refused (rc {rc}); not running outside a popslot memory scope.", file=_sys.stderr)
        _sys.exit(rc if rc > 0 else 3)
_popguard()

import sys, os, re, random, hashlib, argparse, difflib, tempfile, signal, time
# Exact rationals are gmpy2.mpq (GMP, lowest terms like fractions.Fraction, same str of numerator/denominator), not
# fractions.Fraction: Fraction normalizes every op with CPython math.gcd, which is quadratic. On the heavy cosmo GLS rows
# (entries of 0.5M-3.9M bits) that made la4-heavy2 sit 90 min in op gls (probe: cosmo_spread.p2e1000000_32x2 S^-1 X 98.7 s and
# the exact S x == v check > 600 s with Fraction; 0.45 s and 94 s with mpq; whole row 213 s). gmpy2 is the pinned mpmath 1.3.0
# backend in the oracle venv (no new dep); refuse to run without it rather than fall back to the slow path.
try:
    from gmpy2 import mpq as F
except ImportError as e:
    raise SystemExit(f"gen_la_goldens.py: gmpy2 missing ({e}); run with ~/team-env/oracle-2026-10-30/bin/python (mpmath backend gmpy)")
from math import isqrt

SYMPY_PIN = "1.14.0"
MPMATH_PIN = "1.3.0"
import sympy, mpmath
from mpmath import iv
if sympy.__version__ != SYMPY_PIN or mpmath.__version__ != MPMATH_PIN:
    raise SystemExit(f"gen_la_goldens.py: oracle pin mismatch: sympy {sympy.__version__} (need {SYMPY_PIN}), "
                     f"mpmath {mpmath.__version__} (need {MPMATH_PIN}); run with ~/team-env/oracle-2026-10-30/bin/python")
sys.set_int_max_str_digits(0)

GEN = "scripts/gen_la_goldens.py"
SEED = 20261024
BITS = 256          # default enclosure width (relative)
ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
SLICES = os.path.join(ROOT, "tests", "goldens", "slices")
SLICE_OF = {"cholesky": "cholesky", "lu": "lu", "solve": "lu", "logdet": "cholesky", "qr": "qr", "gls": "gls", "eig": "eigen", "eig_general": "eigen"}
SLICE_ID_OP = {"eig_general": "LA5b"}
SLICE_ID = {"cholesky": "LA1", "lu": "LA2", "qr": "LA3", "gls": "LA4", "eigen": "LA5a"}
GATE_MAX_DIM = 12
EXT_BITS = (4096,)        # extra Cholesky/QR/logdet enclosure widths for gate-tier cases, emitted on the heavy tier (la_<op>_enc<b>_heavy.txt)
GATE_MAX_BITS = 10 ** 4 + 64
LOGDET_EXT_BITS = (4096, 16384)  # logdet headroom widths (cosmology x >= 64): every gate-tier case plus every cosmo_* case, heavy tier

def is_heavy(*mats):
    for M in mats:
        if not M: continue
        if isinstance(M[0], list):
            if len(M) > GATE_MAX_DIM or len(M[0]) > GATE_MAX_DIM: return True
            vals = flat(M)
        else:
            if len(M) > GATE_MAX_DIM: return True
            vals = M
        if any(F(v).numerator.bit_length() > GATE_MAX_BITS or F(v).denominator.bit_length() > GATE_MAX_BITS for v in vals): return True
    return False

def err(variant, detail):
    """Refusal tail: ERR <zb BoundError variant> <detail> (Spence/Harper 2026-10-24)."""
    assert variant in ("NotPositiveDefinite", "Singular", "DimensionMismatch", "NotSymmetric")
    return ["ERR", variant, detail]

def ldl_err(kind, who):
    return err("NotSymmetric", who) if kind == "not_symmetric" else err("NotPositiveDefinite", who)

def parse_dims(t):
    if "x" in t:
        m, n = t.split("x"); return int(m), int(n)
    return int(t), int(t)

def check_fields(text, expect):
    """Assert a row's field count against its dims. expect(m, n, outs_tail) -> (n_inputs, allowed output counts)."""
    f = text.split(" ")
    assert all(f), "empty field / double space"
    m, n = parse_dims(f[1])
    nin, nouts = expect(m, n)
    tail = len(f) - 2 - nin
    ok = (f[2 + nin] == "ERR" and tail == 3) or (tail in nouts and "ERR" not in f)
    assert ok, (f[0], len(f), nin, nouts)


# ---------------------------------------------------------------- formatting
def q(x):
    x = F(x)
    return f"{x.numerator}/{x.denominator}"

def dy(m, e):
    """Canonical dyadic m*2^e: m odd (or m == 0 with e == 0)."""
    if m == 0:
        return "0*2^0"
    while m % 2 == 0:
        m //= 2; e += 1
    return f"{m}*2^{e}"

def dyval(m, e):
    return F(m) * (F(2) ** e)

def ilog2(x):
    """floor(log2 |x|) for nonzero rational x, exactly."""
    x = abs(F(x)); n, d = x.numerator, x.denominator
    k = n.bit_length() - d.bit_length()
    if F(2) ** k > x: k -= 1
    if F(2) ** (k + 1) <= x: k += 1
    return k

def floor_q(x): return x.numerator // x.denominator
def ceil_q(x): return -((-x.numerator) // x.denominator)

def grid(x, bits):
    """Exponent e of the grid 2^e used for a value of magnitude |x|: e = floor(log2|x|) - bits + 1."""
    return ilog2(x) - bits + 1

def enclose_rational(x, bits=BITS):
    x = F(x)
    if x == 0: return (0, 0), (0, 0)
    e = grid(x, bits); s = x / (F(2) ** e)
    return (floor_q(s), e), (ceil_q(s), e)

def enclose_csqrt(c, s, bits=BITS):
    """Enclose y = c*sqrt(s), c, s rational, s >= 0, by exact integer isqrt."""
    c, s = F(c), F(s)
    assert s >= 0
    if c == 0 or s == 0: return (0, 0), (0, 0)
    y2 = c * c * s
    e = ilog2(y2) // 2 - bits + 1           # grid exponent ~ floor(log2 |y|) - bits + 1
    t = y2 / (F(4) ** e)                    # (|y| / 2^e)^2
    m = isqrt(floor_q(t))
    hi = m if F(m * m) == t else m + 1
    lo = m
    assert F(lo * lo) <= t <= F(hi * hi)
    if c < 0: lo, hi = -hi, -lo
    L, H = dyval(lo, e), dyval(hi, e)
    # containment check: lo <= y <= hi, exact via squares and sign
    if c > 0: assert 0 <= L and L * L <= y2 <= H * H
    else: assert H <= 0 and H * H <= y2 <= L * L
    return (lo, e), (hi, e)

def pair(lo, hi): return f"{dy(*lo)} {dy(*hi)}"

# ---------------------------------------------------------------- exact matrix helpers
def ident(n): return [[F(int(i == j)) for j in range(n)] for i in range(n)]
def mm(A, B): return [[sum((A[i][k] * B[k][j] for k in range(len(B))), F(0)) for j in range(len(B[0]))] for i in range(len(A))]
def tr(A): return [list(r) for r in zip(*A)]
def flat(A): return [x for r in A for x in r]
def issym(A): return all(A[i][j] == A[j][i] for i in range(len(A)) for j in range(len(A)))

def lu_pp(A):
    """Partial pivoting: pivot = max |U[i][k]| over i >= k, ties -> smallest i; zero column -> no step."""
    n = len(A); U = [list(map(F, r)) for r in A]; L = ident(n); perm = list(range(n))
    for k in range(n):
        p = max(range(k, n), key=lambda i: (abs(U[i][k]), -i))
        if U[p][k] == 0: continue
        if p != k:
            U[k], U[p] = U[p], U[k]; perm[k], perm[p] = perm[p], perm[k]
            for j in range(k): L[k][j], L[p][j] = L[p][j], L[k][j]
        for i in range(k + 1, n):
            l = U[i][k] / U[k][k]
            if l:
                L[i][k] = l
                U[i] = [U[i][j] - l * U[k][j] for j in range(n)]
    PA = [list(map(F, A[perm[i]])) for i in range(n)]
    assert mm(L, U) == PA
    return perm, L, U

def perm_sign(perm):
    s, seen = 1, [False] * len(perm)
    for i in range(len(perm)):
        if not seen[i]:
            j, c = i, 0
            while not seen[j]: seen[j] = True; j = perm[j]; c += 1
            if c % 2 == 0: s = -s
    return s

def det(A):
    perm, L, U = lu_pp(A); d = F(perm_sign(perm))
    for i in range(len(A)): d *= U[i][i]
    return d

def solve(A, b):
    n = len(A); perm, L, U = lu_pp(A)
    if any(U[i][i] == 0 for i in range(n)): return None
    y = []
    for i in range(n): y.append(F(b[perm[i]]) - sum((L[i][j] * y[j] for j in range(i)), F(0)))
    x = [F(0)] * n
    for i in reversed(range(n)): x[i] = (y[i] - sum((U[i][j] * x[j] for j in range(i + 1, n)), F(0))) / U[i][i]
    assert [sum((F(A[i][j]) * x[j] for j in range(n)), F(0)) for i in range(n)] == [F(v) for v in b]
    return x

def ldl(A):
    n = len(A)
    if not issym(A): return "not_symmetric"
    L = ident(n); D = [F(0)] * n
    for j in range(n):
        D[j] = F(A[j][j]) - sum((L[j][k] ** 2 * D[k] for k in range(j)), F(0))
        if D[j] <= 0: return "not_spd"
        for i in range(j + 1, n):
            L[i][j] = (F(A[i][j]) - sum((L[i][k] * L[j][k] * D[k] for k in range(j)), F(0))) / D[j]
    assert mm(mm(L, [[D[i] if i == j else F(0) for j in range(n)] for i in range(n)]), tr(L)) == [list(map(F, r)) for r in A]
    return L, D

def gram_schmidt(A):
    m, n = len(A), len(A[0]); cols = tr([list(map(F, r)) for r in A])
    Qt, Rt, d = [], ident(n), []
    for j in range(n):
        v = list(cols[j])
        for i in range(j):
            Rt[i][j] = sum((a * b for a, b in zip(Qt[i], cols[j])), F(0)) / d[i]
            v = [a - Rt[i][j] * b for a, b in zip(v, Qt[i])]
        dj = sum((a * a for a in v), F(0))
        if dj == 0: return None
        Qt.append(v); d.append(dj)
    QtM = tr(Qt)
    assert mm(QtM, Rt) == [list(map(F, r)) for r in A]
    G = mm(Qt, QtM)
    assert all(G[i][j] == (d[i] if i == j else 0) for i in range(n) for j in range(n))
    return QtM, Rt, d

# ---------------------------------------------------------------- enclosure via mpmath.iv
def raw_to_q(t):
    s, man, e, bc = t
    if man == 0: return F(0)
    v = F(man) * (F(2) ** e)
    return -v if s else v

def iv_from_q(x, bits):
    lo, hi = enclose_rational(x, bits)
    a = iv.mpf(lo[0]) * iv.mpf(2) ** lo[1]
    b = iv.mpf(hi[0]) * iv.mpf(2) ** hi[1]
    return iv.mpf([a.a, b.b])

def iv_log_abs(x, bits):
    """[lo, hi] rational enclosure of ln|x| from mpmath.iv at working precision bits+64."""
    old = iv.prec
    try:
        iv.prec = bits + 64
        r = iv.log(iv_from_q(abs(x), bits + 64))
        a, b = r._mpi_
        return raw_to_q(a), raw_to_q(b)
    finally:
        iv.prec = old

def enclose_log(x, bits):
    a, b = iv_log_abs(x, bits)
    if a == 0 and b == 0: return (0, 0), (0, 0)
    mag = max(abs(a), abs(b))
    e = grid(mag, bits)
    lo = (floor_q(a / F(2) ** e), e); hi = (ceil_q(b / F(2) ** e), e)
    # containment check: independent higher-precision enclosure must lie inside
    a2, b2 = iv_log_abs(x, 2 * bits + 64)
    assert dyval(*lo) <= a2 <= b2 <= dyval(*hi)
    if x == 1 or x == -1: assert dyval(*lo) <= 0 <= dyval(*hi)
    return lo, hi

# ---------------------------------------------------------------- eigenvalues (exact algebraic)
X = sympy.Symbol("x")
def eig_exact(A, bits=BITS):
    n = len(A)
    M = sympy.Matrix([[sympy.Rational(int(a.numerator), int(a.denominator)) for a in r] for r in A])
    P = sympy.Poly(M.charpoly(X).as_expr(), X, domain="QQ")
    coeffs = [F(int(c.p), int(c.q)) for c in P.all_coeffs()]
    _, facs = P.factor_list()
    out = []
    for f, mult in facs:
        if f.degree() == 1:
            a, b = f.all_coeffs()
            rr = -b / a; r = F(int(rr.p), int(rr.q))
            lo, hi = enclose_rational(r, bits)
            for _ in range(mult): out.append((r, q(r), lo, hi))
            continue
        ivs = f.intervals()
        assert sum(k for _, k in ivs) == f.degree(), "non-real root in symmetric charpoly?"
        for (s, t), k in ivs:
            s, t = sympy.Rational(s), sympy.Rational(t)
            while True:
                if s * t > 0:
                    w = min(abs(s), abs(t)); e = grid(F(int(w.p), int(w.q)), bits)
                    if t - s <= sympy.Rational(2) ** e / 4: break
                    s, t = f.refine_root(s, t, eps=sympy.Rational(2) ** e / 8)
                else:
                    s, t = f.refine_root(s, t, eps=(t - s) / 4)
            sF, tF = F(int(s.p), int(s.q)), F(int(t.p), int(t.q))
            lo = (floor_q(sF / F(2) ** e), e); hi = (ceil_q(tF / F(2) ** e), e)
            L, H = dyval(*lo), dyval(*hi)
            assert f.count_roots(sympy.Rational(int(L.numerator), int(L.denominator)), sympy.Rational(int(H.numerator), int(H.denominator))) == 1
            for _ in range(mult): out.append(((sF + tF) / 2, "alg", lo, hi))
    out.sort(key=lambda r: r[0])
    for a, b in zip(out, out[1:]):  # ordering must be unambiguous: distinct eigenvalues get disjoint enclosures
        # two rational eigenvalues closer than the grid may share a cell; the exact p/q tokens disambiguate
        if not ((a[1] != "alg" and b[1] != "alg" and a[0] <= b[0]) or (a[1] == "alg" and a[2:] == b[2:]) or dyval(*a[3]) < dyval(*b[2])): print("ORDER", n, a[1], b[1], a[2:], b[2:]); raise AssertionError
    assert sum(1 for _ in out) == n
    return coeffs, out

# ---------------------------------------------------------------- case corpus
rng = random.Random(SEED)
def rq(lo=-9, hi=9, dmax=9): return F(rng.randint(lo, hi), rng.randint(1, dmax))
def rmat(m, n, rational=True):
    return [[rq() if rational else F(rng.randint(-9, 9)) for _ in range(n)] for _ in range(m)]
def hilbert(n): return [[F(1, i + j + 1) for j in range(n)] for i in range(n)]
def spd(n, rational=True):
    B = rmat(n, n, rational); A = mm(B, tr(B))
    for i in range(n): A[i][i] += 1
    return A
def rsym(n):
    A = [[F(0)] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n): A[i][j] = A[j][i] = F(rng.randint(-9, 9))
    return A
def zero(n): return [[F(0)] * n for _ in range(n)]
H6, T6 = F(2) ** 1000000, F(1, 2 ** 1000000)
B3, T3 = F(2) ** 1000, F(1, 2 ** 1000)
B4 = F(2) ** 10000
EPS = F(1, 2 ** 200)

def square_corpus(kind):
    """Common square cases. kind in lu/solve/logdet."""
    C = []
    C += [("onebyone.pos", [[F(3, 7)]]), ("onebyone.neg", [[F(-5, 2)]]), ("onebyone.zero", [[F(0)]])]
    C += [("zero.%02d" % n, zero(n)) for n in (2, 5)]
    C += [("basic.2", [[F(4), F(3)], [F(6), F(3)]]), ("basic.3", [[F(2), F(-1), F(0)], [F(-1), F(2), F(-1)], [F(0), F(-1), F(2)]])]
    C += [("perm.swap2", [[F(0), F(1)], [F(1), F(0)]]), ("perm.zero_pivot3", [[F(0), F(2), F(1)], [F(1), F(1), F(1)], [F(3), F(0), F(2)]]),
          ("perm.cyclic4", [[F(int((j - i) % 4 == 1)) * (i + 1) for j in range(4)] for i in range(4)]),
          ("perm.tie3", [[F(1), F(2), F(3)], [F(-1), F(5), F(1)], [F(1), F(0), F(7)]])]
    C += [("singular.dup_rows3", [[F(1), F(2), F(3)], [F(1), F(2), F(3)], [F(4), F(5), F(6)]]),
          ("singular.classic3", [[F(1), F(2), F(3)], [F(4), F(5), F(6)], [F(7), F(8), F(9)]])]
    C += [("rankdef.r1_4", [[F(i + 1) * F(j + 1, 3) for j in range(4)] for i in range(4)]),
          ("rankdef.r2_5", mm(rmat(5, 2), rmat(2, 5)))]
    C += [("nearsing.2", [[F(1), F(1)], [F(1), 1 + EPS]]), ("nearsing.3", [[F(1), F(2), F(3)], [F(4), F(5), F(6)], [F(7), F(8), 9 + EPS]])]
    C += [("hilbert.%02d" % n, hilbert(n)) for n in range(5, 13)]
    C += [("scale_huge.1x1", [[H6 + 1]]), ("scale_tiny.1x1", [[T6 * F(3, 5)]]),
          ("scale_mixed.diag2", [[H6, F(0)], [F(0), T6]]),
          ("scale_huge.2", [[B4 + 1, B4], [B4, B4 - 3]]), ("scale_tiny.3", [[a * T3 for a in r] for r in rmat(3, 3)]),
          ("scale_mixed.3", [[B3, F(1), T3], [F(1), T3, B3], [T3, B3, F(1)]])]
    for n in (2, 3, 4, 5, 6, 8):
        C.append(("random_q.%02d" % n, rmat(n, n)))
    for n in (12, 16, 24, 32):
        C.append(("random_z.%02d" % n, rmat(n, n, rational=False)))
    return C

def spd_corpus():
    C = [("onebyone.pos", [[F(9, 4)]]), ("onebyone.irr", [[F(2)]])]
    C += [("basic.2", [[F(4), F(2)], [F(2), F(3)]]), ("basic.tridiag5", [[F(2) if i == j else (F(-1) if abs(i - j) == 1 else F(0)) for j in range(5)] for i in range(5)])]
    C += [("hilbert.%02d" % n, hilbert(n)) for n in range(5, 13)]
    C += [("nearsing.2", [[F(1), F(1)], [F(1), 1 + EPS]])]
    C += [("scale_huge.1x1", [[H6]]), ("scale_tiny.1x1", [[T6 * 2]]), ("scale_mixed.diag2", [[H6, F(0)], [F(0), T6]]),
          ("scale_huge.2", [[B4 * 4, B4], [B4, B4]]), ("scale_tiny.3", [[a * T3 for a in r] for r in spd(3)])]
    for n in (2, 3, 4, 6, 8): C.append(("random_q.%02d" % n, spd(n)))
    for n in (12, 16, 24, 32): C.append(("random_z.%02d" % n, spd(n, rational=False)))
    # refusals
    C += [("refuse.zero2", zero(2)), ("refuse.onebyone_neg", [[F(-1)]]), ("refuse.onebyone_zero", [[F(0)]]),
          ("refuse.indefinite2", [[F(1), F(2)], [F(2), F(1)]]), ("refuse.psd_singular3", [[F(1), F(1), F(0)], [F(1), F(1), F(0)], [F(0), F(0), F(1)]]),
          ("refuse.nonsym2", [[F(2), F(1)], [F(0), F(2)]]), ("refuse.neg_def3", [[-a for a in r] for r in spd(3)]),
          ("refuse.nearsing_neg2", [[F(1), F(1)], [F(1), 1 - EPS]])]
    return C

# ---------------------------------------------------------------- cosmology-scale corpus (LA1 logdet + LA4 gls, 2026-10-24)
# Jeff: zb must exceed cosmology work in range, precision and dimension. Families (case_id family prefix cosmo_*):
#   spread    SPD L D L^T / D + U U^T with eigenvalue-scale spread 10^+-300 and 2^+-(10^3..10^6) inside one matrix
#   cond      near-singular SPD, cond ~ 10^50 .. 10^200 (L D L^T with D down to 10^-200; scaled Hilbert)
#   toeplitz  CMB-style stationary covariances (Kac-Murdock-Szego rho^|i-j|, scaled tridiagonal)
#   blockdiag, lowrank (diag noise + low-rank signal), fisher (J^T C^-1 J), bao / cmb_toeplitz (GLS designs)
# Exact determinants of structured cases come from an exact identity (prod D, Sylvester, block product, KMS closed form, D A D);
# on the gate tier every identity is re-checked against the LU determinant / the generic exact solve.
def P10(k): return F(10) ** k
def P2(k): return F(2) ** k
def unit_lower(n, w=2, band=None):
    L = ident(n)
    for i in range(n):
        for j in range(i):
            if band is None or i - j <= band: L[i][j] = F(rng.randint(-w, w))
    return L
def ldlt(L, D):
    n = len(L)
    return [[sum((L[i][k] * D[k] * L[j][k] for k in range(min(i, j) + 1)), F(0)) for j in range(n)] for i in range(n)]
def spread(n, hi, lo, base):
    """n values base^k, k from hi down to lo in integer steps (exact)."""
    return [F(base) ** (hi - (hi - lo) * i // (n - 1)) for i in range(n)]
def kms(n, rho): return [[F(rho) ** abs(i - j) for j in range(n)] for i in range(n)]
def blockdiag(blocks):
    n = sum(len(b) for b in blocks); A = zero(n); o = 0
    for b in blocks:
        for i in range(len(b)):
            for j in range(len(b)): A[o + i][o + j] = F(b[i][j])
        o += len(b)
    return A
def lowrank_diag(d, U):
    n, k = len(U), len(U[0])
    return [[(F(d[i]) if i == j else F(0)) + sum((U[i][t] * U[j][t] for t in range(k)), F(0)) for j in range(n)] for i in range(n)]
def scale_sym(A, d): return [[d[i] * A[i][j] * d[j] for j in range(len(A))] for i in range(len(A))]
def prod(v):
    r = F(1)
    for x in v: r *= x
    return r
def capacitance(d, U):
    k = len(U[0])
    return [[F(int(s == t)) + sum((U[i][s] * U[i][t] / d[i] for i in range(len(d))), F(0)) for t in range(k)] for s in range(k)]
def det_lowrank(d, U): return prod(d) * det(capacitance(d, U))
def woodbury(d, U):
    """Exact solver for S = diag(d) + U U^T (Sherman-Morrison-Woodbury): no n^3 fill-in at 2^+-10^6 scales."""
    n, k = len(U), len(U[0]); K = capacitance(d, U)
    def sol(v):
        dv = [F(v[i]) / d[i] for i in range(n)]
        w = solve(K, [sum((U[i][t] * dv[i] for i in range(n)), F(0)) for t in range(k)])
        return [dv[i] - sum((U[i][t] * w[t] for t in range(k)), F(0)) / d[i] for i in range(n)]
    return sol
def block_solver(blocks):
    def sol(v):
        out = []; o = 0
        for b in blocks: out += solve(b, v[o:o + len(b)]); o += len(b)
        return out
    return sol
def matvec(A, x): return [sum((A[i][j] * x[j] for j in range(len(x)) if A[i][j]), F(0)) for i in range(len(A))]
def spd_small(n, w=3):
    B = [[F(rng.randint(-w, w)) for _ in range(n)] for _ in range(n)]; A = mm(B, tr(B))
    for i in range(n): A[i][i] += 1
    return A
def imat(m, n, w=3): return [[F(rng.randint(-w, w)) for _ in range(n)] for _ in range(m)]

def logdet_cosmo_corpus():
    """Yields (cid, A, det by exact identity or None) one case at a time (streamed: no corpus list) for the cosmology-scale logdet rows; gate n <= 12, heavy n = 32/64, 2^+-10^5..10^6."""
    def ld(cid, L, D): return (cid, ldlt(L, D), prod(D))
    yield ld("cosmo_spread.e300_08", unit_lower(8), spread(8, 300, -300, 10))
    yield ld("cosmo_spread.e300_12", unit_lower(12, band=2), spread(12, 300, -300, 10))
    yield ld("cosmo_spread.p2e1000_06", unit_lower(6), spread(6, 1000, -1000, 2))
    yield ld("cosmo_spread.indef_e300_05", unit_lower(5), [P10(300), -P10(-300), F(1), P10(-100), P10(100)])
    d = spread(4, 10000, -10000, 2); U = imat(4, 1); U[0][0] = F(1)
    yield ("cosmo_spread.p2e10000_lowrank4", lowrank_diag(d, U), det_lowrank(d, U))
    for k, n in ((50, 6), (100, 8), (200, 10)):
        yield ld(f"cosmo_cond.e{k}_{n:02d}", unit_lower(n, band=1), [F(1)] * (n - 1) + [P10(-k)])
    dh = [P10(-30 * i) for i in range(10)]
    yield ("cosmo_cond.hilbert_scaled10", scale_sym(hilbert(10), dh), det(hilbert(10)) * prod(dh) ** 2)
    rho = 1 - P10(-20)
    yield ("cosmo_toeplitz.kms08_rho1m1e-20", kms(8, rho), (1 - rho * rho) ** 7)
    yield ("cosmo_toeplitz.kms12_rho9_10", kms(12, F(9, 10)), (1 - F(81, 100)) ** 11)
    t12 = [[F(2) if i == j else (F(-1) if abs(i - j) == 1 else F(0)) for j in range(12)] for i in range(12)]
    yield ("cosmo_toeplitz.tridiag12_p2e1000", [[P2(1000) * a for a in r] for r in t12], F(13) * P2(12000))
    bl = [[[P10(s) * a for a in r] for r in spd_small(4)] for s in (300, 0, -300)]
    yield ("cosmo_blockdiag.3x4_e300", blockdiag(bl), prod([det(b) for b in bl]))
    d = spread(12, 300, -300, 10); U = imat(12, 2)
    yield ("cosmo_lowrank.12_r2_e300", lowrank_diag(d, U), det_lowrank(d, U))
    J = imat(12, 6); cd = spread(12, 300, -300, 10)
    Fm = [[sum((J[i][a] * J[i][b] / cd[i] for i in range(12)), F(0)) for b in range(6)] for a in range(6)]
    yield ("cosmo_fisher.06_from12_e300", Fm, None)
    yield ("cosmo_refuse.psd_lowrank6", lowrank_diag([F(0)] * 6, imat(6, 2)), F(0))
    # heavy tier
    d = spread(32, 100000, -100000, 2); U = imat(32, 2)
    yield ("cosmo_spread.p2e100000_lowrank32", lowrank_diag(d, U), det_lowrank(d, U))
    d = spread(64, 1000000, -1000000, 2); U = imat(64, 1)
    yield ("cosmo_spread.p2e1000000_lowrank64", lowrank_diag(d, U), det_lowrank(d, U))
    bl = [[[P2(s) * a for a in r] for r in spd_small(8)] for s in (3000, 2000, 1000, 0, -1000, -2000, -3000, -4000)]
    yield ("cosmo_blockdiag.8x8_p2e3000_64", blockdiag(bl), prod([det(b) for b in bl]))
    rho = 1 - P10(-50)
    yield ("cosmo_toeplitz.kms32_rho1m1e-50", kms(32, rho), (1 - rho * rho) ** 31)
    Dc = [F(1)] * 31 + [P10(-200)]
    yield ("cosmo_cond.e200_32", ldlt(unit_lower(32, band=1), Dc), prod(Dc))

def gls_cosmo_corpus():
    """Yields (cid, X, S, y, exact S^-1 solver or None, S is SPD by construction) cosmology-scale GLS cases."""
    def bao(m, p): return [[F(i + 1, m) ** k for k in range(p)] for i in range(m)]
    def ys(m): return [rq() for _ in range(m)]
    yield ("cosmo_bao.12x3_e300", bao(12, 3), ldlt(unit_lower(12, band=2), spread(12, 300, -300, 10)), ys(12), None, True)
    dk = spread(12, 150, -150, 10)
    yield ("cosmo_cmb_toeplitz.12x2_e300", bao(12, 2), scale_sym(kms(12, F(9, 10)), dk), ys(12), None, True)
    bl = [[[P10(s) * a for a in r] for r in spd_small(4)] for s in (300, 0, -300)]
    yield ("cosmo_blockdiag.12x4_e300", bao(12, 4), blockdiag(bl), ys(12), block_solver(bl), True)
    d = spread(12, 300, -300, 10); U = imat(12, 2)
    yield ("cosmo_lowrank.12x3_e300", bao(12, 3), lowrank_diag(d, U), ys(12), woodbury(d, U), True)
    yield ("cosmo_cond.e200_10x2", bao(10, 2), ldlt(unit_lower(10, band=1), [F(1)] * 9 + [P10(-200)]), ys(10), None, True)
    d = spread(6, 10000, -10000, 2); U = imat(6, 1)
    yield ("cosmo_spread.p2e10000_6x2", bao(6, 2), lowrank_diag(d, U), ys(6), woodbury(d, U), True)
    yield ("cosmo_refuse.sigma_lowrank_psd", bao(6, 2), lowrank_diag([F(0)] * 6, imat(6, 2)), ys(6), None, False)
    yield ("cosmo_refuse.sigma_indef_e300", bao(5, 2), ldlt(unit_lower(5), [P10(300), -P10(-300), F(1), F(1), F(1)]), ys(5), None, False)
    # heavy tier
    d = spread(32, 100000, -100000, 2); U = imat(32, 2)
    yield ("cosmo_lowrank.32x4_p2e100000", bao(32, 4), lowrank_diag(d, U), ys(32), woodbury(d, U), True)
    d = spread(32, 1000000, -1000000, 2); U = imat(32, 1)
    yield ("cosmo_spread.p2e1000000_32x2", bao(32, 2), lowrank_diag(d, U), ys(32), woodbury(d, U), True)
    bl = [[[P2(s) * a for a in r] for r in spd_small(8)] for s in (3000, 2000, 1000, 0, -1000, -2000, -3000, -4000)]
    yield ("cosmo_blockdiag.64x3_p2e3000", bao(64, 3), blockdiag(bl), ys(64), block_solver(bl), True)
    # small gate rows (appended last so earlier rows keep their rng draws): one per family whose big 12-dim / 2^+-10^4 rows
    # moved to the heavy tier (GLS_COSMO_HEAVY) to keep the LA gate set under 1 MB
    yield ("cosmo_bao.6x2_e30", bao(6, 2), ldlt(unit_lower(6, band=2), spread(6, 30, -30, 10)), ys(6), None, True)
    yield ("cosmo_cmb_toeplitz.6x2_e30", bao(6, 2), scale_sym(kms(6, F(9, 10)), spread(6, 15, -15, 10)), ys(6), None, True)
    bl = [[[P10(s) * a for a in r] for r in spd_small(2)] for s in (30, 0, -30)]
    yield ("cosmo_blockdiag.6x2_e30", bao(6, 2), blockdiag(bl), ys(6), block_solver(bl), True)
    d = spread(6, 30, -30, 10); U = imat(6, 2)
    yield ("cosmo_lowrank.6x2_e30", bao(6, 2), lowrank_diag(d, U), ys(6), woodbury(d, U), True)
    d = spread(4, 1000, -1000, 2); U = imat(4, 1)
    yield ("cosmo_spread.p2e1000_4x2", bao(4, 2), lowrank_diag(d, U), ys(4), woodbury(d, U), True)

# Gate-size budget (Spence 2026-10-24: LA gate set < 1 MB): these cosmo GLS rows are small enough for the gate's dims/bits rule but
# are emitted on the heavy tier; the small *_e30 / p2e1000_4x2 rows above keep each family on the gate.
GLS_COSMO_HEAVY = frozenset({"cosmo_bao.12x3_e300", "cosmo_cmb_toeplitz.12x2_e300", "cosmo_blockdiag.12x4_e300",
                             "cosmo_lowrank.12x3_e300", "cosmo_spread.p2e10000_6x2"})

# ---------------------------------------------------------------- per-op emitters
def dims_sq(A): return str(len(A))
def dims_mn(A): return f"{len(A)}x{len(A[0])}"

def row(cid, dims, A, outs, heavy=None):
    return (" ".join([cid, dims] + [q(x) for x in flat(A)] + outs), is_heavy(A) if heavy is None else heavy)

def gen_lu():
    rows = Rows()
    for cid, A in square_corpus("lu") + [("refuse.rect_2x3", [[F(1), F(2), F(3)], [F(4), F(5), F(6)]]), ("refuse.rect_3x1", [[F(1)], [F(0)], [F(-2)]])]:
        begin("lu." + cid)
        if len(A) != len(A[0]):
            rows.append(row("lu." + cid, dims_mn(A), A, err("DimensionMismatch", "rectangular"))); continue
        perm, L, U = lu_pp(A)
        rows.append(row("lu." + cid, dims_sq(A), A, [str(p) for p in perm] + [q(x) for x in flat(L)] + [q(x) for x in flat(U)]))
    hdr = ["op: lu",
           "canonical: PA = LU, partial pivoting. At step k pivot row p = argmax_{i>=k} |U[i][k]|, ties -> smallest i; rows k,p of U and of L[:, :k] swapped. If the pivot column is all zero, no swap and no elimination (L column k stays e_k). L unit lower triangular, U upper triangular (zero/singular rows allowed). Exact rationals.",
           "perm: n integers, 0-based; perm[i] = index of the row of A placed at row i of PA",
           "refusals: non-square A -> 'ERR DimensionMismatch rectangular' (rectangular LU is out of scope); singular A is NOT refused (U has zero pivots).",
           "columns: case_id dims(n | mxn) A[m*n] (perm[n] L[n*n] U[n*n] | ERR DimensionMismatch rectangular)"]
    expect = lambda m, n: (m * n, {n + 2 * n * n})
    return hdr, rows, expect

def gen_solve():
    rows = Rows()
    for cid, A in square_corpus("solve"):
        begin("solve." + cid)
        n = len(A); b = [F(1)] * n if cid.startswith("hilbert") else [F(i + 1, 2) - 1 for i in range(n)]
        x = solve(A, b)
        outs = err("Singular", "singular") if x is None else [q(v) for v in x]
        rows.append((" ".join(["solve." + cid, str(n)] + [q(v) for v in flat(A)] + [q(v) for v in b] + outs), is_heavy(A, b)))
    hdr = ["op: solve",
           "canonical: unique x with A x = b, exact rationals. b = ones for Hilbert cases, b[i] = (i+1)/2 - 1 otherwise. Any singular A (det = 0) -> refusal 'ERR Singular singular'.",
           "columns: case_id dims(n) A[n*n] b[n] (x[n] | ERR Singular singular)"]
    expect = lambda m, n: (n * n + n, {n})
    return hdr, rows, expect

def gen_logdet():
    rows = Rows(); ext = {b: Rows() for b in LOGDET_EXT_BITS}
    for cid, A in square_corpus("logdet"):
        begin("logdet." + cid)
        d = det(A); heavy = is_heavy(A)
        if d == 0: outs = err("Singular", "singular")
        else: outs = [str(1 if d > 0 else -1), q(d), pair(*enclose_log(d, BITS))]
        rows.append(row("logdet." + cid, dims_sq(A), A, outs, heavy))
        if not heavy:
            for b in LOGDET_EXT_BITS:
                ext[b].append(row("logdet." + cid, dims_sq(A), A, outs if d == 0 else [pair(*enclose_log(d, b))], True)[0])
    for cid, A, dk in logdet_cosmo_corpus():
        begin("logdet." + cid)
        heavy = is_heavy(A)
        if heavy and SKIP_HEAVY: continue                  # --skip-heavy: heavy cosmo rows are never computed
        if dk is None or not heavy:
            d = det(A); assert dk is None or d == dk, cid   # gate: identity re-checked against the LU determinant
        else: d = dk                                       # heavy: exact identity (prod D, Sylvester, block product, KMS closed form)
        if d == 0: outs = err("Singular", "singular")
        else: outs = [str(1 if d > 0 else -1), q(d), pair(*enclose_log(d, BITS))]
        rows.append(row("logdet." + cid, dims_sq(A), A, outs, heavy))
        for b in LOGDET_EXT_BITS:
            ext[b].append(row("logdet." + cid, dims_sq(A), A, outs if d == 0 else [pair(*enclose_log(d, b))], True)[0])
    hdr = ["op: logdet",
           "canonical: sign = sign(det A) in {1,-1}; det = exact rational determinant; logabsdet = ln|det A| (irrational unless det = +-1). det = 0 -> refusal 'ERR Singular singular'.",
           "exact form: ln(det) with det the exact rational in the row (sympy: log(Rational(det))).",
           "enclosure: lo hi pair enclosing ln|det|; at bits b the grid is 2^e with e = floor(log2 max(|lo_raw|,|hi_raw|)) - b + 1; computed with mpmath.iv.log at b+64 and checked to contain an independent mpmath.iv enclosure at 2*b+64. Gate and heavy files: b = 256 for every row. Extra widths " + str(LOGDET_EXT_BITS) + " for the gate-tier cases and every cosmo_* case (same case_ids) are heavy-tier files la_logdet_enc<b>_heavy.txt.",
           "cosmo_* rows (cosmology scales): SPD covariances with eigenvalue-scale spread 10^+-300 and 2^+-(10^3..10^6) in one matrix, near-singular SPD with cond ~10^50..10^200, Toeplitz (KMS rho^|i-j|, tridiagonal), block-diagonal, diag + low-rank, Fisher J^T C^-1 J. Heavy rows (n = 32/64, 2^+-10^5..10^6) take det from an exact identity (prod D of L D L^T, Sylvester det(D + U U^T) = det D det(I + U^T D^-1 U), block product, KMS (1-rho^2)^(n-1)); gate rows re-check it against the LU determinant.",
           "columns (gate): case_id dims(n) A[n*n] (sign det lo256 hi256 | ERR Singular singular)",
           "columns (heavy): case_id dims(n) A[n*n] (sign det lo256 hi256 | ERR Singular singular)"]
    expect = lambda m, n: (n * n, {4})
    exts = {f"enc{b}": (hdr[:3] + [f"enc{b}: the b = {b} ln|det| enclosure of every gate-tier logdet case and every cosmo_* case (same case_ids, same order as la_logdet.txt then la_logdet_heavy.txt); sign, det and the 256-bit pair are in those files.",
                               f"columns (heavy): case_id dims(n) A[n*n] (lo{b} hi{b} | ERR Singular singular)"],
                        ext[b], lambda m, n: (n * n, {2})) for b in LOGDET_EXT_BITS}
    return hdr, rows, expect, exts

def gen_cholesky():
    rows = Rows(); ext = {b: Rows() for b in EXT_BITS}
    for cid, A in spd_corpus():
        begin("cholesky." + cid)
        r = ldl(A); n = len(A); heavy = is_heavy(A)
        if isinstance(r, str): outs = ldl_err(r, "a")
        else:
            L, D = r; outs = [q(x) for x in flat(L)] + [q(x) for x in D]
            def gblock(b): return [pair(*enclose_csqrt(L[i][j], D[j], b)) if j <= i else "0*2^0 0*2^0" for i in range(n) for j in range(n)]
            outs += gblock(BITS)
        rows.append(row("cholesky." + cid, dims_sq(A), A, outs, heavy))
        if not heavy:
            for b in EXT_BITS:
                ext[b].append(row("cholesky." + cid, dims_sq(A), A, ldl_err(r, "a") if isinstance(r, str) else gblock(b), True)[0])
    hdr = ["op: cholesky",
           "canonical: A = G G^T with G lower triangular and diag(G) > 0 (unique for SPD A). Exact rational companion: A = L D L^T, L unit lower triangular, D diagonal > 0, computed without pivoting; G[i][j] = L[i][j]*sqrt(D[j]) (this is the exact algebraic form).",
           "refusals: 'ERR NotSymmetric a' if A != A^T; 'ERR NotPositiveDefinite a' if some pivot D[j] <= 0 (indefinite, negative, singular PSD, zero).",
           "enclosure: G[n*n] as lo hi pairs, row-major, upper triangle emitted as exact 0*2^0 0*2^0; at bits b the grid is 2^e with e = floor(log2(G^2))//2 - b + 1, lo/hi from exact integer isqrt; lo == hi iff exactly representable. Gate and heavy files: b = 256 for every row. Extra widths " + str(EXT_BITS) + " for the gate-tier cases (same case_ids) are heavy-tier files la_cholesky_enc<b>_heavy.txt.",
           "columns (gate): case_id dims(n) A[n*n] (L[n*n] D[n] G_lohi256[2*n*n] | ERR variant detail)",
           "columns (heavy): case_id dims(n) A[n*n] (L[n*n] D[n] G_lohi256[2*n*n] | ERR variant detail)"]
    expect = lambda m, n: (n * n, {n * n + n + 2 * n * n})
    exts = {f"enc{b}": (hdr[:3] + [f"enc{b}: the b = {b} G enclosure block of every gate-tier cholesky case (same case_ids, same order as la_cholesky.txt); L, D and the 256-bit block are in la_cholesky.txt.",
                               f"columns (heavy): case_id dims(n) A[n*n] (G_lohi{b}[2*n*n] | ERR variant detail)"],
                        ext[b], lambda m, n: (n * n, {2 * n * n})) for b in EXT_BITS}
    return hdr, rows, expect, exts

def gen_qr():
    rows = Rows(); ext = {b: Rows() for b in EXT_BITS}
    C = [("onebyone.pos", [[F(3)]]), ("onebyone.neg", [[F(-2, 3)]]), ("basic.2", [[F(3), F(1)], [F(4), F(2)]]),
         ("pythag.3", [[F(2), F(1), F(0)], [F(1), F(2), F(1)], [F(2), F(2), F(3)]]),
         ("tall.3x2", [[F(1), F(1)], [F(1), F(2)], [F(1), F(3)]]), ("tall.6x3", rmat(6, 3)), ("tall.20x8", rmat(20, 8, rational=False)),
         ("nearsing.2", [[F(1), F(1)], [F(1), 1 + EPS]])]
    C += [("hilbert.%02d" % n, hilbert(n)) for n in range(5, 13)]
    C += [("scale_huge.1x1", [[H6]]), ("scale_mixed.diag2", [[H6, F(0)], [F(0), T6]]), ("scale_tiny.3", [[a * T3 for a in r] for r in rmat(3, 3)]),
          ("scale_huge.2", [[B4 + 1, B4], [B4, B4 - 3]])]
    for n in (3, 4, 6): C.append(("random_q.%02d" % n, rmat(n, n)))
    for n in (8, 12, 16): C.append(("random_z.%02d" % n, rmat(n, n, rational=False)))
    C += [("refuse.zero2", zero(2)), ("refuse.onebyone_zero", [[F(0)]]), ("refuse.rankdef3", [[F(1), F(2), F(3)], [F(4), F(5), F(6)], [F(7), F(8), F(9)]]),
          ("refuse.dupcol_4x2", [[F(i), F(i)] for i in range(1, 5)]), ("refuse.wide2x3", [[F(1), F(0), F(0)], [F(0), F(1), F(0)]])]
    for cid, A in C:
        begin(OP + "." + cid)
        m, n = len(A), len(A[0])
        g = None if m < n else gram_schmidt(A); heavy = is_heavy(A)
        if g is None: outs = err("DimensionMismatch", "wide") if m < n else err("Singular", "rank_deficient")
        else:
            Qt, Rt, d = g
            outs = [q(x) for x in flat(Qt)] + [q(x) for x in flat(Rt)] + [q(x) for x in d]
            def qrblock(b): return ([pair(*enclose_csqrt(Qt[i][j] / d[j], d[j], b)) for i in range(m) for j in range(n)] +
                                    [pair(*enclose_csqrt(Rt[i][j], d[i], b)) for i in range(n) for j in range(n)])
            outs += qrblock(BITS)
        rows.append(row("qr." + cid, dims_mn(A), A, outs, heavy))
        if not heavy:
            for b in EXT_BITS:
                ext[b].append(row("qr." + cid, dims_mn(A), A, outs[-3:] if g is None else qrblock(b), True)[0])
    hdr = ["op: qr",
           "canonical: thin QR, A (m x n, m >= n, full column rank) = Q R, Q m x n with Q^T Q = I, R n x n upper triangular with diag(R) > 0 (unique).",
           "exact rational companion (exact algebraic form): A = Qt Rt with Qt orthogonal columns (Qt^T Qt = diag(d)), Rt unit upper triangular, d[j] = ||Qt[:,j]||^2; Q[i][j] = Qt[i][j]/sqrt(d[j]), R[i][j] = Rt[i][j]*sqrt(d[i]).",
           "refusals: 'ERR Singular rank_deficient' if some d[j] = 0; 'ERR DimensionMismatch wide' if m < n.",
           "enclosure: per bit width b, Q (m*n) then R (n*n) as lo hi pairs, row-major; grid 2^e with e = floor(log2(x^2))//2 - b + 1, exact integer isqrt. Gate and heavy files: b = 256 for every row. Extra widths " + str(EXT_BITS) + " for the gate-tier cases (same case_ids) are heavy-tier files la_qr_enc<b>_heavy.txt.",
           "columns (gate): case_id dims(mxn) A[m*n] (Qt[m*n] Rt[n*n] d[n] Q_lohi256[2*m*n] R_lohi256[2*n*n] | ERR variant detail)",
           "columns (heavy): case_id dims(mxn) A[m*n] (Qt[m*n] Rt[n*n] d[n] Q_lohi256[2*m*n] R_lohi256[2*n*n] | ERR variant detail)"]
    expect = lambda m, n: (m * n, {m * n + n * n + n + 2 * (m * n + n * n)})
    exts = {f"enc{b}": (hdr[:4] + [f"enc{b}: the b = {b} Q,R enclosure block of every gate-tier qr case (same case_ids, same order as la_qr.txt); Qt, Rt, d and the 256-bit block are in la_qr.txt.",
                               f"columns (heavy): case_id dims(mxn) A[m*n] (Q_lohi{b}[2*m*n] R_lohi{b}[2*n*n] | ERR variant detail)"],
                        ext[b], lambda m, n: (m * n, {2 * (m * n + n * n)})) for b in EXT_BITS}
    return hdr, rows, expect, exts

def gen_gls():
    rows = Rows()
    def ar1(m, rho): return [[rho ** abs(i - j) for j in range(m)] for i in range(m)]
    def diag(v): return [[v[i] if i == j else F(0) for j in range(len(v))] for i in range(len(v))]
    def design(m, p): return [[F(i + 1) ** k for k in range(p)] for i in range(m)]
    C = [("basic.ols_identity", design(4, 2), ident(4), [F(1), F(3), F(2), F(5)]),
         ("basic.wls_diag", design(5, 2), diag([F(1), F(2), F(1, 2), F(4), F(1)]), [F(2), F(1), F(4), F(3), F(6)]),
         ("ar1.6x2", design(6, 2), ar1(6, F(1, 2)), [rq() for _ in range(6)]),
         ("ar1.10x3", design(10, 3), ar1(10, F(9, 10)), [rq() for _ in range(10)]),
         ("exactfit.3x3", design(3, 3), spd(3), [F(1), F(4), F(9)]),
         ("onebyone.1x1", [[F(2)]], [[F(3)]], [F(5)]),
         ("mean.5x1", [[F(1)]] * 5, spd(5), [rq() for _ in range(5)]),
         ("random.8x3", rmat(8, 3), spd(8), [rq() for _ in range(8)]),
         ("random.16x4", rmat(16, 4, rational=False), spd(16, rational=False), [F(rng.randint(-9, 9)) for _ in range(16)]),
         ("hilbert_cov.6x2", design(6, 2), hilbert(6), [F(1)] * 6),
         ("scale_mixed.4x2", design(4, 2), diag([B3, T3, F(1), B3]), [F(1), F(2), F(3), F(4)]),
         ("refuse.sigma_not_spd", design(3, 2), [[F(1), F(2), F(0)], [F(2), F(1), F(0)], [F(0), F(0), F(1)]], [F(1), F(2), F(3)]),
         ("refuse.sigma_singular", design(3, 2), zero(3), [F(1), F(2), F(3)]),
         ("refuse.rankdef_X", [[F(1), F(2)], [F(2), F(4)], [F(3), F(6)]], ident(3), [F(1), F(2), F(3)])]
    for cid, Xm, S, y in C:
        begin("gls." + cid)
        m, p = len(Xm), len(Xm[0])
        r = ldl(S)
        if isinstance(r, str): outs = ldl_err(r, "sigma")
        else:
            SiX = tr([solve(S, col) for col in tr(Xm)])
            Siy = solve(S, y)
            Mx = mm(tr(Xm), SiX)
            if det(Mx) == 0: outs = err("Singular", "rank_deficient")
            else:
                beta = solve(Mx, [sum((Xm[i][k] * Siy[i] for i in range(m)), F(0)) for k in range(p)])
                cov = tr([solve(Mx, [F(int(i == j)) for i in range(p)]) for j in range(p)])
                res = [y[i] - sum((Xm[i][k] * beta[k] for k in range(p)), F(0)) for i in range(m)]
                rs = solve(S, res)
                rss = sum((a * b for a, b in zip(res, rs)), F(0))
                # normal equations check
                assert all(sum((Xm[i][k] * rs[i] for i in range(m)), F(0)) == 0 for k in range(p))
                outs = [q(v) for v in beta] + [q(v) for v in flat(cov)] + [q(rss)]
        rows.append((" ".join(["gls." + cid, f"{m}x{p}"] + [q(v) for v in flat(Xm)] + [q(v) for v in flat(S)] + [q(v) for v in y] + outs), is_heavy(Xm, S, y)))
    for cid, Xm, S, y, sol, spd_known in gls_cosmo_corpus():
        begin("gls." + cid)
        m, p = len(Xm), len(Xm[0]); big = is_heavy(Xm, S, y); heavy = big or cid in GLS_COSMO_HEAVY
        if heavy and SKIP_HEAVY: continue                  # --skip-heavy: heavy cosmo rows are never computed
        r = None if big and spd_known else ldl(S)         # heavy cosmo covariances are SPD by construction (D > 0 plus U U^T, blocks B B^T + I)
        if isinstance(r, str):
            assert not spd_known, cid; outs = ldl_err(r, "sigma")
        else:
            assert spd_known, cid
            gsol = lambda v: solve(S, v)
            def Si(v, sol=sol):
                x = (sol or gsol)(v)
                assert matvec(S, x) == [F(a) for a in v]
                if sol is not None and not big: assert x == gsol(v)
                return x
            SiX = tr([Si(col) for col in tr(Xm)]); Siy = Si(y)
            Mx = mm(tr(Xm), SiX)
            if det(Mx) == 0: outs = err("Singular", "rank_deficient")
            else:
                beta = solve(Mx, [sum((Xm[i][k] * Siy[i] for i in range(m)), F(0)) for k in range(p)])
                cov = tr([solve(Mx, [F(int(i == j)) for i in range(p)]) for j in range(p)])
                res = [y[i] - sum((Xm[i][k] * beta[k] for k in range(p)), F(0)) for i in range(m)]
                rs = Si(res); rss = sum((a * b for a, b in zip(res, rs)), F(0))
                assert all(sum((Xm[i][k] * rs[i] for i in range(m)), F(0)) == 0 for k in range(p))
                assert rss >= 0
                outs = [q(v) for v in beta] + [q(v) for v in flat(cov)] + [q(rss)]
        rows.append((" ".join(["gls." + cid, f"{m}x{p}"] + [q(v) for v in flat(Xm)] + [q(v) for v in flat(S)] + [q(v) for v in y] + outs), heavy))
    hdr = ["op: gls",
           "canonical: generalized least squares, beta = (X^T S^-1 X)^-1 X^T S^-1 y; cov = (X^T S^-1 X)^-1; rss = r^T S^-1 r with r = y - X beta. All exact rationals. Checked: X^T S^-1 r = 0.",
           "refusals: 'ERR NotPositiveDefinite sigma' / 'ERR NotSymmetric sigma' (the covariance S fails LDL^T, see la_cholesky); 'ERR Singular rank_deficient' if X^T S^-1 X is singular.",
           "cosmo_* rows (cosmology scales): BAO-style designs X[i][k] = ((i+1)/m)^k with covariances of eigenvalue-scale spread 10^+-300 and 2^+-(10^4..10^6) (L D L^T, scaled KMS Toeplitz, block-diagonal, diag + low-rank), cond ~10^200, and non-SPD refusals. S^-1 v is exact (Woodbury for diag + low-rank, blockwise for block-diagonal, LU otherwise), checked S x = v; for every row with dims <= 12 and inputs <= 10064 bits also equal to the plain LU solve. Tier: rows over that size, plus the 12-dim e300 rows and cosmo_spread.p2e10000_6x2 (gate-size budget, < 1 MB for the LA gate set), are in la_gls_heavy.txt; each family keeps a small gate row (*_e30, cosmo_spread.p2e1000_4x2, cosmo_cond.e200_10x2, cosmo_refuse.*). rss = chi^2 >= 0 checked.",
           "columns: case_id dims(mxp) X[m*p] S[m*m] y[m] (beta[p] cov[p*p] rss | ERR variant detail)"]
    expect = lambda m, n: (m * n + m * m + m, {n + n * n + 1})
    return hdr, rows, expect

def gen_eig():
    rows = Rows()
    C = [("onebyone.pos", [[F(7, 3)]]), ("zero.03", zero(3)), ("basic.2", [[F(2), F(1)], [F(1), F(2)]]),
         ("irr.2", [[F(1), F(1)], [F(1), F(0)]]), ("tridiag.5", [[F(2) if i == j else (F(-1) if abs(i - j) == 1 else F(0)) for j in range(5)] for i in range(5)]),
         ("repeated.4", [[F(1)] * 4 for _ in range(4)]), ("rankdef.3", [[F(1), F(2), F(3)], [F(2), F(4), F(6)], [F(3), F(6), F(9)]]),
         ("nearsing.2", [[F(1), F(1)], [F(1), 1 + EPS]]), ("diag_dup.3", [[F(2), F(0), F(0)], [F(0), F(-1, 3), F(0)], [F(0), F(0), F(2)]])]
    C += [("hilbert.%02d" % n, hilbert(n)) for n in range(5, 13)]
    C += [("scale_huge.1x1", [[H6]]), ("scale_mixed.diag2", [[H6, F(0)], [F(0), T6]]),
          ("scale_huge.2", [[B3, F(1)], [F(1), B3]]), ("scale_tiny.2", [[T3, T3 * T3], [T3 * T3, 2 * T3]])]
    for n in (3, 4, 6, 8, 12, 16): C.append(("random_z.%02d" % n, rsym(n)))
    C += [("refuse.nonsym2", [[F(0), F(1)], [F(-1), F(0)]])]
    for cid, A in C:
        begin(OP + "." + cid)
        n = len(A)
        if not issym(A): outs = err("NotSymmetric", "a")
        else:
            coeffs, ev = eig_exact(A)
            outs = [q(c) for c in coeffs] + [f"{e[1]} {pair(e[2], e[3])}" for e in ev]
        rows.append(row("eig." + cid, dims_sq(A), A, outs))
    hdr = ["op: eig",
           "canonical: real symmetric A only; eigenvalues listed ascending with algebraic multiplicity. charpoly = det(xI - A), monic, n+1 exact rational coefficients, highest degree first.",
           "exact form: each eigenvalue is the unique root of the charpoly in [lo, hi] at its ascending index (sympy CRootOf of the irreducible factor); token 'p/q' if the eigenvalue is rational, else 'alg'.",
           "enclosure: lo hi per eigenvalue; grid 2^e, e = floor(log2 |lambda|) - 255; from sympy exact rational root isolation/refinement (Poly.intervals/refine_root, exact QQ) then outward rounding; checked with exact Poly.count_roots(lo, hi) == 1 on the irreducible factor. Zero/dyadic rational eigenvalues have lo == hi.",
           "scope: slice LA5a (symmetric). General nonsymmetric eigenvalues are LA5b, la_eig_general.txt in this slice (status: pending (needs complex Ball, L2-1); header only, zero rows); here a nonsymmetric A is refused.",
           "refusals: 'ERR NotSymmetric a'.",
           "columns: case_id dims(n) A[n*n] (charpoly[n+1] {exact lo hi}[n] | ERR NotSymmetric a)"]
    expect = lambda m, n: (n * n, {n + 1 + 3 * n})
    return hdr, rows, expect

def gen_eig_general():
    """LA5b stub: header only, zero rows, until zb has a complex Ball (L2-1)."""
    hdr = ["op: eig_general",
           "status: pending (needs complex Ball, L2-1). This file is a stub: header only, zero rows. Rows land once zb has complex Ball enclosures.",
           "canonical (planned): general real (nonsymmetric) A; complex eigenvalues with algebraic multiplicity; charpoly = det(xI - A), monic, exact rational coefficients, highest degree first.",
           "scope: slice LA5b (general). Symmetric A is covered by la_eig.txt (LA5a)."]
    return hdr, Rows(), lambda m, n: (n * n, set()), {}, "stub"

OPS = [("lu", gen_lu), ("cholesky", gen_cholesky), ("qr", gen_qr), ("solve", gen_solve), ("gls", gen_gls), ("logdet", gen_logdet), ("eig", gen_eig), ("eig_general", gen_eig_general)]

# ---------------------------------------------------------------- streaming output (Paul 2026-10-24: OOM at 58 GB RSS)
# Rows never accumulate in memory: each row is written to an anonymous spool file on disk the moment it is produced (gate and
# heavy apart), files are streamed spool -> destination with an incremental sha256, --check hashes and diffs files chunk by chunk.
# Only case_ids (short strings) are kept, for the duplicate and ext-order checks. Run it under a memory cap (see the docstring).
SKIP_HEAVY = False
SPOOL_DIR = os.path.join(ROOT, "target", "la-goldens-spool")   # /target is gitignored; spool files are unlinked at creation

# Heavy-run progress contract (Jeff, 2026-10-24, after la4-heavy2 sat 90 min at 100% CPU with zero output; selftest
# tools/gen_la_progress_selftest.sh). Whenever heavy rows are computed (not --skip-heavy): each row prints
#   progress: op=<op> case=<case_id> tier=<gate|heavy> row_s=<s> elapsed_s=<s>
# to stderr after it is on disk; with --heavy DIR heavy rows also go to DIR/<slice>/la_<op>.rows<k>.partial, flushed and
# fsynced per row (removed once the final file is written; left behind if the run dies). A per-case wall budget
# (--case-budget, SIGALRM, re-armed at every case start and row) exits 5 naming the op and the case_id being computed.
OP = ""; CASE = ""; PROGRESS = False; PARTIAL_DIR = None; BUDGET = 0.0; T0 = time.monotonic(); T_ROW = T0; N_ROWS_OBJ = 0

class CaseBudget(Exception):
    pass

def _on_alarm(*_):
    raise CaseBudget()

def _arm():
    if BUDGET > 0: signal.setitimer(signal.ITIMER_REAL, BUDGET)

def begin(cid):
    """Marks the case_id being computed (named by the budget failure) and re-arms the per-case budget."""
    global CASE
    CASE = cid; _arm()

class Rows:
    """Row sink: append((text, heavy)) or append(text) (= heavy). Text goes straight to disk."""
    def __init__(self):
        global N_ROWS_OBJ
        os.makedirs(SPOOL_DIR, exist_ok=True)
        heavy_spool, self.partial = tempfile.TemporaryFile(dir=SPOOL_DIR), None
        if PARTIAL_DIR is not None:
            os.makedirs(PARTIAL_DIR, exist_ok=True)
            self.partial = os.path.join(PARTIAL_DIR, f"la_{OP}.rows{N_ROWS_OBJ}.partial"); heavy_spool = open(self.partial, "w+b")
        N_ROWS_OBJ += 1
        self.spool = {False: tempfile.TemporaryFile(dir=SPOOL_DIR), True: heavy_spool}
        self.count = {False: 0, True: 0}; self.ids = []
    def append(self, item):
        global T_ROW
        text, heavy = item if isinstance(item, tuple) else (item, True)
        heavy = bool(heavy); cid = text.split(" ", 1)[0]
        f = self.spool[heavy]; f.write(text.encode("utf-8") + b"\n"); self.count[heavy] += 1
        if heavy and self.partial is not None: f.flush(); os.fsync(f.fileno())
        self.ids.append((cid, heavy))
        if PROGRESS:
            now = time.monotonic()
            print(f"progress: op={OP} case={cid} tier={'heavy' if heavy else 'gate'} row_s={now - T_ROW:.3f} elapsed_s={now - T0:.3f}",
                  file=sys.stderr, flush=True)
            T_ROW = now
        _arm()
    def __bool__(self): return bool(self.ids)
    def lines(self, heavy):
        f = self.spool[heavy]; f.flush(); f.seek(0)
        for l in f: yield l
    def close(self):
        for f in self.spool.values(): f.close()
        if self.partial is not None: os.remove(self.partial)

def header_bytes(name, tier, hdr, nrows):
    h = [f"generator: {GEN}", f"sympy: {SYMPY_PIN}", f"mpmath: {MPMATH_PIN}", f"seed: {SEED} (python random.Random('{SEED}:{name}'))",
         f"slice: {SLICE_ID_OP.get(name, SLICE_ID[SLICE_OF[name]])} ({SLICE_OF[name]})", f"tier: {tier} (gate: all dims <= {GATE_MAX_DIM} and input integers <= {GATE_MAX_BITS} bits; heavy: the rest, checked by sha256 at merge points, file not committed)",
         f"rows: {nrows}"] + [l for l in hdr if not l.startswith("columns (") or l.startswith(f"columns ({tier})")] + [
         "numbers: every exact entry is p/q in lowest terms, q > 0, arbitrary-size signed integers; dyadics are m*2^e with m odd (or 0*2^0); fields are single-space separated; no floats.",
         "case_id: <op>.<family>.<variant>, stable across regenerations; inputs are always listed before outputs; dims gives every matrix size, so a row's field count is fixed by dims; refusal rows end in 'ERR <variant> <detail>' (variant = zb BoundError variant)."]
    return "".join(f"# {l}\n" for l in h).encode("utf-8")

def stream_file(name, tier, hdr, rows, heavy, expect, sink):
    """Header + spooled rows -> sink (a binary file or None); every row's field count is checked on the way. Returns (sha256, bytes)."""
    hs = hashlib.sha256(); n = 0
    def put(b):
        nonlocal n
        hs.update(b); n += len(b)
        if sink is not None: sink.write(b)
    put(header_bytes(name, tier, hdr, rows.count[heavy]))
    for l in rows.lines(heavy):
        check_fields(l[:-1].decode("utf-8"), expect); put(l)
    return hs.hexdigest(), n

def sha_path(p):
    hs = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""): hs.update(c)
    return hs.hexdigest()

def first_diff_line(a, b):
    with open(a, "rb") as fa, open(b, "rb") as fb:
        n = 0
        for la, lb in zip(fa, fb):
            n += 1
            if la != lb: return n
        return n + 1

def sums_header(slice_dir):
    p = os.path.join(slice_dir, "SHA256SUMS")
    if os.path.exists(p):
        return "".join(l for l in open(p, encoding="utf-8") if l.startswith("#"))
    return ""

def read_sums(p):
    return {l.split()[1]: l.split()[0] for l in open(p) if l.strip() and not l.startswith("#")} if os.path.exists(p) else {}

PIN_LINE = re.compile(r"[0-9a-f]{64}  [^/\s]+")

def check_pins(sl, d, names):
    """Regular-gate heavy check, no recompute: HEAVY.SHA256 lines are strict and name exactly the defined heavy files. Returns the failure count."""
    p = os.path.join(d, "HEAVY.SHA256"); bad = 0
    lines = [l.rstrip("\n") for l in open(p, encoding="utf-8") if l.strip() and not l.startswith("#")] if os.path.exists(p) else []
    for l in lines:
        if not PIN_LINE.fullmatch(l): bad += 1; print(f"check: HEAVY.SHA256 malformed line {sl}: {l}")
    pins = {l.split()[1]: l.split()[0] for l in lines if PIN_LINE.fullmatch(l)}
    for k in sorted(names):
        if k in pins: print(f"check: OK pin {sl}/{k} {pins[k][:16]}")
        else: bad += 1; print(f"check: HEAVY.SHA256 pin MISSING {sl}/{k}")
    for k in sorted(set(pins) - set(names)): bad += 1; print(f"check: HEAVY.SHA256 pin STALE {sl}/{k}")
    return bad

HEAVY_HDR = ("# sha256 of heavy-tier LA golds (not committed; rebuilt deterministically by scripts/gen_la_goldens.py --heavy DIR,\n"
             "# verified by scripts/gen_la_goldens.py --check, step 'heavy' at merge points). sha256sum format.\n")

def render(gate_dir, heavy_dir, skip_heavy):
    """Streams every op: gate files to gate_dir/<slice>/ (written as .tmp then renamed), heavy files to heavy_dir/<slice>/ if
    heavy_dir else hash only. Returns {slice: {"gate": {fname: (sha, size)}, "heavy": {fname: (sha, size)}}}."""
    global rng, SKIP_HEAVY, OP, CASE, PROGRESS, PARTIAL_DIR, N_ROWS_OBJ
    SKIP_HEAVY = skip_heavy; PROGRESS = not skip_heavy
    out = {}
    def to(d, sl, k):
        if d is None: return None
        os.makedirs(os.path.join(d, sl), exist_ok=True); return open(os.path.join(d, sl, k + ".tmp"), "wb")
    def done(d, sl, k, f):
        if f is not None: f.close(); os.replace(os.path.join(d, sl, k + ".tmp"), os.path.join(d, sl, k))
    for name, fn in OPS:
        rng = random.Random(f"{SEED}:{name}")   # per-op stream: ops are independent
        OP, CASE, N_ROWS_OBJ = name, name + ".(corpus-setup)", 0
        PARTIAL_DIR = None if heavy_dir is None or skip_heavy else os.path.join(heavy_dir, SLICE_OF[name])
        _arm()
        hdr, rows, expect, *rest = fn(); exts = rest[0] if rest else {}; stub = "stub" in rest
        ids = [i for i, _ in rows.ids]; assert len(ids) == len(set(ids)), "duplicate case_id"
        if stub: assert not rows, "stub file must have zero rows"
        sl = out.setdefault(SLICE_OF[name], {"gate": {}, "heavy": {}, "names": []})
        k = f"la_{name}.txt"; f = to(gate_dir, SLICE_OF[name], k)
        sl["gate"][k] = stream_file(name, "gate", hdr, rows, False, expect, f); done(gate_dir, SLICE_OF[name], k, f)
        if not stub:
            k = f"la_{name}_heavy.txt"; sl["names"].append(k); f = None if skip_heavy else to(heavy_dir, SLICE_OF[name], k)
            r = stream_file(name, "heavy", hdr, rows, True, expect, f); done(heavy_dir, SLICE_OF[name], k, f)
            if not skip_heavy: sl["heavy"][k] = r
        want = [i for i, h in rows.ids if not h or ".cosmo_" in i]
        for tag, (ehdr, erows, eexp) in exts.items():
            assert [i for i, _ in erows.ids] == want, "ext case_ids must equal the gate case_ids (plus cosmo_* heavy ids)"
            k = f"la_{name}_{tag}_heavy.txt"; sl["names"].append(k); f = None if skip_heavy else to(heavy_dir, SLICE_OF[name], k)
            r = stream_file(name, "heavy", ehdr, erows, True, eexp, f); done(heavy_dir, SLICE_OF[name], k, f)
            if not skip_heavy: sl["heavy"][k] = r
            erows.close()
        rows.close()
    if BUDGET > 0: signal.setitimer(signal.ITIMER_REAL, 0)
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", nargs="?", const="", default=None, metavar="DIR")
    ap.add_argument("--heavy", metavar="DIR"); ap.add_argument("--skip-heavy", action="store_true")
    def positive(v):
        x = float(v)
        if not x > 0: raise argparse.ArgumentTypeError(f"must be > 0 seconds, got {v}")
        return x
    ap.add_argument("--case-budget", type=positive, default=1800.0, metavar="SECONDS",
                    help="per-case wall budget; exceeding it exits 5 naming the op and case_id (default 1800)")
    a = ap.parse_args()
    global BUDGET
    BUDGET = a.case_budget; signal.signal(signal.SIGALRM, _on_alarm)
    try:
        run(a)
    except CaseBudget:
        print(f"gen_la_goldens: CASE BUDGET exceeded ({a.case_budget:g} s): op={OP} case={CASE} "
              f"elapsed_s={time.monotonic() - T0:.3f}; no row finished within the budget. Not a pass.", file=sys.stderr, flush=True)
        sys.exit(5)

def run(a):
    if a.check is not None:
        cdir = a.check or tempfile.mkdtemp(dir=(os.makedirs(SPOOL_DIR, exist_ok=True) or SPOOL_DIR))
        print("check: heavy mode=pins (HEAVY.SHA256 pin/sha check only, no recompute)" if a.skip_heavy else "check: heavy mode=recompute", flush=True)
        out = render(cdir, None, a.skip_heavy); bad = 0
        for sl, t in sorted(out.items()):
            d = os.path.join(SLICES, sl); sums = read_sums(os.path.join(d, "SHA256SUMS"))
            for k, (h, _) in sorted(t["gate"].items()):
                p = os.path.join(d, k); g = os.path.join(cdir, sl, k)
                if not os.path.exists(p): bad += 1; print(f"check: MISSING committed {sl}/{k}")
                elif sha_path(p) == h: print(f"check: OK {sl}/{k}")
                else: bad += 1; print(f"check: DIFF {sl}/{k} line {first_diff_line(p, g)}")
                if sums.get(k) != h: bad += 1; print(f"check: SHA256SUMS mismatch {sl}/{k}")
                if not a.check: os.remove(g)
            if a.skip_heavy: bad += check_pins(sl, d, t["names"])
            if t["heavy"]:
                hs = read_sums(os.path.join(d, "HEAVY.SHA256"))
                for k, (h, _) in sorted(t["heavy"].items()):
                    if hs.get(k) == h: print(f"check: OK heavy {sl}/{k} {h[:16]}")
                    else: bad += 1; print(f"check: HEAVY sha256 mismatch {sl}/{k}")
        if not a.check:
            for sl in out: os.rmdir(os.path.join(cdir, sl)) if os.path.isdir(os.path.join(cdir, sl)) else None
            os.rmdir(cdir)
        print("check:", "OK" if not bad else f"FAIL ({bad})"); sys.exit(1 if bad else 0)
    out = render(SLICES, a.heavy, a.skip_heavy)
    for sl, t in sorted(out.items()):
        d = os.path.join(SLICES, sl)
        for k, (_, n) in sorted(t["gate"].items()): print(f"wrote {sl}/{k} ({n} bytes)")
        hdr = sums_header(d)
        open(os.path.join(d, "SHA256SUMS"), "w", encoding="utf-8", newline="\n").write(hdr + "".join(f"{t['gate'][k][0]}  {k}\n" for k in sorted(t["gate"])))
        if a.skip_heavy: print(f"skip-heavy: {sl}/HEAVY.SHA256 left as is (not regenerated)"); continue
        open(os.path.join(d, "HEAVY.SHA256"), "w", encoding="utf-8", newline="\n").write(HEAVY_HDR + "".join(f"{t['heavy'][k][0]}  {k}\n" for k in sorted(t["heavy"])))
        for k, (_, n) in sorted(t["heavy"].items()):
            if a.heavy: print(f"wrote heavy {a.heavy}/{sl}/{k} ({n} bytes)")

if __name__ == "__main__":
    main()
