"""Per-row condition / tightness columns for the exact LA golds (format 2; Spence's spec from Cassie's review, Harper 2026-10-29).

Importable helper, shared by scripts/gen_la_goldens.py (la_logdet, la_solve) and reusable for the LA3 lstsq golds (Rusty).
Pure exact arithmetic on gmpy2.mpq (no floats, no sympy). Inputs are exact rationals (anything gmpy2.mpq accepts).

    la_cond_cols(A, perm, L, U, scale_row, rref_kind, x=None) -> CondCols   (square: logdet, solve)
        A          n x n exact matrix (list of rows), nonsingular
        perm, L, U exact partial-pivoting factors with P A = L U (perm[i] = row of A placed at row i of P A)
        scale_row  exact rational >= 0: scale_row_ln() for ln outputs (logdet), scale_row_solve(x) for solve,
                   scale_row_lstsq(x) for lstsq (both ||x||_inf)
        rref_kind  "logdet" (trace form) or "solve" (forward-error form; needs x)
    make_cond_cols(cond_row, dim, scale_row, rref) -> CondCols
        the shared tail: k_row from (cond_row, dim), containment_only and reason from k_row. Every op builds its CondCols
        here, so cond_tokens and FORMAT_LINES hold for all of them. LA3 (Rusty) adds la_cond_cols_lstsq next to la_cond_cols,
        with Spence's cond_row = (|| |A^+||A| ||*||x|| + || |(A^TA)^-1||A^T||r| ||)/||x|| and dim = m (k_row uses ceil(log2 m)),
        and returns make_cond_cols(cond_row, m, scale_row_lstsq(x), rref).
    cond_tokens(c) -> the six row tokens, in COLUMNS order, formatted as in FORMAT_LINES.

Definitions (G = P^T |L| |U|, so |A| <= G entrywise; |.| entrywise; gamma_n ~ n * 2^-prec, first order):
    cond_row = || |A^-1| G ||_inf                        exact rational >= 1 (componentwise Skeel condition with LU growth)
    k_row    = 4 + ceil(log2 dim) + ceil(log2 cond_row)  exact integer; dim = n (square ops), m rows (lstsq)
    scale_row: 1 (ln outputs), ||x||_inf (solve, lstsq)
    containment_only, reason: one of REASONS (closed): 2 lstsq_x_zero (LA3: x = 0, scale_row undefined; the lstsq builder
                              sets it), else     containment_only, reason: 1 k_row_too_large when prec > k_row + 8 holds for no replay prec in PRECS
                              (k_row + 8 >= max(PRECS)); else 0 tight. Prec-independent: 1 = containment only at every
                              replay prec; 0 = the bound applies at each prec with prec > k_row + 8 (the harness decides per prec).
    rref (diagnostic only, never a gate): the coefficient c with first-order bound c * 2^-prec from |dA| <= gamma_n G:
        logdet: c = n * tr(|A^-1| G)          (d ln|det A| = tr(A^-1 dA))
        solve : c = n * || |A^-1| G |x| ||_inf (dx = -A^-1 dA x)
Replay tightness rule (tests/harness): rad <= 2^(k_row - prec) * max(|v|, scale_row) when prec > k_row + 8, containment only otherwise.
"""
from collections import namedtuple
from gmpy2 import mpq as F

PRECS = (53, 113, 192)
COLUMNS = ("cond_row", "k_row", "scale_row", "containment_only", "reason", "rref")
FORMAT_VERSION = 2
# The one closed containment_only code/reason set (shared with LA3 lstsq; tests/harness/slice_cholesky.rs COND_REASONS):
# 0 tight, 1 k_row_too_large (k_row + 8 >= max(PRECS)), 2 lstsq_x_zero (lstsq x = 0: scale_row undefined, containment only).
REASONS = ((0, "tight"), (1, "k_row_too_large"), (2, "lstsq_x_zero"))
TIGHT, K_ROW_TOO_LARGE = REASONS[:2]   # REASONS[2] is set only by the LA3 lstsq builder
FORMAT_LINES = [
    f"format: {FORMAT_VERSION} (format 1 had no cond columns; readers accept both, keyed on this line)",
    "cond columns: after the value columns of a non-refusal row, in this order: " + " ".join(COLUMNS) + " (refusal rows have none)",
    "cond_row: exact rational p/q = || |A^-1| P^T|L||U| ||_inf (componentwise Skeel condition with growth from the exact L, U of la_lu's pivoting)",
    "k_row: decimal integer = 4 + ceil(log2 n) + ceil(log2 cond_row), exact",
    "scale_row: exact rational p/q; 1/1 for ln outputs (logdet), ||x||_inf for solve",
    f"containment_only reason: closed set {' | '.join(f'{c} {r}' for c, r in REASONS)} (k_row_too_large: k_row + 8 >= {max(PRECS)}: no replay prec in {{{', '.join(map(str, PRECS))}}} has prec > k_row + 8; lstsq_x_zero: lstsq x = 0, scale_row undefined). The code does not depend on prec: 1 or 2 means containment only at EVERY replay prec; 0 does not mean tight at every prec, the replay still applies prec > k_row + 8 per prec",
    "rref: exact rational p/q, diagnostic only (never a gate): c with first-order bound c*2^-prec from |dA| <= gamma_n P^T|L||U|, gamma_n ~ n*2^-prec; logdet c = n*tr(|A^-1| P^T|L||U|), solve c = n*|| |A^-1| P^T|L||U| |x| ||_inf",
    "tightness (replay): rad <= 2^(k_row - prec) * max(|v|, scale_row) when prec > k_row + 8, containment only otherwise",
]

CondCols = namedtuple("CondCols", COLUMNS)

def _q(x):
    x = F(x); return f"{x.numerator}/{x.denominator}"

def ceil_log2(x):
    """ceil(log2 x) for exact rational x > 0."""
    x = F(x); assert x > 0
    n, d = x.numerator, x.denominator
    k = n.bit_length() - d.bit_length()          # 2^(k-1) < x < 2^(k+1)
    while F(2) ** k < x: k += 1
    while F(2) ** (k - 1) >= x: k -= 1
    return k

def scale_row_ln(): return F(1)
def scale_row_solve(x): return max((abs(F(v)) for v in x), default=F(0))
def scale_row_lstsq(x):
    """lstsq wrapper (LA3, Rusty): scale_row = ||x||_inf of the exact least-squares solution x."""
    return scale_row_solve(x)

def _inverse(perm, L, U):
    n = len(L); inv = [[F(0)] * n for _ in range(n)]
    for c in range(n):
        y = []
        for i in range(n): y.append(F(int(perm[i] == c)) - sum((F(L[i][j]) * y[j] for j in range(i)), F(0)))
        x = [F(0)] * n
        for i in reversed(range(n)):
            x[i] = (y[i] - sum((F(U[i][j]) * x[j] for j in range(i + 1, n)), F(0))) / F(U[i][i])
        for i in range(n): inv[i][c] = x[i]
    return inv

def la_cond_cols(A, perm, L, U, scale_row, rref_kind, x=None):
    n = len(A); assert n > 0 and all(F(U[i][i]) != 0 for i in range(n)), "la_cond_cols: singular input"
    LU = [[sum((abs(F(L[i][k])) * abs(F(U[k][j])) for k in range(n)), F(0)) for j in range(n)] for i in range(n)]
    G = [None] * n
    for i in range(n): G[perm[i]] = LU[i]
    for i in range(n):
        for j in range(n): assert abs(F(A[i][j])) <= G[i][j]
    Ai = [[abs(v) for v in r] for r in _inverse(perm, L, U)]
    M = [[sum((Ai[i][k] * G[k][j] for k in range(n)), F(0)) for j in range(n)] for i in range(n)]
    cond = max(sum(r, F(0)) for r in M)
    if rref_kind == "logdet": rref = F(n) * sum((M[i][i] for i in range(n)), F(0))
    elif rref_kind == "solve":
        assert x is not None
        rref = F(n) * max(sum((M[i][j] * abs(F(x[j])) for j in range(n)), F(0)) for i in range(n))
    else: raise ValueError(f"la_cond_cols: rref_kind {rref_kind!r}")
    return make_cond_cols(cond, n, scale_row, rref)

def make_cond_cols(cond_row, dim, scale_row, rref):
    """CondCols from an exact cond_row >= 1 and the op dimension dim >= 1 (n for square ops, m for lstsq); see the module doc."""
    cond_row = F(cond_row); assert cond_row >= 1 and dim >= 1
    k_row = 4 + (dim - 1).bit_length() + ceil_log2(cond_row)
    co, reason = K_ROW_TOO_LARGE if k_row + 8 >= max(PRECS) else TIGHT
    return CondCols(cond_row, k_row, F(scale_row), co, reason, F(rref))

def cond_tokens(c):
    return [_q(c.cond_row), str(c.k_row), _q(c.scale_row), str(c.containment_only), c.reason, _q(c.rref)]

# ---- LA3 lstsq (Rusty 2026-10-29; Spence's definition, CONTRACT 16.LA.5). Separate from la_cond_cols (square LU) above. ----
LSTSQ_FORMAT_LINES = [
    "lstsq cond (16.LA.5): with A m x n full column rank, A^+ = (A^T A)^-1 A^T, x = A^+ b, r = b - A x (all exact):",
    "  cond_row = ( || |A^+||A| ||_inf * ||x||_inf + || |(A^T A)^-1| |A^T| |r| ||_inf ) / ||x||_inf  (exact-fit r = 0: first term only; residual term not capped)",
    "  k_row = 4 + ceil(log2 m) + ceil(log2 cond_row) (m = row count); scale_row = ||x||_inf",
    "  rref (diagnostic, never a gate) = m * || |A^+||A||x| + |(A^T A)^-1||A^T||r| ||_inf",
    "  x = 0 (lstsq_x_zero, in the closed reason set above): cond_row and k_row from the first term || |A^+||A| ||_inf only",
]

def _inv_exact(N):
    """Exact inverse of a nonsingular square rational matrix (Gauss-Jordan, nonzero pivot)."""
    n = len(N); M = [[F(v) for v in r] + [F(int(i == j)) for j in range(n)] for i, r in enumerate(N)]
    for c in range(n):
        p = next(i for i in range(c, n) if M[i][c] != 0); M[c], M[p] = M[p], M[c]
        pv = M[c][c]; M[c] = [v / pv for v in M[c]]
        for i in range(n):
            if i != c and M[i][c] != 0:
                f = M[i][c]; M[i] = [a - f * b for a, b in zip(M[i], M[c])]
    return [r[n:] for r in M]

def _mm(X, Y):
    return [[sum((X[i][k] * Y[k][j] for k in range(len(Y))), F(0)) for j in range(len(Y[0]))] for i in range(len(X))]

def _mv(X, v):
    return [sum((a * b for a, b in zip(r, v)), F(0)) for r in X]

def _abs(X):
    return [[abs(v) for v in r] for r in X]

def la_cond_cols_lstsq(A, b, x):
    """CondCols for an lstsq value row (16.LA.5). A m x n, m >= n >= 1, full column rank; x the exact solution."""
    A = [[F(v) for v in r] for r in A]; b = [F(v) for v in b]; x = [F(v) for v in x]
    m, n = len(A), len(A[0]); assert m >= n >= 1 and len(b) == m and len(x) == n
    At = [[A[i][j] for i in range(m)] for j in range(n)]
    Ni = _inv_exact(_mm(At, A)); Ap = _mm(Ni, At)
    r = [bi - ai for bi, ai in zip(b, _mv(A, x))]
    assert all(v == 0 for v in _mv(At, r)), "la_cond_cols_lstsq: x is not the least-squares solution"
    T1 = _mm(_abs(Ap), _abs(A)); T2 = _mm(_abs(Ni), _abs(At))
    t1 = max(sum(row, F(0)) for row in T1)
    t2v = _mv(T2, [abs(v) for v in r]); t2 = max(t2v)
    xn = scale_row_lstsq(x)
    rref = F(m) * max(u + w for u, w in zip(_mv(T1, [abs(v) for v in x]), t2v))
    if xn == 0:
        cond = t1; reason = "lstsq_x_zero"
    else:
        cond = (t1 * xn + t2) / xn; reason = None
    assert cond >= 1
    k_row = 4 + (m - 1).bit_length() + ceil_log2(cond)
    if reason is None: reason = "k_row_too_large" if k_row + 8 >= max(PRECS) else "tight"
    return CondCols(F(cond), k_row, F(xn), dict((r, c) for c, r in REASONS)[reason], reason, rref)
