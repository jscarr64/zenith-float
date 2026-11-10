#!/home/jsc64/.local/oracle-venv/bin/python
"""LA3 lstsq golds: tests/goldens/slices/qr/la_lstsq.txt, in the la_qr.txt format (Rusty, slice LA3 owner, 2026-10-29).

Usage (oracle venv, under popslot/poprun like gen_la_goldens.py):
    ~/.local/oracle-venv/bin/python scripts/gen_la_lstsq_goldens.py            # writes la_lstsq.txt and its SHA256SUMS line
    ~/.local/oracle-venv/bin/python scripts/gen_la_lstsq_goldens.py --check DIR # writes DIR/qr/la_lstsq.txt, compares sha256 with the committed file
    --small: the 10-row smoke corpus to stdout-dir only (never committed)

Reuses scripts/gen_la_goldens.py (popguard, oracle pins, exact mpq, Rows spool-to-disk streaming,
header_bytes, enclose_rational, err, gram_schmidt, rmat, hilbert). NOTE: gen_la_goldens.py rewrites slices/qr/SHA256SUMS from its own file set; until gen_lstsq() is folded into it, re-add the la_lstsq.txt line by running this script after it.
x = argmin ||A x - b||_2 = (A^T A)^-1 A^T b is RATIONAL, so x is computed exactly (mpq) and enclosed at 256 bits by
enclose_rational (lo == hi iff exact); mpmath at 1200 bits re-checks every enclosure (mp.lu_solve on the normal equations).
Cond columns (format 2, 16.LA.5): scripts/la_cond_cols.py la_cond_cols_lstsq, via the merged harper/zb-la-golds.
Spence 2026-10-29: b.len() != m -> ERR DimensionMismatch b_len (checked before n = 0); n = 0 (m x 0, 0 x 0) -> Ok(empty x).
"""
import os, sys
REPO = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(REPO, "scripts"))
import gen_la_goldens as G          # runs popguard + oracle-venv pin checks on import
import mpmath, hashlib
G.GEN = "scripts/gen_la_lstsq_goldens.py"
F = G.F
HERE = os.path.dirname(os.path.abspath(__file__))
G.SPOOL_DIR = os.path.join(os.environ.get("TMPDIR", os.path.join(REPO, "target")), "la-lstsq-spool")
G.SLICE_OF["lstsq"] = "qr"; G.OP = "lstsq"
G.rng = __import__("random").Random(f"{G.SEED}:lstsq")
from la_cond_cols import la_cond_cols_lstsq, cond_tokens, FORMAT_LINES as COND_FORMAT, LSTSQ_FORMAT_LINES  # 16.LA.5 (merged harper/zb-la-golds)
def cond_cols(A, b, x):
    return cond_tokens(la_cond_cols_lstsq(A, b, x))

# HOOK(Harper): pinned python-flint arb cross-check of every x enclosure goes here (runs after mp_check). Not implemented.
def arb_check(A, b, x, encl):
    return None

def lstsq_exact(A, b):
    m, n = len(A), len(A[0]) if A else 0
    At = G.tr(A); N = G.mm(At, A); c = [sum((At[i][k] * b[k] for k in range(m)), F(0)) for i in range(n)]
    return G.solve(N, c)            # None iff rank deficient

def mp_check(A, b, x, encl):
    mpmath.mp.prec = 1200   # normal equations square cond(A) (Hilbert 12: ~2^110); 1200 bits leaves > 2^-900 accuracy
    At = mpmath.matrix([[mpmath.mpf(v.numerator) / v.denominator for v in r] for r in G.tr(A)])
    Am = At.T; bm = mpmath.matrix([mpmath.mpf(v.numerator) / v.denominator for v in b])
    xm = mpmath.lu_solve(At * Am, At * bm)
    for i, ((lo, hi)) in enumerate(encl):
        l = mpmath.ldexp(lo[0], lo[1]); h = mpmath.ldexp(hi[0], hi[1])
        tol = mpmath.ldexp(1, -200) * (abs(xm[i]) + 1)
        assert l - tol <= xm[i] <= h + tol, ("mpmath cross-check", i)
        assert F(lo[0]) * F(2) ** lo[1] <= x[i] <= F(hi[0]) * F(2) ** hi[1]

def corpus(small):
    C = [("tall.3x2", [[F(1), F(1)], [F(1), F(2)], [F(1), F(3)]], [F(1), F(2), F(2)]),
         ("square.basic2", [[F(3), F(1)], [F(4), F(2)]], [F(1), F(1)]),
         ("exact_fit.3x2", [[F(1), F(0)], [F(0), F(1)], [F(1), F(1)]], [F(1), F(2), F(3)]),
         ("refuse.zero_col_3x2", [[F(1), F(0)], [F(2), F(0)], [F(3), F(0)]], [F(1)] * 3),
         ("refuse.dupcol_4x2", [[F(i), F(i)] for i in range(1, 5)], [F(1)] * 4),
         ("refuse.wide2x3", [[F(1), F(0), F(0)], [F(0), F(1), F(0)]], [F(1), F(1)]),
         ("refuse.blen_3x2_b2", [[F(1), F(1)], [F(1), F(2)], [F(1), F(3)]], [F(1), F(1)]),
         ("refuse.blen_3x2_b4", [[F(1), F(1)], [F(1), F(2)], [F(1), F(3)]], [F(1)] * 4),
         ("x_zero.3x2", [[F(1), F(0)], [F(0), F(1)], [F(0), F(0)]], [F(0), F(0), F(1)]),
         ("empty.3x0", [[], [], []], [F(1)] * 3), ("empty.0x0", [], []),
         ("refuse.blen_3x0_b2", [[], [], []], [F(1)] * 2)]
    if not small:
        C += [("tall.6x3", G.rmat(6, 3), [G.rq() for _ in range(6)]), ("tall.12x4", G.rmat(12, 4), [G.rq() for _ in range(12)])]
        C += [("hilbert.%02d" % n, G.hilbert(n), [F(1)] * n) for n in range(5, 13)]
        C += [("hilbert_tall.%02dx%02d" % (n + 2, n), G.hilbert(n + 2)[:], [F(1)] * (n + 2)) for n in (4, 8)]
        C = [(c, ([r[:n] for r in A] if c.startswith("hilbert_tall") else A), b) for c, A, b in C
             for n in [int(c.split("x")[-1]) if c.startswith("hilbert_tall") else 0]]
    return C

EXACT_CERTIFIED = ["lstsq.refuse.zero_col_3x2", "lstsq.refuse.wide2x3", "lstsq.refuse.blen_3x2_b2", "lstsq.refuse.blen_3x2_b4", "lstsq.refuse.blen_3x0_b2",
                   "lstsq.empty.3x0", "lstsq.empty.0x0"]

def gen(small):
    rows = G.Rows()
    for cid, A, b in corpus(small):
        G.begin("lstsq." + cid)
        m = len(A); n = len(A[0]) if m else 0
        dims = f"{m}x{n}"; head = ["lstsq." + cid, dims, str(len(b))] + [G.q(v) for v in G.flat(A)] + [G.q(v) for v in b]
        if m < n: outs = G.err("DimensionMismatch", "wide")
        elif len(b) != m: outs = G.err("DimensionMismatch", "b_len")
        elif n == 0: outs = []  # Ok(empty x): zero x fields, zero lohi fields (Spence)
        else:
            x = lstsq_exact(A, b)
            if x is None: outs = G.err("Singular", "rank_deficient")
            else:
                encl = [G.enclose_rational(v, G.BITS) for v in x]; mp_check(A, b, x, encl); arb_check(A, b, x, encl)
                outs = [G.q(v) for v in x] + [G.pair(*e) for e in encl] + cond_cols(A, b, x)
        rows.append((" ".join(head + outs), G.is_heavy(A, b)))
    return rows

HDR = ["op: lstsq",
       "canonical: x = argmin ||A x - b||_2 for A m x n, m >= n, full column rank: x = (A^T A)^-1 A^T b, exact rational.",
       "refusals (16.LA.4 format): 'ERR Singular rank_deficient' if rank A < n; 'ERR DimensionMismatch wide' if m < n; 'ERR DimensionMismatch b_len' if len(b) != m (checked before n = 0). n = 0 (m x 0 or 0 x 0): Ok, x empty.",
       "exact-certified (must return the named error, not Resource): " + " ".join(EXACT_CERTIFIED),
       "enclosure: x[n] as lo hi pairs at b = 256, grid 2^e with e = floor(log2|x|) - b + 1 (enclose_rational); re-checked against mpmath (1200 bits, normal equations); python-flint arb hook pending (Harper).",
       "columns (gate): case_id dims(mxn) blen A[m*n] b[blen] (x[n] x_lohi256[2*n] cond_row k_row scale_row containment_only reason rref | ERR variant detail); m x 0 rows have no x, lohi or cond fields"] + COND_FORMAT + LSTSQ_FORMAT_LINES

def write(path, rows):
    hs = hashlib.sha256()
    with open(path, "wb") as f:
        for chunk in [G.header_bytes("lstsq", "gate", HDR, rows.count[False])] + list(rows.lines(False)):
            f.write(chunk); hs.update(chunk)
    return hs.hexdigest()

if __name__ == "__main__":
    small = "--small" in sys.argv
    rows = gen(small)
    assert rows.count[True] == 0, "lstsq corpus has no heavy rows"
    if small or "--check" in sys.argv:
        d = sys.argv[sys.argv.index("--check") + 1] if "--check" in sys.argv else os.environ.get("TMPDIR", "/tmp")
        os.makedirs(os.path.join(d, "qr"), exist_ok=True)
        h = write(os.path.join(d, "qr", "la_lstsq%s.txt" % ("_small" if small else "")), rows)
        if not small:
            want = G.read_sums(os.path.join(G.SLICES, "qr", "SHA256SUMS")).get("la_lstsq.txt")
            print("check:", "OK qr/la_lstsq.txt" if want == h else "DIFF qr/la_lstsq.txt"); rows.close(); sys.exit(0 if want == h else 1)
    else:
        d = os.path.join(G.SLICES, "qr"); h = write(os.path.join(d, "la_lstsq.txt"), rows)
        p = os.path.join(d, "SHA256SUMS"); keep = [l for l in open(p, encoding="utf-8") if not l.rstrip().endswith("  la_lstsq.txt")]
        body = [l for l in keep if not l.startswith("#")] + [f"{h}  la_lstsq.txt\n"]
        open(p, "w", encoding="utf-8", newline="\n").write("".join(l for l in keep if l.startswith("#")) + "".join(sorted(body, key=lambda l: l.split()[1])))
        print(f"wrote qr/la_lstsq.txt ({h[:16]})")
    rows.close(); print("rows:", rows.count)
