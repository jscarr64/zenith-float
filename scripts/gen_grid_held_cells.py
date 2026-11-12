#!/usr/bin/env python3
"""Exact per-part oracle for complex_mul_div_grid_digest dump lines (ExactComplex mul/div).

Decoding (confirmed against the library): value = sign * M * 2^(e - 64*W), M the little-endian words, W
their count; n is the significant-bit count and is NOT the scale (subnormals have n < 64*W).
Every operand is c * t^k with t = 2^EXPONENT_MIN (k = 1 for the subnormal operand, else 0) and c rational. A part is a
Laurent series in t; t is far below any ulp, so the value is the leading term c0 * t^K and the sign of the first
non-zero later term decides a sticky/tie. K = 0: normal range. K = -1: c0 * 2^EXPONENT_MAX (overflow iff the rounded
|c0| >= 1). K <= -2: overflow. K >= 1: below 2^(EXPONENT_MIN-1): judged only with ZF_SUBNORMAL_MODEL=1, else TINY.
Expected per mode at p (None: toward zero at rp(p,64)). Overflow: +-Inf for nearest/away, +-MAX for toward-zero side
(the 12:15 PM CT ruling). Exact rationals only (no mpmath needed). Usage: gen_grid_held_cells.py DUMP OUT_TSV (dump: ZF_DIGEST_DUMP of a run that hashes every cell)
       ZF_SUBNORMAL_MODEL=1 gen_grid_held_cells.py DUMP OUT_TSV (also judge the subnormal range);
       gen_grid_held_cells.py --emit-rust CLASSES_TSV DUMP OUT_DIR
"""
import sys, re, collections
from fractions import Fraction as Fr
EMAX = (1 << 62) - 1
EMIN = -EMAX
import os
SUBNORMAL_MODEL = os.environ.get('ZF_SUBNORMAL_MODEL') == '1'   # gradual underflow, grid 2^(EXPONENT_MIN-P)
R = re.compile(r'e: (-?\d+), s: (Pos|Neg), m: Mantissa \{ m: WordBuf \{ inner: \[([^\]]*)\] \}, n: (\d+) \}')
def dec(s):
    """-> ('NaN',) | ('Inf', sign) | ('V', c, k) value c*t^k"""
    if 'NaN(' in s: return ('NaN',)
    if 'Inf(' in s: return ('Inf', -1 if 'Inf(Neg)' in s else 1)
    m = R.search(s); e, sg, ws = int(m[1]), m[2], m[3]
    W = [w for w in ws.split(',') if w.strip()]
    M = sum(int(w.strip().strip('"')) << (64 * i) for i, w in enumerate(W))
    sgn = -1 if sg == 'Neg' else 1
    if M == 0: return ('V', Fr(0), 0)
    x = e - 64 * len(W)
    if e < -(1 << 61): return ('V', sgn * M * pw(x - EMIN), 1)
    if e > (1 << 61): return ('V', sgn * M * pw(x - EMAX), -1)
    return ('V', sgn * M * pw(x), 0)
def pw(n): return Fr(2) ** n if n >= 0 else Fr(1, 2 ** (-n))
# polynomials in t: dict power -> Fraction
def P(c, k): return {k: c} if c != 0 else {}
def add(a, b, sg=1):
    r = dict(a)
    for k, v in b.items():
        r[k] = r.get(k, 0) + sg * v
        if r[k] == 0: del r[k]
    return r
def mul(a, b):
    r = {}
    for i, x in a.items():
        for j, y in b.items():
            r[i + j] = r.get(i + j, 0) + x * y
    return {k: v for k, v in r.items() if v != 0}
def series_div(n, d, terms=4):
    """n/d as (K, [c0, c1, ...], exact) with exact = the division terminates."""
    if not n: return (0, [Fr(0)], True)
    a = min(n); b = min(d); K = a - b
    rem = dict(n); cs = []
    for i in range(terms):
        k = a + i
        c = rem.get(k, Fr(0)) / d[b]
        cs.append(c)
        if c != 0: rem = add(rem, mul(P(c, k - b), d), -1)
        if not rem: return (K, cs, True)
    return (K, cs, False)
def poly_series(p):
    if not p: return (0, [Fr(0)], True)
    a = min(p); top = max(p)
    return (a, [p.get(a + i, Fr(0)) for i in range(top - a + 1)], True)
def tail_sign(cs, exact):
    for c in cs[1:]:
        if c != 0: return 1 if c > 0 else -1
    return 0 if exact else None   # None: more terms needed (not seen)
def rnd(x, d, P, mode, emin=None):
    """round x + d*eps (eps infinitesimal > 0, d in -1/0/1) to P bits. -> (Fraction, inexact)
    emin: lowest exponent E (value in [2^(E-1), 2^E)); below it the grid stays 2^(emin-P) (gradual underflow)."""
    if x == 0:
        if d == 0: return Fr(0), False
        # an infinitesimal of sign d below the whole grid: 0 or one grid step 2^(emin-P), per mode
        neg = d < 0
        upmag = {'ToZero': False, 'FromZero': True, 'Up': not neg, 'Down': neg, 'ToEven': False}[mode]
        r = pw(emin - P) if upmag else Fr(0)
        return (-r if neg else r), True
    neg = x < 0; a = -x if neg else x; dm = -d if neg else d
    E = 0
    num, den = a.numerator, a.denominator
    E = num.bit_length() - den.bit_length()
    if pw(E - 1) > a: E -= 1
    if a >= pw(E): E += 1
    if emin is not None and E < emin: E = emin
    s = a * pw(P - E); fl = s.numerator // s.denominator; fr = s - fl
    if fr == 0 and dm < 0:
        if fl == 1 << (P - 1) and (emin is None or E - 1 >= emin): E -= 1; fl = (1 << P) - 1
        else: fl -= 1
        frc = (Fr(1), -1)
    else: frc = (fr, dm if fr in (0, Fr(1, 2)) else 0)
    inexact = not (frc[0] == 0 and frc[1] == 0)
    if not inexact: return (-a if neg else a), False
    upmag = {'ToZero': False, 'FromZero': True, 'Up': not neg, 'Down': neg}.get(mode)
    if mode == 'ToEven':
        f, e2 = frc
        if f > Fr(1, 2) or (f == Fr(1, 2) and e2 > 0): upmag = True
        elif f < Fr(1, 2) or (f == Fr(1, 2) and e2 < 0): upmag = False
        else: upmag = fl % 2 == 1
    q = fl + (1 if upmag else 0)
    r = q * pw(E - P)
    return (-r if neg else r), True
def rp(p): return max(64, -(-p // 64) * 64)
def expect(K, cs, exact, p, rm):
    """-> ('V', value, k, inexact) | ('Inf', sign) | ('MAX', sign, P) | ('TINY',)"""
    P_ = rp(p) if rm == 'None' else p
    mode = 'ToZero' if rm == 'None' else rm
    c0 = cs[0]
    if c0 == 0 and K == 0 and exact: return ('V', Fr(0), 0, False)
    d = tail_sign(cs, exact)
    if K >= 1 and SUBNORMAL_MODEL:
        if d is None: return ('UNDECIDED',)
        if K == 1: r, ix = rnd(c0, d, P_, mode, emin=0)
        else: r, ix = rnd(Fr(0), 1 if c0 > 0 else -1, P_, mode, emin=0)
        return ('V', r, 1 if r != 0 else 0, ix)
    if K >= 1: return ('TINY',)
    if d is None: return ('UNDECIDED',)
    sg = 1 if c0 > 0 else -1
    tz = mode == 'ToZero' or (mode == 'Down' and sg > 0) or (mode == 'Up' and sg < 0)
    if K <= -2: return ('MAX', sg, P_) if tz else ('Inf', sg)
    r, ix = rnd(c0, d, P_, mode)
    if K == -1 and abs(r) >= 1: return ('MAX', sg, P_) if tz else ('Inf', sg)
    return ('V', r, K, ix)
def main():
    dump, out = sys.argv[1], sys.argv[2]
    cnt = collections.Counter(); fo = open(out, 'w')
    fo.write("cell (test:line)\tp\tmode\tints\tclass\tlog2|true|-EXPONENT_MAX (K=-1) or K\tlib\texpected\n")
    import math
    rows = {}
    lines = open(dump).read().splitlines()
    for ln, l in enumerate(lines, 1):
        t = l.rstrip('\n').split('\t'); p = int(t[1]); rm = t[2]; k, part = eval(t[3])
        a = [dec(x.rsplit('|', 1)[0]) for x in t[4:8]]
        lib = dec(t[9]); lix = t[10] == 'true'
        if any(x[0] != 'V' for x in a): rows[ln] = ('SKIP-nonfinite-operand', None, None, ''); continue
        A, B, C, D = (P(x[1], x[2]) for x in a)
        if k == 0:
            K, cs, ex = poly_series(add(mul(A, C), mul(B, D), -1) if part == 0 else add(mul(A, D), mul(B, C)))
        else:
            den = add(mul(C, C), mul(D, D))
            if not den: rows[ln] = ('SKIP-pole', None, None, ''); continue
            num = add(mul(A, C), mul(B, D)) if part == 0 else add(mul(B, C), mul(A, D), -1)
            K, cs, ex = series_div(num, den)
        e = expect(K, cs, ex, p, rm)
        mag = ''
        if cs[0] != 0 and K == -1: mag = f"{math.log2(abs(cs[0].numerator)) - math.log2(cs[0].denominator):+.4f}"
        elif cs[0] != 0: mag = f"K={K}"
        if e[0] in ('TINY', 'UNDECIDED'): cls = e[0]
        elif lib[0] == 'NaN': cls = 'WRONG-NaN'
        elif e[0] == 'Inf':
            cls = 'OK-Inf' if lib == ('Inf', e[1]) else 'WRONG-want-Inf'
        elif e[0] == 'MAX':
            if lib[0] == 'Inf': cls = 'WRONG-Inf-want-MAX'
            else:
                mx = e[1] * (1 - pw(-e[2]))
                cls = 'OK-MAX' if lib[0] == 'V' and lib[2] == -1 and lib[1] == mx else 'WRONG-want-MAX'
        else:
            if lib[0] == 'Inf': cls = 'WRONG-Inf-true-finite'
            elif (lib[1], lib[2] if lib[1] != 0 else 0) != (e[1], e[2] if e[1] != 0 else 0): cls = 'WRONG-value'
            else: cls = 'OK-value'
        rows[ln] = (cls, e, lix, mag)
    # inexact is value-wide for ExactComplex (both parts carry re.ix || im.ix): judge the flag per pair
    for ln, l in enumerate(lines, 1):
        t = l.rstrip('\n').split('\t'); p = int(t[1]); rm = t[2]; k, part = eval(t[3])
        cls, e, lix, mag = rows[ln]
        if cls == 'OK-value':
            mate = rows[ln + 1] if part == 0 else rows[ln - 1]
            eix = [e[3] if e[0] == 'V' else True]
            me = mate[1]
            if mate[0] in ('TINY', 'UNDECIDED'): cls = 'OK-value-flag-unjudged'
            else:
                mix = (me[3] if me[0] == 'V' else True) if me is not None else False
                cls = 'OK' if lix == (eix[0] or mix) else 'WRONG-flag-only'
        cnt[(cls, rm)] += 1
        if not cls.startswith('OK') and not cls.startswith('SKIP'):
            fo.write(f"complex_mul_div_grid_digest:{ln}\t{p}\t{rm}\t{[k, part]}\t{cls}\t{mag}\t{t[9][:160]}\t{e[:2] if e[0]!='V' else ('V', float(e[1]) if e[2]==0 else f'{float(e[1])}*2^EMAX', e[3])}\n")
    for kk, v in sorted(cnt.items()): print(kk, v)

def any_model_wrong(line):
    """A subnormal-range (K >= 1) cell: is the library part wrong under every subnormal grid? S = values every model
    can represent: 0 and the P-bit values of magnitude >= 2^(EXPONENT_MIN - 1). A directed result must lie between
    the true value v and the nearest member of S on the rounding side (Down: max S <= v .. v; Up: v .. min S >= v;
    ToZero/None and FromZero by the sign of v). Outside that interval it is wrong in any model."""
    t = line.split('\t'); p = int(t[1]); rm = t[2]; k, part = eval(t[3])
    a = [dec(x.rsplit('|', 1)[0]) for x in t[4:8]]; L = dec(t[9])
    A, B, C, D = (P(x[1], x[2]) for x in a)
    if k == 0:
        K, cs, ex = poly_series(add(mul(A, C), mul(B, D), -1) if part == 0 else add(mul(A, D), mul(B, C)))
    else:
        den = add(mul(C, C), mul(D, D))
        num = add(mul(A, C), mul(B, D)) if part == 0 else add(mul(B, C), mul(A, D), -1)
        K, cs, ex = series_div(num, den)
    c0 = cs[0]; d = tail_sign(cs, ex)
    assert K >= 1 and c0 != 0 and d is not None
    vc, vd = (c0, d) if K == 1 else (Fr(0), 1 if c0 > 0 else -1)      # v = vc*t + vd*eps
    Lt = Fr(0) if L[0] == 'V' and L[1] == 0 else (L[1] if L[0] == 'V' and L[2] == 1 else None)
    if Lt is None: return True
    P_ = rp(p) if rm == 'None' else p
    mode = 'ToZero' if rm == 'None' else rm
    half = Fr(1, 2)
    def s_down(x, dd):
        if abs(x) >= half: return rnd(x, dd, P_, 'Down')[0]
        return Fr(0) if (x > 0 or (x == 0 and dd >= 0)) else -half
    lo, hi = s_down(vc, vd), -s_down(-vc, -vd)
    neg = vc < 0 or (vc == 0 and vd < 0)
    if mode == 'ToZero': mode = 'Up' if neg else 'Down'
    if mode == 'FromZero': mode = 'Down' if neg else 'Up'
    le_v = Lt < vc or (Lt == vc and vd > 0); ge_v = Lt > vc or (Lt == vc and vd < 0)
    if mode == 'Down': return not (lo <= Lt and le_v)
    if mode == 'Up': return not (ge_v and Lt <= hi)
    return False   # nearest: not judged model-independently


def emit_rust(classes_tsv, dump, out_dir):
    """Write the F-25 held-cell lists (1-based cell numbers of complex_mul_div_grid_digest) for `include!`.
    classes_tsv: output of a ZF_SUBNORMAL_MODEL=1 run on `dump`. Normal range: every WRONG-value cell with K = 0.
    Subnormal range: the WRONG-value cells with K >= 1 that are wrong in any subnormal model (any_model_wrong)."""
    lines = open(dump).read().splitlines()
    normal, sub, dep = [], [], []
    for l in list(open(classes_tsv))[1:]:
        t = l.split('\t')
        if t[4] != 'WRONG-value':
            continue
        n = int(t[0].split(':')[1])
        if t[5] == 'K=0': normal.append(n)
        elif any_model_wrong(lines[n - 1]): sub.append(n)
        else: dep.append(n)
    for name, xs, what in (('f25_normal', normal, 'true value in the normal range (K = 0)'),
                           ('f25_subnormal', sub, 'true value below 2^(EXPONENT_MIN - 1) (K >= 1), wrong under any '
                            'subnormal grid (gen_grid_held_cells.py any_model_wrong)')):
        s = ("// Generated by scripts/gen_grid_held_cells.py --emit-rust; do not edit. F-25 cells of\n"
             "// complex_mul_div_grid_digest whose library part is not the correctly rounded value (exact rational\n"
             f"// oracle), {what}. Each number is the cell's 1-based position in the digest loop. {len(xs)} cells.\n[\n")
        for k in range(0, len(xs), 12):
            s += "    " + ", ".join(str(x) for x in sorted(xs)[k:k + 12]) + ",\n"
        open(f"{out_dir}/{name}.in", 'w').write(s + "]\n")
        print(name, len(xs))
    print('model-dependent (kept in the pin):', len(dep), ' '.join(str(x) for x in sorted(dep)))

if __name__ == '__main__':
    if len(sys.argv) == 5 and sys.argv[1] == '--emit-rust': emit_rust(sys.argv[2], sys.argv[3], sys.argv[4])
    else: main()
