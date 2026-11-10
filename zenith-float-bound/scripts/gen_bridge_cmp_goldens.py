#!/home/jsc64/team-env/oracle-2026-10-30/bin/python
"""
Consumer: tests/harness/slice_bridge/harper_s2_goldens.rs (every row asserted exactly). Listed in tests/goldens/slices/bridge/GENERATORS and SHA256SUMS (checked by tools/goldens_regen.sh).
Generates tests/goldens/slices/bridge/bridge_cmp.txt: S2 rows for Mag::cmp_exact_num, Mag::cmp_exact_diff,
Mag::from_exact_num / from_exact_num_lower and Mag::to_exact_num.
TOKENS. ExactNum: 'nan' '+inf' '-inf' '+0' '-0' or '[+-]HEX^E' = +-HEX * 2^E (HEX lowercase, odd).
Mag: '0', 'inf' or 'HEX16:E' = 0.HEX * 2^E (canonical, top bit set), E any i64.
ROWS. 'N R X ORD' (cmp_exact_num(R, X)); 'D R A B ORD' (cmp_exact_diff); 'U X UP LO' (directed conversions);
'T R RES' (to_exact_num; RES an ExactNum token or err:inf|err:max|err:min). ORD in lt eq gt none.
TRUTH. Exact Python integers only: the sign of a sum of signed dyadics is found by clustering terms whose bit ranges
lie within 256 bits of each other (summed exactly), the top nonzero cluster decides (a gap of > 64 bits makes the lower
clusters unable to change the sign). Close-gap rows are cross-checked with mpmath libmp exact mpf_sub/mpf_cmp; up/lo are
cross-checked against mpmath from_man_exp(.., 64, 'f'/'c') plus containment and 1-ulp tightness. Any mismatch exits 3.
STREAMING: rows are written one at a time. SIZE: the 100000-bit values appear only in the few HUGE D rows (low thousands of rows,
well under 1 MB). EDGE: the hand-worked cmp_far_diff rows of tests/harness/slice_bridge/harper_s2_diff_edge.rs (radius above
EXPONENT_MAX: Mag(TOP, E+1), Mag(TOP+1, E+1), Mag(u64::MAX, E+1), Mag(u64::MAX, E)) are emitted with their hand truth asserted.
"""
import os
import resource
import sys

import goldlib
from mpmath import libmp as L

GENERATOR_VERSION = 2
GOLDEN_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tests", "goldens", "slices", "bridge")
GOLDEN_NAME = "bridge_cmp.txt"
EXPONENT_MAX = (1 << 62) - 1
EXPONENT_MIN = -EXPONENT_MAX
I64_MAX, I64_MIN = (1 << 63) - 1, -(1 << 63)
TOP = 1 << 63
SEED_SUITE = 300
RANDOM_PER_KIND = 600


class Bad(Exception):
    pass


def norm(s, m, e):
    """(s, m odd, e) or zero (s, 0, 0)."""
    if m == 0:
        return (s, 0, 0)
    while m % 2 == 0:
        m //= 2
        e += 1
    return (s, m, e)


def top_e(m, e):
    """0.m form exponent of m*2^e."""
    return e + m.bit_length()


def sum_sign(terms):
    """Sign (-1, 0, 1) of sum of (sign, m, e) terms, exact for any i64-scale exponents."""
    ts = sorted([t for t in terms if t[1] != 0], key=lambda t: -top_e(t[1], t[2]))
    clusters = []
    for s, m, e in ts:
        if clusters and top_e(m, e) >= clusters[-1][0] - 256:
            clusters[-1][1].append((s, m, e))
            clusters[-1][0] = min(clusters[-1][0], e)
        else:
            clusters.append([e, [(s, m, e)]])
    for lowe, cl in clusters:
        v = sum((-1 if s else 1) * (m << (e - lowe)) for s, m, e in cl)
        if v != 0:
            return 1 if v > 0 else -1
    return 0


def ordname(c):
    return {-1: "lt", 0: "eq", 1: "gt"}[c]


def mag_val(r):
    """Mag (mant, exp) -> dyadic (0, mant, exp - 64)."""
    return (0, r[0], r[1] - 64)


def cmp_num(r, x):
    if x in ("nan", "+inf", "-inf"):
        return "none"
    if r == "inf":
        return "gt"
    return ordname(sum_sign([mag_val(r), (1, x[1], x[2])]))


def cmp_diff(r, a, b):
    if a in ("nan", "+inf", "-inf") or b in ("nan", "+inf", "-inf"):
        return "none"
    if r == "inf":
        return "gt"
    d = sum_sign([a, (1 - b[0], b[1], b[2])])   # sign of a - b
    if d == 0:
        return cmp_num(r, (0, 0, 0))
    neg = lambda t: (1 - t[0], t[1], t[2])
    absd = [a, neg(b)] if d > 0 else [neg(a), b]   # |a - b| as a sum
    return ordname(sum_sign([mag_val(r)] + [neg(t) for t in absd]))   # sign of r - |a - b|


def mpf(x):
    return L.from_man_exp(-x[1] if x[0] else x[1], x[2]) if x[1] else L.fzero


def cross_diff(r, a, b, got):
    if r == "inf" or isinstance(a, str) or isinstance(b, str):
        return
    if a[1] and b[1] and abs(top_e(a[1], a[2]) - top_e(b[1], b[2])) > 100000:
        return
    if r[0] and abs(r[1]) > 1 << 40:
        return
    d = L.mpf_abs(L.mpf_sub(mpf(a), mpf(b)))
    if abs(r[1]) < 1 << 40 or d == L.fzero:
        w = ordname(L.mpf_cmp(L.from_man_exp(r[0], r[1] - 64) if r[0] else L.fzero, d))
        if w != got:
            raise Bad("mpmath disagrees on D %r %r %r: %s vs %s" % (r, a, b, w, got))


def up_lo(x):
    if x in ("nan", "+inf", "-inf"):
        return "inf", (0, 0)
    s, m, e = x
    if m == 0:
        return (0, 0), (0, 0)
    b = m.bit_length()
    if b <= 64:
        v = (m << (64 - b), e + b)
        up, lo = v, v
    else:
        lo_m = m >> (b - 64)
        lo = (lo_m, e + b)
        up = (lo_m + 1, e + b) if lo_m + 1 < (1 << 64) else (TOP, e + b + 1)
    # containment and tightness, exact
    xv = (0, m, e)
    if sum_sign([mag_val(lo), (1, m, e)]) > 0 or sum_sign([mag_val(up), (1, m, e)]) < 0:
        raise Bad("up/lo containment %r" % (x,))
    if up != lo and not (up[1] == lo[1] and up[0] - lo[0] == 1 or up == (TOP, lo[1] + 1) and lo[0] == (1 << 64) - 1):
        raise Bad("up/lo not one ulp %r" % (x,))
    for mode, want in (("f", lo), ("c", up)):
        v = L.normalize(0, m, e, b, 64, mode)
        if v[1] << (64 - v[3]) != want[0] or v[2] + v[3] != want[1]:
            raise Bad("mpmath disagrees on %s of %r" % (mode, x))
    return up, lo


def to_exact(r):
    if r == "inf":
        return "err:inf"
    if r == (0, 0):
        return "+0"
    mant, exp = r
    if exp > EXPONENT_MAX:
        return "err:max"
    d = EXPONENT_MIN - exp
    if d > 0 and not (d < 64 and (mant & ((1 << d) - 1)) == 0):
        return "err:min"
    return enum(norm(0, mant, exp - 64))


def enum(x):
    if isinstance(x, str):
        return x
    s, m, e = x
    if m == 0:
        return "-0" if s else "+0"
    return "%s%x^%d" % ("-" if s else "+", m, e)


def emag(r):
    if r == "inf":
        return "inf"
    if r == (0, 0):
        return "0"
    if not (r[0] >> 63 == 1 and r[0] < 1 << 64 and I64_MIN <= r[1] <= I64_MAX):
        raise Bad("not a canonical Mag %r" % (r,))
    return "%016x:%d" % r


def xnum(s, m, e):
    x = norm(s, m, e)
    if x[1] and not (EXPONENT_MIN <= top_e(x[1], x[2]) <= EXPONENT_MAX):
        raise Bad("ExactNum out of zf range %r" % (x,))
    return x


def rand_x(rng):
    k = rng.next() % 8
    bits = (64, 1, 63, 65, 128, 200, 1000, 4096)[k]
    m = 0
    while m.bit_length() < bits:
        m = (m << 64) | rng.next()
    m = (m >> max(0, m.bit_length() - bits)) | 1 | (1 << (bits - 1))
    ec = rng.next() % 4
    te = (rng.next() % 4096) - 2048 if ec < 2 else (EXPONENT_MAX - (rng.next() % 512) if ec == 2 else EXPONENT_MIN + (rng.next() % 512) + bits)
    te = min(max(te, EXPONENT_MIN), EXPONENT_MAX)
    return xnum(rng.next() & 1, m, te - bits)


def near_mag(rng, x):
    """A Mag at, just below or just above |x| (ties and last-bit differences)."""
    up, lo = up_lo(x)
    c = rng.next() % 4
    if c == 0 or up == (0, 0):
        return lo if lo != (0, 0) else (TOP, 1)
    if c == 1:
        return up
    if c == 2 and lo[0] > TOP:
        return (lo[0] - 1, lo[1])
    return (up[0] + 1, up[1]) if up[0] + 1 < 1 << 64 else (TOP, up[1] + 1)


def fixed_x():
    out = ["nan", "+inf", "-inf", (0, 0, 0), (1, 0, 0)]
    for s in (0, 1):
        out += [xnum(s, 1, 0), xnum(s, (1 << 64) - 1, EXPONENT_MAX - 64), xnum(s, 1, EXPONENT_MIN - 1),
                                xnum(s, (1 << 65) - 1, EXPONENT_MAX - 65), xnum(s, (1 << 64) + 1, EXPONENT_MIN - 1)]
    return out


def huge_d_rows():
    """The 100000-bit values, in a few D rows only: a = 1 +- 2^-100000 against 1 and against each other, the radius
    1 (gt), 2^-100000 (Mag(TOP, -99999); eq against 1) and 2^-100001 (lt)."""
    one = xnum(0, 1, 0)
    hp, hm = xnum(0, (1 << 100000) + 1, -100000), xnum(0, (1 << 100000) - 1, -100000)
    for a, b in ((hp, one), (one, hm), (hp, hm), (xnum(1, hp[1], hp[2]), xnum(1, 1, 0))):
        for r in ((TOP, 1), (TOP, -99999), (TOP, -100000)):
            yield r, a, b


def edge_d_rows():
    """harper_s2_diff_edge.rs, E = EXPONENT_MAX: a = 2^E (1 - 2^-200), s = 2^(E-130); (r, a, b, hand truth)."""
    e = EXPONENT_MAX
    a = xnum(0, (1 << 200) - 1, e - 200)
    s = xnum(0, 1, e - 130)
    na, ns = xnum(1, a[1], a[2]), xnum(1, s[1], s[2])
    for r, want in (((TOP, e + 1), ("lt", "lt", "gt")), ((TOP + 1, e + 1), ("gt", "gt", "gt")),
                    (((1 << 64) - 1, e + 1), ("gt", "gt", "gt")), (((1 << 64) - 1, e), ("lt", "lt", None))):
        for (x, y), w in zip(((a, ns), (na, s), (a, s)), want):
            if w is not None:
                yield r, x, y, w


FIXED_MAGS = ["inf", (0, 0), (TOP, 1), (TOP | 1, 1), ((1 << 64) - 1, 0), ((1 << 64) - 1, I64_MAX), (TOP, I64_MIN),
              (TOP, EXPONENT_MAX), (TOP, EXPONENT_MAX + 1), (TOP, EXPONENT_MIN), (TOP, EXPONENT_MIN - 1),
              (TOP, EXPONENT_MIN - 64), (TOP | 1, EXPONENT_MIN - 1)]


def rows(limit):
    rng = goldlib.SplitMix64(goldlib.stream_seed(SEED_SUITE, 1))
    n = 0
    fx = fixed_x()
    for x in fx:
        up, lo = up_lo(x)
        yield "U %s %s %s" % (enum(x), emag(up), emag(lo))
        for r in FIXED_MAGS:
            yield "N %s %s %s" % (emag(r), enum(x), cmp_num(r, x))
    for r in FIXED_MAGS:
        yield "T %s %s" % (emag(r), to_exact(r))
    for a in fx:
        for b in fx[:9]:
            for r in ((0, 0), (TOP, 1), "inf"):
                got = cmp_diff(r, a, b)
                cross_diff(r, a, b, got)
                yield "D %s %s %s %s" % (emag(r), enum(a), enum(b), got)
    for r, a, b in huge_d_rows():
        got = cmp_diff(r, a, b)
        cross_diff(r, a, b, got)
        yield "D %s %s %s %s" % (emag(r), enum(a), enum(b), got)
    for r, a, b, want in edge_d_rows():
        got = cmp_diff(r, a, b)
        if got != want:
            raise Bad("edge row %r %r %r: oracle %s, hand truth %s" % (r, a, b, got, want))
        yield "D %s %s %s %s" % (emag(r), enum(a), enum(b), got)
    while n < limit:
        n += 1
        x = rand_x(rng)
        up, lo = up_lo(x)
        yield "U %s %s %s" % (enum(x), emag(up), emag(lo))
        r = near_mag(rng, x)
        yield "N %s %s %s" % (emag(r), enum(x), cmp_num(r, x))
        # cancellation: b = a +- tiny or a mirror; radius near |a - b|
        a = x
        t = rng.next() % 3
        if t == 0:
            b = xnum(a[0], a[1] + (1 << (rng.next() % 64)), a[2]) if a[1] > 1 else a
        elif t == 1:
            b = rand_x(rng)
        else:
            b = xnum(1 - a[0], a[1], a[2])
        dt = sum_sign([a, (1 - b[0], b[1], b[2])])
        r2 = (TOP, 1) if dt == 0 else near_mag(rng, xnum(0, abs(a[1] * (1 << max(a[2] - b[2], 0)) - (-1 if a[0] != b[0] else 1) * b[1] * (1 << max(b[2] - a[2], 0))), min(a[2], b[2]))) if abs(top_e(a[1], a[2]) - top_e(b[1], b[2])) < 100000 else (TOP, top_e(a[1], a[2]) + (rng.next() % 3) - 1)
        got = cmp_diff(r2, a, b)
        cross_diff(r2, a, b, got)
        yield "D %s %s %s %s" % (emag(r2), enum(a), enum(b), got)
        if n % 50 == 0:
            yield "T %s %s" % (emag(r2), to_exact(r2))


def write(path, limit):
    count = 0
    with open(path, "wb") as fh:
        fh.write(("# zenith-bound S2 cmp/conversion rows, scripts/gen_bridge_cmp_goldens.py v%d (see its docstring)\n" % GENERATOR_VERSION).encode())
        for row in rows(limit):
            fh.write((row + "\n").encode("ascii"))
            count += 1
    print("gen_bridge_cmp_goldens: %d rows, max RSS %d KiB" % (count, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))


def main(argv):
    goldlib.require_mpmath()
    args = [a for a in argv if a != "--skip-heavy"]
    if len(args) == 2 and args[0] == "--check":
        os.makedirs(args[1], exist_ok=True)
        write(os.path.join(args[1], GOLDEN_NAME), RANDOM_PER_KIND)
        return 1 if goldlib.check_against(args[1], GOLDEN_DIR, [GOLDEN_NAME]) else 0
    if len(args) == 3 and args[0] == "--sample":
        os.makedirs(args[1], exist_ok=True)
        write(os.path.join(args[1], GOLDEN_NAME), int(args[2]))
        return 0
    if args == []:
        write(os.path.join(GOLDEN_DIR, GOLDEN_NAME), RANDOM_PER_KIND)
        return 0
    raise SystemExit("usage: scripts/gen_bridge_cmp_goldens.py [--check DIR [--skip-heavy] | --sample DIR N]")


if __name__ == "__main__":
    try:
        sys.exit(main(sys.argv[1:]))
    except Bad as e:
        print("gen_bridge_cmp_goldens: %s" % e)
        sys.exit(3)
