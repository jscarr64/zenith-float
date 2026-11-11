#!/usr/bin/env python3
"""
Consumer: tests/harness/slice_bridge/harper_s2_wide_goldens.rs (every row asserted exactly). Listed in
tests/goldens/slices/bridge/GENERATORS and SHA256SUMS.
Generates tests/goldens/slices/bridge/bridge_wide.txt: S2 rows above 1028 bits (16.S2.1, 16.S2.2, 16.S2.3, 16.S2.4) in
the bridge_cmp.txt token and row format ('N R X ORD', 'D R A B ORD', 'U X UP LO'), truth from the same exact-integer
oracle (gen_bridge_cmp_goldens: sum_sign, up_lo with its mpmath normalize cross-check, cross_diff with mpmath mpf_sub/mpf_cmp).
VALUES. Seeded random P-bit mantissas (top and bottom bit set; goldlib.SplitMix64 from stream_seed(301, P)) at
P in {2048, 10_000, 100_000}. Top exponents: cosmology binades of 1e+-61, 1e+-120, 1e+-300 (+-203, +-399, +-997, jittered by
+-3), the unit binade, EXPONENT_MAX - j, EXPONENT_MIN + P + j, and uniform over [EXPONENT_MIN + P, EXPONENT_MAX]. At P 2048
and 10_000 also 10^k for k in +-61, +-120, +-300 rounded down and up to P bits.
PAIRS (D rows), per value a: (0) b = a with a random 64-bit word added at a random bit j < P (cancellation to about
P - j bits); (1) b = a cut to its top P/2 bits (|a - b| is the low half); (2) b = -a mirror; (3) b a fresh random value.
The radius is a Mag at, one ulp below or one ulp above |a - b| (or the top binade of the larger operand when the gap is
>= 100_000 binades). EXCLUDED: D rows whose |a - b| lies above EXPONENT_MAX (the known 16.S2.2 red,
harper_s2_diff_overflow.rs); they are counted, not written.
STREAMING: one row at a time. SIZE: about 1 MB (the 100_000-bit values appear in a handful of rows).
"""
import os
import resource
import sys

import goldlib
import gen_bridge_cmp_goldens as B

GENERATOR_VERSION = 1
GOLDEN_DIR = B.GOLDEN_DIR
GOLDEN_NAME = "bridge_wide.txt"
SEED_SUITE = 301
PLAN = ((2048, 60), (10_000, 16), (100_000, 3))
COSMO = (203, -203, 399, -399, 997, -997)
DECIMAL_K = (61, -61, 120, -120, 300, -300)
SKIPPED = [0]


def wide_mant(rng, p):
    m = 0
    while m.bit_length() < p:
        m = (m << 64) | rng.next()
    return (m >> (m.bit_length() - p)) | 1 | (1 << (p - 1))


def top_exp(rng, p):
    c = rng.next() % 5
    if c == 0:
        return COSMO[rng.next() % len(COSMO)] + (rng.next() % 7) - 3
    if c == 1:
        return (rng.next() % 9) - 4
    if c == 2:
        return B.EXPONENT_MAX - (rng.next() % 256)
    if c == 3:
        return B.EXPONENT_MIN + p + (rng.next() % 256)
    span = B.EXPONENT_MAX - (B.EXPONENT_MIN + p)
    return B.EXPONENT_MIN + p + ((rng.next() << 64 | rng.next()) % span)


def rand_wide(rng, p):
    m = wide_mant(rng, p)
    return B.xnum(rng.next() & 1, m, top_exp(rng, p) - p)


def decimal_x(k, p, up):
    """10^k rounded down (up) to p bits, as an odd-mantissa dyadic."""
    if k >= 0:
        n = 10 ** k
        sh = n.bit_length() - p
        m = n >> sh if sh > 0 else n << -sh
        if up and sh > 0 and (m << sh) != n:
            m += 1
        return B.xnum(0, m, sh)
    d = 10 ** -k
    sh = p + d.bit_length()
    q = (1 << sh) // d
    while q.bit_length() > p:
        q, sh = q >> 1, sh - 1   # floor stays floor: 2^sh / d >> 1 = floor(2^(sh-1) / d)
    if up and q * d != (1 << sh):
        q += 1
    return B.xnum(0, q, -sh)


def abs_diff(a, b):
    """|a - b| as (0, m, e) (exact), or None when the operands are >= 100_000 binades apart."""
    if a[1] and b[1] and abs(B.top_e(a[1], a[2]) - B.top_e(b[1], b[2])) >= 100_000:
        return None
    e = min(a[2], b[2])
    va = (-1 if a[0] else 1) * (a[1] << (a[2] - e))
    vb = (-1 if b[0] else 1) * (b[1] << (b[2] - e))
    return B.norm(0, abs(va - vb), e)


def partner(rng, a, p, kind):
    s, m, e = a
    if kind == 0:
        j = rng.next() % (p - 64)
        return B.xnum(s, m + (rng.next() << j), e)
    if kind == 1:
        h = m.bit_length() // 2
        return B.xnum(s, (m >> h) << h, e)
    if kind == 2:
        return B.xnum(1 - s, m, e)
    return rand_wide(rng, p)


def x_rows(rng, x, p):
    up, lo = B.up_lo(x)
    yield "U %s %s %s" % (B.enum(x), B.emag(up), B.emag(lo))
    for _ in range(2):
        r = B.near_mag(rng, x)
        yield "N %s %s %s" % (B.emag(r), B.enum(x), B.cmp_num(r, x))
    for kind in (0, 1 + rng.next() % 3):
        b = partner(rng, x, p, kind)
        d = abs_diff(x, b)
        if d is not None and d[1] and B.top_e(d[1], d[2]) > B.EXPONENT_MAX:
            SKIPPED[0] += 1
            continue
        if d is None:
            big = max(B.top_e(x[1], x[2]), B.top_e(b[1], b[2]))
            r = (B.TOP, big + (rng.next() % 3) - 1)
        elif d[1] == 0:
            r = (B.TOP, 1)
        else:
            r = B.near_mag(rng, d)
        got = B.cmp_diff(r, x, b)
        B.cross_diff(r, x, b, got)
        yield "D %s %s %s %s" % (B.emag(r), B.enum(x), B.enum(b), got)


def rows():
    for p, count in PLAN:
        rng = goldlib.SplitMix64(goldlib.stream_seed(SEED_SUITE, p))
        if p < 100_000:
            for k in DECIMAL_K:
                for up in (False, True):
                    yield from x_rows(rng, decimal_x(k, p, up), p)
        for _ in range(count):
            yield from x_rows(rng, rand_wide(rng, p), p)


def write(path, limit):
    count = 0
    with open(path, "wb") as fh:
        fh.write(("# zenith-float-bound S2 rows above 1028 bits, scripts/gen_bridge_wide_goldens.py v%d (see its docstring)\n" % GENERATOR_VERSION).encode())
        for row in rows():
            if limit is not None and count >= limit:
                break
            fh.write((row + "\n").encode("ascii"))
            count += 1
    print("gen_bridge_wide_goldens: %d rows, %d skipped (above EXPONENT_MAX), %d bytes, max RSS %d KiB"
          % (count, SKIPPED[0], os.path.getsize(path), resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))


def main(argv):
    goldlib.require_mpmath()
    args = [a for a in argv if a != "--skip-heavy"]
    if len(args) == 2 and args[0] == "--check":
        os.makedirs(args[1], exist_ok=True)
        write(os.path.join(args[1], GOLDEN_NAME), None)
        return 1 if goldlib.check_against(args[1], GOLDEN_DIR, [GOLDEN_NAME]) else 0
    if len(args) == 3 and args[0] == "--sample":
        os.makedirs(args[1], exist_ok=True)
        write(os.path.join(args[1], GOLDEN_NAME), int(args[2]))
        return 0
    if args == []:
        write(os.path.join(GOLDEN_DIR, GOLDEN_NAME), None)
        return 0
    raise SystemExit("usage: scripts/gen_bridge_wide_goldens.py [--check DIR [--skip-heavy] | --sample DIR N]")


if __name__ == "__main__":
    try:
        sys.exit(main(sys.argv[1:]))
    except B.Bad as e:
        print("gen_bridge_wide_goldens: %s" % e)
        sys.exit(3)
