#!/usr/bin/env python3
"""Reference rows for zenith-float-num/tests/sin_huge_exp_red.rs (sin(2^k) past the 2^16-bit
reduction cap in ops/util.rs trig_reduce_extra).

Env: ~/team-env/mpmath-1.3.0 (mpmath 1.3.0, Python 3.12.14), run only via popslot.
Formula: sin(2^k) = sin(2^k mod 2pi) (DLMF 4.14.1 periodicity; pi: OEIS A000796), evaluated with
mp.prec = k + 1024 so the reduction 2^k - n*2pi keeps over 1000 correct bits.
Rows: ToEven = nearest-even at p bits, Up = toward +inf at p bits, None = truncation at rp(p, 64).
Each line: REF <k> <p> <rm> <binary>. Output SHA256 is recorded in the test file header."""
import mpmath, sys
from mpmath import mp, mpf
from mpmath.libmp import normalize, round_nearest, round_ceiling, round_down
print("mpmath", mpmath.__version__, "python", sys.version.split()[0])
def fmt(v, prec, r):
    s0, m, e, bc0 = v._mpf_
    s, man, ex, _ = normalize(s0, m, e, bc0, prec, r)
    b = bin(man)[2:]
    ip = len(b) + ex
    assert ip <= 0
    t = ("0." + "0" * (-ip) + b).rstrip("0")
    return ("-" if s else "") + t
for k in [65600, 70000, 100000]:
    mp.prec = k + 1024
    x = mpf(2) ** k
    v = mp.sin(x - mp.floor(x / (2 * mp.pi)) * 2 * mp.pi)
    mp.prec = 1024
    print("CHECK", k, mp.nstr(+v, 25))
    for p in [53, 100]:
        rp64 = -(-p // 64) * 64
        for name, prec, r in [("ToEven", p, round_nearest), ("Up", p, round_ceiling), ("None", rp64, round_down)]:
            print("REF", k, p, name, fmt(v, prec, r))
