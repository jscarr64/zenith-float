#!/usr/bin/env python3
"""Reference rows for zenith-float-num/tests/xtarget_rem_pi_red.rs.

Env: ~/team-env/mpmath-1.3.0 (mpmath 1.3.0, Python 3.12.14), run only via popslot.
Formula: 10 is in (2pi, 4pi), so rem_pi(10) = 10 - 2*pi (pi: OEIS A000796; DLMF 3.12.1).
Rows: ToEven = nearest-even at p bits, Up = toward +inf at p bits, None = truncation at rp(p, 64) bits.
Output SHA256 is recorded in the test file header."""
import mpmath, sys
from mpmath import mp, mpf
from mpmath.libmp import normalize, round_nearest, round_ceiling, round_down
print("mpmath", mpmath.__version__, "python", sys.version.split()[0])
mp.prec = 4096
v = mpf(10) - 2 * mp.pi
def rnd(p, r):
    s0, m, e, bc0 = v._mpf_
    s, man, ex, bc = normalize(s0, m, e, bc0, p, r)
    return man, ex, bc
for p in [20, 33, 53, 65, 100, 129, 191]:
    rp64 = -(-p // 64) * 64
    for name, prec, r in [("ToEven", p, round_nearest), ("Up", p, round_ceiling), ("None", rp64, round_down)]:
        man, ex, bc = rnd(prec, r)
        b = bin(man)[2:]
        # value = man * 2^ex ; integer part has 2 bits (value in [2,4))
        ip = len(b) + ex
        s = b[:ip] + "." + b[ip:]
        print("REF", p, name, s.rstrip("0").rstrip("."))
