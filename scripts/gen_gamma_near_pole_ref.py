#!/usr/bin/env python3
"""Reference rows for zenith-float-num/tests/i686_gamma_pole_red.rs host rows (k = 256, 320).

Env: ~/team-env/mpmath-1.3.0 (mpmath 1.3.0, Python 3.12.14), run only via popslot.
x = -(1 - 2^-k). mpmath.gamma at mp.prec = 2048, cross-checked against the reflection formula
Gamma(x) = pi / (sin(pi x) Gamma(1 - x)) (DLMF 5.5.3), agreeing to well over 1024 bits.
Rows: ToEven = nearest-even at p = 256. Each line: REF <k> <binary>. SHA256 in the test header."""
import mpmath, sys
from mpmath import mp, mpf
from mpmath.libmp import normalize, round_nearest
print("mpmath", mpmath.__version__, "python", sys.version.split()[0])
mp.prec = 2048
for k in [256, 320]:
    x = -(1 - mpf(2) ** -k)
    g = mp.gamma(x)
    r = mp.pi / (mp.sin(mp.pi * x) * mp.gamma(1 - x))
    assert abs(g - r) < abs(g) * mpf(2) ** -1500
    s0, m, e, bc0 = g._mpf_
    s, man, ex, _ = normalize(s0, m, e, bc0, 256, round_nearest)
    b = bin(man)[2:]
    ip = len(b) + ex
    assert ip > len(b) or ip > 0
    t = b + "0" * (ip - len(b)) if ip >= len(b) else b[:ip] + "." + b[ip:]
    t = t.rstrip("0").rstrip(".") if "." in t else t
    print("CHECK", k, mp.nstr(g, 20))
    print("REF", k, ("-" if s else "") + t)
