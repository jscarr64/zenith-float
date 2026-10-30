#!/usr/bin/env python3
"""Reference rows for zenith-float-num/tests/xtarget_complex_round_red.rs.

Env: ~/team-env/mpmath-1.3.0 (mpmath 1.3.0, Python 3.12.14), run only via popslot. mp.prec = 1024.
  asin(z), z = 1/2 + i/4: principal branch, DLMF 4.23.19 (mpmath.asin).
  I_nu(x), nu = 1/2 + i/4, x = 3/2 + i/4 (call x.bessel_i(&nu)): DLMF 10.25.2 power series (mpmath.besseli).
Rows: ToEven = nearest-even at p bits, Up = toward +inf at p bits, None = truncation at rp(p, 64) bits,
applied to the real and imaginary parts separately.
Each line: REF <fn> <p> <rm> <re-binary> <im-binary>. Output SHA256 is recorded in the test header."""
import mpmath, sys
from mpmath import mp, mpf, mpc
from mpmath.libmp import normalize, round_nearest, round_ceiling, round_down
print("mpmath", mpmath.__version__, "python", sys.version.split()[0])
mp.prec = 1024
z = mpc(mpf(1) / 2, mpf(1) / 4)
nu = mpc(mpf(3) / 2, mpf(1) / 4)
vals = [("asin", mp.asin(z)), ("bessel_i", mp.besseli(z, nu))]
def fmt(v, prec, r):
    s0, m, e, bc0 = v._mpf_
    s, man, ex, _ = normalize(s0, m, e, bc0, prec, r)
    b = bin(man)[2:]
    ip = len(b) + ex
    if ip <= 0:
        t = "0." + "0" * (-ip) + b
    elif ip >= len(b):
        t = b + "0" * (ip - len(b))
    else:
        t = b[:ip] + "." + b[ip:]
    if "." in t:
        t = t.rstrip("0").rstrip(".")
    return ("-" if s else "") + t
print("CHECK", mp.nstr(vals[0][1], 20), mp.nstr(vals[1][1], 20))
for name, v in vals:
    for p in [20, 65]:
        rp64 = -(-p // 64) * 64
        for rn, prec, r in [("ToEven", p, round_nearest), ("Up", p, round_ceiling), ("None", rp64, round_down)]:
            print("REF", name, p, rn, fmt(v.real, prec, r), fmt(v.imag, prec, r))
