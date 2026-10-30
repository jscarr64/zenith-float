#!/usr/bin/env python3
"""Stress reference rows for zenith-float-num/tests/xtarget_rem_pi_stress_red.rs (rem_pi cancellation).

Env: ~/team-env/mpmath-1.3.0 (mpmath 1.3.0, Python 3.12.14), run only via popslot.
Formula: rem_pi(x) = fmod(x, 2*pi) = x - trunc(x / (2*pi)) * 2*pi, in (-2pi, 2pi) (pi: OEIS A000796).
Inputs (all exact binary values):
  - 2^100;
  - 6381956970095103 * 2^797, the IEEE double closest to a multiple of pi/2 (J.-M. Muller, Elementary
    Functions: Algorithms and Implementation, 3rd ed., Ch. 11, worst case for double-precision reduction);
  - N_k = round(2*pi * 10^k) for k = 6, 12, 18 (nearest integer to a multiple of 2pi).
Rows: ToEven = nearest-even at p bits, Up = toward +inf at p bits, None = truncation at rp(p, 64) bits.
Each line: REF <input-binary> <p> <rm> <result-binary>. Output SHA256 is recorded in the test file header."""
import mpmath, sys
from mpmath import mp, mpf
from mpmath.libmp import normalize, round_nearest, round_ceiling, round_down
print("mpmath", mpmath.__version__, "python", sys.version.split()[0])
mp.prec = 4096
twopi = 2 * mp.pi
ins = [2 ** 100, 6381956970095103 * 2 ** 797] + [int(mp.nint(twopi * mpf(10) ** k)) for k in (6, 12, 18)]
def fmt(man, ex):
    b = bin(man)[2:]
    ip = len(b) + ex
    if ip <= 0:
        s = "0." + "0" * (-ip) + b
    elif ip >= len(b):
        s = b + "0" * (ip - len(b))
    else:
        s = b[:ip] + "." + b[ip:]
    return s.rstrip("0").rstrip(".") if "." in s else s
for x in ins:
    v = mpf(x) - mp.floor(mpf(x) / twopi) * twopi
    assert 0 < v < twopi
    for p in [53, 65, 100]:
        rp64 = -(-p // 64) * 64
        for name, prec, r in [("ToEven", p, round_nearest), ("Up", p, round_ceiling), ("None", rp64, round_down)]:
            s0, m, e, bc0 = v._mpf_
            _, man, ex, _ = normalize(s0, m, e, bc0, prec, r)
            print("REF", bin(x)[2:], p, name, fmt(man, ex))
