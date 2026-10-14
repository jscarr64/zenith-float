#!/usr/bin/env python3
"""mpmath goldens for tests/oracle/goldens/specials.

Word-rounded (round_p) values, same canon as the step-2 generator. Inputs are
parsed at the destination precision with ToEven before the oracle call, matching
ExactNum::parse. Audit rows set "exact" so a full-word miss is red. Folded
fixtures keep the bit slack their existing tests already allow ("lost_max").
clausen_cl3 stays on the two-digit smoke (no exact flag); the harness timeout
is the hang check.

  python3 scripts/gen_specials_goldens.py
"""

from __future__ import annotations

import importlib.util
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "zenith-float-num" / "tests" / "oracle" / "goldens" / "specials"
SPEC = importlib.util.spec_from_file_location(
    "oracle_goldens", ROOT / "scripts" / "gen_oracle_goldens.py"
)
g = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(g)

FAST_P = [53, 113]
DEEP_P = [64, 128, 256]
FAST_RM = ["ToEven", "Up"]
ALL_RM = ["ToEven", "ToOdd", "Up", "Down", "ToZero", "FromZero", "None"]

ZIV = {
    "erf", "erfc", "gamma", "ln_gamma", "digamma", "gammainc", "gammainc_upper",
    "ei", "si", "ci", "li", "fresnel_s", "fresnel_c", "ai", "bi", "ai_prime",
    "bi_prime", "bessel_j", "bessel_j_nu", "bessel_y", "bessel_i", "bessel_k",
    "elliptic_k", "elliptic_e_complete", "elliptic_f", "elliptic_e",
    "elliptic_pi", "elliptic_pi_complete", "elliptic_f_phi", "elliptic_e_phi",
    "elliptic_pi_phi", "jacobi_am", "jacobi_sn", "jacobi_cn",
    "jacobi_dn", "jacobi_arcsn", "hypergeom_2f1", "betainc", "normal_cdf",
}

# op -> file. Order is the file list.
FILE_OF = {
    "erf": "erf", "erfc": "erf",
    "gamma": "gamma", "ln_gamma": "gamma", "digamma": "gamma",
    "gammainc": "gamma", "gammainc_upper": "gamma",
    "bessel_j": "bessel", "bessel_j_nu": "bessel", "bessel_y": "bessel",
    "bessel_i": "bessel", "bessel_k": "bessel",
    "elliptic_k": "elliptic", "elliptic_e_complete": "elliptic",
    "elliptic_f": "elliptic", "elliptic_e": "elliptic",
    "elliptic_pi": "elliptic", "elliptic_pi_complete": "elliptic",
    "elliptic_f_phi": "elliptic", "elliptic_e_phi": "elliptic",
    "elliptic_pi_phi": "elliptic",
    "jacobi_am": "jacobi", "jacobi_sn": "jacobi", "jacobi_cn": "jacobi",
    "jacobi_dn": "jacobi", "jacobi_arcsn": "jacobi",
    "ai": "airy", "bi": "airy", "ai_prime": "airy", "bi_prime": "airy",
    "ei": "expint", "si": "expint", "ci": "expint", "li": "expint",
    "fresnel_s": "expint", "fresnel_c": "expint",
    "hypergeom_2f1": "hypergeom_2f1",
    "betainc": "gamma",
    "normal_cdf": "erf",
    "jacobi_cd": "jacobi",
    "jacobi_dc": "jacobi",
    "log": "other",
    "nth_root": "other",
    "hypergeom_1f1": "hypergeom_1f1",
    "hypergeom_0f1": "hypergeom_0f1",
    "lambert_w0": "lambert", "lambert_wm1": "lambert",
    "kelvin_ber": "kelvin", "kelvin_bei": "kelvin", "kelvin_ker": "kelvin",
    "kelvin_kei": "kelvin",
    "scorer_gi": "scorer", "scorer_hi": "scorer",
    "clausen_cl2": "clausen", "clausen_cl3": "clausen",
    "riemann_zeta": "zeta", "hurwitz_zeta": "zeta",
    "barnes_g": "barnes",
    "anger_j": "anger_weber_struve", "weber_e": "anger_weber_struve",
    "struve_h": "anger_weber_struve",
    "hermite_he": "orthopoly", "hermite_h": "orthopoly", "laguerre": "orthopoly",
    "gen_laguerre": "orthopoly", "chebyshev_t": "orthopoly", "chebyshev_u": "orthopoly",
    "gegenbauer": "orthopoly", "jacobi_p": "orthopoly", "legendre_p": "orthopoly",
    "polylog": "polylog", "polygamma": "polylog",
    "meijer_g": "meijer_fox_heun", "fox_h": "meijer_fox_heun",
    "heun_g": "meijer_fox_heun",
}


def mpf_arg(tok: str, p: int):
    import mpmath as mp

    if tok.startswith("n:"):
        return int(tok[2:])
    key = g.input_tuple(tok, p)
    if key[0] == "zero":
        return mp.mpf(0)
    return g.mpf_of(key)


# Functions that round at exactly p bits (2.0.0); the others round at round_p(p).
# The op of the row being evaluated, set by eval_row.
CURRENT_OP = [""]


def round_real(x, p: int, rm: str) -> dict:
    import mpmath as mp

    mp.mp.dps = max(80, g.round_p(p) // 2 + 40)
    if x == 0:
        return g.zero("-" if mp.sign(x) < 0 else "+", g.round_p(p))
    # Round 2: every real special function rounds at exactly p bits.
    sig = p if p >= 2 else None
    return g.round_mpf(mp.mpf(x), p, rm, sig=sig)


def nan() -> dict:
    return g.nan("InvalidArgument")


def call_mp(op: str, args: list, extra: dict):
    """Return an mpmath value, 'nan', or ('c', re, im)."""
    import mpmath as mp

    a = args
    if op == "erf":
        return mp.erf(a[0])
    if op == "erfc":
        return mp.erfc(a[0])
    if op == "gamma":
        # Γ(0) overflows to +∞ in this crate. Negative integers stay NaN.
        if a[0] == 0:
            return mp.inf
        return mp.gamma(a[0])
    if op == "ln_gamma":
        return mp.loggamma(a[0])
    if op == "digamma":
        return mp.digamma(a[0])
    if op == "gammainc":
        # lower γ(s, x); self is s, arg is x
        return mp.gammainc(a[0], 0, a[1], regularized=False)
    if op == "gammainc_upper":
        return mp.gammainc(a[0], a[1], mp.inf, regularized=False)
    if op == "bessel_j":
        return mp.besselj(a[1], a[0])
    if op == "bessel_j_nu":
        return mp.besselj(a[1], a[0])
    if op == "bessel_y":
        return mp.bessely(a[1], a[0])
    if op == "bessel_i":
        return mp.besseli(a[1], a[0])
    if op == "bessel_k":
        return mp.besselk(a[1], a[0])
    if op == "elliptic_k":
        return mp.ellipk(a[0])
    if op == "elliptic_e_complete":
        return mp.ellipe(a[0])
    if op == "elliptic_f":
        # args: x=sin(phi), m. mpmath wants the angle.
        if abs(a[0]) > 1:
            return "nan"
        return mp.ellipf(mp.asin(a[0]), a[1])
    if op == "elliptic_e":
        if abs(a[0]) > 1:
            return "nan"
        return mp.ellipe(mp.asin(a[0]), a[1])
    if op == "elliptic_pi":
        # self=n, x=sin phi, m
        if abs(a[1]) > 1:
            return "nan"
        return mp.ellippi(a[0], mp.asin(a[1]), a[2])
    if op == "elliptic_pi_complete":
        return mp.ellippi(a[0], mp.pi / 2, a[1])
    if op == "jacobi_am":
        # mpmath has no "am". am(u|m) = asin(sn(u|m)) on the principal range.
        if a[0] == 0:
            return mp.mpf(0)
        return mp.asin(mp.ellipfun("sn", a[0], a[1]))
    if op == "jacobi_sn":
        # ellipfun("sn", 0, m) is a tiny nonzero at high dps. sn(0|m) = 0.
        if a[0] == 0:
            return mp.mpf(0)
        return mp.ellipfun("sn", a[0], a[1])
    if op == "jacobi_cn":
        return mp.ellipfun("cn", a[0], a[1])
    if op == "jacobi_dn":
        return mp.ellipfun("dn", a[0], a[1])
    if op == "jacobi_cd":
        return mp.ellipfun("cd", a[0], a[1])
    if op == "jacobi_dc":
        return mp.ellipfun("dc", a[0], a[1])
    if op == "jacobi_arcsn":
        if abs(a[0]) > 1:
            return "nan"
        return mp.ellipf(mp.asin(a[0]), a[1])
    if op == "ai":
        return mp.airyai(a[0])
    if op == "bi":
        return mp.airybi(a[0])
    if op == "ai_prime":
        return mp.airyai(a[0], derivative=1)
    if op == "bi_prime":
        return mp.airybi(a[0], derivative=1)
    if op == "ei":
        return mp.ei(a[0])
    if op == "si":
        return mp.si(a[0])
    if op == "ci":
        return mp.ci(a[0])
    if op == "li":
        return mp.li(a[0])
    if op == "fresnel_s":
        return mp.fresnels(a[0])
    if op == "fresnel_c":
        return mp.fresnelc(a[0])
    if op == "hypergeom_2f1":
        return mp.hyp2f1(a[1], a[2], a[3], a[0])
    if op == "hypergeom_1f1":
        return mp.hyp1f1(a[1], a[2], a[0])
    if op == "hypergeom_0f1":
        return mp.hyp0f1(a[1], a[0])
    if op == "lambert_w0":
        return mp.lambertw(a[0], 0)
    if op == "lambert_wm1":
        return mp.lambertw(a[0], -1)
    if op == "kelvin_ber":
        return mp.ber(0, a[0])
    if op == "kelvin_bei":
        return mp.bei(0, a[0])
    if op == "kelvin_ker":
        return mp.ker(0, a[0])
    if op == "kelvin_kei":
        return mp.kei(0, a[0])
    if op == "scorer_gi":
        return mp.scorergi(a[0])
    if op == "scorer_hi":
        return mp.scorerhi(a[0])
    if op == "clausen_cl2":
        return mp.clsin(2, a[0])
    if op == "clausen_cl3":
        # Cosine Clausen: Cl_3(θ) = Σ cos(kθ)/k^3. clsin is the sine sum.
        if a[0] == 0:
            return mp.zeta(3)
        return mp.clcos(3, a[0])
    if op == "riemann_zeta":
        if a[0] not in (0, -1) and (a[0] < 2 or int(a[0]) != a[0]):
            return "nan"
        return mp.zeta(a[0])
    if op == "hurwitz_zeta":
        return mp.hurwitz(a[0], a[1])
    if op == "barnes_g":
        return mp.barnesg(a[0])
    if op == "anger_j":
        return mp.angerj(a[1], a[0])
    if op == "weber_e":
        return mp.webere(a[1], a[0])
    if op == "struve_h":
        return mp.struveh(a[1], a[0])
    if op == "hermite_h":
        return mp.hermite(a[1], a[0])
    if op == "hermite_he":
        # He_n(x) = 2^{-n/2} H_n(x/sqrt(2))
        n = a[1]
        return mp.power(2, -mp.mpf(n) / 2) * mp.hermite(n, a[0] / mp.sqrt(2))
    if op == "laguerre":
        return mp.laguerre(a[1], 0, a[0])
    if op == "gen_laguerre":
        return mp.laguerre(a[1], a[2], a[0])
    if op == "chebyshev_t":
        return mp.chebyt(a[1], a[0])
    if op == "chebyshev_u":
        return mp.chebyu(a[1], a[0])
    if op == "gegenbauer":
        return mp.gegenbauer(a[1], a[2], a[0])
    if op == "jacobi_p":
        return mp.jacobi(a[1], a[2], a[3], a[0])
    if op == "legendre_p":
        return mp.legenp(a[1], 0, a[0])
    if op == "polylog":
        return mp.polylog(a[1], a[0])
    if op == "polygamma":
        return mp.polygamma(a[1], a[0])
    if op == "meijer_g":
        an = extra.get("an", [])
        ap = extra.get("ap", [])
        bm = extra.get("bm", [])
        bq = extra.get("bq", [])
        return mp.meijerg([an, ap], [bm, bq], a[0])
    if op == "fox_h":
        # Unit scales: the value is the Meijer G of the a,b parameters.
        an = [p[0] for p in extra.get("an", [])]
        ap = [p[0] for p in extra.get("ap", [])]
        bm = [p[0] for p in extra.get("bm", [])]
        bq = [p[0] for p in extra.get("bq", [])]
        return mp.meijerg([an, ap], [bm, bq], a[0])
    if op == "heun_g":
        return mp.mpf(extra["decimal"])
    if op == "betainc":
        return mp.betainc(a[0], a[1], 0, a[2], regularized=True)
    if op == "normal_cdf":
        z = (a[0] - a[1]) / a[2]
        return mp.ncdf(z)
    if op == "log":
        return mp.log(a[0]) / mp.log(a[1])
    if op == "nth_root":
        return mp.power(a[0], mp.mpf(1) / a[1])
    if op.startswith("c."):
        return call_complex(op, a)
    raise ValueError(op)


def call_complex(op: str, a: list):
    import mpmath as mp

    def z(i: int):
        return mp.mpc(a[i], a[i + 1])

    name = op[2:]
    if name == "erf":
        return mp.erf(z(0))
    if name == "gamma":
        return mp.gamma(z(0))
    if name == "ln_gamma":
        # The crate returns ln(Γ(z)). mpmath loggamma differs by 2πi on some
        # branches; the folded fixtures follow ln(Γ).
        return mp.log(mp.gamma(z(0)))
    if name == "ei":
        return mp.ei(z(0))
    if name == "si":
        return mp.si(z(0))
    if name == "ci":
        return mp.ci(z(0))
    if name == "acosh":
        return mp.acosh(z(0))
    if name == "ai":
        return mp.airyai(z(0))
    if name == "bi":
        return mp.airybi(z(0))
    if name == "bessel_j":
        return mp.besselj(z(2), z(0))
    if name == "bessel_y":
        return mp.bessely(z(2), z(0))
    if name == "bessel_i":
        return mp.besseli(z(2), z(0))
    if name == "bessel_k":
        return mp.besselk(z(2), z(0))
    if name == "hyp2f1":
        return mp.hyp2f1(z(0), z(2), z(4), z(6))
    raise ValueError(op)


def pack(value, p: int, rm: str) -> dict:
    import mpmath as mp

    if value == "nan":
        return nan()
    if isinstance(value, tuple) and value[0] == "c":
        return value[1]
    mp.mp.dps = max(80, g.round_p(p) // 2 + 40)
    if isinstance(value, mp.mpc) or (hasattr(value, "imag") and value.imag != 0):
        w = mp.mpc(value)
        return {"k": "c", "re": round_real(w.real, p, rm), "im": round_real(w.imag, p, rm)}
    if mp.isnan(value):
        return nan()
    if mp.isinf(value):
        return g.inf("-" if value < 0 else "+")
    return round_real(value, p, rm)


# Recorded from mpmath (dps 200, nstr 180). phi(x) = arccos(2*sqrt(3)/(1+sqrt(3)+x) - 1).
# m = (2+sqrt(3))/4. The integral identity is checked in main().
HALF_PI = (
    "1.57079632679489661923132169163975144209858469968755291048747229615390820314310449931401741267105853399107404325664115332354692230477529111586267970406424055872514205135096926055278"
)
THREE_PI_OVER_4 = (
    "2.35619449019234492884698253745962716314787704953132936573120844423086230471465674897102611900658780098661106488496172998532038345716293667379401955609636083808771307702645389082917"
)
M_CUBIC = (
    "0.933012701892219323381861585376468091735701313452595157013951744862983254227200009270286546689312143918906535353851675757492486273749738119702913878023593412132023307975576395516994"
)
PHI0 = (
    "1.29953257306887582080339773503650160302735297255117556883708654212600457346289335365431960644782250061318071518211065710940897438485509079286935599598508141086119376863166755133786"
)
PHI1 = (
    "1.64265492275483293872763014546957250537481163728018102972221743087594415906527765906601436244541554006160551657304897484811247961936506448745328548521739496845797918755628708095511"
)
# 1 - 2^{-30}, an exact dyadic just below the complete-integral pole.
M_NEAR_ONE = "0.999999999068677425384521484375"
PHI_OPS = ("elliptic_f_phi", "elliptic_e_phi", "elliptic_pi_phi")


def nonfinite_tok(tok: str) -> bool:
    return tok.lower() in ("inf", "+inf", "-inf", "nan", "+nan", "-nan")


def round_half_even(q):
    """Nearest integer, ties to even. Matches `round(0, ToEven)`."""
    import mpmath as mp

    n = mp.floor(q)
    frac = q - n
    half = mp.mpf("0.5")
    if frac < half:
        return n
    if frac > half:
        return n + 1
    if mp.floor(n / 2) * 2 == n:
        return n
    return n + 1


def phi_not_real(val) -> bool:
    import mpmath as mp

    if val == "nan":
        return True
    if isinstance(val, mp.mpc):
        return val.imag != 0
    if hasattr(val, "imag") and val.imag != 0:
        return True
    return bool(mp.isnan(val))


def eval_phi(op: str, raw_args: list[str], p: int, rm: str) -> dict:
    """Angle-form goldens. mpmath takes φ; poles follow this crate, not ellipf(π/2, 1)."""
    import mpmath as mp

    if any(nonfinite_tok(tok) for tok in raw_args):
        return nan()
    mp.mp.dps = max(100, g.round_p(p) // 2 + 40)
    nums = [mpf_arg(tok, p) for tok in raw_args]
    if op == "elliptic_pi_phi":
        n, phi, m = nums
    else:
        n = None
        phi, m = nums
    # F(φ|0) = E(φ|0) = φ, and Π(0; φ|0) = φ. The harness parses φ at `p`
    # with ToEven, and the kernel rounds that value once. Every mode therefore
    # shares the ToEven encoding of the argument. mp.ellipf(φ, 0) is not that
    # encoding: a directed mode follows the integrator's rounding noise.
    identity = m == 0 and (
        op in ("elliptic_f_phi", "elliptic_e_phi") or (op == "elliptic_pi_phi" and n == 0)
    )
    if identity:
        tok = raw_args[1] if op == "elliptic_pi_phi" else raw_args[0]
        return g.finalize(g.parse_dec(tok), p, "ToEven")
    half = mpf_arg(HALF_PI, p)
    one = mp.mpf(1)
    sign = -1 if phi < 0 else 1
    if phi == half or phi == -half:
        if op == "elliptic_f_phi":
            if m > 1:
                return nan()
            if m == 1:
                return g.inf("-" if sign < 0 else "+")
            val = sign * mp.ellipk(m)
        elif op == "elliptic_e_phi":
            if m > 1:
                return nan()
            val = sign * mp.ellipe(m)
        else:
            # Π(n|m) in this crate is InvalidArgument for n >= 1 or m >= 1,
            # including the pole mpmath reports as +inf.
            if m >= 1 or n >= 1:
                return nan()
            val = sign * mp.ellippi(n, mp.pi / 2, m)
        return pack(val, p, rm)
    k = round_half_even(phi / mp.pi)
    if op == "elliptic_f_phi":
        if m == 1 and k != 0:
            return g.inf("-" if sign < 0 else "+")
        if m > 1 and k != 0:
            return nan()
        val = mp.ellipf(phi, m)
    elif op == "elliptic_e_phi":
        if m > 1 and k != 0:
            return nan()
        val = mp.ellipe(phi, m)
    else:
        if k != 0 and (m >= 1 or n >= 1):
            return nan()
        val = mp.ellippi(n, phi, m)
    if phi_not_real(val):
        return nan()
    if mp.isinf(val):
        return g.inf("-" if val < 0 else "+")
    return pack(val, p, rm)


def eval_row(op: str, raw_args: list[str], p: int, rm: str, extra: dict) -> dict:
    CURRENT_OP[0] = op
    if op in PHI_OPS:
        return eval_phi(op, raw_args, p, rm)
    nums = []
    for tok in raw_args:
        if tok.startswith("n:"):
            nums.append(int(tok[2:]))
        elif tok.startswith("c:"):
            re_s, im_s = tok[2:].split(",", 1)
            nums.append(mpf_arg(re_s, p))
            nums.append(mpf_arg(im_s, p))
        else:
            nums.append(mpf_arg(tok, p))
    # meijer parameter lists are decimals, parsed at p
    parsed_extra = {}
    for key in ("an", "ap", "bm", "bq"):
        if key in extra:
            parsed_extra[key] = []
            for item in extra[key]:
                if isinstance(item, list):
                    parsed_extra[key].append([mpf_arg(str(x), p) for x in item])
                else:
                    parsed_extra[key].append(mpf_arg(str(item), p))
    if "decimal" in extra:
        parsed_extra["decimal"] = extra["decimal"]
    try:
        value = call_mp(op, nums, parsed_extra)
    except (ValueError, ZeroDivisionError):
        value = "nan"
    return pack(value, p, rm)


# (op, args, exact, lost_max or None, note or None, extra or None, only_deep_toeven)
def base_cases() -> list[dict]:
    rows: list[dict] = []

    def add(op, args, exact=False, lost_max=None, note=None, extra=None, grid="family"):
        rows.append({
            "op": op,
            "args": [str(a) for a in args],
            "exact": exact,
            "lost_max": lost_max,
            "note": note,
            "extra": extra or {},
            "grid": grid,
        })

    # Ziv samples: 0 / 1 / inside / a modest large value. Limb-exact.
    add("erf", ["0"])
    add("erf", ["0.5"])
    add("erf", ["-1"])
    add("erfc", ["1"])
    add("gamma", ["5"])
    add("gamma", ["0.5"])
    add("gamma", ["0"], note="pole, +inf")
    add("gamma", ["-1"], note="pole, invalid")
    add("ln_gamma", ["3"])
    add("digamma", ["1"])
    add("gammainc", ["2", "1"])
    add("gammainc_upper", ["2", "1"])
    add("bessel_j", ["0", "n:0"])
    add("bessel_j", ["1", "n:0"])
    add("bessel_j_nu", ["1", "0.5"])
    add("bessel_y", ["1", "0"])
    add("bessel_i", ["1", "0"])
    add("bessel_k", ["1", "0"])
    add("elliptic_k", ["0.5"])
    add("elliptic_e_complete", ["0.5"])
    add("jacobi_sn", ["0", "0.5"])
    add("jacobi_cn", ["0", "0.5"])
    add("jacobi_dn", ["0", "0.5"])
    add("jacobi_am", ["0.4", "0.5"])
    add("ai", ["0"])
    add("bi", ["0"])
    add("ai_prime", ["0"])
    add("bi_prime", ["0"])
    add("ei", ["1"])
    add("si", ["1"])
    add("ci", ["1"])
    add("li", ["2"])
    add("fresnel_s", ["1"])
    add("fresnel_c", ["1"])
    add("hypergeom_2f1", ["0.25", "1", "1", "1"])
    add("normal_cdf", ["0", "0", "1"])

    # Elliptic argument convention. x = sin(phi), m = k^2. Limb-exact (Ziv).
    import mpmath as mp

    mp.mp.dps = 50
    sin05 = mp.nstr(mp.sin(mp.mpf("0.5")), 40, strip_zeros=False)
    sin2 = mp.nstr(mp.sin(2), 40, strip_zeros=False)
    nsin05 = mp.nstr(-mp.sin(mp.mpf("0.5")), 40, strip_zeros=False)
    sin04 = mp.nstr(mp.sin(mp.mpf("0.4")), 40, strip_zeros=False)
    add("elliptic_f", [sin05, "0.3"], note="F(phi=0.5|m=0.3); argument is sin(phi), m=k^2 not k")
    add("elliptic_f", [sin05, "0.25"], note="m=0.25 is k=0.5; m is k^2, not k")
    add("elliptic_f", [nsin05, "0.3"], note="odd in sin(phi); negative amplitude")
    add("elliptic_f", [sin2, "0.5"], note="sin(2) is not F(2|m); phi=2 is past pi/2")
    add("elliptic_f", ["2", "0.5"], note="|x|>1 and phi=2 cannot be passed; InvalidArgument")
    add("elliptic_e", [sin05, "0.3"], note="E uses sin(phi) and m=k^2, same as F")
    add("elliptic_e", ["2", "0.3"], note="|x|>1 is InvalidArgument for E as well")
    add("elliptic_pi", ["0.2", sin04, "0.3"], note="Pi(n=0.2; sin(0.4)|m=0.3); x is sin(phi), m=k^2")
    add("jacobi_arcsn", [sin05, "0.3"], note="arcsn is F(sin phi|m)")

    # Audit misses. exact => limb-identical to the word-rounded value (red today).
    audit = [
        ("lambert_w0", ["0.5"]),
        ("lambert_w0", ["-0.1"]),
        ("lambert_w0", ["2"]),
        ("lambert_wm1", ["-0.1"]),
        ("lambert_wm1", ["-0.2"]),
        ("riemann_zeta", ["2"]),
        ("riemann_zeta", ["-1"]),
        ("riemann_zeta", ["0"]),
        ("scorer_gi", ["0.3"]),
        ("scorer_gi", ["1"]),
        ("scorer_hi", ["0.3"]),
        ("scorer_hi", ["-0.4"]),
        ("kelvin_ber", ["0.5"]),
        ("kelvin_ber", ["1"]),
        ("kelvin_bei", ["0.5"]),
        ("kelvin_ker", ["0.5"]),
        ("kelvin_ker", ["1"]),
        ("kelvin_kei", ["0.5"]),
        ("clausen_cl2", ["0.5"]),
        ("clausen_cl2", ["1"]),
        ("barnes_g", ["2.5"]),
        ("barnes_g", ["1.5"]),
        ("hypergeom_1f1", ["0.3", "0.5", "1.5"]),
        ("hypergeom_0f1", ["0.2", "1.5"]),
        ("hermite_h", ["0.3", "n:4"]),
        ("hermite_he", ["0.3", "n:4"]),
        ("laguerre", ["0.3", "n:4"]),
        ("gen_laguerre", ["0.3", "n:3", "0.5"]),
        ("chebyshev_t", ["0.3", "n:5"]),
        ("chebyshev_u", ["0.3", "n:5"]),
        ("gegenbauer", ["0.3", "n:3", "0.5"]),
        ("jacobi_p", ["0.3", "n:3", "1", "1"]),
        ("polylog", ["0.5", "n:2"]),
        ("polygamma", ["1.5", "n:2"]),
        ("struve_h", ["0.3", "0"]),
        ("anger_j", ["0.3", "0"]),
        ("weber_e", ["0.3", "0"]),
    ]
    for op, args in audit:
        add(op, args, exact=True, note="audit")
    # Clausen Cl3 is the two-digit smoke, not a limb promise. 0.5 must be in both tiers.
    add("clausen_cl3", ["0.5"], exact=False, note="timeout plus two digits")
    add("clausen_cl3", ["0"], exact=False, note="Cl3(0)=zeta(3)")
    add("clausen_cl3", ["1"], exact=False, note="nonzero angle")
    # zeta domain that stays an error, and the two new points are in the audit list.
    add("riemann_zeta", ["0.5"], note="non-integer stays InvalidArgument")
    add("riemann_zeta", ["1"], note="pole stays InvalidArgument")
    add("riemann_zeta", ["4"], exact=True, note="even integer formula")
    add("hurwitz_zeta", ["3", "0.25"], exact=True)
    # Meijer / Fox exp(-1/2). Limb check: 12 stored ulps is still a miss at full word.
    add("meijer_g", ["0.5"], exact=True, extra={"an": [], "ap": [], "bm": ["0"], "bq": []}, note="G exp(-1/2)")
    add(
        "fox_h",
        ["0.5"],
        exact=True,
        extra={"an": [], "ap": [], "bm": [["0", "1"]], "bq": []},
        note="Fox unit-scale exp(-1/2)",
    )
    add(
        "heun_g",
        ["0.25", "2", "0.3", "0.4", "0.5", "0.8", "0.6"],
        lost_max=8,
        extra={"decimal": "1.0532709966795252110568734749696132676374130054028"},
        note="folded heun fixture",
        grid="fixture128",
    )
    # Barnes integers stay exact (factorial product) and should already be green.
    add("barnes_g", ["1"], exact=True, note="integer")
    add("barnes_g", ["2"], exact=True, note="integer")
    add("barnes_g", ["4"], exact=True, note="integer")
    add("legendre_p", ["0.5", "n:2"], exact=True)
    return rows


def phi_cases() -> list[dict]:
    """Angle φ, not sin φ. Appended after the sine-form rows so those lines stay put."""
    rows: list[dict] = []

    def add(op, args, note):
        rows.append({
            "op": op,
            "args": [str(a) for a in args],
            "exact": False,
            "lost_max": None,
            "note": note,
            "extra": {},
            "grid": "phi",
        })

    neg_half = "-" + HALF_PI
    # F(φ|m). ±π/2 is ±K. m = 1 at |φ| >= π/2 is signed infinity. m > 1 is
    # real only while the path stays inside |sin t| < 1/√m.
    f_rows = [
        (["0.5", "0.5"], "F(0.5|0.5); phi inside (0, pi/2)"),
        ([HALF_PI, "0.5"], "F(pi/2|0.5) = K(0.5)"),
        (["2", "0.5"], "F(2|0.5); phi past pi/2"),
        ([THREE_PI_OVER_4, "0.5"], "F(3pi/4|0.5) = 2K - F(pi/4)"),
        (["5", "0.3"], "F(5|0.3)"),
        (["100", "0.5"], "F(100|0.5); large phi, extended reduction"),
        (["-0.7", "0.5"], "odd in phi"),
        (["0", "0.5"], "F(0|m) = 0"),
        (["0.8", "0"], "m = 0 is the angle"),
        (["1", "0.99"], "m close to 1, phi inside (0, pi/2)"),
        (["0.6", M_CUBIC], "m = (2+sqrt(3))/4"),
        (["1.2", "-0.3"], "negative m"),
        (["2", M_NEAR_ONE], "m = 1-2^{-30}, phi past pi/2"),
        ([PHI0, M_CUBIC], "cubic integral phi(0); integral_0^1 dx/sqrt(x^3+1)"),
        ([PHI1, M_CUBIC], "cubic integral phi(1); phi past pi/2"),
        (["0.2", "1.5"], "m > 1, small phi, real"),
        (["1", "1.5"], "m > 1, path not real"),
        (["2", "1"], "F(2|1) = +inf"),
        ([HALF_PI, "1"], "F(pi/2|1) = +inf, same as K(1)"),
        ([HALF_PI, "1.5"], "F(pi/2|m>1) is not real"),
        ([neg_half, "0.5"], "F(-pi/2|0.5) = -K(0.5)"),
        (["inf", "0.5"], "infinity is InvalidArgument"),
    ]
    for args, note in f_rows:
        add("elliptic_f_phi", args, note)
    e_rows = [
        (["0.5", "0.5"], "E(0.5|0.5)"),
        ([HALF_PI, "0.5"], "E(pi/2|0.5) = E(0.5)"),
        (["2", "0.5"], "E(2|0.5)"),
        ([THREE_PI_OVER_4, "0.5"], "E(3pi/4|0.5)"),
        (["5", "0.3"], "E(5|0.3)"),
        (["100", "0.5"], "E(100|0.5)"),
        (["-0.7", "0.5"], "odd in phi"),
        (["0", "0.5"], "E(0|m) = 0"),
        (["0.8", "0"], "m = 0"),
        (["1", "0.99"], "m close to 1"),
        (["2", M_NEAR_ONE], "m = 1-2^{-30}"),
        (["1.2", "-0.3"], "negative m"),
        (["0.2", "1.5"], "m > 1, small phi, real"),
        (["1", "1.5"], "m > 1, path not real"),
        (["2", "1"], "E(2|1) stays real"),
        ([HALF_PI, "1"], "E(pi/2|1) = 1"),
        ([HALF_PI, "1.5"], "E(pi/2|m>1) is not real"),
        ([neg_half, "0.5"], "E(-pi/2|0.5) = -E(0.5)"),
        (["inf", "0.5"], "infinity is InvalidArgument"),
    ]
    for args, note in e_rows:
        add("elliptic_e_phi", args, note)
    # Π(n; φ|m). n on both sides of 1. n = 1 at sin^2 = 1 is InvalidArgument
    # (the complete integral rejects n >= 1; mpmath's +inf is not used).
    p_rows = [
        (["0.2", "0.4", "0.3"], "Pi(0.2; 0.4|0.3)"),
        (["0.2", HALF_PI, "0.3"], "Pi(0.2; pi/2|0.3) = Pi(0.2|0.3)"),
        (["0.2", "2", "0.3"], "Pi past pi/2"),
        (["0.2", "-0.4", "0.3"], "odd in phi"),
        (["0.2", "0", "0.5"], "Pi(n; 0|m) = 0"),
        (["0.5", "0.8", "0"], "m = 0"),
        (["-0.5", "2", "0.5"], "negative n, phi past pi/2"),
        (["2", "0.3", "0.5"], "n > 1, path still real"),
        (["1", "0.4", "0.3"], "n = 1 inside (0, pi/2), real"),
        (["1", HALF_PI, "0.3"], "n = 1 at pi/2 is the pole"),
        (["2", "1", "0.5"], "n > 1, path reaches 1/sqrt(n)"),
        (["0.2", "100", "0.5"], "large phi"),
        (["0.3", "0.5", "-0.25"], "negative m"),
        (["0.2", "0.4", "0.99"], "m close to 1"),
        (["0.2", "1.2", M_NEAR_ONE], "m = 1-2^{-30}, phi past pi/2"),
        (["0.2", "0.2", "1.5"], "m > 1, small phi, real"),
        (["0.2", "1", "1.5"], "m > 1, path not real"),
        (["0.2", neg_half, "0.3"], "Pi(n; -pi/2|m) = -Pi(n|m)"),
        (["0.2", "inf", "0.5"], "infinity is InvalidArgument"),
    ]
    for args, note in p_rows:
        add("elliptic_pi_phi", args, note)
    return rows


def shortcut_cases() -> list[dict]:
    """Identity shortcuts. Each one is the input, or one elementary call, rounded once."""
    import mpmath as mp

    mp.mp.dps = 200
    pi = mp.nstr(mp.pi, 180, strip_zeros=False)
    three_half = mp.nstr(3 * mp.pi / 2, 180, strip_zeros=False)
    rows: list[dict] = []

    def add(op, args, note):
        rows.append({
            "op": op,
            "args": [str(a) for a in args],
            "exact": False,
            "lost_max": None,
            "note": note,
            "extra": {},
            "grid": "phi",
        })

    # m = 0 is the amplitude, including the half-turn and full-turn angles.
    for op in ("elliptic_f_phi", "elliptic_e_phi"):
        add(op, ["-0.8", "0"], "shortcut m=0: the value is the input rounded once")
        add(op, [HALF_PI, "0"], "shortcut m=0 and phi=pi/2: the value is the input rounded once")
        add(op, [pi, "0"], "shortcut m=0 and phi=pi: the value is the input rounded once")
        add(op, [three_half, "0"], "shortcut m=0 and phi=3*pi/2: the value is the input rounded once")
        add(op, ["0", "0"], "shortcut phi=0 and m=0")
    # n = 0 is F. phi = 0 is 0. n = 0 and m = 0 is the amplitude.
    add("elliptic_pi_phi", ["0", "0.7", "0.4"], "shortcut n=0: Pi equals F")
    add("elliptic_pi_phi", ["0", HALF_PI, "0.25"], "shortcut n=0 and phi=pi/2: Pi equals K")
    add("elliptic_pi_phi", ["0", "0.8", "0"], "shortcut n=0 and m=0: Pi equals the input rounded once")
    add("elliptic_pi_phi", ["0", "0", "0.3"], "shortcut phi=0")
    # Sine form. m = 0 is one arcsin. n = 0 is F. x = 0 is 0.
    add("elliptic_f", ["0.5", "0"], "shortcut m=0: F(x|0) is arcsin(x) rounded once")
    add("elliptic_e", ["0.5", "0"], "shortcut m=0: E(x|0) is arcsin(x) rounded once")
    add("elliptic_pi", ["0", "0.5", "0.3"], "shortcut n=0: Pi equals F")
    add("elliptic_f", ["0", "0.3"], "shortcut x=0")
    add("elliptic_e", ["0", "0.3"], "shortcut x=0")
    add("elliptic_pi", ["0.2", "0", "0.3"], "shortcut x=0")
    return rows


def fixture_rows() -> list[dict]:
    rows = []
    root = ROOT / "zenith-float-num" / "tests"

    def add(op, args, p, lost_max, note, extra=None):
        tier_p = p if p in (53, 64, 113, 128, 256) else 256
        rows.append({
            "op": op,
            "args": args,
            "exact": False,
            "lost_max": lost_max,
            "note": note,
            "extra": extra or {},
            "grid": "fixture",
            "p": tier_p,
        })

    text = (root / "special_audit_mpmath.rs").read_text()
    for op, p, argbody, _re, _im in re.findall(
        r'\("([^"]+)",\s*(\d+),\s*&\[([^\]]*)\],\s*"([^"]*)",\s*"([^"]*)"\)',
        text,
    ):
        args = re.findall(r'"([^"]*)"', argbody)
        mapped = []
        for a in args:
            if "," in a and op.startswith("c."):
                mapped.append("c:" + a)
            else:
                mapped.append(a)
        # bessel_j in that file is the real-order form (x, nu), not integer n.
        use = op
        if op == "bessel_j":
            use = "bessel_j_nu"
        if op == "nth_root" and len(mapped) == 2:
            mapped = [mapped[0], "n:" + mapped[1]]
        add(use, mapped, int(p), 2, "fixture special_audit")

    text = (root / "erf_tail_mpmath.rs").read_text()
    for op, p, x, _ref in re.findall(
        r'\(\s*"([^"]+)",\s*(\d+),\s*"([^"]+)",\s*"([^"]+)"\s*,?\s*\)',
        text,
    ):
        # The checked-in tail test asks for min(p-2, 160) bits, and its
        # counter is one bit stricter than an exponent difference. Keep that
        # slack so a 50-digit reference limit is not a new limb demand.
        tier = int(p) if int(p) in (53, 64, 113, 128, 256) else 256
        need = min(tier - 2, 160)
        args = [x, "0", "1"] if op == "normal_cdf" else [x]
        add(op, args, int(p), tier - need + 1, "fixture erf_tail")

    text = (root / "bessel_large_arg_mpmath.rs").read_text()
    for kind, x, nu, _ref in re.findall(
        r"\(\s*'([jyi])',\s*\"([^\"]+)\",\s*\"([^\"]+)\",\s*\"([^\"]+)\"\s*,?\s*\)",
        text,
    ):
        op = {"j": "bessel_j_nu", "y": "bessel_y", "i": "bessel_i"}[kind]
        for p in (128, 256):
            add(op, [x, nu], p, 4, "fixture bessel_large")

    text = (root / "complex_bessel_mpmath.rs").read_text()
    for p, kind, z, nu, _re, _im in re.findall(
        r'\((\d+),\s*\'([a-z])\',\s*"([^"]+)",\s*"([^"]+)",\s*"([^"]+)",\s*"([^"]+)"\)',
        text,
    ):
        op = {"j": "c.bessel_j", "y": "c.bessel_y", "i": "c.bessel_i", "k": "c.bessel_k"}[kind]
        add(op, ["c:" + z, "c:" + nu], int(p), 4, "fixture complex_bessel")

    text = (root / "c2f1_log_mpmath.rs").read_text()
    for p, body, _re, _im in re.findall(
        r'\((\d+),\s*\[([^\]]+)\],\s*"([^"]+)",\s*"([^"]+)"\)', text
    ):
        parts = re.findall(r'"([^"]+)"', body)
        if len(parts) != 4:
            continue
        add("c.hyp2f1", ["c:" + s for s in parts], int(p), 4, "fixture c2f1")

    text = (root / "jacobi_p_mpmath.rs").read_text()
    for n, a, b, x, _ref in re.findall(
        r'\(\s*(\d+),\s*"([^"]+)",\s*"([^"]+)",\s*"([^"]+)",\s*"([^"]+)"\s*,?\s*\)',
        text,
    ):
        for p in (128, 256):
            add("jacobi_p", [x, f"n:{n}", a, b], p, 4, "fixture jacobi_p")

    # catalog_mpmath points at p=128, slack 8.
    import mpmath as mp

    mp.mp.dps = 50
    catalog = [
        ("scorer_gi", ["0.5"]),
        ("scorer_gi", ["2"]),
        ("scorer_gi", ["-2"]),
        ("scorer_hi", ["0.5"]),
        ("struve_h", ["1", "0"]),
        ("barnes_g", ["3.5"]),
        ("polygamma", ["1.5", "n:1"]),
        ("hurwitz_zeta", ["3", "0.25"]),
        ("lambert_w0", ["1.5"]),
        ("hypergeom_1f1", ["0.5", "1", "2"]),
        ("kelvin_ber", ["1.5"]),
        ("kelvin_bei", ["1.5"]),
        ("polylog", ["0.3", "n:2"]),
        ("clausen_cl2", [mp.nstr(mp.pi / 2, 40, strip_zeros=False)]),
    ]
    for op, args in catalog:
        add(op, args, 128, 8, "fixture catalog_mpmath")
    return rows


def expand(cases: list[dict]) -> dict[str, list[dict]]:
    files: dict[str, list[dict]] = {name: [] for name in sorted(set(FILE_OF.values()))}
    files["complex"] = []
    seq = {name: 0 for name in files}
    for case in cases:
        op = case["op"]
        file_name = "complex" if op.startswith("c.") else FILE_OF[op]
        ziv = op in ZIV or op.startswith("c.")
        grid = case["grid"]
        if grid == "fixture":
            slots = [(case["p"], "deep", "ToEven")]
        elif grid == "fixture128":
            slots = [(128, "deep", "ToEven")]
        elif grid == "phi":
            # 256-bit rows are fast so the default suite sees them. Their
            # budget is 10 s: a 256-bit Carlson row can exceed the 2 s fast cap.
            slots = []
            for rm_name in ALL_RM:
                slots.append((256, "fast", rm_name))
                slots.append((512, "deep", rm_name))
        else:
            slots = []
            for p in FAST_P:
                for rm in FAST_RM:
                    slots.append((p, "fast", rm))
            deep_modes = ALL_RM if ziv else ["ToEven"]
            for p in DEEP_P:
                for rm in deep_modes:
                    slots.append((p, "deep", rm))
        for p, tier, rm in slots:
            seq[file_name] += 1
            try:
                expect = eval_row(op, case["args"], p, rm, case["extra"])
            except Exception as exc:  # noqa: BLE001 — surface the case and stop
                raise RuntimeError(f"{op} {case['args']} p={p} {rm}: {exc}") from exc
            row = {
                "id": f"{op}-{p}-{rm}-{seq[file_name]}",
                "tier": tier,
                "op": op,
                "p": p,
                "rm": rm,
                "args": case["args"],
                "budget_s": 10 if grid == "phi" or tier == "deep" else 2,
                "expect": expect,
            }
            if case["exact"]:
                row["exact"] = True
            if case["lost_max"] is not None:
                row["lost_max"] = case["lost_max"]
            if case["note"]:
                row["note"] = case["note"]
            for key in ("an", "ap", "bm", "bq"):
                if key in case["extra"]:
                    row[key] = case["extra"][key]
            if "decimal" in case["extra"]:
                row["heun"] = case["args"][1:]
            files[file_name].append(row)
    return files


def main() -> None:
    import mpmath as mp

    mp.mp.dps = 25
    # Convention pin used by the elliptic notes: sin(2) is about 0.9093, and
    # F(2 | 0.5) is not F(asin(sin(2)) | 0.5).
    f_angle = mp.ellipf(2, mp.mpf("0.5"))
    f_sine = mp.ellipf(mp.asin(mp.sin(2)), mp.mpf("0.5"))
    if abs(f_angle - f_sine) < mp.mpf("0.1"):
        raise SystemExit("elliptic convention pin collapsed")
    z0 = eval_row("riemann_zeta", ["0"], 53, "ToEven", {})
    if z0.get("s") != "-" or z0.get("m") != ["8000000000000000"]:
        raise SystemExit(f"zeta(0) golden is not -1/2: {z0}")
    zm1 = eval_row("riemann_zeta", ["-1"], 53, "ToEven", {})
    if zm1.get("s") != "-":
        raise SystemExit(f"zeta(-1) sign: {zm1}")
    g0 = eval_row("gamma", ["0"], 53, "ToEven", {})
    if g0.get("k") != "inf" or g0.get("s") != "+":
        raise SystemExit(f"gamma(0) golden is not +inf: {g0}")
    gm1 = eval_row("gamma", ["-1"], 53, "ToEven", {})
    if gm1.get("k") != "nan":
        raise SystemExit(f"gamma(-1) golden is not nan: {gm1}")
    big = eval_row("elliptic_f", ["2", "0.5"], 53, "ToEven", {})
    if big.get("k") != "nan":
        raise SystemExit(f"elliptic_f(2) golden is not nan: {big}")
    cl30 = eval_row("clausen_cl3", ["0"], 53, "ToEven", {})
    if cl30.get("m") == ["0000000000000000"]:
        raise SystemExit(f"Cl3(0) golden collapsed to 0: {cl30}")
    sn0 = eval_row("jacobi_sn", ["0", "0.5"], 53, "ToEven", {})
    if sn0.get("m") != ["0000000000000000"]:
        raise SystemExit(f"sn(0) golden is not 0: {sn0}")
    mp.mp.dps = 80
    if abs(mp.mpf(HALF_PI) - mp.pi / 2) > mp.power(10, -60):
        raise SystemExit("recorded pi/2 drifted")
    if abs(mp.mpf(THREE_PI_OVER_4) - 3 * mp.pi / 4) > mp.power(10, -60):
        raise SystemExit("recorded 3pi/4 drifted")
    m_cubic = (2 + mp.sqrt(3)) / 4
    if abs(mp.mpf(M_CUBIC) - m_cubic) > mp.power(10, -60):
        raise SystemExit("recorded cubic m drifted")

    def cubic_phi(x):
        return mp.acos(2 * mp.sqrt(3) / (1 + mp.sqrt(3) + x) - 1)

    phi0 = cubic_phi(0)
    phi1 = cubic_phi(1)
    if abs(mp.mpf(PHI0) - phi0) > mp.power(10, -60):
        raise SystemExit("recorded phi0 drifted")
    if abs(mp.mpf(PHI1) - phi1) > mp.power(10, -60):
        raise SystemExit("recorded phi1 drifted")
    # ∫_0^1 dx/√(x^3+1) = 3^{-1/4} (F(φ1|m) - F(φ0|m)). The printed 20 digits
    # are a truncation; the next digits are 249…, so the residual is ~2.5e-21.
    integral = mp.power(3, -mp.mpf("0.25")) * (
        mp.ellipf(phi1, m_cubic) - mp.ellipf(phi0, m_cubic)
    )
    printed = mp.mpf("0.90960424263889577391")
    if abs(integral - printed) > mp.mpf("1e-20"):
        raise SystemExit(f"cubic integral is not 0.90960424263889577391: {integral}")
    half_inf = eval_row("elliptic_f_phi", [HALF_PI, "1"], 256, "ToEven", {})
    if half_inf.get("k") != "inf" or half_inf.get("s") != "+":
        raise SystemExit(f"F(pi/2|1) golden is not +inf: {half_inf}")
    pole = eval_row("elliptic_pi_phi", ["1", HALF_PI, "0.3"], 256, "ToEven", {})
    if pole.get("k") != "nan":
        raise SystemExit(f"Pi(1; pi/2|0.3) golden is not nan: {pole}")
    cases = base_cases() + fixture_rows() + phi_cases() + shortcut_cases()
    files = expand(cases)
    OUT.mkdir(parents=True, exist_ok=True)
    for name, rows in files.items():
        path = OUT / f"{name}.jsonl"
        path.write_text("".join(json.dumps(row, separators=(",", ":")) + "\n" for row in rows))
        print(f"{name}.jsonl {len(rows)}")


if __name__ == "__main__":
    main()
