#!/usr/bin/env python3
"""Generate step-1 JSONL goldens and cross-check them with MPFR.

The integer rounding below is independent of the library. scripts/mpfr_oracle.c
is the second oracle. A case is written only when the two agree, except for
parser-contract rows (malformed input), which state the 2.0.0 error tag.

  python3 scripts/gen_oracle_goldens.py
  cc -O2 -o /tmp/mpfr_oracle scripts/mpfr_oracle.c -lmpfr -lgmp
"""

from __future__ import annotations

import json
import math
import os
import subprocess
import sys
from pathlib import Path

EMIN = -2147483648
EMAX = 2147483647
WORD = 64
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "zenith-float-num" / "tests" / "oracle" / "goldens"
MODES = ["ToEven", "ToOdd", "Up", "Down", "ToZero", "FromZero", "None"]
FAST_P = [24, 53, 64, 113]
DEEP_P = [128, 256, 1024]


def round_p(p: int) -> int:
    if p <= 0:
        return 0
    return ((p + WORD - 1) // WORD) * WORD


def words_of(mant: int, prec: int) -> list[str]:
    nwords = prec // WORD
    mask = (1 << WORD) - 1
    return [f"{(mant >> (WORD * i)) & mask:016x}" for i in range(nwords)]


def zero(sign: str, prec: int) -> dict:
    return {"k": "f", "s": sign, "e": 0, "m": words_of(0, prec)}


def inf(sign: str) -> dict:
    return {"k": "inf", "s": sign}


def nan(err: str | None = "InvalidArgument") -> dict:
    return {"k": "nan", "err": err}


def floor_log2_ratio(n: int, d: int) -> int:
    if n <= 0 or d <= 0:
        raise ValueError("log2 of a non-positive")
    shift = n.bit_length() - d.bit_length()

    def ge(k: int) -> bool:
        if k >= 0:
            if k > n.bit_length() + 2:
                return False
            return n >= (d << k)
        sh = -k
        if sh > n.bit_length() + 2:
            return False
        return (n << sh) >= d

    if ge(shift):
        return shift + 1 if ge(shift + 1) else shift
    return shift - 1


def cmp_mul_pow2(n: int, exp: int, d: int) -> int:
    """Compare n * 2^exp with d. n and d are non-negative."""
    if n == 0:
        return 0 if d == 0 else -1
    if d == 0:
        return 1
    if exp >= 0:
        if exp >= d.bit_length() + 2:
            return 1
        v = n << exp
        return (v > d) - (v < d)
    sh = -exp
    if d.bit_length() + sh > n.bit_length() + 2:
        return -1
    rhs = d << sh
    return (n > rhs) - (n < rhs)


def round_ratio(numer: int, denom: int, rm: str, positive: bool) -> int:
    q, r = divmod(numer, denom)
    if r == 0:
        return q
    if rm in ("None", "ToZero"):
        return q
    if rm == "FromZero":
        return q + 1
    if rm == "Up":
        return q + 1 if positive else q
    if rm == "Down":
        return q if positive else q + 1
    twice = r * 2
    if twice < denom:
        return q
    if twice > denom:
        return q + 1
    if rm == "ToEven":
        return q if q % 2 == 0 else q + 1
    if rm == "ToOdd":
        return q if q % 2 == 1 else q + 1
    raise ValueError(rm)


def scale(n: int, d: int, exp2: int) -> tuple[int, int]:
    if exp2 >= 0:
        if exp2 > 200000:
            raise OverflowError(exp2)
        return n << exp2, d
    sh = -exp2
    if sh > 200000:
        raise OverflowError(exp2)
    return n, d << sh


def round_abs(sign: str, n: int, d: int, ez: int, prec: int, rm: str) -> dict:
    """Round a positive magnitude that is below 2^EMIN onto the subnormal grid."""
    positive = sign == "+"
    exp2 = ez - (EMIN - prec)
    half = cmp_mul_pow2(n, exp2 + 1, d)
    if half < 0:
        away = rm == "FromZero" or (rm == "Up" and positive) or (rm == "Down" and not positive)
        q = 1 if away else 0
    else:
        numer, denom = scale(n, d, exp2)
        q = round_ratio(numer, denom, rm, positive)
    if q <= 0:
        return zero(sign, prec)
    e = EMIN
    if q >= (1 << prec):
        if q > (1 << prec):
            raise AssertionError("subnormal rounding passed 2^EMIN")
        q = 1 << (prec - 1)
        e = EMIN + 1
        if e > EMAX:
            return inf(sign)
    return {"k": "f", "s": sign, "e": e, "m": words_of(q, prec)}


def round_rational(sign: str, n: int, d: int, ez: int, p: int, rm: str, sig_bits: int | None = None) -> dict:
    prec = round_p(p)
    sig = prec if sig_bits is None else sig_bits
    if sig <= 0 or sig > prec:
        raise ValueError(f"sig {sig} storage {prec}")
    if n == 0:
        return zero(sign, prec)
    g = math.gcd(n, d)
    n //= g
    d //= g
    fl = ez + floor_log2_ratio(n, d)
    if fl > EMAX:
        return inf(sign)
    if fl < EMIN:
        return round_abs(sign, n, d, ez, prec, rm)
    e = fl + 1
    numer, denom = scale(n, d, ez + sig - e)
    positive = sign == "+"
    mant = round_ratio(numer, denom, rm, positive)
    if mant >= (1 << sig):
        mant = 1 << (sig - 1)
        e += 1
    if e > EMAX:
        return inf(sign)
    if mant == 0 or e < EMIN:
        return round_abs(sign, n, d, ez, prec, rm)
    mant <<= prec - sig
    return {"k": "f", "s": sign, "e": e, "m": words_of(mant, prec)}


def parse_dec(text: str) -> tuple:
    s = text.strip().lower()
    sign = "+"
    if s.startswith("+"):
        s = s[1:]
    elif s.startswith("-"):
        sign = "-"
        s = s[1:]
    if s in ("inf",):
        return ("inf", sign)
    if s == "nan":
        return ("nan", None)
    exp = 0
    if "e" in s:
        mant, exp_s = s.split("e", 1)
        if exp_s in ("", "+", "-"):
            raise ValueError(f"bare exponent {text}")
        exp = int(exp_s)
    else:
        mant = s
    if "." in mant:
        ip, fp = mant.split(".", 1)
        frac = len(fp)
        digits = (ip + fp).lstrip("0") or "0"
    else:
        frac = 0
        digits = mant.lstrip("0") or "0"
    n = int(digits)
    if n == 0:
        return ("zero", sign)
    k = exp - frac
    if k >= 0:
        num, den, ez = n * pow(5, k), 1, k
    else:
        num, den, ez = n, pow(5, -k), k
    g = math.gcd(num, den)
    return ("f", sign, num // g, den // g, ez)


def token(tok: str, p: int) -> tuple:
    prec = round_p(p)
    if tok in ("+0", "0"):
        return ("zero", "+")
    if tok == "-0":
        return ("zero", "-")
    if tok in ("+inf", "inf"):
        return ("inf", "+")
    if tok == "-inf":
        return ("inf", "-")
    if tok in ("min", "+min"):
        return ("f", "+", 1, 1, EMIN - prec)
    if tok == "-min":
        return ("f", "-", 1, 1, EMIN - prec)
    if tok in ("3min", "+3min"):
        return ("f", "+", 3, 1, EMIN - prec)
    if tok == "-3min":
        return ("f", "-", 3, 1, EMIN - prec)
    if tok in ("max", "+max"):
        return ("f", "+", (1 << prec) - 1, 1, EMAX - prec)
    if tok == "-max":
        return ("f", "-", (1 << prec) - 1, 1, EMAX - prec)
    return parse_dec(tok)


def negate(v: tuple) -> tuple:
    if v[0] == "f":
        return ("f", "-" if v[1] == "+" else "+", v[2], v[3], v[4])
    if v[0] == "zero":
        return ("zero", "-" if v[1] == "+" else "+")
    if v[0] == "inf":
        return ("inf", "-" if v[1] == "+" else "+")
    return v


def zero_sign(s1: str, s2: str, rm: str) -> str:
    if s1 == s2:
        return s1
    return "-" if rm == "Down" else "+"


def add_raw(a: tuple, b: tuple) -> tuple:
    if a[0] == "nan" or b[0] == "nan":
        return ("nan",)
    if a[0] == "inf" and b[0] == "inf":
        return a if a[1] == b[1] else ("nan",)
    if a[0] == "inf":
        return a
    if b[0] == "inf":
        return b
    if a[0] == "zero" and b[0] == "zero":
        return ("zerosum", a[1], b[1])
    if a[0] == "zero":
        return b
    if b[0] == "zero":
        return a
    s1, n1, d1, e1 = a[1], a[2], a[3], a[4]
    s2, n2, d2, e2 = b[1], b[2], b[3], b[4]
    base = e1 if e1 <= e2 else e2
    sh1, sh2 = e1 - base, e2 - base
    if sh1 > 200000 or sh2 > 200000:
        raise OverflowError("exponent align")
    left = n1 * d2
    right = n2 * d1
    if sh1:
        left <<= sh1
    if sh2:
        right <<= sh2
    if s1 == s2:
        num, sign = left + right, s1
    elif left == right:
        return ("cancel",)
    elif left > right:
        num, sign = left - right, s1
    else:
        num, sign = right - left, s2
    den = d1 * d2
    g = math.gcd(num, den)
    return ("f", sign, num // g, den // g, base)


def mul_raw(a: tuple, b: tuple) -> tuple:
    if a[0] == "nan" or b[0] == "nan":
        return ("nan",)
    if a[0] == "zero" and b[0] == "inf" or b[0] == "zero" and a[0] == "inf":
        return ("nan",)
    if a[0] == "inf" or b[0] == "inf":
        sign = "+" if a[1] == b[1] else "-"
        return ("inf", sign)
    if a[0] == "zero" or b[0] == "zero":
        sign = "+" if a[1] == b[1] else "-"
        return ("zero", sign)
    sign = "+" if a[1] == b[1] else "-"
    num = a[2] * b[2]
    den = a[3] * b[3]
    g = math.gcd(num, den)
    return ("f", sign, num // g, den // g, a[4] + b[4])


def div_raw(a: tuple, b: tuple) -> tuple:
    if a[0] == "nan" or b[0] == "nan":
        return ("nan",)
    if b[0] == "zero":
        if a[0] == "zero":
            return ("nan",)
        if a[0] == "inf":
            return ("nan",)
        sign = "+" if a[1] == b[1] else "-"
        return ("inf", sign)
    if a[0] == "zero":
        sign = "+" if a[1] == b[1] else "-"
        return ("zero", sign)
    if a[0] == "inf" and b[0] == "inf":
        return ("nan",)
    if a[0] == "inf":
        sign = "+" if a[1] == b[1] else "-"
        return ("inf", sign)
    if b[0] == "inf":
        sign = "+" if a[1] == b[1] else "-"
        return ("zero", sign)
    sign = "+" if a[1] == b[1] else "-"
    num = a[2] * b[3]
    den = a[3] * b[2]
    g = math.gcd(num, den)
    return ("f", sign, num // g, den // g, a[4] - b[4])


def isqrt_exact(n: int) -> int | None:
    r = math.isqrt(n)
    return r if r * r == n else None


def sqrt_raw(a: tuple) -> tuple:
    if a[0] == "nan":
        return a
    if a[0] == "zero":
        return a
    if a[0] == "inf":
        return a if a[1] == "+" else ("nan",)
    if a[1] == "-":
        return ("nan",)
    n, d, ez = a[2], a[3], a[4]
    if ez % 2 != 0:
        n *= 2
        ez -= 1
    sn, sd = isqrt_exact(n), isqrt_exact(d)
    if sn is None or sd is None:
        raise ValueError("sqrt oracle only accepts exact squares in this step")
    return ("f", "+", sn, sd, ez // 2)


def icbrt(n: int) -> int | None:
    r = 0
    while (r + 1) ** 3 <= n:
        r += 1
    return r if r * r * r == n else None


def cbrt_raw(a: tuple) -> tuple:
    if a[0] in ("nan", "zero", "inf"):
        return a
    sign = a[1]
    n, d, ez = a[2], a[3], a[4]
    while ez % 3 != 0:
        n *= 2
        ez -= 1
    sn, sd = icbrt(n), icbrt(d)
    if sn is None or sd is None:
        raise ValueError("cbrt oracle only accepts exact cubes in this step")
    return ("f", sign, sn, sd, ez // 3)


def hypot_raw(a: tuple, b: tuple) -> tuple:
    return sqrt_raw(add_raw(mul_raw(a, a), mul_raw(b, b)))


def rem_raw(a: tuple, b: tuple) -> tuple:
    if a[0] == "nan" or b[0] == "nan":
        return ("nan",)
    if b[0] == "zero" or a[0] == "inf":
        return ("nan",)
    if a[0] == "zero":
        return a
    if b[0] == "inf":
        return a
    if a[3] != 1 or b[3] != 1 or a[4] != 0 or b[4] != 0:
        raise ValueError("rem oracle expects integers in this step")
    r = a[2] % b[2]
    if r == 0:
        return ("zero", a[1])
    return ("f", a[1], r, 1, 0)


def finalize(v: tuple, p: int, rm: str) -> dict:
    prec = round_p(p)
    if v[0] == "nan":
        return nan()
    if v[0] == "inf":
        return inf(v[1])
    if v[0] == "zero":
        return zero(v[1], prec)
    if v[0] == "cancel":
        return zero("-" if rm == "Down" else "+", prec)
    if v[0] == "zerosum":
        return zero(zero_sign(v[1], v[2], rm), prec)
    return round_rational(v[1], v[2], v[3], v[4], p, rm)


def expect_of(op: str, p: int, rm: str, args: list[str]) -> dict:
    if op == "parse":
        try:
            v = parse_dec(args[0])
        except ValueError:
            return nan()
        if v[0] == "nan":
            return nan(None)
        return finalize(v, p, rm)
    if op == "chop":
        v = parse_dec(args[0])
        return round_rational(v[1], v[2], v[3], v[4], p, rm, sig_bits=p)
    a = token(args[0], p)
    b = token(args[1], p) if len(args) > 1 else None
    c = token(args[2], p) if len(args) > 2 else None
    if op == "add":
        return finalize(add_raw(a, b), p, rm)
    if op == "sub":
        return finalize(add_raw(a, negate(b)), p, rm)
    if op == "mul":
        return finalize(mul_raw(a, b), p, rm)
    if op == "div":
        return finalize(div_raw(a, b), p, rm)
    if op == "rem":
        return finalize(rem_raw(a, b), p, rm)
    if op == "sqrt":
        return finalize(sqrt_raw(a), p, rm)
    if op == "cbrt":
        return finalize(cbrt_raw(a), p, rm)
    if op == "hypot":
        return finalize(hypot_raw(a, b), p, rm)
    if op == "fma":
        prod = mul_raw(a, b)
        if prod[0] == "nan":
            return nan()
        return finalize(add_raw(prod, c), p, rm)
    raise ValueError(op)


def value_key(st: dict, p: int) -> tuple:
    if st["k"] == "nan":
        return ("nan",)
    if st["k"] == "inf":
        return ("inf", st["s"])
    prec = round_p(p)
    mant = 0
    for i, word in enumerate(st["m"]):
        mant |= int(word, 16) << (WORD * i)
    if mant == 0:
        return ("zero", st["s"])
    return ("f", st["s"], mant, 1, st["e"] - prec)


def norm_dyadic(n: int, d: int, ez: int) -> tuple[int, int, int]:
    g = math.gcd(n, d)
    n //= g
    d //= g
    while d % 2 == 0:
        d //= 2
        ez -= 1
    if n:
        tz = (n & -n).bit_length() - 1
        n >>= tz
        ez += tz
    return n, d, ez


def same(st: dict, p: int, c_line: str) -> bool:
    parts = c_line.split()
    kind = value_key(st, p)
    if parts[0] == "nan":
        return kind[0] == "nan"
    if parts[0] == "inf":
        return kind == ("inf", parts[1])
    if parts[0] == "zero":
        return kind == ("zero", parts[1])
    if parts[0] != "finite" or kind[0] != "f":
        return False
    _k, sign, n, d, ez = kind
    if sign != parts[1]:
        return False
    cz, cez = int(parts[2]), int(parts[3])
    return norm_dyadic(n, d, ez) == norm_dyadic(cz, 1, cez)


def c_request(op: str, p: int, rm: str, args: list[str]) -> str:
    if op == "chop":
        return f"chop {p} {rm} {args[0]} {args[1]}"
    if op == "parse":
        return f"parse {p} {rm} {args[0]}"
    return " ".join([op, str(p), rm, *args])


def build_cases() -> list[dict]:
    rows: list[dict] = []

    def add(tier: str, op: str, p: int, rm: str, args: list[str], **extra: object) -> None:
        rows.append(
            {
                "tier": tier,
                "op": op,
                "p": p,
                "rm": rm,
                "args": args,
                "id": f"{op}-{p}-{rm}-{len(rows)}",
                **extra,
            }
        )

    for tier, precs in (("fast", FAST_P), ("deep", DEEP_P)):
        for p in precs:
            for rm in MODES:
                for src in ("1", "3", "1.5", "18446744073709551617"):
                    add(tier, "chop", p, rm, [src, "256"])
                for pair in (("1", "2"), ("1.5", "-0.5"), ("10", "-3"), ("0.25", "0.25")):
                    for op in ("add", "sub", "mul", "div"):
                        add(tier, op, p, rm, list(pair))
                add(tier, "sqrt", p, rm, ["4"])
                add(tier, "sqrt", p, rm, ["0.25"])
                add(tier, "cbrt", p, rm, ["8"])
                add(tier, "cbrt", p, rm, ["-8"])
                add(tier, "hypot", p, rm, ["3", "4"])
                add(tier, "hypot", p, rm, ["-3", "-4"])
                add(tier, "rem", p, rm, ["10", "3"])
                add(tier, "fma", p, rm, ["2", "3", "4"])
                for a, b in (
                    ("+0", "+0"),
                    ("+0", "-0"),
                    ("-0", "+0"),
                    ("-0", "-0"),
                ):
                    add(tier, "add", p, rm, [a, b])
                    add(tier, "sub", p, rm, [a, b])
                add(tier, "add", p, rm, ["1", "-1"])
                add(tier, "sub", p, rm, ["5", "5"])
                add(tier, "div", p, rm, ["min", "2"])
                add(tier, "div", p, rm, ["-min", "2"])
                add(tier, "div", p, rm, ["min", "4"])
                add(tier, "div", p, rm, ["3min", "4"])
                add(tier, "mul", p, rm, ["min", "min"])
                add(tier, "mul", p, rm, ["max", "2"])
                add(tier, "mul", p, rm, ["-max", "2"])
                add(tier, "add", p, rm, ["max", "max"])
                for a, b in (
                    ("0", "+inf"),
                    ("+inf", "0"),
                    ("+inf", "-inf"),
                    ("-inf", "+inf"),
                    ("+inf", "+inf"),
                    ("-inf", "-inf"),
                ):
                    add(tier, "add", p, rm, [a, b])
                add(tier, "sub", p, rm, ["+inf", "+inf"])
                add(tier, "sub", p, rm, ["+inf", "-inf"])
                add(tier, "div", p, rm, ["+inf", "+inf"])
                add(tier, "div", p, rm, ["-inf", "+inf"])
                add(tier, "div", p, rm, ["0", "0"])
                add(tier, "div", p, rm, ["1", "0"])
                add(tier, "div", p, rm, ["-1", "0"])
                add(tier, "fma", p, rm, ["0", "+inf", "1"])

    junk = [
        "",
        ".",
        "1.2e",
        "1e",
        "1e+",
        "1e-",
        "1E",
        "1.5foo",
        "0x10",
        "++1",
        "nanx",
        "infinity",
        "1.2.3",
        " 1",
        "1 ",
        "e10",
        "+-1",
        "1.2e+",
        "1.2e+3x",
        "+",
        "-",
    ]
    for tier, precs, modes in (
        ("fast", [24, 53, 113], ["ToEven", "Up", "ToZero"]),
        ("deep", [256], MODES),
    ):
        for p in precs:
            for rm in modes:
                for text in junk:
                    add(tier, "parse", p, rm, [text], contract=True)
                for text in ("0", "1", "-1", "0.1", "1.25e-2", "1E-1", "2.5", "-0", "1.", "1e-0"):
                    add(tier, "parse", p, rm, [text], roundtrip=True)
                for text in ("inf", "-inf", "+inf", "nan", "NaN", "+nan"):
                    add(tier, "parse", p, rm, [text])
    return rows


def load_oracle(rows: list[dict]) -> list[str]:
    cc = os.environ.get("CC", "cc")
    binary = Path("/tmp/mpfr_oracle")
    src = ROOT / "scripts" / "mpfr_oracle.c"
    subprocess.check_call(
        [cc, "-O2", "-o", str(binary), str(src), "-lmpfr", "-lgmp"],
        stdout=subprocess.DEVNULL,
    )
    payload = []
    index = []
    for i, row in enumerate(rows):
        if row.get("contract"):
            continue
        payload.append(c_request(row["op"], row["p"], row["rm"], row["args"]))
        index.append(i)
    proc = subprocess.run(
        [str(binary)],
        input="\n".join(payload) + "\n",
        text=True,
        capture_output=True,
        check=True,
    )
    lines = [ln for ln in proc.stdout.splitlines() if ln]
    if len(lines) != len(index):
        raise SystemExit(f"MPFR oracle returned {len(lines)} lines for {len(index)} cases")
    out = [""] * len(rows)
    for i, line in zip(index, lines):
        out[i] = line
    return out


def main() -> None:
    rows = build_cases()
    c_lines = load_oracle(rows)
    grouped: dict[str, list[dict]] = {
        "rounding": [],
        "arith": [],
        "edges": [],
        "parse": [],
    }
    for row, c_line in zip(rows, c_lines):
        if row.get("contract"):
            exp = nan()
        else:
            try:
                exp = expect_of(row["op"], row["p"], row["rm"], row["args"])
            except Exception as exc:
                raise SystemExit(f"{row} {exc}") from exc
        if not row.get("contract"):
            if not same(exp, row["p"], c_line):
                raise SystemExit(
                    f"oracle disagree {row['id']} {row['op']} {row['args']} "
                    f"rm={row['rm']} p={row['p']}\n  python {exp}\n  mpfr   {c_line}"
                )
        item = {
            "id": row["id"],
            "tier": row["tier"],
            "op": row["op"],
            "p": row["p"],
            "rm": row["rm"],
            "args": row["args"],
            "expect": exp,
        }
        if row["op"] == "parse":
            item["rdx"] = 10
        if row.get("roundtrip"):
            item["roundtrip"] = True
        if row["op"] == "chop":
            grouped["rounding"].append(item)
        elif row["op"] == "parse":
            grouped["parse"].append(item)
        elif row["args"][:1] in (["min"], ["-min"], ["3min"], ["max"], ["-max"]) or (
            row["op"] in ("add", "sub", "mul", "div", "fma")
            and any(a in ("+0", "-0", "0", "+inf", "-inf", "min", "-min", "3min", "max", "-max") for a in row["args"])
        ):
            grouped["edges"].append(item)
        else:
            grouped["arith"].append(item)
    OUT.mkdir(parents=True, exist_ok=True)
    for name, items in grouped.items():
        path = OUT / f"{name}.jsonl"
        path.write_text("".join(json.dumps(it, separators=(",", ":")) + "\n" for it in items))
        fast = sum(1 for it in items if it["tier"] == "fast")
        deep = sum(1 for it in items if it["tier"] == "deep")
        print(f"{name}.jsonl fast={fast} deep={deep}")


if __name__ == "__main__":
    main()
