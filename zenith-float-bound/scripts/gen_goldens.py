#!/home/jsc64/team-env/oracle-2026-10-30/bin/python
"""
Generates tests/goldens/mag_ops.txt; run with ~/team-env/oracle-2026-10-30/bin/python (goldlib.require_oracle checks the four pinned versions at import).
"""

import os
import sys
import goldlib  # oracle startup check, before mpmath is imported
import mpmath
from fractions import Fraction
from mpmath.libmp import from_man_exp, mpf_add, mpf_sub, mpf_mul, mpf_div, mpf_sqrt, fzero

if mpmath.__version__ != goldlib.MPMATH_VERSION:
    raise SystemExit("gen_goldens.py needs mpmath " + goldlib.MPMATH_VERSION + "; found " + mpmath.__version__)

I64_MIN = -2**63
I64_MAX = 2**63 - 1
TWO63 = 2**63
MASK64 = 2**64 - 1

RNDUP = 'c'
RNDDOWN = 'f'

def xorshift_mantissas(count, seed=0x5EED):
    state = seed | 1
    result = []
    for _ in range(count):
        x = state
        x ^= x >> 12
        x ^= (x << 25) & MASK64
        x ^= x >> 27
        state = x
        out = (x * 0x2545F4914F6CDD1D) & MASK64
        result.append(out | TWO63)
    return result

ALL_MANTS = [TWO63, TWO63 + 1, TWO63 + 2, 0xC000000000000000, 0xAAAAAAAAAAAAAAAB, MASK64 - 1, MASK64] + xorshift_mantissas(5)
BIN_MANTS = [TWO63, TWO63 + 1, 0xAAAAAAAAAAAAAAAB, MASK64]
EXPS = [-65, -64, -63, 0, 63, 64, 65, I64_MIN, I64_MIN + 1, I64_MAX - 1, I64_MAX]

BIN_SET = [(0, 0)] + [(e, m) for e in EXPS for m in BIN_MANTS]
UNARY_SET = [(0, 0)] + [(e, m) for e in EXPS for m in ALL_MANTS]

WIDE_A = [(e, m) for e in [0, I64_MAX] for m in ALL_MANTS]
WIDE_B = [(e, m) for e in [0, 64, I64_MIN] for m in ALL_MANTS]

SQUARES = []
for r in [0xB504F334, 0xB504F335, 0xC0000000, 0xDEADBEEF, 0xFFFFFFFE, 0xFFFFFFFF]:
    sq = r * r
    for d in [-1, 0, 1]:
        m = sq + d
        if TWO63 <= m <= MASK64:
            for e in [-64, 0, 2, 3, 64]:
                SQUARES.append((e, m))

def to_mpf(mag):
    exp, mant = mag
    return fzero if mant == 0 else from_man_exp(mant, exp - 64)

def from_mpf(x, up):
    sign, man, exp, bc = x
    if man == 0:
        return (0, 0)
    
    # Shift mantissa to 64-bit precision
    mant = man << (64 - bc)
    E = exp + bc
    
    if E > I64_MAX:
        if up:
            return "inf"
        else:
            return (I64_MAX, MASK64)
    elif E < I64_MIN:
        if up:
            return (I64_MIN, TWO63)
        else:
            return (0, 0)
    else:
        return (E, mant)

def mag_lt_or_eq(a, b):
    if a[1] == 0:
        return True
    if b[1] == 0:
        return False
    return a[0] < b[0] or (a[0] == b[0] and a[1] <= b[1])

def add_up(a, b):
    return from_mpf(mpf_add(to_mpf(a), to_mpf(b), 64, RNDUP), True)

def mul_up(a, b):
    return from_mpf(mpf_mul(to_mpf(a), to_mpf(b), 64, RNDUP), True)

def mul_lower(a, b):
    return from_mpf(mpf_mul(to_mpf(a), to_mpf(b), 64, RNDDOWN), False)

def sub_lower(a, b):
    if mag_lt_or_eq(a, b):
        return (0, 0)
    return from_mpf(mpf_sub(to_mpf(a), to_mpf(b), 64, RNDDOWN), False)

def div_up(a, b):
    return from_mpf(mpf_div(to_mpf(a), to_mpf(b), 64, RNDUP), True)

def div_lower(a, b):
    return from_mpf(mpf_div(to_mpf(a), to_mpf(b), 64, RNDDOWN), False)

def sqrt_up(a):
    return from_mpf(mpf_sqrt(to_mpf(a), 64, RNDUP), True)

def sqrt_lower(a):
    return from_mpf(mpf_sqrt(to_mpf(a), 64, RNDDOWN), False)

def value(m, base=0):
    """Exact value of a Mag, scaled by 2^-base so that huge exponents stay computable."""
    exp, mant = m
    if mant == 0:
        return Fraction(0)
    return Fraction(mant) * Fraction(2) ** (exp - 64 - base)

MAXMAG = (I64_MAX, MASK64)
MINPOS = (I64_MIN, TWO63)

def pred(m):
    exp, mant = m
    if mant > TWO63:
        return (exp, mant - 1)
    elif exp == I64_MIN:
        return None
    else:
        return (exp - 1, MASK64)

def succ(m):
    exp, mant = m
    if mant < MASK64:
        return (exp, mant + 1)
    elif exp == I64_MAX:
        return None
    else:
        return (exp + 1, TWO63)

skipped = 0
WINDOW = 3000


def window_base(op, mags):
    """Fractions of 2^(huge exponent) are not computable. Returns the scale base for a row that is cheap enough to check, else None.
    add_up and sub_lower: all finite nonzero exponents within WINDOW of each other (base = the smallest).
    mul, div and sqrt change exponents, so they are checked only when every exponent has absolute value <= WINDOW (base = 0)."""
    exps = [m[0] for m in mags if m != "inf" and m[1] != 0]
    if not exps:
        return 0
    if op in ("add_up", "sub_lower"):
        return min(exps) if max(exps) - min(exps) <= WINDOW else None
    return 0 if max(abs(e) for e in exps) <= WINDOW else None


def exact_value(op, args, base):
    """Exact result of the op in the scale 2^-base (base is 0 except for add_up and sub_lower, see window_base)."""
    if op == "add_up":
        return value(args[0], base) + value(args[1], base)
    if op in ("mul_up", "mul_lower"):
        return value(args[0], base) * value(args[1], base)
    if op == "sub_lower":
        d = value(args[0], base) - value(args[1], base)
        return d if d > 0 else Fraction(0)
    if op in ("div_up", "div_lower"):
        return value(args[0], base) / value(args[1], base)
    if op in ("sqrt_up", "sqrt_lower"):
        return value(args[0], base)
    raise SystemExit("unknown op " + op)


def fail(op, args, result, why):
    raise SystemExit("self-check failed: %s %s = %s: %s" % (op, args, result, why))


def self_check(op, args, result):
    """Independent exact-rational check of one row; True when checked, False when skipped (exponents too far apart, see window_base)."""
    global skipped
    up = op.endswith("_up")
    sqrt = op.startswith("sqrt")
    probe = args + [result]
    if result == "inf":
        probe = args + [MAXMAG]
    elif result == (0, 0) and exact_nonzero(op, args):
        probe = args + [MINPOS]
    base = window_base(op, probe)
    if base is None:
        skipped += 1
        return False
    exact = exact_value(op, args, base)
    # sqrt rows compare squares with the radicand, the other rows compare values
    def scaled(m):
        return value(m, base) ** 2 if sqrt else value(m, base)
    if exact == 0:
        if result != (0, 0):
            fail(op, args, result, "exact value is 0")
    elif up:
        if result == "inf":
            if not exact > value(MAXMAG, base):
                fail(op, args, result, "finite upper bound exists")
        else:
            if result == (0, 0):
                fail(op, args, result, "zero for a positive exact value")
            if scaled(result) < exact:
                fail(op, args, result, "not an upper bound")
            p = pred(result)
            if result != MINPOS and p is not None and scaled(p) >= exact:
                fail(op, args, result, "not the smallest upper bound")
    else:
        if result == (0, 0):
            if not exact < value(MINPOS, base):
                fail(op, args, result, "zero but a positive Mag fits below the exact value")
        else:
            if scaled(result) > exact:
                fail(op, args, result, "not a lower bound")
            q = succ(result)
            if result != MAXMAG and q is not None and scaled(q) <= exact:
                fail(op, args, result, "not the largest lower bound")
    return True


def exact_nonzero(op, args):
    """True unless the exact result is certainly zero (a zero operand, or sub_lower with a <= b), decided without building huge Fractions."""
    if op in ("mul_up", "mul_lower") and (args[0] == (0, 0) or args[1] == (0, 0)):
        return False
    if op in ("sqrt_up", "sqrt_lower", "div_up", "div_lower") and args[0] == (0, 0):
        return False
    if op == "add_up" and args[0] == (0, 0) and args[1] == (0, 0):
        return False
    if op == "sub_lower" and mag_lt_or_eq(args[0], args[1]):
        return False
    return True


def unique(seq):
    """Order-preserving de-duplication, so the output file is deterministic."""
    seen = set()
    out = []
    for item in seq:
        if item not in seen:
            seen.add(item)
            out.append(item)
    return out


def binary_pairs(op):
    pairs = [(a, b) for a in BIN_SET for b in BIN_SET]
    pairs += [(a, b) for a in WIDE_A for b in WIDE_B]
    pairs += [(a, b) for a in WIDE_B for b in WIDE_A]
    pairs = unique(pairs)
    if op in ("div_up", "div_lower"):
        pairs = [(a, b) for (a, b) in pairs if b != (0, 0)]
    return pairs


# L0a-1 remainder (TEST_PLAN 21): rows in their own file, tests/goldens/mag_ops_l0a1.txt, so the S0 per-slice layout can move them whole.
# (1) i64 +-3 neighbours: exponents I64_MIN+2, +3 and I64_MAX-3, -2 (the base set has only the +-0/1 neighbours);
# (2) exponent gaps for add_up and sub_lower: b's exponent below a's by 1, 2, 62, 63, 64, 65, 66, 127, 128, 129, 200 (carry, sticky, and
#     beyond-the-mantissa cases), at mid exponents and next to both i64 ends;
# (3) sqrt odd/even: radicands with every exponent parity next to 0 and next to both i64 ends, with mantissas whose square root is exact
#     only for one parity (2^63 = 2^126 / 2^63 shape) or for neither.
L0A1_NAME = "mag_ops_l0a1.txt"
NB3_EXPS = [I64_MIN + 2, I64_MIN + 3, I64_MAX - 3, I64_MAX - 2]
GAPS = [1, 2, 62, 63, 64, 65, 66, 127, 128, 129, 200]

def l0a1_binary_pairs(op):
    nb = [(e, m) for e in NB3_EXPS for m in BIN_MANTS]
    pairs = [(a, b) for a in nb for b in nb]
    pairs += [(a, b) for a in nb for b in BIN_SET] + [(a, b) for a in BIN_SET for b in nb]
    if op in ("add_up", "sub_lower"):
        for ea in [0, 100, -100, I64_MAX, I64_MAX - 1, I64_MIN + 200, I64_MIN + 300]:
            for g in GAPS:
                eb = ea - g
                if eb < I64_MIN:
                    continue
                for ma in BIN_MANTS:
                    for mb in BIN_MANTS:
                        pairs.append(((ea, ma), (eb, mb)))
                        pairs.append(((eb, mb), (ea, ma)))
    old = set(binary_pairs(op))   # rows already in mag_ops.txt are not repeated here
    pairs = [p for p in unique(pairs) if p not in old]
    if op in ("div_up", "div_lower"):
        pairs = [(a, b) for (a, b) in pairs if b != (0, 0)]
    return pairs

def l0a1_unary_args():
    exps = NB3_EXPS + [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5, I64_MIN, I64_MIN + 1, I64_MAX - 1, I64_MAX]
    mants = [TWO63, TWO63 + 1, 0xC000000000000000, MASK64, 0xB504F333F9DE6484, 0xB504F333F9DE6485]
    old = set(UNARY_SET + SQUARES)
    return [a for a in unique([(e, m) for e in exps for m in mants]) if a not in old]

GOLDEN_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tests", "goldens")
GOLDEN_NAME = "mag_ops.txt"


def main(path=None, name=GOLDEN_NAME):
    if path is None:
        path = os.path.join(GOLDEN_DIR, name)
    l0a1 = name == L0A1_NAME
    global skipped
    skipped = 0
    os.makedirs(os.path.dirname(path), exist_ok=True)
    
    with open(path, "w") as f:
        f.write("# zenith-bound Mag operation goldens%s, generated" % (" (L0a-1: i64 +-3 neighbours, exponent gaps, sqrt parity)" if l0a1 else "") + " by scripts/gen_goldens.py (mpmath.libmp directed rounding to 64 bits)\n")
        f.write("# binary: OP EA MA EB MB R ; unary: OP E M R ; R is 'ER MR' or the word inf\n")
        f.write("# value = mant * 2^(exp - 64); zero is '0 0'\n")
        f.write("# generated in the oracle venv (mpmath 1.3.0, sympy 1.14.0, gmpy2 2.3.2, python-flint 0.9.0); see scripts/README.md\n")
        
        ops = [
            ("add_up", 2),
            ("mul_up", 2),
            ("mul_lower", 2),
            ("sub_lower", 2),
            ("div_up", 2),
            ("div_lower", 2),
            ("sqrt_up", 1),
            ("sqrt_lower", 1)
        ]
        
        line_count = 0
        checked = 0
        
        for op_name, num_args in ops:
            if num_args == 2:
                pairs = l0a1_binary_pairs(op_name) if l0a1 else binary_pairs(op_name)

                for a, b in pairs:
                    if op_name == "add_up":
                        result = add_up(a, b)
                    elif op_name == "mul_up":
                        result = mul_up(a, b)
                    elif op_name == "mul_lower":
                        result = mul_lower(a, b)
                    elif op_name == "sub_lower":
                        result = sub_lower(a, b)
                    elif op_name == "div_up":
                        result = div_up(a, b)
                    elif op_name == "div_lower":
                        result = div_lower(a, b)
                    else:
                        raise ValueError(f"Unknown binary operation: {op_name}")
                    
                    if isinstance(result, str):
                        f.write(f"{op_name} {a[0]} {a[1]} {b[0]} {b[1]} {result}\n")
                    else:
                        f.write(f"{op_name} {a[0]} {a[1]} {b[0]} {b[1]} {result[0]} {result[1]}\n")
                    
                    if self_check(op_name, [a, b], result):
                        checked += 1
                    line_count += 1
            else:  # unary operations
                unique_args = l0a1_unary_args() if l0a1 else unique(UNARY_SET + SQUARES)
                
                for a in unique_args:
                    if op_name == "sqrt_up":
                        result = sqrt_up(a)
                    elif op_name == "sqrt_lower":
                        result = sqrt_lower(a)
                    else:
                        raise ValueError(f"Unknown unary operation: {op_name}")
                    
                    if isinstance(result, str):
                        f.write(f"{op_name} {a[0]} {a[1]} {result}\n")
                    else:
                        f.write(f"{op_name} {a[0]} {a[1]} {result[0]} {result[1]}\n")
                    
                    if self_check(op_name, [a], result):
                        checked += 1
                    line_count += 1
    
    print("rows %d, checked with Fractions %d, skipped (exponents too far apart to check) [window %d] %d" % (line_count, checked, WINDOW, skipped))

if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "--check":
        # Regenerate into DIR (never into tests/goldens) and compare byte for byte with the committed file (scripts/goldlib.py, `--check DIR`).
        os.makedirs(sys.argv[2], exist_ok=True)
        main(os.path.join(sys.argv[2], GOLDEN_NAME))
        main(os.path.join(sys.argv[2], L0A1_NAME), L0A1_NAME)
        raise SystemExit(1 if goldlib.check_against(sys.argv[2], GOLDEN_DIR, [GOLDEN_NAME, L0A1_NAME]) else 0)
    elif len(sys.argv) == 1:
        main()
        main(None, L0A1_NAME)
    else:
        raise SystemExit("usage: ~/team-env/oracle-2026-10-30/bin/python scripts/gen_goldens.py [--check DIR]")
