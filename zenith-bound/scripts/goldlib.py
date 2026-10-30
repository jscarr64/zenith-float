#!/home/jsc64/.local/oracle-venv/bin/python
"""
goldlib: shared helpers for the golden generators (TEST_PLAN section 9). Run self-tests with:  ~/.local/oracle-venv/bin/python scripts/goldlib.py --selftest
Rules: the oracle venv only (require_oracle: mpmath 1.3.0, sympy 1.14.0, gmpy2 2.3.2, python-flint 0.9.0), no Python `random` (its algorithms changed across 3.x), integer-only arithmetic or Fraction,
no dict/set iteration order in output, LF line endings, no timestamps/hostnames/paths in generated files, no network.

`--check DIR` DESIGN (every generator script implements it with check_against below):
  1. `GENERATOR --check DIR` writes the file(s) it would produce into DIR (never into tests/goldens),
  2. then compares each, byte for byte, with the committed file of the same name in tests/goldens,
  3. prints `check: OK <name>` or the first differing line (`check: DIFF <name> line N`, committed vs regenerated) and exits 1 on any
     difference, on a missing committed file, or on a missing regenerated file.
tools/goldens_regen.sh runs this for every generator into a fresh mktemp dir under /tmp, then `sha256sum -c tests/goldens/SHA256SUMS`.
"""
import os
import sys
from fractions import Fraction

MASK64 = (1 << 64) - 1

# ---- oracle startup check (Jeff 9:22 AM CT 2026-10-23: goldens come from the isolated oracle venv) ----------------------------
# Exact pins, the same as zenith-float scripts/oracle-requirements.lock and scripts/oracle_pins.py at 7c97b3c (read, not imported: zf is
# not a dependency). The 'flint' in ~/.local is an unrelated tool, so nothing may be imported from ~/.local outside the venv.
ORACLE_VENV = os.path.expanduser("~/.local/oracle-venv")
ORACLE_PY = os.path.join(ORACLE_VENV, "bin", "python")
USER_LOCAL = os.path.expanduser("~/.local")
PINS = (("mpmath", "1.3.0"), ("sympy", "1.14.0"), ("gmpy2", "2.3.2"), ("flint", "0.9.0"))
MPMATH_VERSION = "1.3.0"


def _inside(path, root):
    return path == root or path.startswith(root + os.sep)


def oracle_problems(found, prefix, base_prefix, path, venv=ORACLE_VENV, user_local=USER_LOCAL):
    """Every reason the running interpreter is not the pinned oracle venv (empty list: ok). Pure; selftest covers each case."""
    probs = []
    for name, want in PINS:
        if found.get(name) != want:
            probs.append("%s is %s, want %s" % (name, found.get(name), want))
    venv_r = os.path.realpath(venv)
    if os.path.realpath(prefix) != venv_r:
        probs.append("sys.prefix is %s, not the oracle venv %s" % (prefix, venv))
    if prefix == base_prefix:
        probs.append("sys.prefix equals sys.base_prefix (not a venv)")
    local_r = os.path.realpath(user_local)
    for entry in path:
        r = os.path.realpath(entry) if entry else ""
        if r and _inside(r, local_r) and not _inside(r, venv_r):
            probs.append("~/.local path on sys.path: %s" % entry)
    return probs


def _found_versions():
    import importlib
    found = {}
    for name, _want in PINS:
        try:
            mod = importlib.import_module(name)
        except ImportError:
            found[name] = "not importable"
            continue
        found[name] = getattr(mod, "__version__", "no __version__")
        if name == "flint" and not hasattr(mod, "arb"):
            found[name] = "not-python-flint"
    return found


def require_oracle():
    """Exit nonzero unless this is the oracle venv with exactly the four pinned versions and no ~/.local on sys.path. Returns mpmath."""
    probs = oracle_problems(_found_versions(), sys.prefix, sys.base_prefix, sys.path)
    if probs:
        raise SystemExit("goldlib: oracle check failed (run with %s or GOLDENS_PYTHON pointing to it): %s" % (ORACLE_PY, "; ".join(probs)))
    import mpmath
    return mpmath


def require_mpmath():
    """Kept for the generators: the full oracle check, then mpmath."""
    return require_oracle()


def require_popslot_scope(name=None):
    """Pop direct-run guard (~/pop-logs/bin/popguard.sh; 11:49 2026-10-24 a bare 58 GB generator run global-OOMed Pop):
    exit 3 outside a capped popslot/poprun memory scope; exit 2 if POPGUARD is unusable on a host that has popguard.sh.
    Hosts without popguard.sh (CI) are not guarded and get a stderr WARNING."""
    import subprocess
    host_guard = os.path.expanduser("~/pop-logs/bin/popguard.sh")
    guard = os.environ.get("POPGUARD") or host_guard
    if not os.access(guard, os.X_OK):
        # POPGUARD pointing at nothing must not switch the guard off (Spence, 2026-10-24)
        if os.path.exists(host_guard):
            sys.stderr.write("goldlib: REFUSING: POPGUARD=%s is missing or not executable, but this host has %s "
                             "(the memory guard cannot be switched off)\n" % (guard, host_guard))
            raise SystemExit(2)
        sys.stderr.write("goldlib: WARNING: no popguard.sh on this host (%s): running WITHOUT the Pop memory guard\n" % guard)
        return
    name = name or os.path.basename(sys.argv[0] or "goldlib")
    if subprocess.run([guard, name]).returncode != 0:
        raise SystemExit(3)


require_popslot_scope()   # at import, before the oracle check: every mpmath generator imports goldlib first
require_oracle()   # at import: every generator imports goldlib before mpmath


# ---- splitmix64: mirrors tests/common/rng.rs and tests/common/seeds.rs -------------------------------------------------------
INCREMENT = 0x9E3779B97F4A7C15
MUL1 = 0xBF58476D1CE4E5B9
MUL2 = 0x94D049BB133111EB
BASE = 0x5A454E4954424F55      # tests/common/seeds.rs BASE
GENERATOR_VERSION = 4          # tests/common/seeds.rs GENERATOR_VERSION


class SplitMix64:
    """One u64 state; next() returns the next output (same draw order as the Rust generator)."""

    def __init__(self, seed):
        self.state = seed & MASK64

    def next(self):
        self.state = (self.state + INCREMENT) & MASK64
        z = self.state
        z = ((z ^ (z >> 30)) * MUL1) & MASK64
        z = ((z ^ (z >> 27)) * MUL2) & MASK64
        return z ^ (z >> 31)


def splitmix64(seed):
    """First output of a generator seeded with `seed` (the Rust `splitmix64`)."""
    return SplitMix64(seed).next()


def stream_seed(suite_id, case_class):
    """splitmix64(BASE + suite_id * 2^32 + case_class), wrapping at 2^64 (the Rust `stream_seed`)."""
    if not (0 <= suite_id <= 0xFFFFFFFF and 0 <= case_class <= 0xFFFFFFFF):
        raise ValueError("suite_id and case_class are u32")
    return splitmix64((BASE + suite_id * (1 << 32) + case_class) & MASK64)


# ---- row formatting and sorted output --------------------------------------------------------------------------------------
def fmt_row(*fields):
    """One golden row: fields joined by single spaces; ints in decimal (as the existing file does), strings as given."""
    return " ".join(str(f) for f in fields)


def fmt_dyadic(sign, exp, mant):
    """A dyadic value as `sign exp mant` (sign 0 or 1, mant an int printed in decimal)."""
    if sign not in (0, 1):
        raise ValueError("sign must be 0 or 1")
    return "%d %d %d" % (sign, exp, mant)


def sorted_rows(rows, key=None):
    """Rows sorted by a documented key (default: the row text itself), duplicates removed; deterministic."""
    return sorted(set(rows), key=key)


def write_lines(path, lines):
    """Write lines with LF endings and a final newline, as bytes (no platform translation)."""
    data = "".join(l + "\n" for l in lines).encode("ascii")
    with open(path, "wb") as f:
        f.write(data)


# ---- --check DIR ------------------------------------------------------------------------------------------------------------
def first_difference(committed, regenerated):
    """None if the two byte strings are identical, else a message naming the first differing line (1-based)."""
    if committed == regenerated:
        return None
    a = committed.split(b"\n")
    b = regenerated.split(b"\n")
    for i in range(max(len(a), len(b))):
        la = a[i] if i < len(a) else None
        lb = b[i] if i < len(b) else None
        if la != lb:
            return "line %d: committed %r, regenerated %r" % (i + 1, la, lb)
    return "contents differ"


def check_against(dirpath, committed_dir, names):
    """Compare dirpath/<name> with committed_dir/<name> for each name; print one result line each; return the number of failures."""
    bad = 0
    for name in names:
        cpath = os.path.join(committed_dir, name)
        rpath = os.path.join(dirpath, name)
        if not os.path.isfile(cpath):
            print("check: MISSING committed %s" % name)
            bad += 1
            continue
        if not os.path.isfile(rpath):
            print("check: MISSING regenerated %s" % name)
            bad += 1
            continue
        with open(cpath, "rb") as f:
            c = f.read()
        with open(rpath, "rb") as f:
            r = f.read()
        d = first_difference(c, r)
        if d is None:
            print("check: OK %s" % name)
        else:
            print("check: DIFF %s %s" % (name, d))
            bad += 1
    return bad


# ---- directed rounding of an exact Fraction, six modes ----------------------------------------------------------------------
MODES = ("ToEven", "ToOdd", "Up", "Down", "ToZero", "FromZero")
# ToEven / ToOdd: round to nearest, an exact tie goes to the neighbour with an even / odd mantissa.
# Up = toward +infinity, Down = toward -infinity, ToZero = toward zero, FromZero = away from zero.
# For prec == 1 the mantissa is always 1, so the tie parity is taken from the unit count of the lower neighbour (documented, not cross-checked).


def _floor_log2(a):
    """floor(log2(a)) for a positive Fraction, exactly."""
    n, d = a.numerator, a.denominator
    e = n.bit_length() - d.bit_length()
    # 2^e <= n/d < 2^(e+1) up to one step; fix exactly
    if e >= 0:
        ok = n >= (d << e)
    else:
        ok = (n << -e) >= d
    return e if ok else e - 1


def round_fraction(x, prec, mode):
    """Round the exact value x (a Fraction) to `prec` significant bits under `mode`.
    Returns (value, inexact): value is an exact Fraction (a dyadic with at most prec significant bits), inexact is True iff value != x.
    There is no exponent range here (the model of an unbounded exponent); zero is exact."""
    if mode not in MODES:
        raise ValueError("unknown rounding mode %r" % (mode,))
    if prec < 1:
        raise ValueError("prec must be >= 1")
    x = Fraction(x)
    if x == 0:
        return Fraction(0), False
    neg = x < 0
    a = -x if neg else x
    e = _floor_log2(a) - prec + 1               # unit = 2^e, so the mantissa q = a / unit satisfies 2^(prec-1) <= q < 2^prec
    q = a / (Fraction(2) ** e) if e >= 0 else a * (1 << -e)
    lo = q.numerator // q.denominator            # lower neighbour in units
    rem = q - lo
    if rem == 0:
        return x, False
    hi = lo + 1
    twice = rem * 2
    if mode == "ToZero":
        pick = lo
    elif mode == "FromZero":
        pick = hi
    elif mode == "Up":
        pick = lo if neg else hi
    elif mode == "Down":
        pick = hi if neg else lo
    else:  # ToEven, ToOdd
        if twice < 1:
            pick = lo
        elif twice > 1:
            pick = hi
        elif mode == "ToEven":
            pick = lo if lo % 2 == 0 else hi
        else:
            pick = lo if lo % 2 == 1 else hi
    v = Fraction(pick) * (Fraction(2) ** e if e >= 0 else Fraction(1, 1 << -e))
    return (-v if neg else v), True


# mpmath's rounding names for the modes it has: 'n' nearest-even, 'c' ceiling (+inf), 'f' floor (-inf), 'd' toward zero, 'u' away from zero.
MPMATH_MODE = {"ToEven": "n", "Up": "c", "Down": "f", "ToZero": "d", "FromZero": "u"}   # no ToOdd in mpmath


def mpmath_round(x, prec, mode):
    """The same rounding through mpmath.libmp.from_rational (only for modes mpmath has). Returns a Fraction."""
    require_mpmath()
    from mpmath.libmp import from_rational
    if mode not in MPMATH_MODE:
        raise ValueError("mpmath has no mode %s" % mode)
    x = Fraction(x)
    t = from_rational(x.numerator, x.denominator, prec, MPMATH_MODE[mode])
    sign, man, exp, _bc = t          # mpmath's tuple: the sign is separate from the unsigned mantissa
    man = -man if sign else man
    return Fraction(man) * (Fraction(2) ** exp if exp >= 0 else Fraction(1, 1 << -exp))


# ---- self-tests -------------------------------------------------------------------------------------------------------------
def _brute_round(x, prec, mode):
    """Independent reference for round_fraction: enumerate the two grid neighbours of x explicitly and choose with exact comparisons."""
    x = Fraction(x)
    if x == 0:
        return Fraction(0)
    s = 1 if x > 0 else -1
    a = abs(x)
    e = 0
    while a >= 2 ** (e + 1):
        e += 1
    while a < 2 ** e:
        e -= 1
    unit = Fraction(2) ** (e - prec + 1)
    lo = (a // unit) * unit
    if lo == a:
        return x
    hi = lo + unit
    neighbours = (lo, hi)
    if mode == "ToZero":
        r = lo
    elif mode == "FromZero":
        r = hi
    elif mode == "Up":
        r = hi if s > 0 else lo
    elif mode == "Down":
        r = lo if s > 0 else hi
    else:
        dl, dh = a - lo, hi - a
        if dl != dh:
            r = lo if dl < dh else hi
        else:
            odd = [n for n in neighbours if (n / unit) % 2 == 1]
            even = [n for n in neighbours if (n / unit) % 2 == 0]
            want = even if mode == "ToEven" else odd
            r = want[0] if len(want) == 1 else lo      # prec == 1: both parities equal, lower unit chosen
    return s * r


def _check(cond, msg):
    if not cond:
        raise SystemExit("goldlib selftest FAIL: " + msg)


def selftest():
    # oracle startup check (pure part): each mismatch must be reported, the good case must be clean
    venv = "/home/u/.local/oracle-venv"
    good = {"mpmath": "1.3.0", "sympy": "1.14.0", "gmpy2": "2.3.2", "flint": "0.9.0"}
    vpath = ["", "/usr/local/lib/python3.12", venv + "/lib/python3.12/site-packages"]
    _check(oracle_problems(good, venv, "/usr/local", vpath, venv, "/home/u/.local") == [], "good oracle reported a problem")
    for name, bad in (("mpmath", "1.2.1"), ("sympy", "1.12"), ("gmpy2", "not importable"), ("flint", "not-python-flint"), ("mpmath", "9.9.9")):
        f = dict(good)
        f[name] = bad
        _check(len(oracle_problems(f, venv, "/usr/local", vpath, venv, "/home/u/.local")) == 1, "%s %s not reported" % (name, bad))
    _check(oracle_problems(good, "/usr", "/usr", vpath, venv, "/home/u/.local") != [], "system prefix not reported")
    _check(oracle_problems(good, venv, venv, vpath, venv, "/home/u/.local") != [], "prefix == base_prefix (not a venv) not reported")
    _check(oracle_problems(good, venv, "/usr/local", vpath + ["/home/u/.local/lib/python3.12/site-packages"], venv, "/home/u/.local") != [], "~/.local on sys.path not reported")
    # splitmix64: the same reference values as tests/harness/*_selftest.rs (computed independently of Rust)
    g = SplitMix64(0)
    _check([g.next() for _ in range(3)] == [0xE220A8397B1DCDAF, 0x6E789E6AA1B965F4, 0x06C45D188009454F], "seed 0 reference outputs")
    g = SplitMix64(1)
    _check([g.next() for _ in range(3)] == [0x910A2DEC89025CC1, 0xBEEB8DA1658EEC67, 0xF893A2EEFB32555E], "seed 1 reference outputs")
    g = SplitMix64(MASK64)
    _check([g.next() for _ in range(3)] == [0xE4D971771B652C20, 0xE99FF867DBF682C9, 0x382FF84CB27281E9], "seed u64::MAX reference outputs")
    _check(stream_seed(0, 0) == 0xA6A3ED836C686CD0, "stream_seed(0,0)")
    _check(stream_seed(1, 1) == 0x8BC77992C3640E5B, "stream_seed(1,1)")
    _check(stream_seed(0xFFFFFFFF, 0xFFFFFFFF) == 0x86CFDEB4D1586C75, "stream_seed(max,max) wraps")
    _check(splitmix64(0xDEADBEEF) == 0x4ADFB90F68C9EB9B, "splitmix64(0xDEADBEEF)")

    # row formatting and sorted output
    _check(fmt_row("add_up", 1, -2, "inf") == "add_up 1 -2 inf", "fmt_row")
    _check(fmt_dyadic(1, -5, 7) == "1 -5 7", "fmt_dyadic")
    _check(sorted_rows(["b", "a", "b"]) == ["a", "b"], "sorted_rows")
    _check(first_difference(b"a\nb\n", b"a\nb\n") is None, "first_difference equal")
    _check("line 2" in first_difference(b"a\nb\n", b"a\nc\n"), "first_difference reports the line")

    # hand-derived rounding cases (prec 3: mantissas 4..7 in units of 2^e)
    F = Fraction
    cases = [
        # (x, prec, mode, expected value)  -- 81/64 = 1.265625; grid at prec 3 is multiples of 1/4 in [1,2)
        (F(81, 64), 3, "Up", F(3, 2)), (F(81, 64), 3, "Down", F(5, 4)), (F(81, 64), 3, "ToZero", F(5, 4)),
        (F(81, 64), 3, "FromZero", F(3, 2)), (F(81, 64), 3, "ToEven", F(5, 4)), (F(81, 64), 3, "ToOdd", F(5, 4)),
        (F(-81, 64), 3, "Up", F(-5, 4)), (F(-81, 64), 3, "Down", F(-3, 2)), (F(-81, 64), 3, "ToZero", F(-5, 4)),
        (F(-81, 64), 3, "FromZero", F(-3, 2)),
        # exact tie 9/8 = 1.125 between 1 (mant 4, even) and 5/4 (mant 5, odd) at prec 3
        (F(9, 8), 3, "ToEven", F(1)), (F(9, 8), 3, "ToOdd", F(5, 4)),
        (F(-9, 8), 3, "ToEven", F(-1)), (F(-9, 8), 3, "ToOdd", F(-5, 4)),
        # exact tie 11/8 = 1.375 between 5/4 (odd) and 3/2 (even)
        (F(11, 8), 3, "ToEven", F(3, 2)), (F(11, 8), 3, "ToOdd", F(5, 4)),
        # carry into the next binade: 31/16 at prec 3 rounds up to 2
        (F(31, 16), 3, "Up", F(2)), (F(31, 16), 3, "ToEven", F(2)), (F(31, 16), 3, "Down", F(7, 4)),
        # exact values do not move in any mode
        (F(5, 4), 3, "Up", F(5, 4)), (F(0), 5, "FromZero", F(0)),
    ]
    for x, p, m, want in cases:
        got, inexact = round_fraction(x, p, m)
        _check(got == want, "round_fraction(%s, %d, %s) = %s, want %s" % (x, p, m, got, want))
        _check(inexact == (got != x), "inexact flag of round_fraction(%s, %d, %s)" % (x, p, m))

    # cross-checks on a deterministic sample: brute-force reference (all 6 modes) and mpmath (the 5 modes it has)
    require_mpmath()
    rng = SplitMix64(stream_seed(0, 2))
    n = 0
    for _ in range(3000):
        den = (rng.next() % 4096) + 1
        num = (rng.next() % (1 << 40)) - (1 << 39)
        shift = (rng.next() % 41) - 20
        x = Fraction(num, den) * (Fraction(2) ** shift if shift >= 0 else Fraction(1, 1 << -shift))
        prec = (rng.next() % 40) + 1
        for mode in MODES:
            got, inexact = round_fraction(x, prec, mode)
            _check(got == _brute_round(x, prec, mode), "brute-force mismatch: x=%s prec=%d %s got %s" % (x, prec, mode, got))
            _check(inexact == (got != x), "inexact flag mismatch: x=%s prec=%d %s" % (x, prec, mode))
            if mode in MPMATH_MODE and (prec >= 2 or mode != "ToEven"):   # prec 1 nearest-even tie parity is not cross-checked
                _check(got == mpmath_round(x, prec, mode), "mpmath mismatch: x=%s prec=%d %s got %s want %s" % (x, prec, mode, got, mpmath_round(x, prec, mode)))
            n += 1
    # ties are rare in the random sample: sweep every tie at prec 2..12 explicitly (both signs)
    for prec in range(2, 13):
        for mant in range(1 << (prec - 1), 1 << prec):
            for sign in (1, -1):
                x = sign * (Fraction(2 * mant + 1, 2)) * Fraction(1, 1 << 3)
                for mode in MODES:
                    got, _ = round_fraction(x, prec, mode)
                    _check(got == _brute_round(x, prec, mode), "tie sweep brute mismatch: x=%s prec=%d %s" % (x, prec, mode))
                    if mode in MPMATH_MODE:
                        _check(got == mpmath_round(x, prec, mode), "tie sweep mpmath mismatch: x=%s prec=%d %s" % (x, prec, mode))
                    n += 1
    # ToOdd and ToEven differ exactly on ties
    for prec in range(2, 9):
        for mant in range(1 << (prec - 1), 1 << prec):
            x = Fraction(2 * mant + 1, 2)
            ev, _ = round_fraction(x, prec, "ToEven")
            od, _ = round_fraction(x, prec, "ToOdd")
            _check(ev != od, "ToEven and ToOdd agree on a tie x=%s prec=%d" % (x, prec))
    print("goldlib selftest: ok (%d rounding comparisons)" % n)


if __name__ == "__main__":
    if sys.argv[1:] == ["--selftest"]:
        selftest()
    else:
        raise SystemExit("usage: ~/.local/oracle-venv/bin/python scripts/goldlib.py --selftest")
