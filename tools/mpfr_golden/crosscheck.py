#!/usr/bin/env python3
"""Cross-check of the MPFR golden streams against MPFR itself (gmpy2) and published constants (r6 item 6).

gen.py computes every answer with python-flint (arb Ziv loop) or exact integers. This script recomputes the
in-range records with gmpy2 (GMP/MPFR, pinned team env oracle-2026-10-30), same precision and rounding mode, MPFR's
widest exponent range, and compares value and ternary sign. Records whose answer or operands lie outside MPFR's
range (huge / tiny, or beyond the range gmpy2 really uses, see _probe_exponent_range) cannot be checked that way and are counted as such. Every distinct record
is checked once (`@k` repeats are the same record).

Also checks the 50-digit published values of pi, ln 2 and Euler's gamma (OEIS A000796, A002162, A001620) against
the const_* records.

Usage: crosscheck.py [--jobs N] [--sample K] [--binary BIN_DIR [--case-map TSV]] TEXT_DIR      exit 0 when no record differs.
  TEXT_DIR: gen mode's readable streams (tests/mpfr/golden-full). --binary: first checks the committed zfgold2 streams in
  BIN_DIR (tests/mpfr/golden) against that text record for record (SHA256SUMS, query digests, stored answers, record
  kinds, repeats); --case-map writes the per-stream record counts.
  --sample K: every K-th distinct record of each file (default 1 = all).
"""
import hashlib
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'zenith-float-bound', 'scripts'))
import goldlib  # noqa: E402  (popguard memory scope + pinned env, fail closed)
import gmpy2  # noqa: E402
from gmpy2 import mpfr, mpq  # noqa: E402

RND = {'N': gmpy2.RoundToNearest, 'Z': gmpy2.RoundToZero, 'U': gmpy2.RoundUp, 'D': gmpy2.RoundDown,
       'A': gmpy2.RoundAwayZero}


def _probe_exponent_range():
    """The exponent range gmpy2 really computes in. gmpy2 2.3.2 accepts a context emax/emin up to 2^62 - 1 but leaves
    MPFR at its default range (2^30 - 1), so 2^(2^30) silently became inf and 2^-(2^31) zero; found by bisection."""
    def ok(k):
        v = gmpy2.mul_2exp(mpfr(1, 2), k)
        return gmpy2.is_finite(v) and v != 0
    lo, hi = 0, gmpy2.get_emax_max()
    while lo < hi:
        mid = (lo + hi + 1) // 2
        lo, hi = (mid, hi) if ok(mid) else (lo, mid - 1)
    top = lo
    lo, hi = gmpy2.get_emin_min(), 0
    while lo < hi:
        mid = (lo + hi) // 2  # floor; mid < hi
        lo, hi = (lo, mid) if ok(mid) else (mid + 1, hi)
    return lo + 1, top + 1  # MPFR exponents: 2^k has exponent k + 1


EMIN, EMAX = _probe_exponent_range()

PUBLISHED = {  # 50 significant decimal digits
    'const_pi': '3.1415926535897932384626433832795028841971693993751',
    'const_log2': '0.69314718055994530941723212145817656807550013436026',
    'const_euler': '0.57721566490153286060651209008240243104215933593992',
}


def parse_val(s):
    if s == 'nan':
        return mpfr('nan')
    neg = s[0] == '-'
    b = s[1:]
    if b == 'inf':
        return mpfr('-inf') if neg else mpfr('inf')
    if b in ('huge', 'tiny'):
        return None
    if b == '0':
        return mpfr('-0') if neg else mpfr('0')
    m, e = b.split('p')
    m, e = int(m, 16), int(e)
    if not (EMIN <= e + m.bit_length() <= EMAX):
        return None
    # exact: the sign goes on the integer, so no step rounds (negating outside this context would round to 53 bits)
    with gmpy2.context(precision=max(2, m.bit_length()), emax=EMAX, emin=EMIN):
        return gmpy2.mul_2exp(mpfr(-m if neg else m, max(2, m.bit_length())), e)


def args_of(rest):
    out = []
    while rest:
        if rest[0] == 's':
            ln, _, tail = rest[1:].partition(':')
            ln = int(ln)
            out.append(tail[:ln])
            rest = tail[ln + 1:]
            continue
        tok, _, rest = rest.partition(' ')
        out.append(int(tok[1:]) if tok[0] in 'iu' else parse_val(tok))
    return out


F1 = {k: getattr(gmpy2, k) for k in ('sqrt', 'cbrt', 'exp', 'exp2', 'exp10', 'expm1', 'log', 'log2', 'log10', 'log1p',
                                     'sin', 'cos', 'tan', 'asin', 'acos', 'atan', 'sinh', 'cosh', 'tanh', 'asinh',
                                     'acosh', 'atanh', 'gamma', 'digamma', 'erf', 'erfc', 'eint', 'ai')}
F1['lngamma'] = gmpy2.lngamma
F2 = {'add': lambda x, y: x + y, 'sub': lambda x, y: x - y, 'mul': lambda x, y: x * y, 'div': lambda x, y: x / y,
      'atan2': gmpy2.atan2, 'hypot': gmpy2.hypot, 'pow': lambda x, y: x ** y, 'remainder': gmpy2.remainder,
      'gamma_inc': gmpy2.gamma_inc}


def compute(op, a):
    if op in F1:
        return (F1[op](a[0]),)
    if op in F2:
        return (F2[op](a[0], a[1]),)
    if op == 'sqr':
        return (a[0] * a[0],)
    if op == 'fma':
        return (gmpy2.fma(a[0], a[1], a[2]),)
    if op in ('pow_ui', 'pow_si'):
        # through mpfr_pow on an exact 64-bit y: gmpy2's mpfr ** int gives (-0)^1 = +0, MPFR's pow_ui/pow_si give -0
        return (a[0] ** mpfr(a[1], 64),)
    if op == 'rootn_ui':
        return (gmpy2.rootn(a[0], a[1]),)
    if op == 'add_ui':
        return (a[0] + a[1],)
    if op == 'div_ui':
        return (a[0] / a[1],)
    if op == 'ui_div':
        return (a[0] / a[1],)
    if op == 'jn':
        return (gmpy2.jn(a[0], a[1]),)
    if op == 'yn':
        return (gmpy2.yn(a[0], a[1]),)
    if op == 'const_pi':
        return (gmpy2.const_pi(),)
    if op == 'const_log2':
        return (gmpy2.const_log2(),)
    if op == 'const_euler':
        return (gmpy2.const_euler(),)
    if op == 'sin_cos':
        return gmpy2.sin_cos(a[0])
    if op == 'sinh_cosh':
        return gmpy2.sinh_cosh(a[0])
    if op == 'parse10':
        return (mpfr(a[0]),)
    raise KeyError(op)


def same(x, y):
    if gmpy2.is_nan(x) or gmpy2.is_nan(y):
        return gmpy2.is_nan(x) and gmpy2.is_nan(y)
    return x == y and gmpy2.is_signed(x) == gmpy2.is_signed(y)


def check_record(rec):
    """'ok' | 'skip:<why>' | 'DIFF:<text>'."""
    q, a = rec.split(' => ')
    parts = q.split(' ', 3)
    op, prec, rnd = parts[0], int(parts[1]), parts[2]
    if op.startswith('q_'):
        x, y = (mpq(*map(int, t[1:].split('/'))) for t in parts[3].split(' '))
        r = {'q_add': x + y, 'q_sub': x - y, 'q_mul': x * y, 'q_div': x / y}[op]
        return 'ok' if 'q%d/%d' % (r.numerator, r.denominator) == a else 'DIFF:%s' % rec[:300]
    args = args_of(parts[3] if len(parts) > 3 else '')
    if any(v is None for v in args):
        return 'skip:operand out of MPFR range'
    toks = a.split(' ')
    want = [(parse_val(toks[k]), int(toks[k + 1])) for k in range(0, len(toks), 2)]
    if any(v is None for v, _ in want):
        return 'skip:answer out of MPFR range'
    try:
        with gmpy2.context(precision=prec, round=RND[rnd], emax=EMAX, emin=EMIN, subnormalize=False,
                           trap_underflow=False, trap_overflow=False, trap_invalid=False, trap_divzero=False) as c:
            got = compute(op, args)
            if c.overflow or c.underflow:
                return 'skip:MPFR range reached'
        with gmpy2.context(precision=prec, round=gmpy2.RoundUp, emax=EMAX, emin=EMIN, subnormalize=False,
                           trap_underflow=False, trap_overflow=False, trap_invalid=False, trap_divzero=False):
            up = compute(op, args)
        with gmpy2.context(precision=prec, round=gmpy2.RoundDown, emax=EMAX, emin=EMIN, subnormalize=False,
                           trap_underflow=False, trap_overflow=False, trap_invalid=False, trap_divzero=False):
            down = compute(op, args)
    except KeyError:
        return 'skip:op %s not in gmpy2 map' % op
    for (w, t), g, u, d in zip(want, got, up, down):
        if not same(w, g):
            return 'DIFF:value %s' % rec[:300]
        if gmpy2.is_nan(g):
            continue
        # exactness by value: an exact zero sum is +0 rounded up and -0 rounded down, both exact
        tern = 0 if u == d else (1 if same(g, u) else -1)
        if (t > 0) - (t < 0) != tern:
            return 'DIFF:ternary %d vs MPFR %d: %s' % (t, tern, rec[:300])
    return 'ok'


def records(path, sample):
    seen = 0
    with open(path, encoding='ascii') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('@'):
                continue
            seen += 1
            if seen % sample == 0:
                yield line


def published(dirn):
    bad = []
    n = 0
    for name in sorted(os.listdir(dirn)):
        if not name.endswith('.gold'):
            continue
        for rec in records(os.path.join(dirn, name), 1):
            op = rec.split(' ', 1)[0]
            if op not in PUBLISHED:
                continue
            q, a = rec.split(' => ')
            prec = int(q.split(' ')[1])
            if prec < 168 or not a.startswith('+'):
                continue  # 50 digits pin the first ~166 bits
            v = parse_val(a.split(' ')[0])
            with gmpy2.context(precision=400):
                d = abs(mpfr(v) - mpfr(PUBLISHED[op])) / mpfr(PUBLISHED[op])
            n += 1
            if d > mpfr(10) ** -49:
                bad.append('%s: %s' % (name, rec[:120]))
    return n, bad


EXACT_OPS = ('add', 'sub', 'mul', 'div', 'remainder')  # gold.rs stores only their digest; exact.rs computes them


def _varint(b, at):
    v, sh = 0, 0
    while True:
        x = b[at]
        at += 1
        v |= (x & 0x7f) << sh
        if not x & 0x80:
            return v, at
        sh += 7


def _unzigzag(u):
    return (u >> 1) ^ -(u & 1)


def _answer(b, at):
    n, at = _varint(b, at)
    toks = []
    for _ in range(n):
        tag = b[at]
        at += 1
        if tag == 0:
            neg = b[at]
            at += 1
            e, at = _varint(b, at)
            ln, at = _varint(b, at)
            m = int.from_bytes(b[at:at + ln], 'little')
            at += ln
            toks.append('%s%xp%d' % ('-' if neg else '+', m, _unzigzag(e)))
        elif tag == 1:
            u, at = _varint(b, at)
            toks.append('%d' % _unzigzag(u))
        elif tag == 2:
            ln, at = _varint(b, at)
            toks.append(b[at:at + ln].decode('ascii'))
            at += ln
        else:
            raise ValueError('answer tag %d' % tag)
    return ' '.join(toks), at


def binary_streams(bindir):
    """{name: bytes} of every zfgold2 stream in bindir (parts joined), each file checked against SHA256SUMS."""
    sums = {}
    with open(os.path.join(bindir, 'SHA256SUMS'), encoding='ascii') as f:
        for line in f:
            h, n = line.rstrip('\n').split('  ', 1)
            sums[n] = h
    files = sorted(n for n in os.listdir(bindir) if re.fullmatch(r'.+\.gold(\.[0-9]+)?', n))
    if sorted(sums) != files:
        raise SystemExit('crosscheck: SHA256SUMS lists %s, the dir has %s' % (sorted(set(sums) ^ set(files)), 'that set'))
    parts = {}
    for n in files:
        with open(os.path.join(bindir, n), 'rb') as f:
            data = f.read()
        if hashlib.sha256(data).hexdigest() != sums[n]:
            raise SystemExit('crosscheck: %s SHA256 differs from SHA256SUMS' % n)
        m = re.fullmatch(r'(.+\.gold)(?:\.([0-9]+))?', n)
        parts.setdefault(m.group(1), []).append((int(m.group(2) or 0), data))
    out = {}
    for name, ps in parts.items():
        ps.sort()
        ks = [k for k, _ in ps]
        if ks != [0] and ks != list(range(1, len(ks) + 1)):
            raise SystemExit('crosscheck: %s parts %s are not 1..n' % (name, ks))
        out[name] = b''.join(d for _, d in ps)
    return out


def compare_binary(textdir, bindir, case_map=None):
    """Every binary stream against its text form, record for record: same record count and kinds, the query digests,
    the stored answers, the repeat targets. Returns (problems, total records)."""
    bins = binary_streams(bindir)
    texts = sorted(n for n in os.listdir(textdir) if n.endswith('.gold'))
    bad = ['stream set differs: %s' % sorted(set(bins) ^ set(texts))] if sorted(bins) != texts else []
    rows, total = [], 0
    for name in texts:
        b = bins.get(name)
        if b is None:
            continue
        if not b.startswith(b'zfgold2\n'):
            bad.append('%s: not zfgold2' % name)
            continue
        at, n, stored, exact, rep = 8, 0, 0, 0, 0
        with open(os.path.join(textdir, name), encoding='ascii') as f:
            for line in f:
                line = line.rstrip('\n')
                n += 1
                if at >= len(b):
                    bad.append('%s: binary ends before text record %d' % (name, n))
                    break
                tag = b[at]
                at += 1
                if line.startswith('@'):
                    k, at = _varint(b, at) if tag == 3 else (None, at)
                    rep += 1
                    if tag != 3 or '@%d' % k != line:
                        bad.append('%s record %d: text %s, binary tag %d' % (name, n, line, tag))
                        break
                    continue
                q, a = line.split(' => ')
                dg = hashlib.sha256(q.encode('ascii')).digest()[:16]
                if b[at:at + 16] != dg:
                    bad.append('%s record %d: digest differs' % (name, n))
                    break
                at += 16
                qt = q.split(' ')  # exact.rs: `<op> <prec> <rnd> <x> <y>`, op in EXACT_OPS, two float operands
                is_exact = len(qt) == 5 and qt[0] in EXACT_OPS and not any(t[0] in 'ius' for t in qt[3:])
                if tag == 1 and not is_exact:
                    got, at = _answer(b, at)
                    stored += 1
                    if got != a:
                        bad.append('%s record %d: stored answer %s, text %s' % (name, n, got[:80], a[:80]))
                        break
                elif tag == 2 and is_exact:
                    exact += 1
                else:
                    bad.append('%s record %d: tag %d for %s' % (name, n, tag, q[:80]))
                    break
        if at != len(b):
            bad.append('%s: %d binary bytes after the last text record' % (name, len(b) - at))
        total += n
        rows.append((name, n, stored, exact, rep, len(b)))
    if case_map:
        with open(case_map, 'w') as f:
            f.write('stream\trecords\tstored\texact\trepeats\tbytes\n')
            for r in rows:
                f.write('%s\t%d\t%d\t%d\t%d\t%d\n' % r)
    return bad, total


def main(argv):
    jobs, sample = os.cpu_count() or 1, 1
    if '--jobs' in argv:
        jobs = int(argv[argv.index('--jobs') + 1])
    if '--sample' in argv:
        sample = int(argv[argv.index('--sample') + 1])
    dirn = argv[-1]
    if '--binary' in argv:
        # the committed binary streams against gen mode's text form (golden-full), before the MPFR check of that text
        bad, total = compare_binary(dirn, argv[argv.index('--binary') + 1],
                                    argv[argv.index('--case-map') + 1] if '--case-map' in argv else None)
        print('crosscheck: binary vs text: %d records, %d problems' % (total, len(bad)), flush=True)
        for x in bad[:20]:
            print('  ' + x)
        if bad:
            return 1
    names = sorted(n for n in os.listdir(dirn) if n.endswith('.gold'))
    tally, diffs = {}, []
    with ProcessPoolExecutor(max_workers=jobs) as ex:
        for name in names:
            recs = list(records(os.path.join(dirn, name), sample))
            for r in ex.map(check_record, recs, chunksize=64):
                # DIFFs tallied per kind and op, so a script bug shows up as one op rather than a total
                k = ' '.join(r.split(' ', 2)[:2]) if r.startswith('DIFF:value') else (
                    'DIFF:ternary ' + r.split(': ', 1)[1].split(' ', 1)[0] if r.startswith('DIFF:ternary') else r)
                tally[k] = tally.get(k, 0) + 1
                if r.startswith('DIFF'):
                    diffs.append('%s: %s' % (name, r))
            print('crosscheck: %s %d records' % (name, len(recs)), flush=True)
    n_pub, bad_pub = published(dirn)
    for k in sorted(tally):
        print('crosscheck: %-40s %d' % (k, tally[k]))
    print('crosscheck: published 50-digit constants checked %d, wrong %d' % (n_pub, len(bad_pub)))
    for d in diffs[:50] + bad_pub[:20]:
        print('  ' + d)
    print('crosscheck: team-env=%s gmpy2 %s MPFR %s' % (goldlib.TEAM_ENV, gmpy2.version(), gmpy2.mpfr_version()))
    return 1 if diffs or bad_pub else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
