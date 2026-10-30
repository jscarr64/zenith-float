# Security policy

## What counts as a security issue

zenith-bound's only promise is **containment**: the true result of every guaranteed operation lies inside the returned ball, or the call returns an error or `Indeterminate`. Downstream users (Accumath, automated tools, signed result traces) rely on that promise without re-checking it.

So these are treated as **security-grade** vulnerabilities, not ordinary bugs:

- **Containment failures**: any input where a function marked *Guaranteed* returns a ball that does not contain the true value. This includes a radius that is too small, a wrong midpoint with a zero radius, or a finite ball returned where the contract requires `Domain`, `Overflow` or `Indeterminate`.
- Decoder acceptance of non-canonical or malformed `zf-ball-enc` bytes, and any way a trace can pass the `zenith-bound-check` checker without being valid.
- Panics, unbounded allocation or non-termination reachable from untrusted input (decoding, decimal parsing, huge exponent gaps).
- Results that differ by platform, build profile or thread count (determinism breaks).

Tightness problems (balls that are correct but wider than they need to be) are ordinary performance bugs. Report them in public issues.

## Supported versions

zenith-bound is pre-release (0.0.0). Nothing is guaranteed yet, and no release is supported. Once 0.1 ships, the latest 0.x release gets fixes.

## Reporting

Please report privately through GitHub's **"Report a vulnerability"** (Security Advisories) on `jscarr64/zenith-bound`. Don't open a public issue. Include the inputs (preferably as `zf-ball-enc` hex), `prec`, rounding mode, crate versions, target, and the oracle value showing the miss.

We aim to acknowledge reports within 7 days. A confirmed containment failure gets a fix, a regression set, a RustSec advisory and a yanked affected release.
