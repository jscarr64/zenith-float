# Provenance

Audit of this tree against [astro-float](https://github.com/stencillogic/astro-float) (MIT, Copyright (c) 2022 stencillogic).

Method: rename `ExactNumNumber` to `BigFloatNumber`, `ExactNum` to `BigFloat`, and the `zenith_float` paths to `astro_float`, then compare stripped lines with Python `difflib.SequenceMatcher` against the same relative path in astro-float 0.9.6 (commit `c0d9a26bd5d66871212f66b54475177b46b9bb66`, 2026-08-07; crates `astro-float` 0.9.6, `astro-float-num` 0.3.7, `astro-float-macro` 0.4.6). The same files were also compared to 0.9.5 (`cd4e409`, num 0.3.6), 0.9.4, 0.9.3, 0.9.1, 0.9.0, 0.7.1, and 0.7.0. 0.9.6 is the best match except `zenith-float-num/src/num.rs`, which is 0.916 against 0.9.5 / num 0.3.6 and 0.914 against 0.9.6.

Classes: **verbatim** (line ratio ≥ 0.98 after that rename), **renamed** (≥ 0.90: the same source with the `BigFloat` names changed to `ExactNum`), **adapted** (≥ 0.45, or a short module file that is the astro-float file plus later additions). License for every row is MIT. The notice is in `NOTICE`.

| File | Origin | License | How much |
| --- | --- | --- | --- |
| `zenith-float-num/src/common/int.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | verbatim (1.000) |
| `zenith-float-num/src/mantissa/mod.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | verbatim (1.000) |
| `tests/tests/expr.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | verbatim (1.000) |
| `doc/README.md` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | verbatim (1.000) |
| `zenith-float-num/src/ops/sinh.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | verbatim (0.982) |
| `zenith-float-num/src/ops/consts/ln2.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | verbatim (0.981) |
| `zenith-float-num/src/ops/consts/pi.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | renamed (0.977) |
| `zenith-float-num/src/ops/consts/ln10.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | renamed (0.973) |
| `zenith-float-num/src/ops/consts/e.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | renamed (0.970) |
| `zenith-float-num/src/mantissa/mantissa.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | renamed (0.953) |
| `zenith-float-num/src/ops/tan.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | renamed (0.953) |
| `zenith-float-num/src/common/consts.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | renamed (0.951) |
| `zenith-float-num/src/mantissa/sqrt.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | renamed (0.951) |
| `zenith-float-num/src/ops/atan.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | renamed (0.950) |
| `zenith-float-num/src/parser.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | renamed (0.940) |
| `zenith-float-num/src/mantissa/cbrt.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | renamed (0.934) |
| `zenith-float-num/src/mantissa/toom2.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | renamed (0.932) |
| `zenith-float-num/src/mantissa/toom3.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | renamed (0.931) |
| `zenith-float-num/src/ops/log.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | renamed (0.928) |
| `zenith-float-num/src/ops/cbrt.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | renamed (0.925) |
| `zenith-float-num/src/conv.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | renamed (0.925) |
| `zenith-float-num/src/ops/sqrt.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | renamed (0.924) |
| `zenith-float-num/src/mantissa/fft.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | renamed (0.922) |
| `zenith-float-num/src/num.rs` | astro-float-num 0.3.6 (`cd4e409`) slightly over 0.3.7; same path | MIT, Copyright (c) 2022 stencillogic | renamed (0.917) |
| `zenith-float-num/src/ops/tanh.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | renamed (0.915) |
| `zenith-float-num/src/ops/asin.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | renamed (0.912) |
| `zenith-float-num/src/ops/atanh.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | renamed (0.908) |
| `zenith-float-num/src/ops/acos.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | renamed (0.903) |
| `zenith-float-num/src/mantissa/conv.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | adapted (0.899) |
| `zenith-float-num/src/macro_util.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | adapted (0.897) |
| `zenith-float-num/src/ops/acosh.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | adapted (0.895) |
| `zenith-float-num/tests/mpfr/compare_const_test.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | adapted (0.893) |
| `zenith-float-num/src/ops/asinh.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | adapted (0.889) |
| `zenith-float-num/src/ops/cosh.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | adapted (0.889) |
| `zenith-float-num/tests/integration_tests.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | adapted (0.889) |
| `zenith-float-num/src/mantissa/util.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | adapted (0.887) |
| `zenith-float-num/src/ctx.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | adapted (0.875) |
| `zenith-float-num/src/ops/util.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | adapted (0.873) |
| `zenith-float-num/src/ops/pow.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | adapted (0.851) |
| `zenith-float-num/src/ops/tests.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | adapted (0.849) |
| `zenith-float-num/src/mantissa/div.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | adapted (0.838) |
| `zenith-float-num/src/common/mod.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | adapted (0.833) |
| `zenith-float-num/src/strop.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | adapted (0.833) |
| `zenith-float-num/src/common/util.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | adapted (0.832) |
| `src/lib.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | adapted (0.824) |
| `zenith-float-num/src/ops/sin.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | adapted (0.822) |
| `zenith-float-num/src/defs.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | adapted (0.787) |
| `zenith-float-num/tests/mpfr/compare_special_test.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | adapted (0.783) |
| `zenith-float-num/src/ops/series.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | adapted (0.778) |
| `zenith-float-num/src/ops/cos.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | adapted (0.775) |
| `zenith-float-num/src/ext.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | adapted (0.754) |
| `zenith-float-num/tests/mpfr/common.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | adapted (0.748) |
| `zenith-float-num/src/ops/mod.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | adapted (0.746) |
| `zenith-float-num/tests/mpfr/compare_ops_test.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | adapted (0.702) |
| `zenith-float-num/src/mantissa/mul.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | adapted (0.693) |
| `zenith-float-macro/src/util.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | adapted (0.673) |
| `tests/mod.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | adapted (0.656) |
| `zenith-float-num/tests/mpfr/mod.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | adapted (0.571) |
| `zenith-float-num/src/common/buf.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | adapted (0.525) |
| `zenith-float-macro/src/lib.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | adapted (0.511) |
| `zenith-float-num/src/lib.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | adapted (0.495) |
| `zenith-float-num/src/ops/consts/mod.rs` | astro-float 0.9.6 (`c0d9a26`), same path | MIT, Copyright (c) 2022 stencillogic | adapted (0.493) |

## Reviewed and not treated as a substantial copy

| File | What was checked | Result |
| --- | --- | --- |
| `zenith-float-num/src/for_3rd/de.rs`, `ser.rs`, `codec.rs` | astro-float serde visitors | line ratio about 0.12–0.16; rewritten for `ExactNum` arrays and `@p=` |
| `zenith-float-num/src/for_3rd/mod.rs` | astro-float `for_3rd/mod.rs` | six lines of `mod` declarations; the bodies are not the astro-float serde impl |
| `zenith-float-macro/src/cplx.rs` | astro-float-macro `lib.rs` | three shared error strings; the complex parser is new (line ratio 0.29, no copied block longer than a match arm) |
| `zenith-float-num/src/ops/sin_cos.rs`, `sinh_cosh.rs`, `log1p.rs`, `expm1.rs`, `atan2.rs`, `fma.rs`, `hypot.rs`, `recip.rs`, `nroot.rs`, `special.rs`, `airy.rs`, `jacobi.rs`, `ieee.rs`, `sign.rs` | same-named or nearest astro-float ops | no astro-float counterpart, or only shared error-doc boilerplate |
| `zenith-float-num/src/gauss_kronrod.rs`, `quadrature.rs` | QUADPACK / SLATEC Fortran names (`xgk`, `wgk`, `dqk15`) | design is cited (Laurie 1997, Piessens et al. 1983); the Fortran is not in this tree |
| `zenith-float-num/src/hash.rs` | FIPS 180-4 SHA-256 / SHA-512 | standard constants and a short implementation; no third-party source matched |
| `zenith-float-num/src/heun_fox.rs`, `catalog.rs`, `complex_*.rs`, `ieee_soft/`, `ode.rs`, `roots.rs`, `poly.rs`, `ball.rs`, `rational.rs`, `integer.rs`, `modular.rs`, `dsp.rs`, `orthopoly.rs`, `chebyshev.rs`, `random_dist.rs`, `dist.rs`, `csvfmt.rs`, `binfmt.rs`, `radix_float.rs` | astro-float tree | not present in astro-float |
| mpmath / SymPy strings under `tests/` | numeric golds | decimal oracles, not ported source |
| `rustfmt.toml` | astro-float `rustfmt.toml` | two generic formatter keys, identical; not a substantial copy |
| `zenith-float-compare/` | astro-float 0.9.6 as a Cargo dependency | benchmark backend, not a source copy. Not published |

The rest of the repository (benches, `nostd-tests/`, `golds/`, `scripts/`, user guides other than `doc/README.md`) had no substantial match to astro-float.

