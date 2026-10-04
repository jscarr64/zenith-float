# expr!

> **Please upgrade to zenith-float 1.0.9.** Versions 1.0.0 through 1.0.4
> return wrong values from `erf`, `erfc`, `normal_cdf`, `hypergeom_2f1`,
> `betainc`, the gamma family, Bessel and several other special functions for
> some inputs, with no error reported, and their quadrature, root-finding and
> ODE routines return `None` for intervals such as `[2, 3]`. 1.0.5 fixes these
> with no API changes:
> `cargo update -p zenith-float -p zenith-float-num -p zenith-float-macro`.
> Versions 1.0.0 through 1.0.4 have been yanked from crates.io.
> Details: [1.0.5 release notes](https://github.com/jscarr64/zenith-float/releases/tag/v1.0.5). 1.0.6 is an additive SoftFloat pack. 1.0.7 is a compile-only fix: the 1.85-required unsafe wrappers around `_addcarry_u64`/`_subborrow_u64` stay, and `unused_unsafe` is allowed so the crate builds on rustc 1.93+ (those intrinsics became safe in 1.93, not 1.87). No math change. 1.0.8 corrects the README upgrade target; no math. 1.0.9 sizes the Bessel \(I\) series guard to cancellation instead of \(1.5\lvert z\rvert+16\).

Procedural macro crate for [zenith-float](https://github.com/jscarr64/zenith-float). Depend on `zenith-float` and use `expr!`; you do not need to depend on this crate directly.
