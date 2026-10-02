# expr!

> **Please upgrade to zenith-float 1.0.6.** Versions 1.0.0 through 1.0.4
> return wrong values from `erf`, `erfc`, `normal_cdf`, `hypergeom_2f1`,
> `betainc`, the gamma family, Bessel and several other special functions for
> some inputs, with no error reported, and their quadrature, root-finding and
> ODE routines return `None` for intervals such as `[2, 3]`. 1.0.5 fixes these
> with no API changes:
> `cargo update -p zenith-float -p zenith-float-num -p zenith-float-macro`.
> Versions 1.0.0 through 1.0.4 have been yanked from crates.io.
> Details: [1.0.5 release notes](https://github.com/jscarr64/zenith-float/releases/tag/v1.0.5).

Procedural macro crate for [zenith-float](https://github.com/jscarr64/zenith-float). Depend on `zenith-float` and use `expr!`; you do not need to depend on this crate directly.
