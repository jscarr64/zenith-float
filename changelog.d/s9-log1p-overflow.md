### Fixed
- `log1p(x)` for huge x no longer overflows to +Inf: it uses ln(x) + 1/x once 1/x is below the working precision.
