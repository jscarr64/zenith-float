### Fixed
- `erf` near ±1 and `erfc` near 1 and 2 round correctly under directed modes: when the tail is below half an ulp they round `base ± tiny` once instead of returning the exact constant, and `erfc(tiny)` no longer exhausts the Ziv retry budget.
