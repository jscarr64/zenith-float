### Changed
- S9 golds: digamma rows at p >= 16384 (10 rows) run only with `MPFR_TIER=full`, which also raises the budgets to 120 s per row and 900 s per table. ψ(½) at p=65536 spent 44 s (release) computing Euler γ before S3; with S3's binary-splitting γ it is 0.49 s.
