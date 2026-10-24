### Changed
- S9 golds: digamma rows at p >= 16384 (10 rows) run only with `MPFR_TIER=full`, which also raises the budgets to 120 s per row and 900 s per table. ψ(½) at p=65536 spends 44 s (release) computing Euler γ until Rusty's S3 rigorous Euler γ lands.
