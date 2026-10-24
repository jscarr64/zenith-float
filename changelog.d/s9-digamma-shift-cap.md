### Fixed
- digamma: the asymptotic Bernoulli term cap rises from 256 to 1024 (`p/8`). With 256 terms the second Ziv pass near the root x0 (p=8448) needed z ≥ 2^22, so the shift loop ran into STIRLING_SHIFT_MAX (2^20 divisions, 16.2 s release) and gave a wrong result; it now does 4095 shifts and the row takes 1.0 s.
