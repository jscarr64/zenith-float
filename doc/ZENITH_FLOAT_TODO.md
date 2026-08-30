# zenith-float — Special Functions To-Do

Inventory of specials this crate owns. Walk top to bottom. Do not mark a row done without a gold on the `ExactNum` method callers will use.

Rules for every item:

- Software limb arithmetic only — no hardware float in any calculation
- Series / AGM / CF with working precision `p_wrk = p + WORD_BIT_SIZE`; `MAX_PREC_RETRY` bounds retries (`PrecisionRetryExhausted`, not a domain error)
- Per-op precision at caller-chosen bits
- `expr!` leaf where the arity fits
- Domain guards → `NaN` (`InvalidArgument`); overflow → `±Inf`; do not invent values outside the stated domain
- Hang-gates that a host copied (`n≤48`, \(\lvertν\rvert\le 32\), `n≤1024`) are **not** in this crate
- MPFR oracle under `mpfr-tests` where MPFR has the function

---

## Done

| Function(s) | Date | Notes |
| --- | --- | --- |
| `erf`, `erfc`, `gamma`, `ln_gamma`, integer `J_n`, `Consts::euler_gamma` | 2026-08-29 | |
| `ei`, `si`, `ci`, `li`, `fresnel_s`, `fresnel_c` | 2026-08-29 | PV `Ei` for \(x<0\); `li` for \(x>0\), \(x\neq 1\) |
| `digamma`, `gammainc`, `gammainc_upper`, `bessel_j_nu`, `bessel_y`, `bessel_i`, `bessel_k` | 2026-08-29 | Digamma reflection; \(K\) unrestricted in \(\lvertν\rvert\) except math domain |
| complete + incomplete elliptic \(K,E,\Pi\) (Carlson, \(m=k^2\), \(x=\sin\varphi\)) | 2026-08-29 | \(K(1)=+\infty\); \(K(m>1)=m^{-1/2}K(1/m)\) |
| `legendre_p`, `assoc_legendre_p` | 2026-08-29 | \(p+O(n)\) recurrence; no \(n\) hang-gate |
| `hypergeom_2f1`, `betainc` | 2026-08-29 | Real continuation for \(z\le -1\) when defined |
| `Error::PrecisionRetryExhausted` | 2026-08-29 | Split from `InvalidArgument` |
| `Ieee32` / `Ieee64` | 2026-08-29 | Software binary32/binary64; bits golds for add/mul/div/sqrt/FMA |
| `Ieee32Array` / `Ieee64Array` / `ExactNumArray` | 2026-08-29 | 1-D elementwise, `dot`, specials via `ExactNum` |
| 2-D arrays + software matmul | 2026-08-30 | Row-major `from_shape` / `get2` / `matmul`; 1-D is shape `(1, n)` |

---

## Cross-cutting leftovers

- More array ufuncs (remaining `ExactNum` specials)
- Integer SIMD for IEEE arrays (still not an FPU)
- MPFR has `eint`, `jn`, `yn`, `digamma` — oracles in `compare_special_fn_test.rs`. No MPFR `si` / `ci` / `li` / Fresnel; those stay identity/series golds
- Certified interval `sin`/`exp` (`Ball` is `add`/`mul` only)
- Complex specials (permanent: no complex kernel)

Hung searches and invented closed forms stay `NaN` / `InvalidArgument`.
