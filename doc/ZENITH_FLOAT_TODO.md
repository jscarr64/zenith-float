# zenith-float — Special Functions To-Do

What Accumath already golds on SoftFloat that this crate must own. Walk top to bottom. Do not mark a row done without a gold on the implementation Accumath will call.

Rules for every item:

- Software limb arithmetic only — no hardware float in any calculation
- Series / AGM / CF with working precision `p_wrk = p + WORD_BIT_SIZE`; `MAX_PREC_RETRY` bounds retries
- Per-op precision at caller-chosen bits (not global `SOFT_PREC`)
- `expr!` leaf where the arity fits
- Domain guards → `NaN` (`InvalidArgument`); do not invent values outside the stated domain
- MPFR oracle under `mpfr-tests` where MPFR has the function
- Accumath `eval.gold` / specials identities must pass after SoftFloat is pointed here

---

## Done

| Function(s) | Date | Notes |
| --- | --- | --- |
| `erf`, `erfc`, `gamma`, `ln_gamma`, integer `J_n` (`n≤1024`), `Consts::euler_gamma` | 2026-08-29 | SoftFloat wired |
| `ei`, `si`, `ci`, `li`, `fresnel_s`, `fresnel_c` | 2026-08-29 | Full \(f,g\) / factorial large-\(x\). MPFR oracles still open |
| `digamma`, `gammainc`, `bessel_j_nu`, `bessel_y`, `bessel_i`, `bessel_k` | 2026-08-29 | SoftFloat wired |
| complete + incomplete elliptic \(K,E,\Pi\) (Carlson, \(m=k^2\), \(x=\sin\varphi\)) | 2026-08-29 | SoftFloat wired |
| `legendre_p`, `assoc_legendre_p` | 2026-08-29 | Cap 48; Condon–Shortley. \(Y_l^m\) stays Accumath-composed |
| `hypergeom_2f1`, `betainc` | 2026-08-29 | Same domain as Accumath; regularized \(I_x\) |

---

## 1. Digamma and incomplete gamma

Unblocks Accumath \(\psi\) and \(\gamma(s,x)\) golds, and later regularized \(P/Q\).

| Function | Definition | Domain | Status |
| --- | --- | --- | --- |
| `digamma(z)` — \(\psi(z)=\Gamma'(z)/\Gamma(z)\) | Recurrence to large \(z\), then \(\ln z-1/(2z)-\sum B_{2n}/(2n z^{2n})\) | \(z>0\); poles / non-positive → `NaN` | ✅ 2026-08-29 |
| `gammainc(s, x)` — lower \(\gamma(s,x)=\int_0^x t^{s-1}e^{-t}\,dt\) | Series \(x^s e^{-x}\sum x^k/(s)_{k+1}\); \(\Gamma(s)\) when \(e^{-x}\) underflows | \(s>0\), \(x\ge 0\) | ✅ 2026-08-29 |

**Golds:** \(\psi(1)=-\gamma\), \(\psi(2)=-\gamma+1\), \(\psi(1/2)=-\gamma-2\ln 2\); \(\gamma(s,0)=0\); \(\gamma(1,1)=1-e^{-1}\).

**`expr!` leaves:** `digamma`, `gammainc`

---

## 2. Full Bessel family

CAPABILITIES §23 / Additions “will not add Y/I/K” is a **backlog, not a permanent reject**. Accumath already ships these.

| Function | Notes | Status |
| --- | --- | --- |
| non-integer `J_ν` | Integer path exists (`n≤1024`). Series in \((x/2)^{ν}\) / \(\Gamma(ν+1)\) | ✅ 2026-08-29 |
| `Y_ν` | Log + second series; half-integer closed forms Accumath already golds | ✅ 2026-08-29 |
| `I_ν` | \(I_{-n}=I_n\); series | ✅ 2026-08-29 |
| `K_ν` | \(x>0\); \(K_{-ν}=K_ν\); Accumath cap \(\lvertν\rvert\le 32\) | ✅ 2026-08-29 |

**Golds:** \(J_0(0)=1\); \(J_{1/2}(\pi/2)\); \(Y_{1/2}\) closed forms; \(I_0(0)=1\); \(K_{1/2}(1)=\sqrt{\pi/2}\,e^{-1}\); Wronskian \(I_0 K_1+I_1 K_0=1\) at \(x=1\); \(x\le 0\) → `NaN` for \(K\).

**`expr!` leaves:** `bessel_j` (already; extend order), `bessel_y`, `bessel_i`, `bessel_k`

---

## 3. Elliptic integrals

Accumath uses Carlson \(R_F,R_C,R_D,R_J\) with parameter \(m=k^2\) and incomplete \(x=\sin\varphi\). Pick **one** convention and gold it; SoftFloat must match.

| Function | Accumath name | Status |
| --- | --- | --- |
| complete \(K,E,\Pi\) | `elliptic_k`, `elliptic_e_complete`, `elliptic_pi_complete` | ✅ 2026-08-29 |
| incomplete \(F,E,\Pi\) | `elliptic_f`, `elliptic_e`, `elliptic_pi` | ✅ 2026-08-29 |

**Golds:** \(K(0)=E(0)=\pi/2\); \(E(1)=1\); complete = incomplete at \(\varphi=\pi/2\); existing `eval.gold` elliptic rows.

**`expr!` leaves:** `elliptic_k`, `elliptic_e`, `elliptic_pi`, `elliptic_f`, `elliptic_e_inc`, `elliptic_pi_inc`

---

## 4. Legendre

| Function | Notes | Status |
| --- | --- | --- |
| `P_n(x)` | Three-term recurrence; named cap (Accumath `LEGENDRE_N_MAX=48`) | ✅ 2026-08-29 |
| `P_n^m(x)` | Condon–Shortley \((-1)^m\); \(m>n\to\) `NaN` | ✅ 2026-08-29 |

Spherical \(Y_l^m\) may stay Accumath-composed once \(P_n^m\) is here.

**`expr!` leaves:** `legendre_p`, `legendre_p_assoc`

---

## 5. Gaussian \({}_2F_1\) and incomplete beta

| Function | Domain | Status |
| --- | --- | --- |
| `hypergeom_2f1(a,b,c,z)` | Series / Gauss / Pfaff as Accumath; no invented \(z>1\) | ✅ 2026-08-29 |
| `betainc(a,b,x)` | Regularized \(I_x(a,b)\) via \({}_2F_1\); \(a>0\), \(b>0\), \(x\in[0,1]\) | ✅ 2026-08-29 |

**Golds:** \({}_2F_1(\cdot;0)=1\); terminating / Gauss / \(\ln 2\) identities already in Accumath.

**`expr!` leaves:** `hypergeom_2f1`; `betainc` if implemented here

---

## Cross-cutting leftovers

- MPFR 1-ULP oracles for `ei` / `si` / `ci` / `li` / Fresnel
- `scripts/ci_full.sh` green before marking a row done
- `LIBRARY.md` §15, `EXPR.md`, `ZENITH_FLOAT_CAPABILITIES.md` §12, this file

Hung searches and invented closed forms stay `NaN` / `InvalidArgument`.
