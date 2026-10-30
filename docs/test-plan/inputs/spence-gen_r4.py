# Spence 2026-10-30: generates r4-required-results.csv from Harper edge CSV (ruling-needed rows). Read-only on inputs.
import csv,os
SRC=os.path.expanduser("~/harper-work/precision-matrix-edges-2026-10-30.csv")
OUT=os.path.expanduser("~/spence-work/pm-review/r4-required-results.csv")
G={
"G1-accessor":"re im mid rad num den coeff coeffs degree shape radix precision emax emin err inexact inexact_accumulated classify as_raw_parts significant_bits cached_bit_len cache_info inner value source deref deref_mut",
"G2-construct-copy":"clone default nan zero one i with with_radix try_new filled filled_2d reshape transpose conj inv_sign soft_neg map_exact invert",
"G3-exact-int-modular":"gcd mod_pow mod_inv miller_rabin pollard_rho soft_cmp compare contains commensurable_shift",
"G4-exact-arith-cr":"sqr scalb two_sum two_product fused_sum fused_dot polyval eval clenshaw ilogb frexp",
"G5-ieee-soft":"soft_add soft_sub soft_mul soft_div",
"G6-special-real":"anger_j barnes_g clausen_cl2 clausen_cl3 kelvin_bei kelvin_bei_nu kelvin_ber kelvin_ber_nu kelvin_kei kelvin_kei_nu kelvin_ker kelvin_ker_nu scorer_gi scorer_hi weber_e fox_h meijer_g phi",
"G7-complex-branch":"arg nth_root ci si li",
"G8-dist":"chi_squared_cdf normal_cdf normal_pdf poisson_pmf student_t_pdf",
"G9-poly-algebra":"compose derivative integral companion_matrix",
"G10-transform-window":"dct idct dst idst fft ifft fft_real ifft_real hann_window hamming_window blackman_window kaiser_window rectangular_window matmul",
"G11-decomposition":"lu_decomp qr_decomp svd_decomp eigen_decomp",
"G12-iterative-numeric":"brent illinois gauss_kronrod_interval integrate_adaptive_gk integrate_adaptive_gk_infinite kronrod_pair wynn_epsilon rk4 rk45_adaptive",
"G13-random":"random_fill random_gaussian random_normal random_seed random_uniform reseed_random seeded_random",
"G14-meta-error-model":"accuracy_of doc_line compute_added_err compute_added_err_near_one complex_cancel_bits ode_min_step",
}
fam2g={f:g for g,s in G.items() for f in s.split()}
CR="correctly rounded value at p in rm (ternary-consistent); oracle MPFR or mpmath Ziv pair p+256/p+512"
RMC="correctly rounded at p in this mode; None = RNDZ of exact value at rp(p,64), bytes identical i686/x86_64 (R6)"
IA="ExactNum::nan(Some(Error::InvalidArgument)) / Err(Error::InvalidArgument)"
OVF="exact/correct bits, else Error::ExponentOverflow(sign) per overflow_modes (Up/FromZero toward +-Inf side, Down/ToZero clamp per overflow_modes)"
UNF="correct bits incl. subnormal per subnormal_ops, else signed zero of the exact sign per mode"
NANP="NaN propagates with first-operand Error kind (nan_error_kinds.rs Rule R)"
def base(g):
  # class -> (required_result, variant/value, source)
  if g in("G1-accessor","G2-construct-copy"):
    ex=("exact: returned value/bits equal the stored field or the defined constructor value; independent of target","exact bits/bool/usize","CONTRACT/doc of the fn + to_bytes identity")
    d={c:ex for c in CL}
    d["prec-boundary"]=("exact at every axis point; p=0 and p>MAX_PRECISION refused where the fn takes p",IA,"tests/precision_limits.rs")
    d["nan"]=("exact: NaN preserved with its Error kind (err() returns that kind)","NaN(kind)","nan_error_kinds.rs")
    d["signed-zero"]=("exact: sign of zero preserved (conj/inv_sign flip it exactly)","-0 / +0","IEEE 754-2019 5.5.1 (copy/negate/abs are sign-exact)")
    d["out-of-domain"]=("refusal for invalid shape/index/radix/p; else exact",IA+" or documented None","fn doc")
    for m in RM: d[m]=("mode-independent: identical bytes in all 7 modes; if the signature has no RoundingMode, row = signature proof checked by second agent","exact","std 1; design 2 verify-signature rule")
    return d
  if g=="G3-exact-int-modular":
    ex=("exact integer/bool/Ordering (no rounding); oracle python-flint/gmpy2 in pinned env","exact","gmpy2 / python-flint (team-env oracle-2026-10-30)")
    d={c:ex for c in CL}
    d["zero"]=("exact: gcd(0,0)=0; mod_pow/mod_inv with modulus 0 refused; miller_rabin(0)=false; pollard_rho(0)=None","DivisionByZero or None per doc","fn docs modular.rs/integer.rs")
    d["out-of-domain"]=("mod_inv non-invertible -> None; pollard_rho prime or cap hit -> None (must never return a non-factor); compare/contains with NaN -> unordered/false","None / false","modular.rs:88,226 doc")
    d["nan"]=("compare/soft_cmp/contains: NaN unordered (None/false); integer fns: n/a by type (proof: ExactInt has no NaN)","None / false","IEEE 754-2019 5.11")
    d["signed-zero"]=("compare: -0 == +0 (Equal); integer fns n/a by type (signature proof)","Ordering::Equal","IEEE 754-2019 5.11")
    for m in RM: d[m]=("mode-independent exact result; signature proof where no rm","exact","std 1")
    d["huge-exp"]=("exact for huge integers; allocation failure -> MemoryAllocation","Error::MemoryAllocation","defs.rs:155 Error")
    return d
  d={c:(CR,"value",src(g)) for c in CL}
  d["nan"]=("NaN result",NANP,"nan_error_kinds.rs")
  d["inf"]=("IEEE-like limit value (e.g. +Inf, 0, or NaN(InvalidArgument) for Inf-Inf / 0*Inf)","+-Inf / 0 / NaN(InvalidArgument)","IEEE 754-2019 7.2, 9.2 limits; DLMF asymptotics for specials")
  d["huge-exp"]=(OVF,"ExponentOverflow(Sign)","overflow_modes.rs")
  d["tiny-exp"]=(UNF,"value / signed 0","subnormal_ops; IEEE 754-2019 7.5")
  d["subnormal-like"]=("unnormalized input either normalized exactly then correctly rounded, or refused by from_raw_parts checks",IA,"from_raw_parts_checks.rs")
  d["zero"]=("exact limit value at 0 (DLMF), pole -> DivisionByZero NaN","value or NaN(DivisionByZero)",src(g))
  d["signed-zero"]=("sign of zero result per MPFR/IEEE sign rules for the underlying op","+-0","IEEE 754-2019 6.3; MPFR doc 4.3")
  d["out-of-domain"]=("typed refusal; pole -> DivisionByZero; outside real domain -> InvalidArgument; non-convergence -> PrecisionRetryExhausted (never a wrong number)","NaN(InvalidArgument)/NaN(DivisionByZero)/PrecisionRetryExhausted","defs.rs:155 Error; DLMF")
  d["tie"]=("exact halfway at p: ToEven even, ToOdd odd, directed modes unaffected; for transcendental fns a true tie is impossible except exact points (e.g. f(0)), which must be tested at those points","value","IEEE 754-2019 4.3; MPFR doc")
  d["prec-boundary"]=("correctly rounded at every axis point incl. p=1,2, word edges, MAX_PRECISION; p=0 / MAX+1 refused",IA,"precision_limits.rs")
  for m in RM: d[m]=(RMC,"value","defs.rs:306 RoundingMode; R6")
  if g=="G5-ieee-soft":
    for c in CL: d[c]=("IEEE 754-2019 result of the soft format, bit-exact (to-nearest-even; fn is fixed RNE)","bits","IEEE 754-2019 5.4, 6, 7")
    d["nan"]=("quiet NaN out; sNaN -> qNaN + invalid",""+"qNaN","IEEE 754-2019 6.2, 7.2")
    d["inf"]=("Inf-Inf, 0*Inf, Inf/Inf -> qNaN; else signed Inf","qNaN / +-Inf","IEEE 754-2019 7.2")
    d["zero"]=("x/0 -> signed Inf (divideByZero); 0/0 -> qNaN","+-Inf / qNaN","IEEE 754-2019 7.3")
    d["signed-zero"]=("(+0)+(-0)=+0 in RNE; x-x=+0; sign of product/quotient = xor","+-0","IEEE 754-2019 6.3")
    d["huge-exp"]=("overflow -> +-Inf under RNE","+-Inf","IEEE 754-2019 7.4")
    d["tiny-exp"]=("gradual underflow, subnormal correctly rounded","subnormal/+-0","IEEE 754-2019 7.5")
    for m in RM: d[m]=("fixed RNE: mode-independent; signature proof (no rm argument)","bits","fn doc ieee_soft/mod.rs:107")
  if g=="G4-exact-arith-cr":
    d["cancellation"]=("exact cancellation gives +0 (-0 only in Down), else correctly rounded; two_sum/two_product: hi=RN(a op b), lo exact error (Knuth/Dekker), lo exactly representable","value pair","IEEE 754-2019 6.3; Muller et al. Handbook FP ch.4")
    d["huge-exp"]=(OVF+"; scalb with i64 exponent extremes: overflow/underflow, never wrap","ExponentOverflow(Sign)","overflow_modes.rs")
  if g=="G9-poly-algebra":
    for c in ("normal","negative","cancellation"): d[c]=("each output coefficient correctly rounded at p (derivative/integral/compose are exact rational maps of inputs)","coeffs","exact rational oracle (python-flint fmpq_poly)")
    d["zero"]=("zero polynomial: degree/derivative per doc (empty or [0]); integral constant 0","exact","poly.rs doc")
  if g=="G10-transform-window":
    for c in("normal","negative","cancellation"): d[c]=("each output element correctly rounded at p of the exact transform/window value (Ziv per element)","element vector","mpmath fft/dct definitions (FFTW conventions for normalisation)")
    d["out-of-domain"]=("empty input / length 0 or 1 window / dimension mismatch -> refusal",IA+" or None","dsp.rs docs")
  if g=="G11-decomposition":
    for c in("normal","negative","cancellation","prec-boundary"): d[c]=("each factor entry correctly rounded at p of the exact factor under the pinned normalisation (LU: partial pivoting, first max |pivot|; QR: diag(R)>=0; SVD: Sigma descending, first nonzero of each u_i positive; eig: ascending, unit-norm, sign-normalised); OR refusal","entries or None","LAPACK getrf/geqrf/gesvd conventions; UNCLEAR, see summary")
    d["out-of-domain"]=("singular (LU zero pivot), non-convergence (SVD_ITER_MAX), NaN/Inf entry, empty -> None","None (Option)","array.rs:1310-1472 docs")
    d["nan"]=d["inf"]=("NaN/Inf entry -> None","None","array.rs:1470 doc")
  if g=="G12-iterative-numeric":
    d["normal"]=("result correct: roots/integrals correctly rounded at p (Ziv-certified), rk4 = correctly rounded exact fixed-step RK4 recurrence; rk45_adaptive = true solution within requested tol, else refusal","value","mpmath findroot/quad/odefun; UNCLEAR for ode")
    d["out-of-domain"]=("no sign change (brent/illinois), divergent integral, step underflow (ode_min_step), p<GK_MIN_PREC=32 -> refusal","PrecisionRetryExhausted / InvalidArgument / None per signature","gauss_kronrod.rs:52; roots.rs doc")
  if g=="G13-random":
    for c in CL: d[c]=("deterministic: for a fixed seed the stream bytes are identical on i686/x86_64 and equal the committed digest; values lie in documented range/distribution (KS test digest-pinned)","digest","test_rng.rs DEFAULT_RANDOM_SEED; golds SHA256SUMS")
    d["prec-boundary"]=("every axis p: exactly p significant bits, uniform in [0,1) or documented range; p=0 refused",IA,"precision_limits.rs")
    for m in RM: d[m]=("mode-independent where no rm; else values rounded per mode and identical bytes across targets","digest","std 1")
  if g=="G14-meta-error-model":
    for c in CL: d[c]=("exact deterministic integer/string from the documented formula; error bounds must be upper bounds (never under-estimate), checked against exact computation","exact usize/str","macro_util.rs / accuracy.rs docs")
    d["huge-exp"]=("exponent extremes: no wrap (checked arithmetic); result saturates to documented cap or refusal",IA,"Definition of done 4 (no saturating_ in new code)")
  if g=="G7-complex-branch":
    d["signed-zero"]=("arg(+-0+0i)=+-0, arg(+-0-0i)... per atan2: atan2(+-0,+0)=+-0, atan2(+-0,-0)=+-pi; nth_root/li/ci/si on cut: side chosen by sign of imaginary zero (conj-symmetric) per C Annex G clog/cpow","+-0 / +-pi","C11 Annex G.6.3.2, IEEE 754-2019 9.2.1 atan2; UNCLEAR vs LIBRARY.md:399 (ln(-1)=i*pi only)")
    d["out-of-domain"]=("Ci(0), li(1) pole -> NaN(DivisionByZero) or documented NaN; cut points: principal value with +0 imag side",IA+" / DivisionByZero","LIBRARY.md:399-410; DLMF 6.2, 6.11")
  if g=="G8-dist":
    d["out-of-domain"]=("sigma<=0, df<=0, lambda<0, non-integer k -> refusal; x outside support -> exact 0 (pdf/pmf) or 0/1 (cdf)",IA+" or exact 0/1","DLMF 7.18, 8.11; scipy.stats support conventions")
    d["inf"]=("cdf(-Inf)=0, cdf(+Inf)=1, pdf(+-Inf)=+0","exact 0/1","limits")
  return d
CL="normal prec-boundary zero signed-zero negative huge-exp tiny-exp subnormal-like nan inf out-of-domain tie cancellation rm-Up rm-Down rm-ToZero rm-FromZero rm-ToEven rm-None".split()
RM=[c for c in CL if c.startswith("rm-")]
def src(g): return {"G6-special-real":"DLMF (ch 11 Kelvin/Anger/Weber/Scorer, 5.17 Barnes G, 25.12 Clausen, 16.17 Meijer G); mpmath","G8-dist":"DLMF 7, 8; mpmath","G7-complex-branch":"DLMF 6; C11 Annex G"}.get(g,"MPFR doc / mpmath Ziv")
rows=[x for x in csv.DictReader(open(SRC)) if x["status"]=="ruling-needed"]
miss=set(x["family"] for x in rows)-set(fam2g); assert not miss, miss
cache={}
with open(OUT,"w",newline="") as f:
  w=csv.writer(f); w.writerow(["family","edge_class","required_result","error_variant_or_value","convention_source"])
  n=0
  for x in rows:
    g=fam2g[x["family"]]; d=cache.setdefault(g,base(g))
    rr,v,s=d[x["class"]]; w.writerow([x["family"],x["class"],rr,v,s]); n+=1
print("rows",n,"families",len(set(x["family"] for x in rows)),"groups",len(cache))
