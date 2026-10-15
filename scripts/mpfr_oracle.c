/* Independent MPFR oracle for zenith-float step-1 goldens.
   Reads operations on stdin and prints the rounded value as z * 2^ez.
   The exact operation is done in MPFR with a wide exponent range, then
   rounded with integer arithmetic so subnormals at EMIN stay exact.
   For the five modes MPFR implements, that rounding is checked against
   mpfr_subnormalize. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gmp.h>
#include <mpfr.h>

#define EMIN (-2147483647 - 1)
#define EMAX 2147483647

static int round_p(int p) {
    if (p <= 0) {
        return 0;
    }
    return ((p + 63) / 64) * 64;
}

static int positive_of(int sign) { return sign >= 0; }

static void q_round(mpz_t q, mpz_t numer, mpz_t denom, const char *rm, int positive) {
    mpz_t r, twice;
    mpz_init(r);
    mpz_tdiv_qr(q, r, numer, denom);
    if (mpz_sgn(r) == 0 || strcmp(rm, "None") == 0 || strcmp(rm, "ToZero") == 0) {
        mpz_clear(r);
        return;
    }
    if (strcmp(rm, "FromZero") == 0) {
        mpz_add_ui(q, q, 1);
        mpz_clear(r);
        return;
    }
    if (strcmp(rm, "Up") == 0) {
        if (positive) {
            mpz_add_ui(q, q, 1);
        }
        mpz_clear(r);
        return;
    }
    if (strcmp(rm, "Down") == 0) {
        if (!positive) {
            mpz_add_ui(q, q, 1);
        }
        mpz_clear(r);
        return;
    }
    mpz_init(twice);
    mpz_mul_2exp(twice, r, 1);
    int cmp = mpz_cmp(twice, denom);
    if (cmp > 0) {
        mpz_add_ui(q, q, 1);
    } else if (cmp == 0) {
        int odd = mpz_odd_p(q) ? 1 : 0;
        if ((strcmp(rm, "ToEven") == 0 && odd) || (strcmp(rm, "ToOdd") == 0 && !odd)) {
            mpz_add_ui(q, q, 1);
        }
    }
    mpz_clear(twice);
    mpz_clear(r);
}

static void emit_sign_z(int sign, mpz_t z, long ez) {
    if (mpz_sgn(z) == 0) {
        printf("zero %s\n", sign < 0 ? "-" : "+");
        return;
    }
    printf("finite %s ", sign < 0 ? "-" : "+");
    mpz_out_str(stdout, 10, z);
    printf(" %ld\n", ez);
}

/* Overflow: largest finite when the mode rounds toward zero for this sign, else Inf. */
static void emit_overflow(int sign, int prec, const char *rm) {
    int toward_zero = (strcmp(rm, "ToZero") == 0) || (strcmp(rm, "Down") == 0 && sign >= 0) ||
                      (strcmp(rm, "Up") == 0 && sign < 0);
    if (toward_zero) {
        mpz_t z;
        mpz_init(z);
        mpz_set_ui(z, 1);
        mpz_mul_2exp(z, z, (mp_bitcnt_t)prec);
        mpz_sub_ui(z, z, 1);
        emit_sign_z(sign, z, (long)EMAX - (long)prec);
        mpz_clear(z);
    } else {
        printf("inf %s\n", sign < 0 ? "-" : "+");
    }
}

/* Round sign * z * 2^ez to prec bits. z is non-negative. */
static void emit_rounded(mpz_t z, long ez, int sign, int prec, const char *rm) {
    if (mpz_sgn(z) == 0) {
        printf("zero %s\n", sign < 0 ? "-" : "+");
        return;
    }
    int positive = positive_of(sign);
    long fl = ez + (long)mpz_sizeinbase(z, 2) - 1;
    mpz_t q, numer, denom, limit;
    mpz_init(q);
    mpz_init(numer);
    mpz_init(denom);
    mpz_init(limit);

    if (fl > EMAX) {
        emit_overflow(sign, prec, rm);
        goto done;
    }

    if (fl < EMIN) {
        long exp2 = ez - (EMIN - (long)prec);
        int away = strcmp(rm, "FromZero") == 0 || (strcmp(rm, "Up") == 0 && positive) ||
                   (strcmp(rm, "Down") == 0 && !positive);
        if (exp2 < 0 && (unsigned long)(-exp2) > mpz_sizeinbase(z, 2) + 2) {
            mpz_set_ui(q, away ? 1 : 0);
        } else if (exp2 >= 0) {
            mpz_mul_2exp(q, z, (mp_bitcnt_t)exp2);
        } else {
            mpz_set(numer, z);
            mpz_set_ui(denom, 1);
            mpz_mul_2exp(denom, denom, (mp_bitcnt_t)(-exp2));
            q_round(q, numer, denom, rm, positive);
        }
        if (mpz_sgn(q) <= 0) {
            printf("zero %s\n", sign < 0 ? "-" : "+");
            goto done;
        }
        mpz_set_ui(limit, 1);
        mpz_mul_2exp(limit, limit, (mp_bitcnt_t)prec);
        if (mpz_cmp(q, limit) > 0) {
            fprintf(stderr, "subnormal q too large\n");
            printf("err range\n");
            goto done;
        }
        goto sub_print;
    }

    {
        long e = fl + 1;
        long exp2 = ez + (long)prec - e;
        if (exp2 >= 0) {
            mpz_mul_2exp(numer, z, (mp_bitcnt_t)exp2);
            mpz_set_ui(denom, 1);
        } else {
            mpz_set(numer, z);
            mpz_set_ui(denom, 1);
            mpz_mul_2exp(denom, denom, (mp_bitcnt_t)(-exp2));
        }
        q_round(q, numer, denom, rm, positive);
        mpz_set_ui(limit, 1);
        mpz_mul_2exp(limit, limit, (mp_bitcnt_t)prec);
        if (mpz_cmp(q, limit) >= 0) {
            mpz_tdiv_q_2exp(q, q, 1);
            e += 1;
        }
        if (e > EMAX) {
            emit_overflow(sign, prec, rm);
            goto done;
        }
        emit_sign_z(sign, q, e - (long)prec);
        goto done;
    }

sub_print:
    mpz_set_ui(limit, 1);
    mpz_mul_2exp(limit, limit, (mp_bitcnt_t)prec);
    if (mpz_cmp(q, limit) == 0) {
        mpz_tdiv_q_2exp(q, q, 1);
        emit_sign_z(sign, q, (long)EMIN + 1 - (long)prec);
    } else {
        emit_sign_z(sign, q, (long)EMIN - (long)prec);
    }

done:
    mpz_clear(q);
    mpz_clear(numer);
    mpz_clear(denom);
    mpz_clear(limit);
}

static int set_tok(mpfr_t x, const char *tok, mpfr_prec_t prec) {
    if (strcmp(tok, "+0") == 0 || strcmp(tok, "0") == 0) {
        mpfr_set_zero(x, 1);
        return 1;
    }
    if (strcmp(tok, "-0") == 0) {
        mpfr_set_zero(x, -1);
        return 1;
    }
    if (strcmp(tok, "+inf") == 0 || strcmp(tok, "inf") == 0) {
        mpfr_set_inf(x, 1);
        return 1;
    }
    if (strcmp(tok, "-inf") == 0) {
        mpfr_set_inf(x, -1);
        return 1;
    }
    if (strcmp(tok, "min") == 0 || strcmp(tok, "+min") == 0 || strcmp(tok, "-min") == 0) {
        mpfr_set_ui_2exp(x, 1, (mpfr_exp_t)EMIN - (mpfr_exp_t)prec, MPFR_RNDN);
        if (tok[0] == '-') {
            mpfr_neg(x, x, MPFR_RNDN);
        }
        return 1;
    }
    if (strcmp(tok, "3min") == 0 || strcmp(tok, "+3min") == 0 || strcmp(tok, "-3min") == 0) {
        mpfr_set_ui_2exp(x, 3, (mpfr_exp_t)EMIN - (mpfr_exp_t)prec, MPFR_RNDN);
        if (tok[0] == '-') {
            mpfr_neg(x, x, MPFR_RNDN);
        }
        return 1;
    }
    if (strcmp(tok, "max") == 0 || strcmp(tok, "+max") == 0 || strcmp(tok, "-max") == 0) {
        mpfr_set_ui(x, 1, MPFR_RNDN);
        mpfr_mul_2ui(x, x, (unsigned long)prec, MPFR_RNDN);
        mpfr_sub_ui(x, x, 1, MPFR_RNDN);
        mpfr_mul_2si(x, x, (mpfr_exp_t)EMAX - (mpfr_exp_t)prec, MPFR_RNDN);
        if (tok[0] == '-') {
            mpfr_neg(x, x, MPFR_RNDN);
        }
        return 1;
    }
    return mpfr_set_str(x, tok, 10, MPFR_RNDN) == 0;
}

static mpfr_rnd_t rnd_of(const char *rm) {
    if (strcmp(rm, "ToZero") == 0 || strcmp(rm, "None") == 0) {
        return MPFR_RNDZ;
    }
    if (strcmp(rm, "Up") == 0) {
        return MPFR_RNDU;
    }
    if (strcmp(rm, "Down") == 0) {
        return MPFR_RNDD;
    }
    if (strcmp(rm, "FromZero") == 0) {
        return MPFR_RNDA;
    }
    return MPFR_RNDN;
}

static int do_op(mpfr_t dst, const char *op, mpfr_t a, mpfr_t b, mpfr_t c, mpfr_rnd_t rnd) {
    if (strcmp(op, "add") == 0) {
        return mpfr_add(dst, a, b, rnd);
    }
    if (strcmp(op, "sub") == 0) {
        return mpfr_sub(dst, a, b, rnd);
    }
    if (strcmp(op, "mul") == 0) {
        return mpfr_mul(dst, a, b, rnd);
    }
    if (strcmp(op, "div") == 0) {
        return mpfr_div(dst, a, b, rnd);
    }
    if (strcmp(op, "rem") == 0) {
        return mpfr_fmod(dst, a, b, rnd);
    }
    if (strcmp(op, "sqrt") == 0) {
        return mpfr_sqrt(dst, a, rnd);
    }
    if (strcmp(op, "cbrt") == 0) {
        return mpfr_cbrt(dst, a, rnd);
    }
    if (strcmp(op, "hypot") == 0) {
        return mpfr_hypot(dst, a, b, rnd);
    }
    if (strcmp(op, "fma") == 0) {
        return mpfr_fma(dst, a, b, c, rnd);
    }
    if (strcmp(op, "exp2") == 0) {
        return mpfr_exp2(dst, a, rnd);
    }
    if (strcmp(op, "exp10") == 0) {
        return mpfr_exp10(dst, a, rnd);
    }
    if (strcmp(op, "pow") == 0) {
        return mpfr_pow(dst, a, b, rnd);
    }
    if (strcmp(op, "sqrt2") == 0) {
        return mpfr_sqrt_ui(dst, 2, rnd);
    }
    if (strcmp(op, "phi") == 0) {
        int rc = mpfr_sqrt_ui(dst, 5, rnd);
        if (rc < 0 && !mpfr_number_p(dst)) {
            return rc;
        }
        mpfr_add_ui(dst, dst, 1, rnd);
        mpfr_div_2ui(dst, dst, 1, rnd);
        return 0;
    }
    return -2;
}

/* Round an operand to the stored width before a transcendental, matching parse. */
static void round_operand(mpfr_t x, int prec) {
    if (prec < 2) {
        prec = 2;
    }
    mpfr_prec_round(x, (mpfr_prec_t)prec, MPFR_RNDN);
}

static int set_step2_tok(mpfr_t x, const char *tok, int p, int prec) {
    if (strcmp(tok, "+pow2") == 0 || strcmp(tok, "-pow2") == 0) {
        mpfr_set_ui_2exp(x, 1, -(mpfr_exp_t)p, MPFR_RNDN);
        if (tok[0] == '-') {
            mpfr_neg(x, x, MPFR_RNDN);
        }
        return 1;
    }
    if (!set_tok(x, tok, prec)) {
        return 0;
    }
    round_operand(x, prec);
    return 1;
}

static void emit_component(mpfr_t v, int prec, const char *rm) {
    if (mpfr_nan_p(v)) {
        printf("nan + 0 0");
        return;
    }
    if (mpfr_inf_p(v)) {
        printf("inf %s 0 0", mpfr_signbit(v) ? "-" : "+");
        return;
    }
    if (mpfr_zero_p(v)) {
        printf("zero %s 0 0", mpfr_signbit(v) ? "-" : "+");
        return;
    }
    mpz_t z;
    mpz_init(z);
    int sign = mpfr_signbit(v) ? -1 : 1;
    mpfr_exp_t ez = mpfr_get_z_2exp(z, v);
    if (mpz_sgn(z) < 0) {
        mpz_neg(z, z);
    }
    /* emit_rounded prints a full line; build the same fields without a newline. */
    int positive = sign >= 0;
    long fl = (long)ez + (long)mpz_sizeinbase(z, 2) - 1;
    mpz_t q, numer, denom, limit;
    mpz_init(q);
    mpz_init(numer);
    mpz_init(denom);
    mpz_init(limit);
    if (fl > EMAX) {
        printf("inf %s 0 0", sign < 0 ? "-" : "+");
        goto done;
    }
    if (fl < EMIN) {
        printf("zero %s 0 0", sign < 0 ? "-" : "+");
        goto done;
    }
    {
        long e = fl + 1;
        long exp2 = (long)ez + (long)prec - e;
        if (exp2 >= 0) {
            mpz_mul_2exp(numer, z, (mp_bitcnt_t)exp2);
            mpz_set_ui(denom, 1);
        } else {
            mpz_set(numer, z);
            mpz_set_ui(denom, 1);
            mpz_mul_2exp(denom, denom, (mp_bitcnt_t)(-exp2));
        }
        q_round(q, numer, denom, rm, positive);
        mpz_set_ui(limit, 1);
        mpz_mul_2exp(limit, limit, (mp_bitcnt_t)prec);
        if (mpz_cmp(q, limit) >= 0) {
            mpz_tdiv_q_2exp(q, q, 1);
            e += 1;
        }
        if (mpz_sgn(q) == 0) {
            printf("zero %s 0 0", sign < 0 ? "-" : "+");
        } else {
            printf("finite %s ", sign < 0 ? "-" : "+");
            mpz_out_str(stdout, 10, q);
            printf(" %ld", e - (long)prec);
        }
    }
done:
    mpz_clear(q);
    mpz_clear(numer);
    mpz_clear(denom);
    mpz_clear(limit);
    mpz_clear(z);
}

/* Principal asin on the real cut, and ln, rounded per component.
   The inputs are already rounded to `prec` bits. The transcendental is
   evaluated much wider than that, then each component is rounded once. */
static int emit_complex(const char *op, mpfr_t re, mpfr_t im, int prec, const char *rm) {
    mpfr_prec_t work = (mpfr_prec_t)prec * 4 + 512;
    if (work < 8192) {
        work = 8192;
    }
    mpfr_t wre, wim, out_re, out_im, t, u;
    mpfr_init2(wre, work);
    mpfr_init2(wim, work);
    mpfr_init2(out_re, work);
    mpfr_init2(out_im, work);
    mpfr_init2(t, work);
    mpfr_init2(u, work);
    mpfr_set(wre, re, MPFR_RNDN);
    mpfr_set(wim, im, MPFR_RNDN);
    if (strcmp(op, "cln") == 0) {
        mpfr_sqr(t, wre, MPFR_RNDN);
        mpfr_sqr(u, wim, MPFR_RNDN);
        mpfr_add(t, t, u, MPFR_RNDN);
        mpfr_log(t, t, MPFR_RNDN);
        mpfr_div_2ui(out_re, t, 1, MPFR_RNDN);
        mpfr_atan2(out_im, wim, wre, MPFR_RNDN);
    } else if (strcmp(op, "casin") == 0) {
        int below = mpfr_signbit(wim) ? 1 : 0;
        mpfr_abs(t, wre, MPFR_RNDN);
        if (!(mpfr_cmp_ui(t, 1) > 0)) {
            printf("err domain\n");
            mpfr_clears(wre, wim, out_re, out_im, t, u, (mpfr_ptr)0);
            return 0;
        }
        /* acosh(|x|) = ln(|x| + sqrt(x^2 - 1)). */
        mpfr_sqr(u, t, MPFR_RNDN);
        mpfr_sub_ui(u, u, 1, MPFR_RNDN);
        mpfr_sqrt(u, u, MPFR_RNDN);
        mpfr_add(u, t, u, MPFR_RNDN);
        mpfr_log(u, u, MPFR_RNDN);
        mpfr_const_pi(out_re, MPFR_RNDN);
        mpfr_div_2ui(out_re, out_re, 1, MPFR_RNDN);
        if (mpfr_signbit(wre)) {
            mpfr_neg(out_re, out_re, MPFR_RNDN);
        }
        /* Above the cut matches mpmath: imag is -acosh for x>1 and +acosh for x<-1.
           Below the cut is the conjugate. */
        if (mpfr_signbit(wre)) {
            mpfr_set(out_im, u, MPFR_RNDN);
        } else {
            mpfr_neg(out_im, u, MPFR_RNDN);
        }
        if (below) {
            mpfr_neg(out_im, out_im, MPFR_RNDN);
        }
    } else {
        mpfr_clears(wre, wim, out_re, out_im, t, u, (mpfr_ptr)0);
        return -2;
    }
    printf("c ");
    emit_component(out_re, prec, rm);
    printf(" ");
    emit_component(out_im, prec, rm);
    printf("\n");
    mpfr_clears(wre, wim, out_re, out_im, t, u, (mpfr_ptr)0);
    return 0;
}

int main(void) {
    char line[8192];
    mpfr_set_emin(MPFR_EMIN_MIN);
    mpfr_set_emax(MPFR_EMAX_MAX);

    while (fgets(line, sizeof line, stdin) != NULL) {
        char op[32], rm[32], a_s[512], b_s[512], c_s[512];
        int p = 0;
        a_s[0] = b_s[0] = c_s[0] = 0;
        if (sscanf(line, "%31s %d %31s %511s %511s %511s", op, &p, rm, a_s, b_s, c_s) < 3) {
            printf("err parse\n");
            continue;
        }
        int prec = round_p(p);
        if (prec < 2) {
            prec = 2;
        }
        int sig = prec;
        /* add, sub, mul, div, and chop round at exactly p bits; operand tokens are the zenith inputs at storage precision. */
        int exact_bits = strcmp(op, "chop") == 0 || strcmp(op, "add") == 0 ||
                         strcmp(op, "sub") == 0 || strcmp(op, "mul") == 0 || strcmp(op, "div") == 0;
        if (exact_bits && p >= 2) {
            sig = p;
        }
        mpfr_prec_t wide = (mpfr_prec_t)prec * 4 + 512;
        if (wide < 8192) {
            wide = 8192;
        }
        mpfr_rnd_t rnd = rnd_of(rm);
        mpfr_t a, b, c, exact;
        mpfr_init2(a, wide);
        mpfr_init2(b, wide);
        mpfr_init2(c, wide);
        mpfr_init2(exact, wide);
        int step2 = strcmp(op, "cln") == 0 || strcmp(op, "casin") == 0 || strcmp(op, "exp2") == 0 ||
                    strcmp(op, "exp10") == 0 || strcmp(op, "pow") == 0;
        if (!step2 && a_s[0] && !set_tok(a, a_s, prec)) {
            printf("err operand\n");
            mpfr_clears(a, b, c, exact, (mpfr_ptr)0);
            continue;
        }
        if (!step2 && b_s[0] && strcmp(op, "chop") != 0 && strcmp(op, "parse") != 0 &&
            !set_tok(b, b_s, prec)) {
            printf("err operand\n");
            mpfr_clears(a, b, c, exact, (mpfr_ptr)0);
            continue;
        }
        if (c_s[0] && strcmp(op, "fma") == 0 && !set_tok(c, c_s, prec)) {
            printf("err operand\n");
            mpfr_clears(a, b, c, exact, (mpfr_ptr)0);
            continue;
        }

        int rc = 0;
        if (strcmp(op, "cln") == 0 || strcmp(op, "casin") == 0) {
            if (!set_step2_tok(a, a_s, p, prec) || !set_step2_tok(b, b_s, p, prec)) {
                printf("err operand\n");
            } else if (emit_complex(op, a, b, prec, rm) == -2) {
                printf("err op\n");
            }
            mpfr_clears(a, b, c, exact, (mpfr_ptr)0);
            continue;
        }
        if (strcmp(op, "parse") == 0) {
            rc = mpfr_set_str(exact, a_s, 10, MPFR_RNDN);
        } else if (strcmp(op, "chop") == 0) {
            int src = atoi(b_s);
            int sp = round_p(src);
            if (sp < prec) {
                sp = prec;
            }
            mpfr_set_prec(exact, (mpfr_prec_t)sp);
            rc = mpfr_set_str(exact, a_s, 10, MPFR_RNDN);
        } else if (strcmp(op, "exp2") == 0 || strcmp(op, "exp10") == 0 || strcmp(op, "pow") == 0) {
            if (a_s[0] && !set_step2_tok(a, a_s, p, prec)) {
                printf("err operand\n");
                mpfr_clears(a, b, c, exact, (mpfr_ptr)0);
                continue;
            }
            if (strcmp(op, "pow") == 0 && !set_step2_tok(b, b_s, p, prec)) {
                printf("err operand\n");
                mpfr_clears(a, b, c, exact, (mpfr_ptr)0);
                continue;
            }
            rc = do_op(exact, op, a, b, c, MPFR_RNDN);
        } else {
            rc = do_op(exact, op, a, b, c, MPFR_RNDN);
        }
        if (rc == -2) {
            printf("err op\n");
            mpfr_clears(a, b, c, exact, (mpfr_ptr)0);
            continue;
        }
        if (mpfr_nan_p(exact)) {
            printf("nan\n");
            mpfr_clears(a, b, c, exact, (mpfr_ptr)0);
            continue;
        }
        if (mpfr_inf_p(exact)) {
            printf("inf %s\n", mpfr_signbit(exact) ? "-" : "+");
            mpfr_clears(a, b, c, exact, (mpfr_ptr)0);
            continue;
        }
        if (mpfr_zero_p(exact)) {
            /* Recompute zeros at the library exponent range so the sign matches the mode. */
            mpfr_set_zero(a, 1);
            mpfr_set_zero(b, 1);
            mpfr_set_zero(c, 1);
            mpfr_set_zero(exact, 1);
            mpfr_set_emin(EMIN);
            mpfr_set_emax(EMAX);
            mpfr_set_prec(a, (mpfr_prec_t)prec);
            mpfr_set_prec(b, (mpfr_prec_t)prec);
            mpfr_set_prec(c, (mpfr_prec_t)prec);
            mpfr_set_prec(exact, (mpfr_prec_t)prec);
            if (a_s[0]) {
                set_tok(a, a_s, prec);
            }
            if (b_s[0] && strcmp(op, "parse") != 0 && strcmp(op, "chop") != 0) {
                set_tok(b, b_s, prec);
            }
            if (c_s[0] && strcmp(op, "fma") == 0) {
                set_tok(c, c_s, prec);
            }
            if (strcmp(op, "parse") == 0) {
                mpfr_set_str(exact, a_s, 10, rnd);
            } else if (strcmp(op, "chop") == 0) {
                mpfr_set_str(exact, a_s, 10, rnd);
            } else {
                do_op(exact, op, a, b, c, rnd);
            }
            mpfr_subnormalize(exact, 0, rnd);
            if (mpfr_nan_p(exact)) {
                printf("nan\n");
            } else if (mpfr_inf_p(exact)) {
                printf("inf %s\n", mpfr_signbit(exact) ? "-" : "+");
            } else if (mpfr_zero_p(exact)) {
                printf("zero %s\n", mpfr_signbit(exact) ? "-" : "+");
            } else {
                mpz_t z;
                mpz_init(z);
                mpfr_exp_t ez = mpfr_get_z_2exp(z, exact);
                if (mpz_sgn(z) < 0) {
                    mpz_neg(z, z);
                }
                emit_sign_z(mpfr_signbit(exact) ? -1 : 1, z, (long)ez);
                mpz_clear(z);
            }
            mpfr_set_emin(MPFR_EMIN_MIN);
            mpfr_set_emax(MPFR_EMAX_MAX);
            mpfr_clears(a, b, c, exact, (mpfr_ptr)0);
            continue;
        }

        mpz_t z;
        mpz_init(z);
        int sign = mpfr_signbit(exact) ? -1 : 1;
        mpfr_exp_t ez = mpfr_get_z_2exp(z, exact);
        if (mpz_sgn(z) < 0) {
            mpz_neg(z, z);
        }
        emit_rounded(z, (long)ez, sign, sig, rm);
        mpz_clear(z);
        mpfr_clears(a, b, c, exact, (mpfr_ptr)0);
    }
    return 0;
}
