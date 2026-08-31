//! `cexpr!` — complex expressions with the same working-precision / cancellation loop as `expr!`.

use proc_macro2::TokenStream;
use quote::quote;
use syn::spanned::Spanned;
use syn::{BinOp, Error, Expr, ExprBinary, ExprCall, ExprLit, ExprPath, ExprUnary, Lit, UnOp};
use zenith_float_num::Consts;

use crate::util::{check_arg_num, str_to_exact_num_expr};
use crate::{MacroInput, SPEC_ADD_ERR};

fn wrap_real(inner: TokenStream) -> TokenStream {
    quote!(zenith_float::ExactComplex::from_real(#inner, p_wrk))
}

fn traverse_binary(
    expr: &ExprBinary,
    err: &mut Vec<usize>,
    cc: &mut Consts,
) -> Result<TokenStream, Error> {
    let left_expr = traverse_expr(&expr.left, err, cc)?;
    let right_expr = traverse_expr(&expr.right, err, cc)?;
    let errs_id = err.len();

    let ts = match expr.op {
        BinOp::Add(_) => {
            err.push(2);
            err.push(2);
            let im_id = errs_id + 1;
            quote!({
                let arg1 = #left_expr;
                let arg2 = #right_expr;
                let ret = zenith_float::ExactComplex::add(
                    &arg1,
                    &arg2,
                    p_wrk,
                    zenith_float::RoundingMode::None,
                );
                let (re_err, im_err) =
                    zenith_float::macro_util::complex_cancel_bits(&arg1, &arg2, &ret, true);
                let mut retry = false;
                if let Some(newerr) = re_err {
                    if errs[#errs_id] < newerr {
                        errs[#errs_id] = newerr;
                        retry = true;
                    }
                }
                if let Some(newerr) = im_err {
                    if errs[#im_id] < newerr {
                        errs[#im_id] = newerr;
                        retry = true;
                    }
                }
                if retry {
                    continue;
                }
                ret
            })
        }
        BinOp::Sub(_) => {
            err.push(2);
            err.push(2);
            let im_id = errs_id + 1;
            quote!({
                let arg1 = #left_expr;
                let arg2 = #right_expr;
                let ret = zenith_float::ExactComplex::sub(
                    &arg1,
                    &arg2,
                    p_wrk,
                    zenith_float::RoundingMode::None,
                );
                let (re_err, im_err) =
                    zenith_float::macro_util::complex_cancel_bits(&arg1, &arg2, &ret, false);
                let mut retry = false;
                if let Some(newerr) = re_err {
                    if errs[#errs_id] < newerr {
                        errs[#errs_id] = newerr;
                        retry = true;
                    }
                }
                if let Some(newerr) = im_err {
                    if errs[#im_id] < newerr {
                        errs[#im_id] = newerr;
                        retry = true;
                    }
                }
                if retry {
                    continue;
                }
                ret
            })
        }
        BinOp::Mul(_) => {
            err.push(3);
            quote!(zenith_float::ExactComplex::mul(
                &(#left_expr),
                &(#right_expr),
                p_wrk,
                zenith_float::RoundingMode::None
            ))
        }
        BinOp::Div(_) => {
            err.push(3);
            quote!(zenith_float::ExactComplex::div(
                &(#left_expr),
                &(#right_expr),
                p_wrk,
                zenith_float::RoundingMode::None
            ))
        }
        _ => {
            return Err(Error::new(
                expr.span(),
                "unexpected binary operator. Only \"+\", \"-\", \"*\", and \"/\" are allowed in cexpr!.",
            ))
        }
    };
    Ok(ts)
}

fn one_arg(
    fun: TokenStream,
    expr: &ExprCall,
    initial_err: usize,
    err: &mut Vec<usize>,
    cc: &mut Consts,
    use_cc: bool,
) -> Result<TokenStream, Error> {
    check_arg_num(1, expr)?;
    let arg = traverse_expr(&expr.args[0], err, cc)?;
    err.push(initial_err);
    Ok(if use_cc {
        quote!(#fun(&(#arg), p_wrk, zenith_float::RoundingMode::None, cc))
    } else {
        quote!(#fun(&(#arg), p_wrk, zenith_float::RoundingMode::None))
    })
}

fn one_arg_real(
    fun: TokenStream,
    expr: &ExprCall,
    initial_err: usize,
    err: &mut Vec<usize>,
    cc: &mut Consts,
    use_cc: bool,
) -> Result<TokenStream, Error> {
    check_arg_num(1, expr)?;
    let arg = traverse_expr(&expr.args[0], err, cc)?;
    err.push(initial_err);
    let call = if use_cc {
        quote!(#fun(&(#arg), p_wrk, zenith_float::RoundingMode::None, cc))
    } else {
        quote!(#fun(&(#arg), p_wrk, zenith_float::RoundingMode::None))
    };
    Ok(quote!(zenith_float::ExactComplex::from_real(#call, p_wrk)))
}

fn three_arg(
    fun: TokenStream,
    expr: &ExprCall,
    initial_err: usize,
    err: &mut Vec<usize>,
    cc: &mut Consts,
) -> Result<TokenStream, Error> {
    check_arg_num(3, expr)?;
    let arg1 = traverse_expr(&expr.args[0], err, cc)?;
    let arg2 = traverse_expr(&expr.args[1], err, cc)?;
    let arg3 = traverse_expr(&expr.args[2], err, cc)?;
    err.push(initial_err);
    Ok(quote!(#fun(
        &(#arg1),
        &(#arg2),
        &(#arg3),
        p_wrk,
        zenith_float::RoundingMode::None
    )))
}

fn three_arg_cc(
    fun: TokenStream,
    expr: &ExprCall,
    initial_err: usize,
    err: &mut Vec<usize>,
    cc: &mut Consts,
) -> Result<TokenStream, Error> {
    check_arg_num(3, expr)?;
    let arg1 = traverse_expr(&expr.args[0], err, cc)?;
    let arg2 = traverse_expr(&expr.args[1], err, cc)?;
    let arg3 = traverse_expr(&expr.args[2], err, cc)?;
    err.push(initial_err);
    Ok(quote!(#fun(
        &(#arg1),
        &(#arg2),
        &(#arg3),
        p_wrk,
        zenith_float::RoundingMode::None,
        cc
    )))
}

fn two_arg(
    fun: TokenStream,
    expr: &ExprCall,
    initial_err: usize,
    err: &mut Vec<usize>,
    cc: &mut Consts,
    use_cc: bool,
) -> Result<TokenStream, Error> {
    check_arg_num(2, expr)?;
    let arg1 = traverse_expr(&expr.args[0], err, cc)?;
    let arg2 = traverse_expr(&expr.args[1], err, cc)?;
    err.push(initial_err);
    Ok(if use_cc {
        quote!(#fun(&(#arg1), &(#arg2), p_wrk, zenith_float::RoundingMode::None, cc))
    } else {
        quote!(#fun(&(#arg1), &(#arg2), p_wrk, zenith_float::RoundingMode::None))
    })
}

fn traverse_call(
    expr: &ExprCall,
    err: &mut Vec<usize>,
    cc: &mut Consts,
) -> Result<TokenStream, Error> {
    let errmes = "unexpected function name. Only \"recip\", \"sqrt\", \"cbrt\", \"root\", \"ln\", \"log2\", \"log10\", \"log\", \"log1p\", \"exp\", \"exp2\", \"exp10\", \"expm1\", \"pow\", \"sin\", \"cos\", \"tan\", \"asin\", \"acos\", \"atan\", \"hypot\", \"fma\", \"mul_add\", \"sinh\", \"cosh\", \"tanh\", \"asinh\", \"acosh\", \"atanh\", \"abs\", \"arg\", \"conj\", \"ldexp\", \"scalb\", \"logb\", \"erf\", \"erfc\", \"gamma\", \"ln_gamma\", \"digamma\", \"ei\", \"si\", \"ci\", \"li\", \"fresnel_s\", \"fresnel_c\", \"ai\", \"bi\", \"bessel_j_nu\", \"bessel_y\", \"bessel_i\", \"bessel_k\", \"elliptic_k\", \"elliptic_e\", \"elliptic_e_inc\", \"elliptic_f\", \"elliptic_pi\", \"elliptic_pi_inc\" are allowed in cexpr!.";
    let Expr::Path(fun) = expr.func.as_ref() else {
        return Err(Error::new(expr.span(), errmes));
    };
    let Some(fname) = fun.path.get_ident() else {
        return Err(Error::new(expr.span(), errmes));
    };
    match fname.to_string().as_str() {
        "recip" => one_arg(
            quote!(zenith_float::ExactComplex::reciprocal),
            expr,
            2,
            err,
            cc,
            false,
        ),
        "sqrt" => one_arg(
            quote!(zenith_float::ExactComplex::sqrt),
            expr,
            1,
            err,
            cc,
            true,
        ),
        "cbrt" => one_arg(
            quote!(zenith_float::ExactComplex::cbrt),
            expr,
            1,
            err,
            cc,
            true,
        ),
        "root" => {
            check_arg_num(2, expr)?;
            let arg = traverse_expr(&expr.args[0], err, cc)?;
            let n = &expr.args[1];
            err.push(1);
            Ok(quote!(zenith_float::ExactComplex::nth_root(
                &(#arg),
                #n as usize,
                p_wrk,
                zenith_float::RoundingMode::None,
                cc
            )))
        }
        "ln" => one_arg(
            quote!(zenith_float::ExactComplex::ln),
            expr,
            SPEC_ADD_ERR,
            err,
            cc,
            true,
        ),
        "log2" => one_arg(
            quote!(zenith_float::ExactComplex::log2),
            expr,
            SPEC_ADD_ERR,
            err,
            cc,
            true,
        ),
        "log10" => one_arg(
            quote!(zenith_float::ExactComplex::log10),
            expr,
            SPEC_ADD_ERR,
            err,
            cc,
            true,
        ),
        "log" => two_arg(
            quote!(zenith_float::ExactComplex::log),
            expr,
            SPEC_ADD_ERR,
            err,
            cc,
            true,
        ),
        "log1p" => one_arg(
            quote!(zenith_float::ExactComplex::log1p),
            expr,
            SPEC_ADD_ERR,
            err,
            cc,
            true,
        ),
        "exp" => one_arg(
            quote!(zenith_float::ExactComplex::exp),
            expr,
            zenith_float_num::EXPONENT_BIT_SIZE + 1,
            err,
            cc,
            true,
        ),
        "exp2" => one_arg(
            quote!(zenith_float::ExactComplex::exp2),
            expr,
            zenith_float_num::EXPONENT_BIT_SIZE + 1,
            err,
            cc,
            true,
        ),
        "exp10" => one_arg(
            quote!(zenith_float::ExactComplex::exp10),
            expr,
            zenith_float_num::EXPONENT_BIT_SIZE + 1,
            err,
            cc,
            true,
        ),
        "expm1" => one_arg(
            quote!(zenith_float::ExactComplex::expm1),
            expr,
            zenith_float_num::EXPONENT_BIT_SIZE + 1,
            err,
            cc,
            true,
        ),
        "pow" => two_arg(
            quote!(zenith_float::ExactComplex::pow),
            expr,
            zenith_float_num::EXPONENT_BIT_SIZE + SPEC_ADD_ERR,
            err,
            cc,
            true,
        ),
        "hypot" => two_arg(
            quote!(zenith_float::ExactComplex::hypot),
            expr,
            1,
            err,
            cc,
            true,
        ),
        "fma" => three_arg(quote!(zenith_float::ExactComplex::fma), expr, 2, err, cc),
        "mul_add" => three_arg(
            quote!(zenith_float::ExactComplex::mul_add),
            expr,
            2,
            err,
            cc,
        ),
        "sin" => one_arg(
            quote!(zenith_float::ExactComplex::sin),
            expr,
            SPEC_ADD_ERR,
            err,
            cc,
            true,
        ),
        "cos" => one_arg(
            quote!(zenith_float::ExactComplex::cos),
            expr,
            SPEC_ADD_ERR,
            err,
            cc,
            true,
        ),
        "tan" => one_arg(
            quote!(zenith_float::ExactComplex::tan),
            expr,
            SPEC_ADD_ERR,
            err,
            cc,
            true,
        ),
        "asin" => one_arg(
            quote!(zenith_float::ExactComplex::asin),
            expr,
            SPEC_ADD_ERR / 2,
            err,
            cc,
            true,
        ),
        "acos" => one_arg(
            quote!(zenith_float::ExactComplex::acos),
            expr,
            SPEC_ADD_ERR / 2,
            err,
            cc,
            true,
        ),
        "atan" => one_arg(
            quote!(zenith_float::ExactComplex::atan),
            expr,
            2,
            err,
            cc,
            true,
        ),
        "sinh" => one_arg(
            quote!(zenith_float::ExactComplex::sinh),
            expr,
            zenith_float_num::EXPONENT_BIT_SIZE + 1,
            err,
            cc,
            true,
        ),
        "cosh" => one_arg(
            quote!(zenith_float::ExactComplex::cosh),
            expr,
            zenith_float_num::EXPONENT_BIT_SIZE + 1,
            err,
            cc,
            true,
        ),
        "tanh" => one_arg(
            quote!(zenith_float::ExactComplex::tanh),
            expr,
            2,
            err,
            cc,
            true,
        ),
        "asinh" => one_arg(
            quote!(zenith_float::ExactComplex::asinh),
            expr,
            2,
            err,
            cc,
            true,
        ),
        "acosh" => one_arg(
            quote!(zenith_float::ExactComplex::acosh),
            expr,
            SPEC_ADD_ERR,
            err,
            cc,
            true,
        ),
        "atanh" => one_arg(
            quote!(zenith_float::ExactComplex::atanh),
            expr,
            SPEC_ADD_ERR,
            err,
            cc,
            true,
        ),
        "abs" => one_arg_real(
            quote!(zenith_float::ExactComplex::abs),
            expr,
            1,
            err,
            cc,
            false,
        ),
        "arg" => one_arg_real(
            quote!(zenith_float::ExactComplex::arg),
            expr,
            SPEC_ADD_ERR,
            err,
            cc,
            true,
        ),
        "conj" => {
            check_arg_num(1, expr)?;
            let arg = traverse_expr(&expr.args[0], err, cc)?;
            err.push(0);
            Ok(quote!(zenith_float::ExactComplex::conj(&(#arg))))
        }
        "ldexp" | "scalb" => {
            check_arg_num(2, expr)?;
            let arg = traverse_expr(&expr.args[0], err, cc)?;
            let n = &expr.args[1];
            err.push(1);
            let fun = if fname == "scalb" {
                quote!(zenith_float::ExactComplex::scalb)
            } else {
                quote!(zenith_float::ExactComplex::ldexp)
            };
            Ok(quote!(#fun(
                &(#arg),
                #n as zenith_float::Exponent,
                p_wrk,
                zenith_float::RoundingMode::None
            )))
        }
        "logb" => one_arg_real(
            quote!(zenith_float::ExactComplex::logb),
            expr,
            1,
            err,
            cc,
            false,
        ),
        "erf" => one_arg(
            quote!(zenith_float::ExactComplex::erf),
            expr,
            2,
            err,
            cc,
            true,
        ),
        "erfc" => one_arg(
            quote!(zenith_float::ExactComplex::erfc),
            expr,
            2,
            err,
            cc,
            true,
        ),
        "gamma" => one_arg(
            quote!(zenith_float::ExactComplex::gamma),
            expr,
            zenith_float_num::EXPONENT_BIT_SIZE + 1,
            err,
            cc,
            true,
        ),
        "ln_gamma" => one_arg(
            quote!(zenith_float::ExactComplex::ln_gamma),
            expr,
            zenith_float_num::EXPONENT_BIT_SIZE + 1,
            err,
            cc,
            true,
        ),
        "digamma" => one_arg(
            quote!(zenith_float::ExactComplex::digamma),
            expr,
            zenith_float_num::EXPONENT_BIT_SIZE + 1,
            err,
            cc,
            true,
        ),
        "ei" => one_arg(
            quote!(zenith_float::ExactComplex::ei),
            expr,
            2,
            err,
            cc,
            true,
        ),
        "si" => one_arg(
            quote!(zenith_float::ExactComplex::si),
            expr,
            2,
            err,
            cc,
            true,
        ),
        "ci" => one_arg(
            quote!(zenith_float::ExactComplex::ci),
            expr,
            2,
            err,
            cc,
            true,
        ),
        "li" => one_arg(
            quote!(zenith_float::ExactComplex::li),
            expr,
            2,
            err,
            cc,
            true,
        ),
        "fresnel_s" => one_arg(
            quote!(zenith_float::ExactComplex::fresnel_s),
            expr,
            2,
            err,
            cc,
            true,
        ),
        "fresnel_c" => one_arg(
            quote!(zenith_float::ExactComplex::fresnel_c),
            expr,
            2,
            err,
            cc,
            true,
        ),
        "ai" => one_arg(
            quote!(zenith_float::ExactComplex::ai),
            expr,
            2,
            err,
            cc,
            true,
        ),
        "bi" => one_arg(
            quote!(zenith_float::ExactComplex::bi),
            expr,
            2,
            err,
            cc,
            true,
        ),
        "bessel_j_nu" => two_arg(
            quote!(zenith_float::ExactComplex::bessel_j_nu),
            expr,
            2,
            err,
            cc,
            true,
        ),
        "bessel_y" => two_arg(
            quote!(zenith_float::ExactComplex::bessel_y),
            expr,
            2,
            err,
            cc,
            true,
        ),
        "bessel_i" => two_arg(
            quote!(zenith_float::ExactComplex::bessel_i),
            expr,
            2,
            err,
            cc,
            true,
        ),
        "bessel_k" => two_arg(
            quote!(zenith_float::ExactComplex::bessel_k),
            expr,
            2,
            err,
            cc,
            true,
        ),
        "elliptic_k" => one_arg(
            quote!(zenith_float::ExactComplex::elliptic_k),
            expr,
            2,
            err,
            cc,
            true,
        ),
        "elliptic_e" => one_arg(
            quote!(zenith_float::ExactComplex::elliptic_e_complete),
            expr,
            2,
            err,
            cc,
            true,
        ),
        "elliptic_f" => two_arg(
            quote!(zenith_float::ExactComplex::elliptic_f),
            expr,
            2,
            err,
            cc,
            true,
        ),
        "elliptic_e_inc" => two_arg(
            quote!(zenith_float::ExactComplex::elliptic_e),
            expr,
            2,
            err,
            cc,
            true,
        ),
        "elliptic_pi" => two_arg(
            quote!(zenith_float::ExactComplex::elliptic_pi_complete),
            expr,
            2,
            err,
            cc,
            true,
        ),
        "elliptic_pi_inc" => three_arg_cc(
            quote!(zenith_float::ExactComplex::elliptic_pi),
            expr,
            2,
            err,
            cc,
        ),
        _ => Err(Error::new(expr.span(), errmes)),
    }
}

fn traverse_path(expr: &ExprPath) -> Result<TokenStream, Error> {
    Ok(if expr.path.is_ident("I") {
        quote!(zenith_float::ExactComplex::i(p_wrk))
    } else if expr.path.is_ident("pi") {
        wrap_real(quote!({ cc.pi(p_wrk, zenith_float::RoundingMode::None) }))
    } else if expr.path.is_ident("e") {
        wrap_real(quote!({ cc.e(p_wrk, zenith_float::RoundingMode::None) }))
    } else if expr.path.is_ident("ln_2") {
        wrap_real(quote!({ cc.ln_2(p_wrk, zenith_float::RoundingMode::None) }))
    } else if expr.path.is_ident("ln_10") {
        wrap_real(quote!({
            cc.ln_10(p_wrk, zenith_float::RoundingMode::None)
        }))
    } else if expr.path.is_ident("sqrt2") {
        wrap_real(quote!({
            cc.sqrt2(p_wrk, zenith_float::RoundingMode::None)
        }))
    } else if expr.path.is_ident("phi") {
        wrap_real(quote!({ cc.phi(p_wrk, zenith_float::RoundingMode::None) }))
    } else if expr.path.is_ident("euler_gamma") {
        wrap_real(quote!({
            cc.euler_gamma(p_wrk, zenith_float::RoundingMode::None)
        }))
    } else {
        quote!({
            let mut arg = zenith_float::ExactComplex::from_ext(
                (#expr).clone(),
                p_wrk,
                zenith_float::RoundingMode::ToEven,
                cc,
            );
            arg.set_inexact(false);
            arg = zenith_float::macro_util::check_complex_exponent_range(arg, emin, emax);
            arg
        })
    })
}

fn traverse_unary(
    expr: &ExprUnary,
    err: &mut Vec<usize>,
    cc: &mut Consts,
) -> Result<TokenStream, Error> {
    let op_expr = traverse_expr(&expr.expr, err, cc)?;
    match expr.op {
        UnOp::Neg(_) => Ok(quote!({
            let z = #op_expr;
            zenith_float::ExactComplex::new(z.re().neg(), z.im().neg())
        })),
        _ => Err(Error::new(
            expr.span(),
            "unexpected unary operator. Only \"-\" is allowed.",
        )),
    }
}

fn traverse_lit(expr: &ExprLit, cc: &mut Consts) -> Result<TokenStream, Error> {
    let inner = match &expr.lit {
        Lit::Str(v) => str_to_exact_num_expr(&v.value(), expr.span(), cc)?,
        Lit::Int(v) => str_to_exact_num_expr(v.base10_digits(), expr.span(), cc)?,
        Lit::Float(v) => str_to_exact_num_expr(v.base10_digits(), expr.span(), cc)?,
        _ => return Err(Error::new(
            expr.span(),
            "unexpected literal. Only string, integer, or floating point literals are supported.",
        )),
    };
    Ok(wrap_real(inner))
}

fn traverse_expr(expr: &Expr, err: &mut Vec<usize>, cc: &mut Consts) -> Result<TokenStream, Error> {
    match expr {
        Expr::Binary(e) => traverse_binary(e, err, cc),
        Expr::Call(e) => traverse_call(e, err, cc),
        Expr::Group(e) => traverse_expr(&e.expr, err, cc),
        Expr::Lit(e) => traverse_lit(e, cc),
        Expr::Paren(e) => traverse_expr(&e.expr, err, cc),
        Expr::Path(e) => traverse_path(e),
        Expr::Unary(e) => traverse_unary(e, err, cc),
        _ => Err(Error::new(expr.span(), "unexpected expression in cexpr!.")),
    }
}

pub fn cexpr(input: proc_macro::TokenStream) -> proc_macro::TokenStream {
    let pmi = syn::parse_macro_input!(input as MacroInput);
    let MacroInput { expr, ctx } = pmi;
    let mut err = Vec::new();
    let mut cc = Consts::new().expect("Failed to initialize constant cache.");
    let expr = traverse_expr(&expr, &mut err, &mut cc).unwrap_or_else(|e| e.to_compile_error());
    let err_sz = err.len();
    quote!({
        use zenith_float::FromExt;
        use zenith_float::ctx::Contextable;

        let mut ctx = &mut (#ctx);
        let p: usize = ctx.precision();
        let rm = ctx.rounding_mode();
        let emin = ctx.emin();
        let emax = ctx.emax();
        let cc = ctx.consts();

        let mut p_rnd = p + zenith_float::WORD_BIT_SIZE;
        let mut errs: [usize; #err_sz] = [#(#err, )*];

        loop {
            let p_wrk = p_rnd.saturating_add(errs.iter().sum());
            let mut ret: zenith_float::ExactComplex = #expr;
            if let Err(err) = ret.set_precision(p, rm) {
                ret = zenith_float::ExactComplex::new(
                    zenith_float::ExactNum::nan(Some(err)),
                    zenith_float::ExactNum::nan(Some(err)),
                );
            }
            break zenith_float::macro_util::check_complex_exponent_range(ret, emin, emax);
        }
    })
    .into()
}
