//! Utility functions.

use proc_macro2::Span;
use proc_macro2::TokenStream;
use quote::quote;
use syn::spanned::Spanned;
use syn::Error;
use syn::ExprCall;
use zenith_float_num::Consts;
use zenith_float_num::ExactNum;
use zenith_float_num::Radix;
use zenith_float_num::RoundingMode;

pub fn str_to_exact_num_literal(s: &str, span: Span) -> Result<TokenStream, Error> {
    let mut cc = Consts::new().map_err(|e| Error::new(span, format!("{e}")))?;
    let f = ExactNum::parse(s, Radix::Dec, usize::MAX, RoundingMode::ToEven, &mut cc);
    if let Some(err) = f.err() {
        return Err(Error::new(
            span,
            format!("failed to parse ExactNum from {s}: {err}"),
        ));
    }

    if f.inexact() {
        return Err(Error::new(
            span,
            format!("literal {s} is inexact at compile time"),
        ));
    }

    if let Some((m, n, sign, e, inexact)) = f.as_raw_parts() {
        let stoken = if sign.is_positive() {
            quote!(zenith_float::Sign::Pos)
        } else {
            quote!(zenith_float::Sign::Neg)
        };
        Ok(quote!(zenith_float::ExactNum::from_raw_parts(&[#(#m),*], #n, #stoken, #e, #inexact)))
    } else {
        Ok(quote!(zenith_float::ExactNum::nan()))
    }
}

pub fn str_to_exact_num_expr(s: &str, span: Span, cc: &mut Consts) -> Result<TokenStream, Error> {
    let f = ExactNum::parse(s, Radix::Dec, usize::MAX, RoundingMode::ToEven, cc);
    if let Some(err) = f.err() {
        return Err(Error::new(
            span,
            format!("failed to parse ExactNum from {}: {}", s, err),
        ));
    }

    let q = if f.inexact() {
        quote!(zenith_float::macro_util::check_exponent_range(zenith_float::ExactNum::parse(#s, zenith_float::Radix::Dec, p_wrk, zenith_float::RoundingMode::ToEven, cc), emin, emax))
    } else if let Some((m, n, s, e, inexact)) = f.as_raw_parts() {
        let stoken = if s.is_positive() {
            quote!(zenith_float::Sign::Pos)
        } else {
            quote!(zenith_float::Sign::Neg)
        };
        quote!(zenith_float::macro_util::check_exponent_range(zenith_float::ExactNum::from_raw_parts(&[#(#m),*], #n, #stoken, #e, #inexact), emin, emax))
    } else {
        quote!(zenith_float::ExactNum::nan())
    };

    Ok(q)
}

pub fn check_arg_num(narg: usize, expr: &ExprCall) -> Result<(), Error> {
    if expr.args.len() != narg {
        return Err(Error::new(
            expr.func.span(),
            if narg == 1 { "expected 1 argument." } else { "expected 2 arguments." },
        ));
    }
    Ok(())
}
