// Derived from astro-float (https://github.com/stencillogic/astro-float),
// Copyright (c) 2022 stencillogic, MIT License.
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

    raw_parts_expr(&f, span)
}

/// Tokens rebuilding `f` through `ExactNum::from_raw_parts`. The parts are checked at compile
/// time (a rejected literal is a compile error); a runtime rejection yields a NaN carrying the
/// returned error, never a bare NaN.
fn raw_parts_expr(f: &ExactNum, span: Span) -> Result<TokenStream, Error> {
    let raw = f.as_raw_parts();
    if let Err(err) = ExactNum::from_raw_parts(&raw) {
        return Err(Error::new(
            span,
            format!("literal raw parts rejected: {err}"),
        ));
    }
    // A proc macro runs on the host: its words are the host's width, not the target's. Emit the
    // mantissa as 32-bit and as 64-bit words (least significant first) and let the expanding
    // crate's `target_pointer_width` pick, the same `cfg` that picks `zenith_float::Word`.
    let mut m32: Vec<u32> = Vec::new();
    for w in raw.words {
        let x = u64::from(*w);
        let mut k = 0;
        while k < zenith_float_num::Word::BITS {
            m32.push(u32::try_from((x >> k) & 0xffff_ffff).unwrap_or(0));
            k += 32;
        }
    }
    if m32.len() % 2 == 1 {
        // A 32-bit host with an odd word count: one zero low word keeps the value.
        m32.insert(0, 0);
    }
    let m64: Vec<u64> = m32
        .chunks(2)
        .map(|c| match c {
            [lo, hi] => u64::from(*lo) | (u64::from(*hi) << 32),
            _ => 0,
        })
        .collect();
    let (p32, p64) = (m32.len() * 32, m64.len() * 64);
    let (negative, e, inexact, op_inexact) =
        (raw.negative, raw.exponent, raw.inexact, raw.op_inexact);
    let kind = match raw.kind {
        zenith_float_num::RawKind::Finite => quote!(zenith_float::RawKind::Finite),
        zenith_float_num::RawKind::Zero => quote!(zenith_float::RawKind::Zero),
        zenith_float_num::RawKind::Inf => quote!(zenith_float::RawKind::Inf),
        zenith_float_num::RawKind::NaN => quote!(zenith_float::RawKind::NaN),
    };
    let ts = quote!(match zenith_float::ExactNum::from_raw_parts(&zenith_float::RawParts {
        kind: #kind,
        words: {
            #[cfg(target_pointer_width = "32")]
            let w: &[zenith_float::Word] = &[#(#m32),*];
            #[cfg(not(target_pointer_width = "32"))]
            let w: &[zenith_float::Word] = &[#(#m64),*];
            w
        },
        precision: {
            #[cfg(target_pointer_width = "32")]
            let p: usize = #p32;
            #[cfg(not(target_pointer_width = "32"))]
            let p: usize = #p64;
            p
        },
        negative: #negative,
        exponent: #e,
        inexact: #inexact,
        op_inexact: #op_inexact,
    }) {
        Ok(v) => v,
        Err(err) => zenith_float::ExactNum::nan(Some(err)),
    });
    Ok(ts)
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
    } else {
        let r = raw_parts_expr(&f, span)?;
        quote!(zenith_float::macro_util::check_exponent_range(#r, emin, emax))
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

/// The constant cache for a macro expansion, or a `compile_error!` naming the failure.
pub(crate) fn consts_or_compile_error(
    r: Result<zenith_float_num::Consts, zenith_float_num::Error>,
) -> Result<zenith_float_num::Consts, proc_macro2::TokenStream> {
    r.map_err(|e| {
        syn::Error::new(
            proc_macro2::Span::call_site(),
            format!("zenith-float: failed to initialize the constant cache: {e:?}"),
        )
        .to_compile_error()
    })
}

#[cfg(test)]
mod tests {
    use super::consts_or_compile_error;

    #[test]
    fn consts_failure_is_compile_error_not_panic() {
        let r = consts_or_compile_error(Err(zenith_float_num::Error::MemoryAllocation));
        let ts = match r {
            Ok(_) => String::new(),
            Err(ts) => ts.to_string(),
        };
        assert!(ts.contains("compile_error"), "{ts}");
        assert!(ts.contains("MemoryAllocation"), "{ts}");
    }

    #[test]
    fn consts_success_passes_through() {
        assert!(consts_or_compile_error(zenith_float_num::Consts::new()).is_ok());
    }
}
