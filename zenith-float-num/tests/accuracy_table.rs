//! Every public numeric function has an `Accuracy` rustdoc line and a table entry.
//!
//! The class is fixed here before the numeric edits, so a later change cannot
//! drop a Ziv function to `Unspecified` without this test failing.

use std::collections::BTreeSet;

use zenith_float_num::{accuracy_of, Accuracy, ACCURACY};

const ELEMENTARY: &[&str] = &[
    "abs",
    "add",
    "add_full_prec",
    "sub",
    "sub_full_prec",
    "mul",
    "mul_full_prec",
    "div",
    "rem",
    "fma",
    "mul_add",
    "hypot",
    "set_precision",
    "sqrt",
    "cbrt",
    "nth_root",
    "reciprocal",
    "exp",
    "ln",
    "log",
    "log2",
    "log10",
    "log1p",
    "expm1",
    "exp2",
    "exp10",
    "pow",
    "powi",
    "powsi",
    "sin",
    "cos",
    "tan",
    "asin",
    "acos",
    "atan",
    "atan2",
    "sinh",
    "cosh",
    "tanh",
    "asinh",
    "acosh",
    "atanh",
    "sin_cos",
    "sinh_cosh",
    "floor",
    "ceil",
    "int",
    "fract",
    "round",
    "ldexp",
    "scalb",
    "logb",
    "copysign",
    "next_after",
    "two_sum",
    "two_product",
    "fused_sum",
    "fused_dot",
    "polyval",
];

const SKIP: &[&str] = &[
    "contains",
    "new",
    "nan",
    "is_inf_pos",
    "is_inf_neg",
    "is_inf",
    "is_nan",
    "is_int",
    "err",
    "cmp",
    "abs_cmp",
    "inv_sign",
    "is_positive",
    "is_negative",
    "is_subnormal",
    "is_zero",
    "clamp",
    "max",
    "min",
    "signum",
    "parse",
    "format",
    "with_radix",
    "random_normal",
    "classify",
    "exponent",
    "precision",
    "max_value",
    "min_value",
    "min_positive",
    "min_positive_normal",
    "from_word",
    "as_raw_parts",
    "from_raw_parts",
    "from_words",
    "sign",
    "set_exponent",
    "mantissa_max_bit_len",
    "is_inline",
    "set_sign",
    "mantissa_digits",
    "convert_from_radix",
    "convert_to_radix",
    "inexact",
    "inexact_accumulated",
    "set_inexact",
    "try_set_precision",
    "frexp",
    "ilogb",
    "zero",
    "one",
    "i",
    "re",
    "im",
    "conj",
    "from_real",
    "into_parts",
    "mid",
];

fn src(rel: &str) -> Result<String, String> {
    let path = format!("{}/src/{rel}", env!("CARGO_MANIFEST_DIR"));
    let text = std::fs::read_to_string(&path).map_err(|err| format!("read {path}: {err}"))?;
    Ok(match text.find("\n#[cfg(test)]") {
        Some(cut) => text[..cut].to_string(),
        None => text,
    })
}

fn pub_fn_names(text: &str) -> Vec<String> {
    let mut names = Vec::new();
    for line in text.lines() {
        let trimmed = line.trim();
        let Some(rest) = trimmed.strip_prefix("pub fn ") else {
            continue;
        };
        let name: String = rest
            .chars()
            .take_while(|c| c.is_ascii_alphanumeric() || *c == '_')
            .collect();
        if name.is_empty() || name.starts_with("test_") || SKIP.contains(&name.as_str()) {
            continue;
        }
        names.push(name);
    }
    names
}

fn macro_names(text: &str) -> Vec<String> {
    let mut names = Vec::new();
    let bytes = text.as_bytes();
    let mut i = 0;
    while i + 12 < text.len() {
        if !text.is_char_boundary(i) || !text[i..].starts_with("gen_wrapper_") {
            i += 1;
            continue;
        }
        let Some(q) = text[i..].find('"') else {
            break;
        };
        let start = i + q;
        let mut j = start + 1;
        while j < bytes.len() {
            if bytes[j] == b'\\' {
                j += 2;
                continue;
            }
            if bytes[j] == b'"' {
                break;
            }
            j += 1;
        }
        let doc = text.get(start..j).unwrap_or("");
        let after = text.get(j + 1..).unwrap_or("");
        let after_trim = after.trim_start_matches([' ', '\n', '\r', ',']);
        let name: String = after_trim
            .chars()
            .take_while(|c| c.is_ascii_alphanumeric() || *c == '_')
            .collect();
        let keep = ELEMENTARY.contains(&name.as_str())
            || doc.contains("Ziv correct-rounding")
            || doc.contains("not a Ziv leaf");
        if keep && !name.is_empty() {
            names.push(name);
        }
        i = j + 1;
    }
    names
}

fn names_near(text: &str, marker: &str) -> Vec<String> {
    let mut names = Vec::new();
    let mut from = 0;
    while let Some(rel) = text[from..].find(marker) {
        let abs = from + rel;
        let mut end = text.len().min(abs + 2500);
        while end > abs && !text.is_char_boundary(end) {
            end -= 1;
        }
        let window = &text[abs..end];
        if let Some(pos) = window.find("pub fn ") {
            let rest = &window[pos + "pub fn ".len()..];
            let name: String = rest
                .chars()
                .take_while(|c| c.is_ascii_alphanumeric() || *c == '_')
                .collect();
            if !name.is_empty() {
                names.push(name);
            }
        }
        from = abs + marker.len();
    }
    names
}

fn discovered() -> Result<BTreeSet<String>, String> {
    let mut set = BTreeSet::new();
    for name in pub_fn_names(&src("catalog.rs")?)
        .into_iter()
        .chain(pub_fn_names(&src("orthopoly.rs")?))
        .chain(pub_fn_names(&src("heun_fox.rs")?))
    {
        set.insert(format!("ExactNum::{name}"));
    }
    let ext = src("ext.rs")?;
    for name in pub_fn_names(&ext) {
        if ELEMENTARY.contains(&name.as_str()) {
            set.insert(format!("ExactNum::{name}"));
        }
    }
    for name in names_near(&ext, "Ziv correct-rounding")
        .into_iter()
        .chain(names_near(&ext, "not a Ziv leaf"))
        .chain(macro_names(&ext))
    {
        if !SKIP.contains(&name.as_str()) {
            set.insert(format!("ExactNum::{name}"));
        }
    }
    for rel in [
        "complex.rs",
        "complex_airy.rs",
        "complex_bessel.rs",
        "complex_ei.rs",
        "complex_elliptic.rs",
        "complex_hypergeom.rs",
        "complex_special.rs",
    ] {
        for name in pub_fn_names(&src(rel)?) {
            set.insert(format!("ExactComplex::{name}"));
        }
    }
    set.insert("Ball::rad".to_string());
    set.insert("ComplexBall::rad".to_string());
    set.insert("GkResult::abs_err_est".to_string());
    Ok(set)
}

fn phrase_near(blob: &str, path: &str, phrase: &str) -> bool {
    let name = path.rsplit("::").next().unwrap_or(path);
    let needles = [
        format!("pub fn {name}"),
        format!("\n        {name},"),
        format!(", {name},"),
        format!("pub {name}:"),
    ];
    for needle in needles {
        let mut from = 0;
        while let Some(rel) = blob[from..].find(&needle) {
            let abs = from + rel;
            let mut begin = abs.saturating_sub(2500);
            while begin > 0 && !blob.is_char_boundary(begin) {
                begin -= 1;
            }
            if blob[begin..abs].contains(phrase) {
                return true;
            }
            from = abs + needle.len();
        }
    }
    false
}

#[test]
fn accuracy_labels_cover_public_numeric_functions() -> Result<(), Box<dyn std::error::Error>> {
    let mut table = BTreeSet::new();
    for (path, _) in ACCURACY {
        assert!(table.insert((*path).to_string()), "duplicate {path}");
    }
    let found = discovered()?;
    let missing: Vec<_> = found.difference(&table).collect();
    let extra: Vec<_> = table.difference(&found).collect();
    assert!(
        missing.is_empty() && extra.is_empty(),
        "missing from table: {missing:?}\nextra in table: {extra:?}"
    );

    let mut blob = String::new();
    for rel in [
        "ext.rs",
        "catalog.rs",
        "orthopoly.rs",
        "heun_fox.rs",
        "complex.rs",
        "complex_airy.rs",
        "complex_bessel.rs",
        "complex_ei.rs",
        "complex_elliptic.rs",
        "complex_hypergeom.rs",
        "complex_special.rs",
        "ball.rs",
        "gauss_kronrod.rs",
        "accuracy.rs",
    ] {
        blob.push_str(&src(rel)?);
        blob.push('\n');
    }
    for (path, class) in ACCURACY {
        let phrase = class.doc_line();
        assert!(
            phrase_near(&blob, path, &phrase),
            "{path} rustdoc missing {phrase}"
        );
        if let Accuracy::WorkingPrecision { lost_bits } = class {
            assert!(
                phrase.contains(&lost_bits.to_string()),
                "{path} rustdoc must state {lost_bits}"
            );
        }
    }

    assert_eq!(
        accuracy_of("ExactNum::add"),
        Some(Accuracy::CorrectlyRounded)
    );
    assert_eq!(
        accuracy_of("ExactNum::sqrt"),
        Some(Accuracy::CorrectlyRounded)
    );
    assert_eq!(
        accuracy_of("ExactNum::erf"),
        Some(Accuracy::CorrectlyRounded)
    );
    assert_eq!(
        accuracy_of("ExactNum::hypergeom_2f1"),
        Some(Accuracy::CorrectlyRounded)
    );
    assert_eq!(
        accuracy_of("ExactNum::elliptic_f"),
        Some(Accuracy::CorrectlyRounded)
    );
    assert_eq!(
        accuracy_of("ExactNum::clausen_cl3"),
        Some(Accuracy::Unspecified)
    );
    assert_eq!(
        accuracy_of("ExactNum::weber_e"),
        Some(Accuracy::Unspecified)
    );
    assert_eq!(
        accuracy_of("ExactNum::barnes_g"),
        Some(Accuracy::Unspecified)
    );
    assert_eq!(
        accuracy_of("ExactNum::riemann_zeta"),
        Some(Accuracy::Unspecified)
    );
    assert_eq!(
        accuracy_of("ExactNum::lambert_w0"),
        Some(Accuracy::Unspecified)
    );
    assert_eq!(
        accuracy_of("ExactNum::hermite_h"),
        Some(Accuracy::WorkingPrecision { lost_bits: 32 })
    );
    assert_eq!(
        accuracy_of("ExactNum::legendre_p"),
        Some(Accuracy::WorkingPrecision { lost_bits: 32 })
    );
    assert_eq!(accuracy_of("Ball::rad"), Some(Accuracy::Enclosure));
    assert_eq!(
        accuracy_of("GkResult::abs_err_est"),
        Some(Accuracy::Enclosure)
    );
    assert!(ACCURACY
        .iter()
        .all(|(_, class)| *class != Accuracy::Faithful));
    Ok(())
}
