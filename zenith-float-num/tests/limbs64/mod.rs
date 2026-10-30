//! Gold limbs are 64-bit hex text (16 digits, low limb first); `Word` is 64 bits on 64-bit
//! builds and 32 bits on i686. These two helpers split gold limbs into Words and regroup
//! Words into gold limb text, so one gold file checks both widths. Tests only: the library's
//! limb layout is unchanged.

use zenith_float_num::Word;

/// Hex digits in one gold limb.
const LIMB_HEX: usize = 16;

/// Hex digits in one Word.
fn word_hex() -> Result<usize, String> {
    usize::try_from(Word::BITS / 4).map_err(|e| format!("Word width: {e}"))
}

/// Words per gold limb (1 on 64-bit, 2 on i686).
fn words_per_limb() -> Result<usize, String> {
    let w = word_hex()?;
    match LIMB_HEX.checked_rem(w) {
        Some(0) => LIMB_HEX
            .checked_div(w)
            .ok_or_else(|| format!("Word width {w}")),
        _ => Err(format!(
            "Word width {w} hex digits does not divide a {LIMB_HEX}-digit limb"
        )),
    }
}

/// One gold limb (exactly 16 hex digits) as Words, low Word first.
pub fn words_of_limb(s: &str) -> Result<Vec<Word>, String> {
    if s.len() != LIMB_HEX || !s.bytes().all(|b| b.is_ascii_hexdigit()) {
        return Err(format!("limb {s}: not {LIMB_HEX} hex digits"));
    }
    let w = word_hex()?;
    let mut out = Vec::with_capacity(words_per_limb()?);
    let mut end = LIMB_HEX;
    while let Some(start) = end.checked_sub(w) {
        let chunk = s
            .get(start..end)
            .ok_or_else(|| format!("limb {s}: no digits {start}..{end}"))?;
        out.push(Word::from_str_radix(chunk, 16).map_err(|e| format!("limb {s}: {e}"))?);
        end = start;
    }
    Ok(out)
}

/// Words (low first) as gold limb text (low limb first). The mantissa is top-aligned, so an
/// odd Word count on i686 gets zero low Words to fill the lowest limb.
pub fn limbs_text(words: &[Word]) -> Result<Vec<String>, String> {
    let w = word_hex()?;
    let per = words_per_limb()?;
    let rem = words.len().checked_rem(per).ok_or("no Words per limb")?;
    let pad = match rem {
        0 => 0,
        r => per.checked_sub(r).ok_or("pad underflow")?,
    };
    let mut padded: Vec<Word> = vec![0; pad];
    padded.extend_from_slice(words);
    Ok(padded
        .chunks(per)
        .map(|limb| limb.iter().rev().map(|x| format!("{x:0w$x}")).collect())
        .collect())
}
