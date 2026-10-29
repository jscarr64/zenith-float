//! Gold limbs are 64-bit hex strings, least significant first. These helpers turn them into
//! this target's words and print any mantissa by value (most significant word first, trailing
//! zeros dropped), so the same golds compare equal on 64- and 32-bit targets.

use zenith_float_num::Word;

/// One 64-bit gold limb as this target's words, least significant first.
pub fn words_of(s: &str) -> Result<impl Iterator<Item = Word>, String> {
    let x = u64::from_str_radix(s, 16).map_err(|e| format!("limb {s}: {e}"))?;
    let mask = u64::from(Word::MAX);
    Ok((0..u64::BITS / Word::BITS)
        .map(move |i| Word::try_from((x >> (i * Word::BITS)) & mask).unwrap_or(0)))
}

/// The mantissa's value as hex: the same text for 64- and 32-bit words.
pub fn hex_of(m: &[Word]) -> String {
    let width = usize::try_from(Word::BITS / 4).unwrap_or(16);
    let s: String = m.iter().rev().map(|w| format!("{w:0width$x}")).collect();
    s.trim_end_matches('0').to_string()
}
