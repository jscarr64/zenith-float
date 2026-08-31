//! Pure-Rust CSV for [`Ieee64Array`] and [`ExactNumArray`].
//!
//! Binary64 cells are unsigned integer bit patterns (same as serde IEEE arrays).
//! Empty / `nan` → [`Ieee64::NAN`]. No hardware IEEE arithmetic.

#[cfg(feature = "std")]
use crate::defs::DEFAULT_P;
use crate::ieee_soft::Ieee64Array;
#[cfg(feature = "std")]
use crate::Consts;
use crate::Error;
#[cfg(feature = "std")]
use crate::ExactNum;
#[cfg(feature = "std")]
use crate::ExactNumArray;
use crate::Ieee64;
#[cfg(feature = "std")]
use crate::Radix;
#[cfg(feature = "std")]
use crate::RoundingMode;
use alloc::string::String;
use alloc::string::ToString;
use alloc::vec::Vec;

/// Maximum data rows accepted by [`Ieee64Array::from_csv_str`].
pub const CSV_MAX_ROWS: usize = 1_048_576;

/// Maximum columns accepted in one CSV row.
pub const CSV_MAX_COLS: usize = 4_096;

#[cfg(feature = "std")]
const P_MARK: &str = "@p=";

impl Ieee64Array {
    /// Encode row-major binary64 bit patterns as unsigned decimals.
    pub fn to_csv_string(&self) -> Result<String, Error> {
        if self.shape().0 > CSV_MAX_ROWS {
            return Err(Error::InvalidArgument);
        }
        encode_bits_csv(self.shape(), self.as_bits())
    }

    /// Parse CSV of bit-pattern cells. Empty / `nan` → [`Ieee64::NAN`].
    pub fn from_csv_str(text: &str) -> Result<Self, Error> {
        let (rows, cols, bits) = decode_bits_csv(text)?;
        Self::from_parts(rows, cols, bits)
    }

    /// Write [`Self::to_csv_string`] to `path`.
    #[cfg(feature = "std")]
    pub fn to_csv<P: AsRef<std::path::Path>>(&self, path: P) -> Result<(), Error> {
        let s = self.to_csv_string()?;
        std::fs::write(path, s.as_bytes()).map_err(|_| Error::InvalidArgument)
    }

    /// Read [`Self::from_csv_str`] from `path`.
    #[cfg(feature = "std")]
    pub fn from_csv<P: AsRef<std::path::Path>>(path: P) -> Result<Self, Error> {
        let s = std::fs::read_to_string(path).map_err(|_| Error::InvalidArgument)?;
        Self::from_csv_str(&s)
    }
}

#[cfg(feature = "std")]
impl ExactNumArray {
    /// Encode each cell as `Display@p=<bits>`.
    pub fn to_csv_string(&self) -> Result<String, Error> {
        if self.shape().0 > CSV_MAX_ROWS {
            return Err(Error::InvalidArgument);
        }
        let (rows, cols) = self.shape();
        let mut out = String::new();
        for i in 0..rows {
            if i > 0 {
                out.push('\n');
            }
            for j in 0..cols {
                if j > 0 {
                    out.push(',');
                }
                let v = self.get2(i, j).ok_or(Error::InvalidArgument)?;
                out.push_str(&encode_exact(v));
            }
        }
        if rows > 0 {
            out.push('\n');
        }
        Ok(out)
    }

    /// Parse CSV of `Display` / `Display@p=` cells. Empty → NaN.
    pub fn from_csv_str(text: &str) -> Result<Self, Error> {
        let lines = data_lines(text)?;
        if lines.is_empty() {
            return Ok(Self::new(DEFAULT_P));
        }
        let cols = row_width(lines[0])?;
        let rows = lines.len();
        let n = rows.checked_mul(cols).ok_or(Error::InvalidArgument)?;
        let mut vals = Vec::new();
        vals.try_reserve_exact(n)?;
        let mut p = None;
        for line in &lines {
            let fields = split_row(line, cols)?;
            for cell in fields {
                let v = decode_exact_cell(cell)?;
                if p.is_none() {
                    p = Some(v.precision().unwrap_or(DEFAULT_P));
                }
                vals.push(v);
            }
        }
        let p = p.unwrap_or(DEFAULT_P);
        Self::from_parts(p, rows, cols, vals)
    }

    /// Write [`Self::to_csv_string`] to `path`.
    pub fn to_csv<P: AsRef<std::path::Path>>(&self, path: P) -> Result<(), Error> {
        let s = self.to_csv_string()?;
        std::fs::write(path, s.as_bytes()).map_err(|_| Error::InvalidArgument)
    }

    /// Read [`Self::from_csv_str`] from `path`.
    pub fn from_csv<P: AsRef<std::path::Path>>(path: P) -> Result<Self, Error> {
        let s = std::fs::read_to_string(path).map_err(|_| Error::InvalidArgument)?;
        Self::from_csv_str(&s)
    }
}

fn encode_bits_csv(shape: (usize, usize), bits: &[u64]) -> Result<String, Error> {
    let (rows, cols) = shape;
    let n = rows.checked_mul(cols).ok_or(Error::InvalidArgument)?;
    if n != bits.len() {
        return Err(Error::InvalidArgument);
    }
    if cols > CSV_MAX_COLS {
        return Err(Error::InvalidArgument);
    }
    let mut out = String::new();
    for i in 0..rows {
        if i > 0 {
            out.push('\n');
        }
        for j in 0..cols {
            if j > 0 {
                out.push(',');
            }
            let idx = i
                .checked_mul(cols)
                .and_then(|b| b.checked_add(j))
                .ok_or(Error::InvalidArgument)?;
            out.push_str(&bits[idx].to_string());
        }
    }
    if rows > 0 {
        out.push('\n');
    }
    Ok(out)
}

fn decode_bits_csv(text: &str) -> Result<(usize, usize, Vec<u64>), Error> {
    let lines = data_lines(text)?;
    if lines.is_empty() {
        return Ok((0, 0, Vec::new()));
    }
    let cols = row_width(lines[0])?;
    let rows = lines.len();
    let n = rows.checked_mul(cols).ok_or(Error::InvalidArgument)?;
    let mut bits = Vec::new();
    bits.try_reserve_exact(n)?;
    for line in &lines {
        for cell in split_row(line, cols)? {
            bits.push(parse_ieee64_cell(cell)?);
        }
    }
    Ok((rows, cols, bits))
}

fn data_lines(text: &str) -> Result<Vec<&str>, Error> {
    let mut lines: Vec<&str> = text
        .split('\n')
        .map(|l| l.strip_suffix('\r').unwrap_or(l))
        .collect();
    while lines.last().is_some_and(|l| l.is_empty()) {
        lines.pop();
    }
    if lines.len() > CSV_MAX_ROWS {
        return Err(Error::InvalidArgument);
    }
    Ok(lines)
}

fn row_width(line: &str) -> Result<usize, Error> {
    let n = line.split(',').count();
    if n == 0 || n > CSV_MAX_COLS {
        return Err(Error::InvalidArgument);
    }
    Ok(n)
}

fn split_row<'a>(line: &'a str, cols: usize) -> Result<Vec<&'a str>, Error> {
    let mut fields: Vec<&str> = line.split(',').collect();
    if fields.len() > cols {
        return Err(Error::InvalidArgument);
    }
    if fields.len() > CSV_MAX_COLS {
        return Err(Error::InvalidArgument);
    }
    while fields.len() < cols {
        fields.push("");
    }
    Ok(fields)
}

fn parse_ieee64_cell(s: &str) -> Result<u64, Error> {
    let t = s.trim();
    if t.is_empty() || t.eq_ignore_ascii_case("nan") {
        return Ok(Ieee64::NAN.to_bits());
    }
    if t.eq_ignore_ascii_case("inf") || t.eq_ignore_ascii_case("+inf") {
        return Ok(Ieee64::INFINITY.to_bits());
    }
    if t.eq_ignore_ascii_case("-inf") {
        return Ok(Ieee64::NEG_INFINITY.to_bits());
    }
    t.parse::<u64>().map_err(|_| Error::InvalidArgument)
}

#[cfg(feature = "std")]
fn encode_exact(n: &ExactNum) -> String {
    let p = n.precision().unwrap_or(DEFAULT_P);
    let mut s = n.to_string();
    s.push_str(P_MARK);
    s.push_str(&p.to_string());
    s
}

#[cfg(feature = "std")]
fn decode_exact_cell(s: &str) -> Result<ExactNum, Error> {
    let t = s.trim();
    if t.is_empty() {
        return Ok(crate::NAN.clone());
    }
    let (body, p) = match t.rfind(P_MARK) {
        Some(i) => {
            let p = t[i + P_MARK.len()..]
                .parse::<usize>()
                .map_err(|_| Error::InvalidArgument)?;
            (&t[..i], p.max(1))
        }
        None => (t, DEFAULT_P),
    };
    let mut cc = Consts::new()?;
    let n = ExactNum::parse(body, Radix::Dec, p, RoundingMode::ToEven, &mut cc);
    if n.is_nan() && !body.eq_ignore_ascii_case("nan") && !body.eq_ignore_ascii_case("err") {
        return Err(Error::InvalidArgument);
    }
    Ok(n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn csv_ieee64_100x3_roundtrip() {
        let mut vals = Vec::new();
        for i in 0..300 {
            vals.push(Ieee64::from_i32(i));
        }
        let a = Ieee64Array::from_shape(100, 3, &vals).unwrap();
        let text = a.to_csv_string().unwrap();
        let b = Ieee64Array::from_csv_str(&text).unwrap();
        assert_eq!(b.shape(), (100, 3));
        assert_eq!(a.as_bits(), b.as_bits());
    }

    #[test]
    fn csv_missing_cell_is_nan() {
        let b = Ieee64Array::from_csv_str("1,,3\n").unwrap();
        assert_eq!(b.shape(), (1, 3));
        assert_eq!(b.get(0).unwrap().to_bits(), 1);
        assert_eq!(b.get(1).unwrap().to_bits(), Ieee64::NAN.to_bits());
        assert_eq!(b.get(2).unwrap().to_bits(), 3);
        let short = Ieee64Array::from_csv_str("5,6\n7\n").unwrap();
        assert_eq!(short.shape(), (2, 2));
        assert_eq!(short.get2(1, 1).unwrap().to_bits(), Ieee64::NAN.to_bits());
    }

    #[test]
    fn csv_exact_num_roundtrip_and_too_many_cols() {
        let p = 64;
        let n = |k: u8| ExactNum::from_u8(k, p);
        let a = ExactNumArray::from_shape(p, 2, 2, &[n(1), n(2), n(3), n(4)]).unwrap();
        let text = a.to_csv_string().unwrap();
        assert!(text.contains("@p=64"));
        let b = ExactNumArray::from_csv_str(&text).unwrap();
        assert_eq!(b.shape(), (2, 2));
        assert_eq!(b.precision(), p);
        for i in 0..4 {
            assert_eq!(a.get(i).unwrap().cmp(b.get(i).unwrap()), Some(0));
        }
        assert!(Ieee64Array::from_csv_str("1,2\n3,4,5\n").is_err());
    }

    #[cfg(feature = "std")]
    #[test]
    fn csv_path_roundtrip() {
        let a = Ieee64Array::from_shape(
            2,
            2,
            &[
                Ieee64::from_i32(1),
                Ieee64::from_i32(2),
                Ieee64::from_i32(3),
                Ieee64::from_i32(4),
            ],
        )
        .unwrap();
        let path = std::env::temp_dir().join("zenith_csv_path_gold.csv");
        a.to_csv(&path).unwrap();
        let b = Ieee64Array::from_csv(&path).unwrap();
        let _ = std::fs::remove_file(&path);
        assert_eq!(a.as_bits(), b.as_bits());
        assert!(Ieee64Array::from_csv("/no/such/zenith_csv_missing.csv").is_err());
    }
}
