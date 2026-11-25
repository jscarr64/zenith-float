//! Golden streams for the MPFR comparison tests (Paul calls 2026-11-24: r6 item 6 option B; 9:13 PM CT format).
//!
//! Every value-producing MPFR call of a test is one record in that test's golden stream. The answers come from
//! `tools/mpfr_golden/gen.py` (python-flint, Ziv loop to the correctly rounded result with an unbounded exponent, run in
//! the pinned team env `oracle-2026-10-30`); `tools/mpfr_golden/gen.sh` regenerates every stream and `golden/SHA256SUMS`.
//!
//! Format (binary, `zfgold2`): the magic line `zfgold2\n`, then one entry per record:
//! - `1` STORED: the query's digest (first 16 bytes of its SHA-256), then the answer: a token count and per token either
//!   `0` value (sign byte, zigzag exponent, mantissa bytes little-endian: `[+-]<hex>p<exp>`), `1` integer (zigzag) or
//!   `2` text (length, bytes). All counts and integers are LEB128 varints.
//! - `2` EXACT: the query's digest only. The op is one of `exact::OPS` (add, sub, mul, div, remainder); its answer is
//!   computed by `exact.rs` with exact big-integer arithmetic. In gen mode every such answer is checked against gen.py's
//!   before the stream is written, so a stored and a computed answer cannot differ.
//! - `3` REPEAT: a varint k, "record k again" (1-based; same query, so the same answer).
//!
//! A stream over 90 MB is split into `<name>.gold.1`, `<name>.gold.2`, ... read in sequence; otherwise it is
//! `<name>.gold`. A test run checks, for each stream:
//! - every file's SHA-256 matches its own `golden/SHA256SUMS` line (checked before the first record is used; the parts
//!   are exactly the ones listed there, numbered 1..n);
//! - every query the test makes has the next record's digest (same op, precision, rounding mode and operands, so the case
//!   list cannot drift), and the record kind matches the op (exact ops are computed, the rest stored);
//! - the test consumed every record (checked when the test thread ends; a short stream aborts the run).
//!
//! With the `mpfr-golden-gen` feature the same calls send each query to `gen.py` instead and write the stream at the end
//! of the test, plus its readable text form (`<query> => <answer>` lines, `@k` repeats) in `golden-full/` (not
//! committed): `tools/mpfr_golden/crosscheck.py` checks that text against MPFR (gmpy2) and against the binary stream.

use crate::mpfr::exact;
use std::cell::RefCell;
use std::io::Write;

pub const GOLDEN_DIR: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/mpfr/golden");
/// Gen mode's readable copy of every stream (not committed; read by crosscheck.py).
#[cfg(feature = "mpfr-golden-gen")]
pub const FULL_DIR: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/mpfr/golden-full");

const MAGIC: &[u8] = b"zfgold2\n";
/// Streams above this many bytes are split into numbered parts of this size.
#[cfg(feature = "mpfr-golden-gen")]
const PART_BYTES: usize = 90_000_000;
const DIGEST: usize = 16;
const T_STORED: u8 = 1;
const T_EXACT: u8 = 2;
const T_REPEAT: u8 = 3;
const A_VALUE: u8 = 0;
const A_INT: u8 = 1;
const A_TEXT: u8 = 2;

struct Stream {
    name: String,
    /// Check mode: the stream's bytes (all parts, in order).
    #[cfg(not(feature = "mpfr-golden-gen"))]
    bytes: Vec<u8>,
    /// Check mode: per record, the offset of its digest and of its stored answer (`None`: an exact op); repeats resolved.
    #[cfg(not(feature = "mpfr-golden-gen"))]
    recs: Vec<(usize, Option<usize>)>,
    #[cfg(not(feature = "mpfr-golden-gen"))]
    total: usize,
    /// Gen mode: the binary stream being written.
    #[cfg(feature = "mpfr-golden-gen")]
    out: Vec<u8>,
    /// Gen mode: its readable text form.
    #[cfg(feature = "mpfr-golden-gen")]
    full: String,
    /// Gen mode: query -> its 1-based record number, for the repeats.
    #[cfg(feature = "mpfr-golden-gen")]
    seen: std::collections::HashMap<String, usize>,
    used: usize,
    panicked: bool,
    #[cfg(feature = "mpfr-golden-gen")]
    child: Option<(
        std::process::Child,
        std::io::BufReader<std::process::ChildStdout>,
    )>,
}

thread_local! {
    static STREAM: RefCell<Option<Stream>> = const { RefCell::new(None) };
}

static HOOK: std::sync::Once = std::sync::Once::new();

fn install_panic_mark() {
    HOOK.call_once(|| {
        let prev = std::panic::take_hook();
        std::panic::set_hook(Box::new(move |info| {
            let _ = STREAM.try_with(|s| {
                if let Ok(mut b) = s.try_borrow_mut() {
                    if let Some(st) = b.as_mut() {
                        st.panicked = true;
                    }
                }
            });
            prev(info);
        }));
    });
}

/// The value of a harness step that must succeed. A failure fails the test right here (an assert, so the
/// message names the step), and the default value it returns is never used.
#[track_caller]
pub fn must<T: Default>(r: Result<T, String>) -> T {
    assert_eq!(r.as_ref().err(), None, "mpfr test harness");
    r.unwrap_or_default()
}

/// `module__test` for the current libtest thread (`mpfr::module::test`).
fn stream_name() -> Result<String, String> {
    let t = std::thread::current();
    let n = t.name().unwrap_or("");
    let m = n.strip_prefix("mpfr::").ok_or_else(|| {
        format!("golden stream: thread name {n:?} is not an mpfr test (run with at least 2 test threads)")
    })?;
    Ok(m.replace("::", "__"))
}

/// The query's digest: the first 16 bytes of its SHA-256.
fn digest(query: &str) -> [u8; DIGEST] {
    let h = sha256(query.as_bytes());
    let mut d = [0u8; DIGEST];
    for (o, x) in d.iter_mut().zip(h) {
        *o = x;
    }
    d
}

fn put_varint(out: &mut Vec<u8>, v: u128) {
    let mut v = v;
    loop {
        let b = (v & 0x7f) as u8;
        v >>= 7;
        if v == 0 {
            out.push(b);
            return;
        }
        out.push(b | 0x80);
    }
}

fn get_byte(b: &[u8], at: &mut usize) -> Result<u8, String> {
    let x = *b
        .get(*at)
        .ok_or_else(|| format!("golden: truncated at byte {at}"))?;
    *at += 1;
    Ok(x)
}

fn get_varint(b: &[u8], at: &mut usize) -> Result<u128, String> {
    let mut v = 0u128;
    let mut sh = 0u32;
    loop {
        let x = get_byte(b, at)?;
        if sh > 126 {
            return Err(format!("golden: varint too long at byte {at}"));
        }
        v |= u128::from(x & 0x7f) << sh;
        if x & 0x80 == 0 {
            return Ok(v);
        }
        sh += 7;
    }
}

fn get_usize(b: &[u8], at: &mut usize) -> Result<usize, String> {
    let v = get_varint(b, at)?;
    usize::try_from(v).map_err(|e| format!("golden: count {v}: {e}"))
}

fn zigzag(v: i128) -> u128 {
    ((v << 1) ^ (v >> 127)) as u128
}

fn unzigzag(u: u128) -> i128 {
    ((u >> 1) as i128) ^ -((u & 1) as i128)
}

/// A value token `[+-]<hex>p<exp>` in gen.py's canonical form (lowercase hex without leading zeros, plain decimal
/// exponent): `(negative, mantissa, exponent)`.
fn value_token(t: &str) -> Option<(bool, crate::mpfr::mp::Big, i128)> {
    use crate::mpfr::mp::Nat;
    let neg = match t.as_bytes().first() {
        Some(b'-') => true,
        Some(b'+') => false,
        _ => return None,
    };
    let (m, e) = t[1..].split_once('p')?;
    if m.is_empty()
        || m.starts_with('0')
        || !m.bytes().all(|c| matches!(c, b'0'..=b'9' | b'a'..=b'f'))
    {
        return None;
    }
    let e: i128 = e.parse().ok()?;
    let big = crate::mpfr::mp::Big::from_hex(m);
    (format!("{}p{e}", big.to_hex()) == t[1..]).then_some((neg, big, e))
}

/// Appends answer text `a` in the binary answer encoding (module doc) and checks that it decodes back to `a`.
fn encode_answer(a: &str, out: &mut Vec<u8>) -> Result<(), String> {
    let start = out.len();
    let toks: Vec<&str> = a.split(' ').collect();
    put_varint(out, toks.len() as u128);
    for t in toks {
        if let Some((neg, m, e)) = value_token(t) {
            out.push(A_VALUE);
            out.push(u8::from(neg));
            put_varint(out, zigzag(e));
            let mut bytes: Vec<u8> = m.limbs().iter().flat_map(|w| w.to_le_bytes()).collect();
            while bytes.last() == Some(&0) {
                bytes.pop();
            }
            put_varint(out, bytes.len() as u128);
            out.append(&mut bytes);
        } else if let Some(i) = t.parse::<i64>().ok().filter(|i| i.to_string() == t) {
            out.push(A_INT);
            put_varint(out, zigzag(i128::from(i)));
        } else {
            out.push(A_TEXT);
            put_varint(out, t.len() as u128);
            out.extend(t.bytes());
        }
    }
    let mut at = start;
    let back = decode_answer(out, &mut at)?;
    if back != a || at != out.len() {
        return Err(format!(
            "golden gen: answer {a} does not round-trip (got {back})"
        ));
    }
    Ok(())
}

/// Decodes one answer at `at` (moved past it) back to gen.py's text.
fn decode_answer(b: &[u8], at: &mut usize) -> Result<String, String> {
    use crate::mpfr::mp::Nat;
    let n = get_usize(b, at)?;
    let mut s = String::new();
    for k in 0..n {
        if k > 0 {
            s.push(' ');
        }
        match get_byte(b, at)? {
            A_VALUE => {
                let neg = get_byte(b, at)? != 0;
                let e = unzigzag(get_varint(b, at)?);
                let len = get_usize(b, at)?;
                let bytes = b
                    .get(*at..*at + len)
                    .ok_or_else(|| format!("golden: truncated mantissa at byte {at}"))?;
                *at += len;
                let words: Vec<u64> = bytes
                    .chunks(8)
                    .map(|c| {
                        let mut w = [0u8; 8];
                        for (o, x) in w.iter_mut().zip(c) {
                            *o = *x;
                        }
                        u64::from_le_bytes(w)
                    })
                    .collect();
                let m = crate::mpfr::mp::Big::from_vec(words);
                s.push(if neg { '-' } else { '+' });
                s.push_str(&m.to_hex());
                s.push('p');
                s.push_str(&e.to_string());
            }
            A_INT => s.push_str(&unzigzag(get_varint(b, at)?).to_string()),
            A_TEXT => {
                let len = get_usize(b, at)?;
                let t = b
                    .get(*at..*at + len)
                    .ok_or_else(|| format!("golden: truncated text at byte {at}"))?;
                *at += len;
                s.push_str(std::str::from_utf8(t).map_err(|e| format!("golden: text token: {e}"))?);
            }
            x => return Err(format!("golden: unknown answer token tag {x} at byte {at}")),
        }
    }
    Ok(s)
}

/// The files of stream `name` with their SHA256SUMS hashes: `<name>.gold`, or the parts `<name>.gold.1..n` in order.
#[cfg(not(feature = "mpfr-golden-gen"))]
fn stream_files(name: &str) -> Result<Vec<(String, String)>, String> {
    let sums = std::fs::read_to_string(format!("{GOLDEN_DIR}/SHA256SUMS"))
        .map_err(|e| format!("golden: cannot read SHA256SUMS: {e}"))?;
    files_in_sums(&sums, name)
}

/// `stream_files` on the text of a SHA256SUMS file.
fn files_in_sums(sums: &str, name: &str) -> Result<Vec<(String, String)>, String> {
    let whole = format!("{name}.gold");
    let part = format!("{name}.gold.");
    let mut single = None;
    let mut parts: Vec<(usize, String, String)> = Vec::new();
    for l in sums.lines() {
        let Some((h, f)) = l.split_once("  ") else {
            continue;
        };
        if f == whole {
            single = Some((f.to_owned(), h.to_owned()));
        } else if let Some(k) = f.strip_prefix(&part).and_then(|k| k.parse::<usize>().ok()) {
            parts.push((k, f.to_owned(), h.to_owned()));
        }
    }
    parts.sort();
    match (single, parts.is_empty()) {
        (Some(s), true) => Ok(vec![s]),
        (None, false) => {
            for (i, (k, f, _)) in parts.iter().enumerate() {
                if *k != i + 1 {
                    return Err(format!(
                        "golden: {f}: parts of {name} are not numbered 1..n"
                    ));
                }
            }
            Ok(parts.into_iter().map(|(_, f, h)| (f, h)).collect())
        }
        (Some(_), false) => Err(format!(
            "golden: {name} has both a whole file and parts in SHA256SUMS"
        )),
        (None, true) => Err(format!("golden: {whole} has no SHA256SUMS line")),
    }
}

#[cfg(not(feature = "mpfr-golden-gen"))]
fn open(name: String) -> Result<Stream, String> {
    let mut bytes = Vec::new();
    for (f, want) in stream_files(&name)? {
        let path = format!("{GOLDEN_DIR}/{f}");
        let mut b = std::fs::read(&path).map_err(|e| format!("golden: cannot read {path}: {e}"))?;
        let got = sha256_hex(&b);
        if got != want {
            return Err(format!(
                "golden: {path} SHA256 {got} differs from SHA256SUMS {want}"
            ));
        }
        bytes.append(&mut b);
    }
    if !bytes.starts_with(MAGIC) {
        return Err(format!("golden: {name} is not a zfgold2 stream"));
    }
    let mut recs: Vec<(usize, Option<usize>)> = Vec::new();
    let mut at = MAGIC.len();
    while at < bytes.len() {
        let n = recs.len() + 1;
        let r = match get_byte(&bytes, &mut at)? {
            T_STORED => {
                let dg = at;
                at += DIGEST;
                let ans = at;
                decode_answer(&bytes, &mut at)?;
                (dg, Some(ans))
            }
            T_EXACT => {
                let dg = at;
                at += DIGEST;
                (dg, None)
            }
            T_REPEAT => {
                let k = get_usize(&bytes, &mut at)?;
                *recs.get(k.wrapping_sub(1)).ok_or_else(|| {
                    format!("golden: {name} record {n}: repeat of {k} is not an earlier record")
                })?
            }
            t => return Err(format!("golden: {name} record {n}: unknown tag {t}")),
        };
        if at > bytes.len() {
            return Err(format!("golden: {name} record {n} is truncated"));
        }
        recs.push(r);
    }
    let total = recs.len();
    Ok(Stream {
        name,
        bytes,
        recs,
        used: 0,
        total,
        panicked: false,
    })
}

#[cfg(feature = "mpfr-golden-gen")]
fn open(name: String) -> Result<Stream, String> {
    use std::process::{Command, Stdio};
    let script = concat!(env!("CARGO_MANIFEST_DIR"), "/../tools/mpfr_golden/gen.py");
    let mut child = Command::new("python3")
        .arg(script)
        .arg("--serve")
        .arg(&name)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .map_err(|e| format!("golden gen: cannot start gen.py: {e}"))?;
    let out = child.stdout.take().ok_or("golden gen: no stdout")?;
    Ok(Stream {
        name,
        out: MAGIC.to_owned(),
        full: String::new(),
        seen: std::collections::HashMap::new(),
        used: 0,
        panicked: false,
        child: Some((child, std::io::BufReader::new(out))),
    })
}

/// Answer of the next record; `query` must be the recorded query.
pub fn answer(query: &str) -> String {
    install_panic_mark();
    must(STREAM.with(|s| {
        let mut b = s.borrow_mut();
        if b.is_none() {
            *b = Some(open(stream_name()?)?);
        }
        let st = b.as_mut().ok_or("golden stream: not open")?;
        next(st, query)
    }))
}

#[cfg(not(feature = "mpfr-golden-gen"))]
fn next(st: &mut Stream, query: &str) -> Result<String, String> {
    let (dg, ans) = *st.recs.get(st.used).ok_or_else(|| {
        format!(
            "golden {}: record {} missing (stream has {}), query {}",
            st.name,
            st.used + 1,
            st.total,
            short(query)
        )
    })?;
    st.used += 1;
    if st.bytes.get(dg..dg + DIGEST) != Some(&digest(query)[..]) {
        return Err(format!(
            "golden {}: record {} query differs (digest)\n made: {}",
            st.name,
            st.used,
            short(query)
        ));
    }
    match (ans, exact::answer(query)?) {
        (Some(at), None) => {
            let mut at = at;
            decode_answer(&st.bytes, &mut at)
        }
        (None, Some(a)) => Ok(a),
        (Some(_), Some(_)) => Err(format!(
            "golden {}: record {} is stored, but its query is an exact op: {}",
            st.name,
            st.used,
            short(query)
        )),
        (None, None) => Err(format!(
            "golden {}: record {} is an exact-op record, but its query is not: {}",
            st.name,
            st.used,
            short(query)
        )),
    }
}

#[cfg(feature = "mpfr-golden-gen")]
fn next(st: &mut Stream, query: &str) -> Result<String, String> {
    use std::io::BufRead;
    let (child, out) = st.child.as_mut().ok_or("golden gen: closed")?;
    let inp = child.stdin.as_mut().ok_or("golden gen: no stdin")?;
    writeln!(inp, "{query}").map_err(|e| format!("golden gen: write: {e}"))?;
    inp.flush().map_err(|e| format!("golden gen: flush: {e}"))?;
    let mut a = String::new();
    out.read_line(&mut a)
        .map_err(|e| format!("golden gen: read: {e}"))?;
    let a = a.trim_end().to_owned();
    if a.is_empty() || a.starts_with("ERROR") {
        return Err(format!("golden gen {}: {query}: {a}", st.name));
    }
    let mine = exact::answer(query)?;
    if let Some(m) = &mine {
        if *m != a {
            return Err(format!(
                "golden gen {}: exact.rs differs from gen.py\n query  {query}\n gen.py {a}\n exact  {m}",
                st.name
            ));
        }
    }
    st.used += 1;
    match st.seen.get(query) {
        Some(&k) => {
            st.out.push(T_REPEAT);
            put_varint(&mut st.out, k as u128);
            st.full.push_str(&format!("@{k}\n"));
        }
        None => {
            st.out.push(if mine.is_some() { T_EXACT } else { T_STORED });
            st.out.extend(digest(query));
            if mine.is_none() {
                encode_answer(&a, &mut st.out)?;
            }
            st.full.push_str(&format!("{query} => {a}\n"));
            st.seen.insert(query.to_owned(), st.used);
        }
    }
    Ok(a)
}

#[cfg(not(feature = "mpfr-golden-gen"))]
fn short(s: &str) -> String {
    if s.len() > 300 {
        format!("{}...({} bytes)", &s[..300], s.len())
    } else {
        s.to_owned()
    }
}

/// Gen mode: writes stream `name` (whole, or in parts over `PART_BYTES`) and its text form; removes its old files.
#[cfg(feature = "mpfr-golden-gen")]
fn write_stream(name: &str, out: &[u8], full: &str) -> Result<(), String> {
    let whole = format!("{name}.gold");
    let part = format!("{name}.gold.");
    let dir =
        std::fs::read_dir(GOLDEN_DIR).map_err(|e| format!("golden gen: {GOLDEN_DIR}: {e}"))?;
    for ent in dir {
        let ent = ent.map_err(|e| format!("golden gen: {GOLDEN_DIR}: {e}"))?;
        let f = ent.file_name();
        let f = f.to_string_lossy();
        if f == whole
            || f.strip_prefix(&part)
                .is_some_and(|k| k.bytes().all(|c| c.is_ascii_digit()))
        {
            std::fs::remove_file(ent.path()).map_err(|e| format!("golden gen: remove {f}: {e}"))?;
        }
    }
    let put = |f: &str, b: &[u8]| {
        let path = format!("{GOLDEN_DIR}/{f}");
        std::fs::write(&path, b).map_err(|e| format!("golden gen: cannot write {path}: {e}"))
    };
    if out.len() <= PART_BYTES {
        put(&whole, out)?;
    } else {
        for (k, c) in out.chunks(PART_BYTES).enumerate() {
            put(&format!("{part}{}", k + 1), c)?;
        }
    }
    std::fs::create_dir_all(FULL_DIR).map_err(|e| format!("golden gen: {FULL_DIR}: {e}"))?;
    let path = format!("{FULL_DIR}/{whole}");
    std::fs::write(&path, full).map_err(|e| format!("golden gen: cannot write {path}: {e}"))
}

impl Drop for Stream {
    fn drop(&mut self) {
        if self.panicked {
            return;
        }
        #[cfg(feature = "mpfr-golden-gen")]
        {
            if let Some((mut child, _)) = self.child.take() {
                drop(child.stdin.take());
                let _ = child.wait();
            }
            if let Err(e) = write_stream(&self.name, &self.out, &self.full) {
                let _ = writeln!(std::io::stderr(), "{e}");
                std::process::abort();
            }
        }
        #[cfg(not(feature = "mpfr-golden-gen"))]
        if self.used != self.total {
            let _ = writeln!(
                std::io::stderr(),
                "golden {}: the test used {} of {} records (cases dropped or reordered)",
                self.name,
                self.used,
                self.total
            );
            std::process::abort();
        }
    }
}

/// SHA-256 (FIPS 180-4) as hex, for the golden fingerprints; no crate needed.
pub fn sha256_hex(data: &[u8]) -> String {
    sha256(data).iter().map(|x| format!("{x:02x}")).collect()
}

/// SHA-256 (FIPS 180-4).
pub fn sha256(data: &[u8]) -> [u8; 32] {
    const K: [u32; 64] = [
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4,
        0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe,
        0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f,
        0x4a7484aa, 0x5cb0a9dc, 0x76f988da, 0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
        0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc,
        0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b,
        0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070, 0x19a4c116,
        0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7,
        0xc67178f2,
    ];
    let mut h: [u32; 8] = [
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab,
        0x5be0cd19,
    ];
    // the whole 64-byte blocks of `data` in place, then the padded tail (one or two blocks)
    let (blocks, rest) = data.as_chunks::<64>();
    let mut tail = [0u8; 128];
    tail[..rest.len()].copy_from_slice(rest);
    tail[rest.len()] = 0x80;
    let tail_len = if rest.len() < 56 { 64 } else { 128 };
    let bit_len = (data.len() as u64).wrapping_mul(8);
    tail[tail_len - 8..tail_len].copy_from_slice(&bit_len.to_be_bytes());
    for chunk in blocks.iter().chain(tail[..tail_len].as_chunks::<64>().0) {
        let mut w = [0u32; 64];
        for (i, b) in chunk.as_chunks::<4>().0.iter().enumerate() {
            w[i] = u32::from_be_bytes(*b);
        }
        for i in 16..64 {
            let s0 = w[i - 15].rotate_right(7) ^ w[i - 15].rotate_right(18) ^ (w[i - 15] >> 3);
            let s1 = w[i - 2].rotate_right(17) ^ w[i - 2].rotate_right(19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16]
                .wrapping_add(s0)
                .wrapping_add(w[i - 7])
                .wrapping_add(s1);
        }
        let [mut a, mut b, mut c, mut d, mut e, mut f, mut g, mut hh] = h;
        for (ki, wi) in K.iter().zip(w.iter()) {
            let s1 = e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25);
            let ch = (e & f) ^ (!e & g);
            let t1 = hh
                .wrapping_add(s1)
                .wrapping_add(ch)
                .wrapping_add(*ki)
                .wrapping_add(*wi);
            let s0 = a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22);
            let maj = (a & b) ^ (a & c) ^ (b & c);
            let t2 = s0.wrapping_add(maj);
            hh = g;
            g = f;
            f = e;
            e = d.wrapping_add(t1);
            d = c;
            c = b;
            b = a;
            a = t1.wrapping_add(t2);
        }
        for (x, y) in h.iter_mut().zip([a, b, c, d, e, f, g, hh]) {
            *x = x.wrapping_add(y);
        }
    }
    let mut out = [0u8; 32];
    for (o, w) in out.as_chunks_mut::<4>().0.iter_mut().zip(h) {
        *o = w.to_be_bytes();
    }
    out
}

#[test]
fn sha256_known_answers() {
    assert_eq!(
        sha256_hex(b""),
        "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
    );
    // 56 bytes: the length no longer fits the first padded block (two tail blocks)
    assert_eq!(
        sha256_hex(b"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq"),
        "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1"
    );
    assert_eq!(
        sha256_hex(b"abc"),
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    );
    let long = vec![b'a'; 1_000_000];
    assert_eq!(
        sha256_hex(&long),
        "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0"
    );
}

#[test]
fn answer_encoding_round_trips() {
    for a in [
        "+1p0 -1",
        "-2abp-11 1",
        "nan 0",
        "+inf 0",
        "-0 0",
        "+huge 1",
        "-tiny -1",
        "q3/4",
        "+0p5 0",
        "+1fffffffffffffffffffffffffffffffffp-4611686018427387904 0 -5bb74c48bcba9345p4611686018427387903 1",
    ] {
        let mut out = Vec::new();
        must(encode_answer(a, &mut out));
        let mut at = 0;
        assert_eq!(must(decode_answer(&out, &mut at)), a);
        assert_eq!(at, out.len(), "{a}: bytes left over");
    }
}

#[test]
fn split_streams_are_read_as_their_listed_parts_in_order() {
    let sums = "aa  x.gold.2\nbb  x.gold.1\ncc  y.gold\ndd  x.gold.extra\n";
    let got = must(files_in_sums(sums, "x"));
    let want = [("x.gold.1", "bb"), ("x.gold.2", "aa")];
    assert_eq!(got.len(), want.len());
    for ((f, h), (wf, wh)) in got.iter().zip(want) {
        assert_eq!((f.as_str(), h.as_str()), (wf, wh));
    }
    assert_eq!(must(files_in_sums(sums, "y")).len(), 1);
    // a gap, a whole file next to parts, a stream with no line: refused
    assert!(files_in_sums("aa  x.gold.1\nbb  x.gold.3\n", "x").is_err());
    assert!(files_in_sums("aa  x.gold\nbb  x.gold.1\n", "x").is_err());
    assert!(files_in_sums("aa  y.gold\n", "x").is_err());
}
