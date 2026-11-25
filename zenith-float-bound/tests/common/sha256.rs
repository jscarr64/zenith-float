//! SHA-256 (FIPS 180-4) for the committed heavy-tier goldens: a heavy test checks the file it reads against its slice's
//! HEAVY.SHA256 before using it (r6 item b, Paul). Integer-only, no crate (none may be added).

use std::path::Path;

const K: [u32; 64] = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
];
const H0: [u32; 8] = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
];

fn compress(h: &mut [u32; 8], block: &[u8; 64]) {
    let mut w: Vec<u32> = block
        .as_chunks::<4>()
        .0
        .iter()
        .map(|b| u32::from_be_bytes(*b))
        .collect();
    for i in 16..64 {
        let at = |k: usize| w.get(i - k).map_or(0, |x| *x);
        let (a15, a2) = (at(15), at(2));
        let s0 = a15.rotate_right(7) ^ a15.rotate_right(18) ^ (a15 >> 3);
        let s1 = a2.rotate_right(17) ^ a2.rotate_right(19) ^ (a2 >> 10);
        let v = at(16).wrapping_add(s0).wrapping_add(at(7)).wrapping_add(s1);
        w.push(v);
    }
    let [mut a, mut b, mut c, mut d, mut e, mut f, mut g, mut hh] = *h;
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

/// SHA-256 of `data` as 64 lowercase hex digits (sha256sum's form).
pub fn sha256_hex(data: &[u8]) -> String {
    let (blocks, rest) = data.as_chunks::<64>();
    let mut h = H0;
    for b in blocks {
        compress(&mut h, b);
    }
    let mut tail = Vec::with_capacity(128);
    tail.extend_from_slice(rest);
    tail.push(0x80);
    while tail.len() % 64 != 56 {
        tail.push(0);
    }
    let bits = u64::try_from(data.len())
        .unwrap_or(u64::MAX)
        .wrapping_mul(8);
    tail.extend_from_slice(&bits.to_be_bytes());
    for b in tail.as_chunks::<64>().0 {
        compress(&mut h, b);
    }
    h.iter().map(|x| format!("{x:08x}")).collect()
}

/// `dir/name` as text, after checking its SHA-256 against the `<hex>  name` line of `dir/HEAVY.SHA256`.
pub fn read_pinned(dir: &Path, name: &str) -> Result<String, String> {
    let pins = dir.join("HEAVY.SHA256");
    let pins = std::fs::read_to_string(&pins)
        .map_err(|e| format!("cannot read {}: {e}", pins.display()))?;
    let want = pins
        .lines()
        .find_map(|l| {
            l.split_once("  ")
                .filter(|(_, f)| *f == name)
                .map(|(h, _)| h)
        })
        .ok_or_else(|| format!("{}: no HEAVY.SHA256 line for {name}", dir.display()))?;
    let path = dir.join(name);
    let bytes = std::fs::read(&path).map_err(|e| format!("cannot read {}: {e}", path.display()))?;
    let got = sha256_hex(&bytes);
    if got != want {
        return Err(format!(
            "{}: sha256 {got} differs from its HEAVY.SHA256 pin {want}",
            path.display()
        ));
    }
    String::from_utf8(bytes).map_err(|e| format!("{}: {e}", path.display()))
}
