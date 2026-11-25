//! COVERS: common::sha256::sha256_hex, common::sha256::read_pinned
//! Oracle: the FIPS 180-4 / NIST example digests (empty, "abc", the 448-bit two-block message, one million 'a'), and
//! sha256sum for the committed heavy files (their HEAVY.SHA256 lines).

use crate::common::goldens::golden_dir;
use crate::common::sha256::{read_pinned, sha256_hex};

#[test]
fn sha256_matches_the_fips_examples() -> Result<(), String> {
    let long = vec![b'a'; 1_000_000];
    let cases: [(&[u8], &str); 4] = [
        (
            b"",
            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
        ),
        (
            b"abc",
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad",
        ),
        (
            b"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq",
            "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1",
        ),
        (
            &long,
            "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0",
        ),
    ];
    for (data, want) in cases {
        let got = sha256_hex(data);
        if got != want {
            return Err(format!(
                "sha256 of {} bytes: {got}, want {want}",
                data.len()
            ));
        }
    }
    Ok(())
}

#[test]
fn read_pinned_refuses_a_file_without_its_pin() -> Result<(), String> {
    // la_qr.txt is a gate file: it has no HEAVY.SHA256 line, so read_pinned must refuse it rather than read it.
    match read_pinned(&golden_dir().join("slices").join("qr"), "la_qr.txt") {
        Err(e) if e.contains("no HEAVY.SHA256 line") => Ok(()),
        other => Err(format!("read_pinned of an unpinned file: {other:?}")),
    }
}
