### zenith-float-bound heavy S6 matrix and LA3 QR goldens are committed; their tests always run

- `tests/goldens/slices/matrix/matrix_rows_heavy.txt` (58.8 MB) and
  `tests/goldens/slices/qr/la_qr_enc4096_heavy.txt` (3.7 MB) are in the tree, byte for byte the files
  pinned in each slice's `HEAVY.SHA256`. Both are left out of any package (`exclude` in Cargo.toml).
- `s6_heavy_rows_replay_through_zb` and `la3_qr_hilbert_256_contains_by_heavy_tier` lose their
  `#[ignore]`, their `ZB_MATRIX_HEAVY` / `ZB_QR_HEAVY_DIR` variables and their SKIP path: every
  `cargo test` replays all 30 heavy matrix rows and re-decides the Hilbert QR rows against enc4096.
- The full gate's `matrix-heavy` step checks both pins (`heavy_check.sh --pin` now matches the pin
  line for the file's own name) and runs both tests; it fails when a heavy file is missing or unpinned.
