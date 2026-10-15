# Source from a CI job. The Rust 1.99 line is required (1.99.0, 1.99.1, ...); another minor release is not a substitute.
if ! command -v rustup >/dev/null 2>&1; then
  echo "error: rustup is required to install Rust 1.99; refusing to use another compiler" >&2
  return 1 2>/dev/null || exit 1
fi
if ! rustup toolchain install 1.99 --component clippy,rustfmt; then
  echo "error: Rust 1.99 could not be installed; refusing to substitute another version" >&2
  return 1 2>/dev/null || exit 1
fi
got="$(rustup run 1.99 rustc --version)"
case "$got" in
  "rustc 1.99."*) ;;
  *)
    echo "error: channel 1.99 is ${got}, not the 1.99 line; refusing to substitute another version" >&2
    return 1 2>/dev/null || exit 1
    ;;
esac
export RUSTUP_TOOLCHAIN=1.99
