# Source from a CI job. Rust 1.99.0 is required; another release is not a substitute.
if ! command -v rustup >/dev/null 2>&1; then
  echo "error: rustup is required to install Rust 1.99.0; refusing to use another compiler" >&2
  return 1 2>/dev/null || exit 1
fi
if ! rustup toolchain install 1.99.0 --component clippy,rustfmt; then
  echo "error: Rust 1.99.0 could not be installed; refusing to substitute another version" >&2
  return 1 2>/dev/null || exit 1
fi
if ! rustup toolchain install 1.99 --component clippy,rustfmt; then
  echo "error: Rust 1.99.0 could not be installed; refusing to substitute another version" >&2
  return 1 2>/dev/null || exit 1
fi
got="$(rustup run 1.99 rustc --version)"
case "$got" in
  *"1.99.0"*) ;;
  *)
    echo "error: channel 1.99 is ${got}, not 1.99.0; refusing to substitute another version" >&2
    return 1 2>/dev/null || exit 1
    ;;
esac
export RUSTUP_TOOLCHAIN=1.99
