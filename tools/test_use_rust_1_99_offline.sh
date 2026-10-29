#!/usr/bin/env bash
# use_rust_1_99.sh must not need the network when Rust 1.99 (with clippy and rustfmt) is already installed.
# Stub rustup: "toolchain install" fails the way an offline rustup does; "run 1.99 ..." works.
set -u
root="$(cd "$(dirname "$0")/.." && pwd)"
tmp="$(mktemp -d "${TMPDIR:-$HOME}/use199.XXXXXX")"
trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp/bin"
cat > "$tmp/bin/rustup" <<'STUB'
#!/usr/bin/env bash
case "$1" in
  toolchain) echo "error: could not download file (dns error)" >&2; exit 1 ;;
  run)
    shift 2
    case "$1" in
      rustc) echo "rustc 1.99.1 (stub)" ;;
      cargo) echo "clippy 0.1.99 (stub)" ;;
      rustfmt) echo "rustfmt 1.8.0 (stub)" ;;
      *) exit 1 ;;
    esac ;;
  *) exit 1 ;;
esac
STUB
chmod +x "$tmp/bin/rustup"
if PATH="$tmp/bin:$PATH" bash -c ". '$root/scripts/use_rust_1_99.sh' && [ \"\$RUSTUP_TOOLCHAIN\" = 1.99 ]"; then
  echo "ok: installed 1.99 used offline"
else
  echo "FAIL: use_rust_1_99.sh needed the network with 1.99 already installed" >&2
  exit 1
fi
