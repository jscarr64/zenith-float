#!/usr/bin/env bash
# Selftest for tools/no_suppress_scan.sh (Jeff 2026-11-24 1:08/1:17 PM CT: zero #[allow] and zero unsafe in zenith-float;
# no allowlist, no exemptions). Each case is a fixture git repo; a planted suppression must FAIL by name, the same text in
# comments / strings / forbid / deny / identifiers must PASS. Unsafe passes nowhere: the pending inventory
# (tools/unsafe_pending.tsv) is closed at 0 since r6 item 6, any entry fails.
# usage: tools/no_suppress_selftest.sh
set -uo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
scan="$root/tools/no_suppress_scan.sh"
fail=0; n=0
bad() { echo "FAIL: $*"; fail=1; }
scratch="${ZF_SCRATCH:-$HOME/.cache/zenith-float}/no-suppress-selftest"  # job scratch, not /tmp (Spence m6); like fold_pr1_selftest
mkdir -p "$scratch" || exit 2
w=$(mktemp -d "$scratch/run.XXXXXX") || exit 2
trap 'rm -rf "$w"' EXIT
# fixture: a clean two-crate workspace
mk() {
  local r="$w/r$n"; rm -rf "$r"; mkdir -p "$r/zenith-float-num/src" "$r/src" "$r/tests" "$r/tools" "$r/scripts"
  printf '[workspace]\nmembers = ["zenith-float-num"]\n\n[package]\nname = "zenith-float"\nversion = "0.0.0"\n' > "$r/Cargo.toml"
  printf '[package]\nname = "zenith-float-num"\nversion = "0.0.0"\n\n[lints.rust]\nunsafe_code = "forbid"\n' > "$r/zenith-float-num/Cargo.toml"
  printf '#![forbid(unsafe_code)]\n#![deny(unused)]\n/// doc\npub fn f() -> u32 { 1 }\n' > "$r/zenith-float-num/src/lib.rs"
  printf '#![forbid(unsafe_code)]\npub use zenith_float_num::f;\n' > "$r/src/lib.rs"
  printf '#[test]\nfn t() { assert_eq!(1, 1); }\n' > "$r/tests/t.rs"
  printf 'disallowed-macros = []\n' > "$r/clippy.toml"
  printf '# path<TAB>count<TAB>note\n' > "$r/tools/unsafe_pending.tsv"
  echo "$r"
}
commit() { git -C "$1" init -q && git -C "$1" add -A && git -C "$1" -c user.name=t -c user.email=t@t commit -qm f; }
# case WANT LABEL SETUP-SHELL [EXPECTED-OUTPUT-REGEX]
case_() {
  local want=$1 label=$2 setup=$3 rx=${4:-}
  n=$((n + 1)); local r; r=$(mk)
  ( cd "$r" && eval "$setup" ) || { bad "$label: setup failed"; return; }
  commit "$r"
  "$scan" --root "$r" > "$w/out$n" 2>&1; local rc=$?
  if [[ ! -x "$scan" ]]; then bad "$label: $scan missing"; return; fi
  if [[ $want == pass ]]; then
    [[ $rc -eq 0 ]] || bad "$label: want PASS, rc=$rc: $(tr '\n' '|' < "$w/out$n" | cut -c1-300)"
  else
    [[ $rc -eq 1 ]] || bad "$label: want FAIL rc=1, rc=$rc: $(tr '\n' '|' < "$w/out$n" | cut -c1-300)"
    [[ -z $rx ]] || grep -qE "$rx" "$w/out$n" || bad "$label: output lacks /$rx/: $(tr '\n' '|' < "$w/out$n" | cut -c1-300)"
  fi
}
L=zenith-float-num/src/lib.rs
case_ pass "clean fixture" ':'
case_ fail "#[allow] in library" "printf '#[allow(dead_code)]\nconst X: u32 = 1;\n' >> $L" 'zenith-float-num/src/lib.rs:5 \[allow\]'
case_ fail "#![allow] inner attribute" "printf '#![allow(clippy::needless_range_loop)]\n' > src/x.rs" 'src/x.rs:1 \[allow\]'
case_ fail "#[expect]" "printf '#[expect(dead_code)]\nconst Y: u32 = 2;\n' >> $L" '\[expect\]'
case_ fail "spaced attribute # [ allow (" "printf '# [ allow ( unused ) ]\nconst Z: u32 = 3;\n' >> $L" '\[allow\]'
case_ fail "cfg_attr(.., allow)" "printf '#[cfg_attr(test, allow(unused))]\nconst W: u32 = 4;\n' >> $L" '\[cfg_attr-allow\]'
case_ fail "#[allow] in a test file" "printf '#[allow(dead_code)]\nfn h() {}\n' >> tests/t.rs" 'tests/t.rs:3 \[allow\]'
case_ fail "#[allow] in nostd-tests" "mkdir -p nostd-tests/src && printf '#![forbid(unsafe_code)]\n#[allow(unused)]\nfn g() {}\n' > nostd-tests/src/lib.rs" 'nostd-tests/src/lib.rs:2 \[allow\]'
case_ fail "unsafe block in library" "printf 'pub fn u() { unsafe {} }\n' >> $L" 'zenith-float-num/src/lib.rs:5 \[unsafe\]'
case_ fail "unsafe in a test file" "printf 'fn v() { unsafe {} }\n' >> tests/t.rs" 'tests/t.rs:3 \[unsafe\]'
# r6 item 6: the pending inventory is closed at 0. A listed file fails, and listing it does not exempt its unsafe.
case_ fail "pending inventory entry (closed at 0)" "printf 'fn v() { unsafe {} }\n' > tests/m.rs && printf 'tests/m.rs\t1\tmpfr FFI\n' >> tools/unsafe_pending.tsv" 'tests/m.rs \[pending\] the inventory is closed at 0'
case_ fail "a pending entry does not exempt its unsafe" "printf 'fn v() { unsafe {} }\n' > tests/m.rs && printf 'tests/m.rs\t1\tmpfr FFI\n' >> tools/unsafe_pending.tsv" 'tests/m.rs:1 \[unsafe\]'
case_ fail "pending entry for a library file" "printf 'pub fn u() { unsafe {} }\n' >> $L && printf '$L\t1\tno\n' >> tools/unsafe_pending.tsv" 'closed at 0'
case_ fail "Cargo [lints] allow" "printf '[lints.rust]\ndead_code = \"allow\"\n' >> Cargo.toml" 'Cargo.toml.*\[cargo-lint\] lints.rust.dead_code = allow'
case_ fail "Cargo [workspace.lints] warn (table form)" "printf '[workspace.lints.clippy]\nx = { level = \"warn\", priority = 1 }\n' >> Cargo.toml" 'workspace.lints.clippy.x = warn'
case_ fail "clippy.toml threshold key" "printf 'too-many-arguments-threshold = 9\n' >> clippy.toml" 'clippy.toml:2 \[clippy-config\] too-many-arguments-threshold'
case_ fail "clippy.toml max-fn-params-bools" "printf 'max-fn-params-bools = 8\n' >> clippy.toml" '\[clippy-config\] max-fn-params-bools'
case_ fail "clippy.toml boolean allow-unwrap-in-tests (Spence R2)" "printf 'allow-unwrap-in-tests = true\n' >> clippy.toml" 'clippy.toml:2 \[clippy-config\] allow-unwrap-in-tests'
case_ fail "clippy.toml any non-disallowed key, even at its default" "printf 'allow-expect-in-tests = false\n' >> clippy.toml" '\[clippy-config\] allow-expect-in-tests'
case_ fail "nested .clippy.toml key" "mkdir -p zenith-float-num && printf 'msrv = \"1.80\"\n' > zenith-float-num/.clippy.toml" 'zenith-float-num/.clippy.toml:1 \[clippy-config\] msrv'
case_ pass "clippy.toml disallowed-* lists only" "printf 'disallowed-methods = [\n    { path = \"std::process::exit\", reason = \"no\" },\n]\ndisallowed-types = []\n' >> clippy.toml"
case_ fail "RUSTFLAGS -A in a script" "printf 'RUSTFLAGS=\"-A warnings\" cargo build\n' > scripts/b.sh" 'scripts/b.sh:1 \[lint-flag\]'
case_ fail "cargo clippy -- -A lint" "printf 'cargo clippy -- -A clippy::type_complexity\n' > scripts/c.sh" '\[lint-flag\]'
case_ fail "--cap-lints" "printf 'RUSTFLAGS=--cap-lints=allow cargo test\n' > scripts/d.sh" '\[lint-flag\]'
case_ fail ".cargo/config rustflags -A" "mkdir -p .cargo && printf '[build]\nrustflags = [\"-A\", \"dead_code\"]\n' > .cargo/config.toml" '\.cargo/config.toml.*\[lint-flag\]'
case_ fail "crate root without forbid(unsafe_code)" "printf 'pub use zenith_float_num::f;\n' > src/lib.rs" 'src/lib.rs \[no-forbid\]'
case_ pass "root covered by Cargo [lints.rust] unsafe_code = forbid" "printf '#![deny(unused)]\n/// d\npub fn f() -> u32 { 1 }\n' > $L"
case_ pass "tokens in comments, strings, raw strings, forbid/deny, identifiers" \
  "printf '// #[allow(dead_code)] unsafe\n/* #[allow(x)] unsafe { } */\n//// #[expect(y)] unsafe (four slashes: not a doc comment)\n/// a doc string \"unsafe #[allow(w)]\" is a literal\npub const S: &str = \"#[allow(z)] unsafe\";\npub const R: &str = r#\"#[allow(q)] unsafe\"#;\n#[deny(unused)]\npub const UNSAFE_CODE: u32 = 0;\npub fn is_unsafe(unsafe_code: u32) -> u32 { unsafe_code }\n' >> $L && printf 'git add -A\n' > scripts/g.sh"
# Spence m1 on 0a1cb6c: doc comments are doctest source, so their text is scanned as code (string literals still stripped)
case_ fail "/// #[expect] unsafe in a doc comment" "printf '/// #[expect(y)] unsafe\npub const D: u32 = 0;\n' >> $L" 'zenith-float-num/src/lib.rs:5 \[expect\]'
case_ fail "unsafe in a doctest code block" "printf '/// \x60\x60\x60\n/// let v = unsafe { core::mem::zeroed::<u8>() };\n/// \x60\x60\x60\npub const E: u32 = 0;\n' >> $L" 'zenith-float-num/src/lib.rs:6 \[unsafe\]'
case_ fail "hidden doctest line # #[allow]" "printf '/// \x60\x60\x60\n/// # #[allow(dead_code)]\n/// # fn h() {}\n/// \x60\x60\x60\npub const F: u32 = 0;\n' >> $L" 'zenith-float-num/src/lib.rs:6 \[allow\]'
case_ fail "//! inner doc comment" "printf '//! #[allow(unused)] in crate docs\n' > src/y.rs" 'src/y.rs:1 \[allow\]'
case_ fail "/** block doc comment */" "printf '/**\n * unsafe { }\n */\npub const G: u32 = 0;\n' >> $L" 'zenith-float-num/src/lib.rs:6 \[unsafe\]'
case_ fail "every finding is listed, then a total" "printf '#[allow(dead_code)]\nconst A: u32 = 1;\n#[allow(unused)]\nconst B: u32 = 2;\n' >> $L" 'no_suppress_scan: FAIL 2 finding'
[[ $fail -eq 0 ]] && { echo "no_suppress_selftest: ok ($n cases)"; exit 0; }
echo "no_suppress_selftest: FAIL"; exit 1
