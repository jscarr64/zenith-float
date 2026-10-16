#!/usr/bin/env bash
# Self-test for tools/no_clone.py: each banned form is caught (test code included), Copy derives are not, the allowlist works.
set -u
here="$(cd "$(dirname "$0")" && pwd)"; t="$(mktemp -d)"; trap 'rm -rf "$t"' EXIT
mkdir -p "$t/src" "$t/tools"
run() { NO_CLONE_ROOT="$t" python3 "$here/no_clone.py" "$@"; }
fail() { echo "SELFTEST FAIL: $*"; exit 1; }
printf 'fn clean() { let a = 1; }\n#[derive(Clone, Copy)]\nstruct S;\n' > "$t/src/a.rs"
run >/dev/null || fail "clean tree reported violations"
# test code is counted: cfg(test) items, cfg(all(test, ..)) items, ops/tests.rs and the test dirs
printf 'fn clean() {}\n#[cfg(test)]\nmod tests { fn t() { let x = v.clone(); let y = s.to_vec(); } }\n' > "$t/src/a.rs"; run >/dev/null && fail "cfg(test) clones not counted"
printf 'fn a() {}\n#[cfg(all(test, feature = "std"))]\n#[test]\nfn t() { let y = x.clone(); }\n' > "$t/src/a.rs"; run >/dev/null && fail "cfg(all(test, ..)) clone not counted"
printf 'fn a() {}\n' > "$t/src/a.rs"; mkdir -p "$t/src/ops" "$t/zenith-float-num/tests"
printf 'fn t() { let y = x.clone(); }\n' > "$t/src/ops/tests.rs"; run >/dev/null && fail "ops/tests.rs clone not counted"; rm "$t/src/ops/tests.rs"
printf 'fn t() { let y = x.clone(); }\n' > "$t/zenith-float-num/tests/it.rs"; run >/dev/null && fail "integration-test clone not counted"; rm "$t/zenith-float-num/tests/it.rs"
printf 'fn clean() {}\n' > "$t/src/a.rs"
for body in 'let y = x.clone();' 'let y = s.to_vec();' 'let y = it.cloned();' 'let y = Clone::clone(&x);' 'y.clone_from(&x);'; do
  printf 'fn f() { %s }\n' "$body" > "$t/src/a.rs"; run >/dev/null && fail "missed: $body"
done
printf '#[derive(Clone)]\nstruct Big;\n' > "$t/src/a.rs"; run >/dev/null && fail "missed derive(Clone)"
printf 'impl Clone for Big { fn clone(&self) -> Self { Big } }\n' > "$t/src/a.rs"; run >/dev/null && fail "missed impl Clone"
printf 'fn f() { let y = x.clone(); }\n' > "$t/src/a.rs"
printf '## Exception: constant caches\n' > "$t/tools/NO_CLONE_RULES.md"
printf 'src/a.rs\tf\tproof\tapproved\ttools/NO_CLONE_RULES.md\t1\n' > "$t/tools/clone_allowlist.txt"; run >/dev/null || fail "allowlist ignored"
printf 'src/a.rs\tf\tproof\tapproved\n' > "$t/tools/clone_allowlist.txt"; run >/dev/null && fail "4-field (undocumented) allowlist entry accepted"
printf 'src/a.rs\tg\tproof\tapproved\ttools/NO_CLONE_RULES.md\t1\n' > "$t/tools/clone_allowlist.txt"; run >/dev/null && fail "stale / wrong-fn allowlist entry accepted"
echo "no_clone selftest: ok"
