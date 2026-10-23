# zenith-float no-clone rules (tools/)

This file documents what `tools/no_clone.sh`, `tools/rules_scan.sh` and `tools/full_gate.sh` enforce. Every baseline and
the exception allowlist only go down; all of them are CODEOWNERS-protected (`@jscarr64`).

## Counters and ratchets

| What | Counted by | Baseline file | Scope |
|---|---|---|---|
| `.clone(` occurrences | `rules_scan.py` `clone_count()` | `tools/clone_count_baseline.txt` | library dirs **and test code** (see below) |
| hidden copy helpers | `no_clone.py --helpers` / `helper_copy_count()` | `tools/helper_copy_count_baseline.txt` | library dirs and test code |
| exceptions: excepted `.clone()` sites (allowlist, empty) + the Ziv working-copy helper (1) | `no_clone.py` allowlist, `rules_scan.py` `ziv_helper_check()` | `tools/exception_count_baseline.txt` | none (baseline 1: the one helper) |
| `unwrap`/`expect`/`panic!`/`unreachable!`/`todo!`/`unimplemented!` in **test code** | `rules_scan.py` `test_panic_sites()` (`--list-test-panics`) | `tools/test_panic_count_baseline.txt` | `#[cfg(test)]` items, test files and dirs, the compare crate |
| `add_ref(` / `sub_ref(` uses | `rules_scan.py` | `tools/ref_sum_count_baseline.txt` | **library code only** (`REF_SCOPE`) |
| `mul_ref(` / `div_ref(` uses | `rules_scan.py` | `tools/mul_ref_count_baseline.txt` | **library code only** (`REF_SCOPE`) |

### Test code is counted

`#[cfg(test)]` items, `ops/tests.rs`, and the test dirs `zenith-float-num/tests`, `zenith-float-num/benches`, `tests`,
`nostd-tests` (`no_clone.TEST_DIRS`) are scanned like library code by `no_clone.py` and by the clone count. `tools/` and
`target/` are never scanned. The ref ratchets keep their old meaning (library code only), so test uses of
`add_ref` etc. do not count against them; `rules_scan.sh --count-only` prints the with-tests figure as a note.
The diff scan (`rules_scan.sh A..B`) applies the clone rules (`.clone(`, `.to_vec(`, `.cloned(`) to added test lines too
(`#[cfg(test)]` items, `tests.rs`, `tests/`, `benches/`, `examples/` files under the library dirs). Test code stays exempt
from the other rules (`unsafe`, `allow`, `unwrap`, `expect`, `panic!`, `f32`/`f64`), whose meaning is unchanged.

### Test code may not crash either (Jeff)

`unwrap(`, `unwrap_err(`, `expect(`, `expect_err(`, `panic!`, `unreachable!`, `todo!` and `unimplemented!` are banned in test
code as in library code (`rules_scan.PANIC_RULES`). `assert!`, `assert_eq!`, `assert_ne!` and `debug_assert!` are fine: a
test returns `Result<(), E>` and uses `?`, or asserts. The ratchet counts every such site in test code (the compare crate
is counted whole) and may only go down (`tools/test_panic_count_baseline.txt`, CODEOWNERS `@jscarr64`).
`rules_scan.sh --list-test-panics` prints each site. The per-line diff scan flags added test lines in `--staged`,
`--working` and `--diff-file` modes (and in test-dir files outside the library dirs); a RANGE scan leaves test code to
the ratchet, because the legacy sites of old commits would fail every range. An unwrap that is truly necessary must be
reported to Jeff with file:line and the reason a `Result` return cannot work; there is no allowlist for it.

### Helper-copy counter: how it counts

The clone scan only sees the spelling `.clone()`. These helpers are a `.clone()` plus `set_inexact(false)` or a
precision change and hide the copy: `exact_copy_ext`, `exact_copy`, `exact_operand`, `copy_p`, `svd_copy_prec`,
`lift_num`. The counter counts every identifier occurrence of those names in code (comments and string literals are
ignored) after removing definitions: `fn NAME` and `macro_rules! NAME` are cut out of the line, and `use` lines are
skipped. It also counts `Rc::clone` / `Arc::clone` in `zenith-float-compare/src` (handle shares, which the clone scan
does not match because it looks for `.clone(` and `Clone::clone`). `no_clone.py --helpers --list` prints every site.

## Exception: Ziv working copy

**Approved** by Jeff (relayed by Paul), under four conditions, all enforced by `rules_scan.sh`:

1. **One helper, one file, CODEOWNERS.** `ExactNumNumber::ziv_pass_operand`, defined in
   `zenith-float-num/src/num/ziv_operand.rs` and nowhere else; the file has a `.github/CODEOWNERS` entry (`@jscarr64`).
2. **A non-clone-shaped name.** `ziv_pass_operand`, not `clone`/`copy`; the helper does not call `.clone()` or `Mantissa::clone`.
3. **One allocation plus one limb copy per retry pass, at that pass's precision.** A Ziv loop reads the same exact input
   at a growing precision on every pass and its kernel consumes the argument, so each pass needs a fresh owned value
   (shape A of `ziv-loop-shapes.md`). The helper allocates the new mantissa once and copies the input's words into it
   once; it never rounds (precision may only grow), and it carries the value, sign, exponent and both inexact flags over.
4. **The exception baseline only goes down, counted by the existing counters.** The definition counts 1 in
   `tools/exception_count_baseline.txt` (0 -> 1, the only raise, approved). Every *use* of `ziv_pass_operand` counts in the
   helper-copy counter (`no_clone.HELPERS`; uses inside its own file, its unit tests, are not counted).
   **The scan rejects a second such helper**: any other `fn` whose name has `ziv` together with pass/operand/copy/clone/
   dup/widen/working/extend, a second definition of `ziv_pass_operand`, the definition in another file, or a missing
   CODEOWNERS entry fails, whatever the baselines say.

The module is compiled in every build since the ln/log2 conversion; its uses are limited to `tools/ziv_shape_a_sites.txt`.
A loop converted to the helper replaces a `.clone()?` (clone count down) with a use (helper-copy count up by one);
each such conversion is reported to Jeff before it is committed.

## Exception: constant caches

**Withdrawn.** This exception was approved by Jeff in chat on 2026-10-15 (relayed by Paul: "Please add the other 4 as
well") and then withdrawn the same day by his ruling R3: borrowed constant arguments (an ordinary `&ExactNum` fn
argument, not a ref form) and retained integer series state replace it, and the constant-cache exception is dropped.
`tools/clone_allowlist.txt` is empty again and `tools/exception_count_baseline.txt` is 0. The five former sites
(`pi.rs`, `e.rs`, `ln2.rs`, `ln10.rs` `for_prec`, and `CachedFBig::round` in `consts/mod.rs`) are real `.clone()` calls
and count as ordinary violations in `tools/clone_count_baseline.txt` until they are converted.

This section is kept (the scan requires a doc with this heading for any allowlist entry, and the selftests exercise the
mechanism) as the worked example of how an exception is documented: a named exception needs a concrete failing example
with real compiler output, Jeff's approval, and an entry
`path<TAB>fn<TAB>reason<TAB>approval<TAB>doc<TAB>sites`; `sites` is the exact number of `.clone()` sites allowed in that
fn. A listed entry that matches no `.clone()` site is **stale** and fails the scan, sites beyond `sites` fail as
ordinary violations, only `.clone()` can be excepted, and the total number of excepted sites can never exceed
`tools/exception_count_baseline.txt` (it can only be lowered, `rules_scan.sh --update-baseline`).

The (withdrawn) entries were:

| File | fn | Why it was proposed |
|---|---|---|
| `zenith-float-num/src/ops/consts/pi.rs` | `for_prec` | hands out an owned value of pi at `k` bits and keeps the cached series value |
| `zenith-float-num/src/ops/consts/e.rs` | `for_prec` | same, for e |
| `zenith-float-num/src/ops/consts/ln2.rs` | `for_prec` | same, for ln 2 |
| `zenith-float-num/src/ops/consts/ln10.rs` | `for_prec` | same, for ln 10 |
| `zenith-float-num/src/ops/consts/mod.rs` | `CachedFBig::round` | returns the cached value rounded to `p` bits, the cache keeps the wide value |

### R3: what replaces it

- A constant needed as an argument is passed borrowed: an ordinary `&ExactNum` fn argument (not an `add_ref`-style ref
  form), so no copy is made to hand it to a kernel.
- Integer series state (the retained `p`/`q`/`t` terms of a constant's series) may be kept between calls.
- A getter that must return a value rounded to the caller's precision while keeping a wider value is rebuilt from the
  retained integer series state, not copied from a cached `ExactNum`.

### Why the copy looked unavoidable (historical, withdrawn)

The getter must return a value rounded to the caller's precision while the cache keeps the wider value for later
callers, so there are two live owners by construction. Getter (pi; e, ln 2 and ln 10 have the same shape):

```rust
    /// Return value of PI with precision `k`.
    pub(crate) fn for_prec(&mut self, k: usize, rm: RoundingMode) -> Result<ExactNumNumber, Error> {
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = round_p(k)
            .checked_add(p_inc)
            .ok_or(Error::InvalidArgument)?;

        loop {
            let kext = (k + 46 + WORD_BIT_SIZE) / 47;

            if self.b > kext as u64 {
                let mut ret = self.val.clone()?;

                if ret.try_set_precision(k, rm, p_wrk)? {
                    return Ok(ret);
                }

                bump_prec_retry(&mut p_wrk, &mut p_inc, k)?;
            }

            let mut pk;
            let mut qk;
```

Public wrapper, which the example below calls (`Consts::pi` takes `&mut self` and returns an owned `ExactNum`):

```rust
    pub fn pi(&mut self, p: usize, rm: RoundingMode) -> ExactNum {
        match self.pi_num(p, rm) {
            Ok(v) => v.into(),
            Err(e) => ExactNum::nan(Some(e)),
        }
    }
```

`CachedFBig::round`:

```rust
    /// Round the cached value to `p` bits.
    pub fn round(&self, p: usize, rm: RoundingMode) -> ExactNum {
        let mut v = self.inner.clone();
        let _ = v.set_precision(p, rm);
        v
    }
}
```

### Example: two owned values from one cache

```rust
let a = cc.pi(128, rm);
let b = cc.pi(128, rm);
// a and b are two owned ExactNum values, and `cc` still holds its own cached pi.
```

If the getter returned `&ExactNum` (a reference into the cache) instead of an owned value, this does not compile.
Reproduced with a throwaway crate (`/tmp/zf-borrowck-demo`, cargo 1.99.1, `Consts::pi_ref(&mut self, p) -> &ExactNum`
delegating to `PiCache::for_prec_ref(&mut self, k) -> &ExactNum`), real compiler output:

```text
error[E0499]: cannot borrow `cc` as mutable more than once at a time
  --> src/bin/bad_ref.rs:17:13
   |
16 |     let a = cc.pi_ref(128);
   |             -- first mutable borrow occurs here
17 |     let b = cc.pi_ref(128);
   |             ^^ second mutable borrow occurs here
18 |     println!("{:?} {:?}", a, b);
   |                           - first borrow later used here

For more information about this error, try `rustc --explain E0499`.
error: could not compile `borrowck-demo` (bin "bad_ref") due to 1 previous error
```

(Growing the cache needs `&mut self`, so every getter call borrows `cc` mutably for as long as the returned reference lives.) For `CachedFBig::round(&self, p, rm)` a reference is
impossible for a different reason: the rounded value does not exist in the cache, so it has to be a new owned value:

```text
error[E0515]: cannot return reference to local variable `v`
 --> src/bin/bad_round.rs:9:9
  |
9 |         &v
  |         ^^ returns a reference to data owned by the current function

For more information about this error, try `rustc --explain E0515`.
error: could not compile `borrowck-demo` (bin "bad_round") due to 1 previous error
```

The owned form compiles and runs in the same crate (`cargo run --bin borrowck-demo` prints two values).

### Not covered by the withdrawn exception

`ExtraCache::for_prec` in `zenith-float-num/src/ops/consts/mod.rs` (two `.clone()` sites, the shared cache for
phi/euler-style constants) and `#[derive(Clone)]` on `CachedFBig` are NOT listed and still count as violations.

## Helper-copy ruling (Jeff, 2026-10-16 8:25 AM CT, relayed by Paul)

"helper_copy may rise only in a commit that lowers clones by at least as many; combined clone+helper_copy ratchet
that only goes down; a pinned, CODEOWNERS-protected (@jscarr64) list of allowed shape-A Ziv sites; loop conversions
may resume after Stage A is gate-green."

Site-list sign-off (Jeff, 2026-10-23 7:29 AM CT, @jscarr64, relayed by Paul: "Go Ahead with the list sign off"):
`tools/ziv_shape_a_sites.txt` is signed off as exactly six sites: ln (log.rs:133), log2 (log.rs:266), sinh
(sinh.rs:49), sin_cos (sin_cos.rs:66), sin (sin.rs:68), cos (cos.rs:71). A seventh site needs a new sign-off.

How `tools/rules_scan.py` enforces it:

| rule | check | file |
|---|---|---|
| combined count only goes down | `ruling_check()` ratchet on clone + helper_copy | `tools/clone_helper_combined_baseline.txt` (806) |
| helper_copy rises only with an equal-or-larger clone drop, per commit | `per_commit_check()` in range mode, commits after the rule base | `tools/helper_copy_rule_base.txt` |
| `ziv_pass_operand` used only at pinned sites | `ruling_check()`, by file and enclosing `fn` | `tools/ziv_shape_a_sites.txt` |
| the pinned list is owned by @jscarr64 | `codeowners_cover()` | `.github/CODEOWNERS` |
