# zenith-float no-clone rules (tools/)

This file documents what `tools/no_clone.sh`, `tools/rules_scan.sh` and `tools/full_gate.sh` enforce. Every baseline and
the exception allowlist only go down; all of them are CODEOWNERS-protected (`@jscarr64`).

## Counters and ratchets

| What | Counted by | Baseline file | Scope |
|---|---|---|---|
| `.clone(` occurrences | `rules_scan.py` `clone_count()` | `tools/clone_count_baseline.txt` | library dirs **and test code** (see below) |
| hidden copy helpers | `no_clone.py --helpers` / `helper_copy_count()` | `tools/helper_copy_count_baseline.txt` | library dirs and test code |
| excepted `.clone()` sites | `no_clone.py` allowlist | `tools/exception_count_baseline.txt` | the five sites in `tools/clone_allowlist.txt` |
| `add_ref(` / `sub_ref(` uses | `rules_scan.py` | `tools/ref_sum_count_baseline.txt` | **library code only** (`REF_SCOPE`) |
| `mul_ref(` / `div_ref(` uses | `rules_scan.py` | `tools/mul_ref_count_baseline.txt` | **library code only** (`REF_SCOPE`) |

### Test code is counted

`#[cfg(test)]` items, `ops/tests.rs`, and the test dirs `zenith-float-num/tests`, `zenith-float-num/benches`, `tests`,
`nostd-tests` (`no_clone.TEST_DIRS`) are scanned like library code by `no_clone.py` and by the clone count. `tools/` and
`target/` are never scanned. The ref ratchets keep their old meaning (library code only), so test uses of
`add_ref` etc. do not count against them; `rules_scan.sh --count-only` prints the with-tests figure as a note.
The diff scan (`rules_scan.sh A..B`) still skips added lines inside `#[cfg(test)]` items for the non-clone rules
(`unwrap`, `panic!`, `f64`, ...); new test clones are caught by the clone-count ratchet.

### Helper-copy counter: how it counts

The clone scan only sees the spelling `.clone()`. These helpers are a `.clone()` plus `set_inexact(false)` or a
precision change and hide the copy: `exact_copy_ext`, `exact_copy`, `exact_operand`, `copy_p`, `svd_copy_prec`,
`lift_num`. The counter counts every identifier occurrence of those names in code (comments and string literals are
ignored) after removing definitions: `fn NAME` and `macro_rules! NAME` are cut out of the line, and `use` lines are
skipped. It also counts `Rc::clone` / `Arc::clone` in `zenith-float-compare/src` (handle shares, which the clone scan
does not match because it looks for `.clone(` and `Clone::clone`). `no_clone.py --helpers --list` prints every site.

## Exception: constant caches

`tools/clone_allowlist.txt` is empty except for one named exception: the constant caches. Each entry is
`path<TAB>fn<TAB>reason<TAB>approval<TAB>doc<TAB>sites`; `sites` is the exact number of `.clone()` sites allowed in that
fn. A listed entry that matches no `.clone()` site is **stale** and fails the scan, sites beyond `sites` fail as
ordinary violations, only `.clone()` can be excepted, and the total number of excepted sites can never exceed
`tools/exception_count_baseline.txt` (it can only be lowered, `rules_scan.sh --update-baseline`).
Approved by Jeff in chat on 2026-10-15 (relayed by Paul: "Please add the other 4 as well").

| File | fn | Why |
|---|---|---|
| `zenith-float-num/src/ops/consts/pi.rs` | `for_prec` | hands out an owned value of pi at `k` bits and keeps the cached series value |
| `zenith-float-num/src/ops/consts/e.rs` | `for_prec` | same, for e |
| `zenith-float-num/src/ops/consts/ln2.rs` | `for_prec` | same, for ln 2 |
| `zenith-float-num/src/ops/consts/ln10.rs` | `for_prec` | same, for ln 10 |
| `zenith-float-num/src/ops/consts/mod.rs` | `CachedFBig::round` | returns the cached value rounded to `p` bits, the cache keeps the wide value |

### Why a copy is unavoidable here

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

### Not covered by this exception

`ExtraCache::for_prec` in `zenith-float-num/src/ops/consts/mod.rs` (two `.clone()` sites, the shared cache for
phi/euler-style constants) and `#[derive(Clone)]` on `CachedFBig` are NOT listed and still count as violations.
