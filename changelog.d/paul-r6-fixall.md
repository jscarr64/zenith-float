### Changed (PR #40 r6 fix-all, paul, 2026-11-24; Jeff 4:40/4:42 PM CT, Paul 9:41/9:53/10:00 PM CT)

Covers this round's commits that carried no changelog.d entry of their own (P17e; no history rewrite):
4bfa874, b5393f5 (item 3 clones 1/n, 2/n), ef223c6 (exp bit-burst inexact flag), 22d1a11 (rem_pi |x| < 2π),
dc9fcfb (zero-word tests count words), 12bacc4 (heavy tests check HEAVY.SHA256), 38f2234 (pubfn_coverage),
39dc8a9 (bound_layer_scan --relative), 6bd7acd (unlisted-clone check), 194e8e4 (item 3 batch 3) and the batch 4 commit.

- zenith-float-bound: pubfn_coverage covers CBall::mul / CBall::div and harper_s2_far_halves.rs has a direct
  cmp_exact_diff test; bound_layer_scan diffs with `git diff --relative`, which exposed (and fixed) asserts in the
  cball.rs R1 test.
- tools/rules_scan.py: every `.clone(` line needs its own tools/rules_allowlist.txt entry (`unlisted-clone` check,
  `--list-unlisted-clones`, tools/unlisted_clone_baseline.txt ratchet, now 0). Entries that share a line key each
  cover the line nearest the path:LINE they name.
- Clones (P16): 464 -> 291 over the round; from 407 at 6bd7acd, 116 removed by borrowing, moving or restructuring
  (Option seeds read through the borrow, rows that own the state, mem::replace / Option::replace, into_parts,
  add_ref, rounded_to where its error behaviour matches) and 142 kept, each with a one-line reason
  ([Jeff 4:42 PM CT] entries). Every change bit-identical: digest identical to base, MPFR check mode green.
- tests/mpfr: mp.rs drops derive(Clone) on Big / Float / Integer; assert_float_close_ref borrows the MPFR value.
