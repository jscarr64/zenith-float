### Fixed (gate tooling)

- The full gate no longer leaves Python bytecode caches in the tree. The official gate4d on 5227f77d passed every
  step but its footer said `RESULT: NOT A GATE (tree dirty at end)` because of
  `zenith-float-bound/scripts/__pycache__/{goldlib,gen_la_goldens}.cpython-312.pyc` and
  `tools/__pycache__/no_clone.cpython-312.pyc`. `tools/oracle_env_selftest.sh` `run_forged` clears the environment with
  `env -i`, so it now sets `PYTHONDONTWRITEBYTECODE=1` itself; every other python call reachable from
  `tools/full_gate.sh` (rules_scan.sh, no_clone.sh, rules_scan_selftest.sh, no_clone_selftest.sh, fresh_test_run.sh,
  test_rounding_mode_matrices.sh, scripts/zenith_prepublish.sh, check_leaves.sh, check_precision_docs.sh) runs
  `python3 -B` (including the per-commit `--print-counts` child that `tools/rules_scan.py` starts with
  `sys.executable`, which re-imported `tools/no_clone.py`), and `tools/full_gate.sh` exports `PYTHONDONTWRITEBYTECODE=1` as a backstop. No check changed what it tests.
