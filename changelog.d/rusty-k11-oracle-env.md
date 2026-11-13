### Oracle scripts: no machine paths; the env is checked at runtime and fails closed (K11)

- All 15 `zenith-float-bound/scripts` oracle scripts use `#!/usr/bin/env python3` (they named one
  machine's interpreters, 4 of them the retired `~/.local/oracle-venv`).
- `gen_la_goldens.py` (and `gen_la_lstsq_goldens.py`, which imports it) now runs
  `goldlib.require_oracle` before any oracle import, so it refuses anything but a verified team env
  with the four pins; before, it checked versions only and ran in `~/.local/oracle-venv`.
- `gen_la_cond_selftest.sh` defaults to the team env, not `~/.local/oracle-venv`.
- `tools/mutate.sh` refuses a MUTATE_TMP / MUTATE_TARGET / MUTATE_PERSIST inside any git checkout
  (was: a hard-coded home path or this repo).
- New gate step `oracle-env` (`tools/oracle_env_selftest.sh`): no machine path in tracked scripts,
  every generator refuses a forged wrong-pin env and a pinned non-team env, control passes.
