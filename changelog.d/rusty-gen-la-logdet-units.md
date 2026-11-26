### gen_la_goldens: heavy cosmo logdet cases can run as separate processes (byte-identical)

- `scripts/gen_la_goldens.py --logdet-unit CASE:PART --row-cache DIR` computes one `cosmo_*` logdet
  case (PART `all`, `main`, `enc4096` or `enc16384`) with the serial code and rng stream and writes
  its row text to DIR; a normal run given `--row-cache DIR` uses those rows (each must start with
  the case's own input fields, else exit 7) and computes the rest. Output bytes are unchanged:
  15 of 15 rows from 5 cases match the serial generator's bytes, and `--check` stays OK through the
  cache (`tools/gen_la_rowcache_selftest.sh`, 18 cases). No limit or case changed; the per-case
  budget stays a `--case-budget` parameter.
