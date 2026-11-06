Line endings (vera, 2026-11-03, vera/zf-td41): a new `.gitattributes` pins
text files to LF on checkout (`* text=auto eol=lf`) and pins golden data
byte-exact (`*.jsonl -text`, `*.tsv -text`, `SHA256SUMS -text`,
`zenith-float-num/tests/oracle/goldens/** -text`, `golds/hex/** -text`), so a
Windows checkout with `core.autocrlf=true` cannot change golden bytes or
leave a `\r` in a parsed field. No result changes: the CSV parsers already
stripped one trailing `\r` per line, so no parser code changed; the new tests
are guards (`zenith-float-num/tests/crlf_inputs.rs`,
`tools/gitattributes_selftest.sh`, and a CRLF case in
`tools/golden_sums_selftest.sh`).
