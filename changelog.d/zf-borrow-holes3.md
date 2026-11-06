- rules_scan BORROW: data-flow detection of operand reads flowing into owned sinks (push/insert/Box::new/map..collect,
  index and deref assignment), bans `.to_string()` and `format!`, scans zenith-float-num/build.rs, exceptions keyed
  by file + impl + fn + line anchor with exact allowances ('lower allowance ... from N to M' on over-allowance),
  std-receiver method resolution (drops the ieee_soft get/map fallback exceptions); documented lexical known limit.
- div_basic: with d == 1 the divisor is borrowed instead of copied into buf2 (values unchanged).
