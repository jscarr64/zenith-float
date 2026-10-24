### Fixed
- `ln_gamma` is ln|Γ|: negative non-integer arguments work (reflection on the fractional part), poles at 0 and the negative integers return +Inf exactly, and overflow is directional (+MAX under Down/ToZero).
