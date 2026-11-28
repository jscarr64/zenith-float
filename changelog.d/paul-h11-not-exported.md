### Changed (bound layer, HOLD H11)

- `CBall::mul` and `CBall::div` are no longer in the public API. They were refusing stubs (`Resource`, "not
  implemented: waits zf borrowing ops") and are held with their value rows (HOLD H11) until they are implemented on
  the borrowing once-rounded zf ops. The two refusal tests went with them; `gold_eval` reports `mul`/`div` golden rows
  as held. Held slices H01-H08 were already not compiled and are not exported.
