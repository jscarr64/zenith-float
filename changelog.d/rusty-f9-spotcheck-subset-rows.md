### Mode-subset loops: the named modes are asserted (Spence spot-check at cc85211)

- `round_reports_its_own_rounding`: Up and FromZero now assert the inexact 1 that the label claimed for a positive
  subnormal rounded to an integer. ToOdd moves to the zero group, because it is nearest with ties to odd; the old
  label wrongly said it rounds up.
- `wide_primitives` 1/-2^k: the label now names None as NONE-GAP F-B. None gives -Inf today where the ruling is -MAX,
  and that row is queued for ZF-M27.
