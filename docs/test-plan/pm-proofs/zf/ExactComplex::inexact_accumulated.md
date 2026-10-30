# NA_PROOF: ExactComplex::inexact_accumulated tie row

Claim: `ExactComplex::inexact_accumulated` (zenith-float-num/src/complex.rs:127) takes no RoundingMode and performs no rounding (Spence G1-G3: accessor, construct-copy or exact integer/modular), so an exact half-ulp tie at p cannot arise.

Evidence: TO BE WRITTEN by the proving agent (cite the signature and body at zenith-float-num/src/complex.rs:127).

Second-agent sign-off: UNSIGNED

<!-- pm-coverage FAILs this row until the line above reads `Second-agent sign-off: <agent> <date> <commit>` -->
