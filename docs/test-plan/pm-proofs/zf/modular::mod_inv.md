# NA_PROOF: modular::mod_inv tie row

Claim: `modular::mod_inv` (zenith-float-num/src/modular.rs:89) takes no RoundingMode and performs no rounding (Spence G1-G3: accessor, construct-copy or exact integer/modular), so an exact half-ulp tie at p cannot arise.

Evidence: TO BE WRITTEN by the proving agent (cite the signature and body at zenith-float-num/src/modular.rs:89).

Second-agent sign-off: UNSIGNED

<!-- pm-coverage FAILs this row until the line above reads `Second-agent sign-off: <agent> <date> <commit>` -->
