# NA_PROOF: ComplexBall::rad tie row

Claim: `ComplexBall::rad` (zenith-float-num/src/ball.rs:270) takes no RoundingMode and performs no rounding (Spence G1-G3: accessor, construct-copy or exact integer/modular), so an exact half-ulp tie at p cannot arise.

Evidence: TO BE WRITTEN by the proving agent (cite the signature and body at zenith-float-num/src/ball.rs:270).

Second-agent sign-off: UNSIGNED

<!-- pm-coverage FAILs this row until the line above reads `Second-agent sign-off: <agent> <date> <commit>` -->
