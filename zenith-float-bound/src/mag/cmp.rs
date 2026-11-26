//! Total order on `Mag`: zero least, `Infinity` greatest, finite by the canonical key.

use super::Inner;
use super::Mag;
use core::cmp::Ordering;

impl PartialEq for Mag {
    fn eq(&self, other: &Mag) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for Mag {}

impl PartialOrd for Mag {
    fn partial_cmp(&self, other: &Mag) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Mag {
    /// Zero is least and equals only zero. `Infinity` is greatest and equals only
    /// `Infinity`. A canonical finite value orders by `(exp, mant)` (CONTRACT 14.5).
    fn cmp(&self, other: &Mag) -> Ordering {
        match (&self.inner, &other.inner) {
            (Inner::Infinity, Inner::Infinity) => Ordering::Equal,
            (Inner::Infinity, _) => Ordering::Greater,
            (_, Inner::Infinity) => Ordering::Less,
            (Inner::Finite { exp: ea, mant: ma }, Inner::Finite { exp: eb, mant: mb }) => {
                match (*ma == 0, *mb == 0) {
                    (true, true) => Ordering::Equal,
                    (true, false) => Ordering::Less,
                    (false, true) => Ordering::Greater,
                    (false, false) => (ea, ma).cmp(&(eb, mb)),
                }
            }
        }
    }
}
