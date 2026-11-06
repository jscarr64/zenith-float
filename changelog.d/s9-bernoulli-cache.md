### Changed
- digamma: the O(K^3) Bernoulli tangent recurrence is replaced by a shared, thread-safe cache of exact rational B_2k (std RwLock, borrowed slices, no Clone/unsafe), built by an O(K^2) boustrophedon recurrence and grown on demand, with per-precision rounded slots. Verified exactly against sympy bernoulli for k=1,2,3,6,10,30,64,100.
