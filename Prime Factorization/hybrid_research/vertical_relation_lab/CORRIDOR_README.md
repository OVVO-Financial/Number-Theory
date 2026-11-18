# Newton-tangent / exact-lattice-corridor factorization proof of concept

This is an experimental, **deterministic, balanced-factor-only** N-only region-exclusion method. It is a stand-alone C++17/GMP research prototype, **not** a faster general-purpose factorizer or a GNFS replacement.

## Build and run

```sh
make
./corridor_gmp 978508015703 10000000 corridor
./corridor_gmp 978508015703 10000000 trial
./corridor_gmp 187 100000 count 14 19
python3 verify_geometry.py
python3 test_corridor.py
python3 benchmark.py
```

The benchmark's `trial` baseline is a proper wheel-10 iterator (only 1,3,7,9 mod 10) using GMP exact divisibility; both methods receive the identical `q` domain `[ceil(sqrt(N)), floor(sqrt(2*N))]`. This is the larger-factor range under the assumption `p <= q < 2p`, not an unconditional factor range.

The search reports `factor_found`, `budget_exhausted`, `no_factor_in_balanced_interval`, or `small_prime`. The exhausted interval is NOT a primality certificate for general inputs.

## Exact mathematical certificate

For `N>0` and integer `Q<=q<=Q+H=B`, consider the decreasing convex factor boundary `p(q)=N/q`, with

    p'(q)=-N/q^2,     p''(q)=2N/q^3 > 0.

Its tangent at `Q` is a lower bound and its chord from `Q` to `B` is an upper bound:

    lower(q) = N*(2*Q-q) / Q^2
    upper(q) = N*(Q+B-q) / (Q*B)

(the upper numerator `Q+B-q` is equivalent to `2Q+H-q`). Therefore, if `N=p*q` and the larger integer factor `q` lies in this interval,

    lower(q) <= p <= upper(q).

The difference `upper(q)-lower(q)` is at most `N*H^2 / Q^3`. Thus a short interval produces a narrow rational linear corridor guaranteed to contain any integer factor pair in the region.

Because a semiprime coprime to 10 has both factors in `1,3,7,9 mod 10`, for every admissible `q mod 10`, `p mod 10` is determined uniquely by `N * inv(q) mod 10`.

For one lane `(rp,rq)`, with `q` running through that residue class, the exact number of integer points in the corridor is

    sum_{q == rq (mod 10), Q<=q<=B}
      [floor((upper(q)-rp)/10) - ceil((lower(q)-rp)/10) + 1].

Both summands are floors of affine rational functions of the progression index. The Euclidean floor-sum algorithm evaluates their sum in O(log magnitude) divisions, **without visiting individual q positions**. This is a count of possible lattice points inside the envelope, not a count of factor pairs. It is an N-only certificate that can exclude an entire region when the count is zero.

The recursive search bisects every nonempty region and recomputes a tighter exact corridor. At a singleton integer q, the corridor is exactly p=N/q, so a nonempty lane signifies an actual factor. Returned factors are independently checked by GMP divisibility.

The method is named for Newton linearization, but it is **not** ordinary Newton-Raphson iterated on the integer sawtooth residual `E(q)=q*floor(N/q)-N`. Such iterations do not locate hidden integer factors; the new ingredient is exact convexity bounds and arithmetic lattice counting.

## Limits

- Completeness holds for a factor pair `p<=q<2p` when the nontrivial factors avoid 2 and 5, with explicit small-factor prechecks. A wider unbalanced factor may be reported as unresolved.
- Every floor sum is an exact arbitrary-precision calculation. The cost of evaluating thousands of region counts is substantial.
- No subexponential bound or GNFS advantage is claimed. The number of regions remains large for RSA-size inputs.
- The broader `|pq-N|<=B_residual` problem can use tangents to `(N-B_residual)/q` and chords to `(N+B_residual)/q`. The modulo-10 factor constraints apply only to zero residual and must NOT be reused unchanged for general nonzero residual relations.

## Reproduced results

For N=978508015703 (factor pair 752867 * 1299709), the balanced larger-factor interval has 409738 q positions. Partitioning it into tiles with `H=100` gives 4057 tiles, 4005 of which have exactly zero admissible integer lattice points inside their rational envelopes. The remaining 52 tiles jointly contain 52 potential integer `(p,q)` pairs; an independent brute-force check found exactly one actual factor pair. The recursive C++ program visited 2928 regions before finding the factor.

For other representative inputs, seven-run median timings measured in this environment are recorded in `benchmarks.json`:

| input class | corridor | wheel-10 trial |
| --- | ---: | ---: |
| 40-bit close | 0.218 ms | **0.024 ms** |
| 40-bit balanced | 35.258 ms | **3.928 ms** |
| 46-bit balanced | 94.166 ms | **19.834 ms** |
| 54-bit balanced | 250.402 ms | **120.011 ms** |

The optimized **ordinary wheel-10 baseline is faster on every one of these tested cases**. A separate single 61-bit example took 1.27 s for corridor and 1.19 s for wheel trial. The method reduces candidate visitation, but GMP floor-sum overhead overwhelms that reduction at these sizes. A native-u128 implementation and a better region scheduler are potential optimizations, **not measured improvements**.

## Independent checks

- `verify_geometry.py`: 1983 independent random rational-envelope checks plus 7932 independent lane checks; several hundred checks also cross-compare Python and C++ counters.
- `test_corridor.py`: 12892 C++ factorization control cases, plus an explicit unbalanced case that remains unresolved.
- `benchmarks.json`: all seven timed repetitions for the four matched method pairs.
