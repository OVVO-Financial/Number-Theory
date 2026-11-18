# Dyadic / Gaussian / vertical-coordinate synthesis (CPU/GMP research)

This research component extends the prior `Number-Theory` hybrid package with **an N-only universal dyadic lane cover**, a proved lower bound on every residue-only three-witness cover, a modular-Newton/squaring conjugacy, an off-strip Gaussian norm-one `p±1` conditional factor path using the `2ab` numerator, and a separate norm-residual GCD channel. The complete proofs are in `DYADIC_VERTICAL_SYNTHESIS.md`.

**These are NOT new general-purpose factoring exponent bounds.** In particular, the Gaussian group method is a standard form of group-order factoring, and the normalizing transform of the precise vertical `Re(z²)=N` line is constant modulo N. This package records positive and negative benchmarks separately.

## Build

Requires C++17 and libgmpxx/libgmp. Python 3.10+ is used for tests and offline number generation, with SymPy for the corpus/census and SciPy for the **optional finite ILP** optimum search. No CUDA device is required.

```sh
make all
./hybrid_native 978508015703 hybrid 256 400000
./hybrid_native 4295491591 hybrid 256 400000
./dual_residual 243866858267 12000 256 both
./newton2_primary 4295491591 48 4
./gaussian_norm_one 978508015703 2048 4
```

The CLI inputs are `N` only, plus work budgets. None accepts hidden `p` or `q`. The `hybrid_native` command does **not** claim primality or factorization when it exhausts a budget: `factor=0`, `method=unresolved`.

## Independent tests

```sh
make test
# optional exact ILP: python opt_lanes.py
# timed one-shot benchmarks: python benchmark_hybrid_native.py
```

`lane_cover.py` constructs an exact factor-independent square-witness covering set. `opt_lanes.py` uses SciPy HiGHS to minimize the number of selected three-witness multiplier classes for finite `k`; this offline optimization is not a factoring algorithm. `test_lane_selector.py` checks *every* odd `N mod 2^k` and odd factor residue at k=3..10.

`test_newton_conjugacy.py`, `test_vertical_group_bridge.py`, `test_gaussian.py`, and `test_norm_channel.py` independently verify the modular identities, prime-field group-order bound, Gaussian norm-one factorizations, extra norm GCD channel, and all returned factors. The unmodified corpus `smooth_order_benchmark.json` includes known `p,q` only for **oracle checking**; C++ workers never read those factors.

## What the programs do

- `hybrid_native`: In-process staged kernel: power-of-two Newton GCD, Gaussian norm-one p±1 stage, Pollard p−1, Lucas p+1, then bounded Brent rho. This is a research composition of existing group-order and rho algorithms. All factors are verified.
- `gaussian_norm_one`: Norm-one Gaussian arithmetic on `(a²−b²+i2ab)/(a²+b²)`; optional `near` seeds let us test whether near-strip geometry improves group-order success.
- `newton2_primary`: Exact `t→t²` modular conjugacy with GCD extraction at t=±1. Useful for a factor whose order has a sufficiently small 2-primary part and no odd component for that seed.
- `smooth_order`: Independent standard p−1/Lucas p+1 and Brent rho controls.
- `dual_residual`: Batched product of the CTM error `E`, the complex norm `H`, or their product; corrects for exact E=0 factor hits.

## Raw data and caveats

- `hybrid_native_benchmark.json` - matched in-process results, 120 semiprimes, including 80 planted smooth-order cases.
- `geometry_vs_fixed.json` - near-strip seeds versus fixed Gaussian seeds.
- `lane_census.json`, `optimum_census.json`, `sample_optimum.json` - N-only lane-cover counts and finite ILP results.
- `dual_residual_census.json` - corrected exact-Fermat-hit accounting.
- `norm_channel_checks.json`, `rsa_bounded_checks.json` - independent geometry data and unsuccessful RSA-scale bounded checks.

All timing numbers are short, environment-specific runs, not hardware scaling predictions. Side-by-side comparisons charge every stage. The stronger Gaussian norm GCD channel helps some inputs but is not a general global factor locator. Production GNFS is not benchmarked here.
