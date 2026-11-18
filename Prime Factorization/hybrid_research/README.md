# Complex-space factorization: unified code and certified research program

This directory contains the factorization algorithms and benchmark results developed in the complex-space / naive-grouping / Fermat / iterated-average research program. They are separated into independently testable engines, **not represented as one benchmark-certified faster-than-GNFS algorithm**.

## Components

- `vf-factor-lab/` — unified VF / trial-product / Fermat wheel-20, shared bounded progress and batched product-tree polynomial GCD; C++17/GMP single- and multiworker executables; Python reference and correctness tests. `Makefile.gmp` and `Makefile.parallel` build independent C++ tools.
- `naive_grouping_gmp/` — corrected deterministic grouping seed worker with exact doubling recurrence, optional GCD extension, conditional balanced-band mode, tests and the exact successful-group / odd-chain counting theorem.
- `vertical_relation_lab/` — GMP squared-green-line Fermat sieve, rational Newton convex lattice envelope, exact floor-sum regional certificates, **regionwise prime-power valuation mass** and smooth-relation exclusion without scanning individual lattice points, original two-sided CTM/grouping relation experiment, and new geometry-selected **one-sided** quadratic-sieve relation modes. Includes raw benchmark evidence and independent mathematical tests.
- `julia/complex_space_factorization.cpp` (existing repository file, updated by this PR) — earlier unified CTM, VF, wheel, IA and Lehman engine. A factor-5 precheck is now independent of its optional stage-0 trial division; no semiprime ending in 5 reaches the active Fermat digit wheel.
- `julia/complex_space_factorization_cuda.cu` (existing repository file) — experimental CUDA route. GPU device performance is **not certified by the included CPU benchmarks**.

## Mathematics and performance

Read `vertical_relation_lab/math/FORMALIZATION.md` and `MONOTONE_ONE_SIDED.md` first. The exact squared-plane identity for `(a,b)` reads

    E=a²-b²-N, Y=2ab, H=a²+b²,
    H²-N²-Y²=E(2N+E).

The *real* green-line hypotenuse `sqrt(N²+Y²)` is increasing and convex, but nearest-integer residuals are not monotone. Rational Newton bounds and modular Fermat restrictions certify empty local lattice regions **without enumerating lattice points**. That does not imply cheap global localization.

A one-sided square-congruence is already available:

    x² == x²-N (mod N).

CTM staircases, grouping bands or the vertical `Y` coordinate can choose `x`; only `x²-N` must be smooth. This is exactly the quadratic-sieve relation, so any algorithmic gain must be demonstrated by selecting better `x` candidates, faster than an equally budgeted conventional QS baseline. The included benchmark is a **small-input relation-generation experiment**, not a production GNFS comparison. It currently shows no consistent overall advantage for the geometry-selected streams.

**Prime-power region theorem.** The new `src/prime_power_region.py` and `math/PRIME_POWER_MASS.md` count valuations of `x²-N` over whole intervals by modular-square-root lifting and arithmetic-progression counts. This yields exact upper bounds on B-smooth values, sometimes certified zero, without testing each x. Finite benchmark data also shows that these bounds can be weak or cost more than a direct scan. They do not identify the factor-bearing region or improve the general asymptotic factorization bound.

## Build and test

```
cd vf-factor-lab
make -f Makefile.gmp selftest
make -f Makefile.parallel

cd ../naive_grouping_gmp
make && make check

cd ../vertical_relation_lab
make
python3 tests/test_vertical.py
python3 tests/test_newton.py
python3 tests/test_relations.py
python3 tests/test_one_sided.py
python3 tests/test_prime_power_region.py
python3 tests/benchmark_prime_power_region.py
python3 tests/benchmark_one_sided.py
```

All paths are relative to this research directory (the C++ preexisting repo engine is built independently). Python reference components require the dependencies listed in `vf-factor-lab/requirements.txt`. The C++/GMP standalones do not require Python or FLINT at runtime.

## Research status and cautions

- A sound boundary or candidate selector does not establish a winning factor-bearing region.
- A CTM residual can be smooth while the two-sided congruence remains unusable. Exact GF(2) dependencies and checked GCD recovery are mandatory.
- The one-sided selectors reduce this to standard QS with a different order. Existing evidence does not establish an advantage over QS or GNFS.
- The uploaded experimental `gnfs_gmp.cpp` has an inconsistent algebraic norm for `alpha²=-c` (`a²+c*b²` is required) and is not merged as a correct GNFS implementation.
- Stated resource savings are limited to their specified hypotheses, input sizes and timed environments; conditional `p<q<2p` filters are not unconditional.

## New: exact global factor-band localization and direct Y=2ab Fermat wheel

- `global_vertical_locator/` — exact nested integer-square-root inverse from any vertical Y-interval to possible integer larger-factor values; full mod-m Fermat paired-residue wheel in `(Y,H)`; paired direct `wheel`/`wheel840` and polynomial `batch`/`batch840` factor-region search; independent executable/CRT oracle tests and full-domain benchmarks on actual 39–49-bit semiprimes.
- `global_vertical_locator/math/VERTICAL_Y_WHEEL_AND_GLOBAL_BATCH.md` — the new **sharpness theorem**: every q-residue coprime to `2m` maps to an allowed Fermat Y/H class, so additional direct Y-wheel sparsity is a lattice-coordinate benefit, not a new q-candidate exclusion beyond the same small-prime wheel. A shared polynomial-product GCD, however, can certify entire candidate blocks without individual divisibility tests. This is the established Pollard–Strassen idea, with fourth-root-type arithmetic scaling, **not** GNFS scaling.

Build this new module from `Prime Factorization/hybrid_research/global_vertical_locator` using `make && make test`. Its `Makefile` includes the shared `vf-factor-lab/cpp/gmp_products.hpp` source directly; it does not require separate compiled libraries. The five-case benchmark additionally requires Python `sympy` for creating the deterministic semiprime test inputs; neither factorization executable receives the known factors.

## Additional geometric witnesses, monotone staircase and global N-only Y cutoff

`vertical_witness_bridge/GEOMETRIC_WITNESSES.md` formalizes the distinction between the original squared-integer complex factor sublattice and **all** integer Pythagorean points `H²-Y²=N²`. For distinct odd prime factors `N=pq` the latter has five nonnegative points, **three** yielding a proper factor by `gcd(N,H-Y)` (the original factor point plus two non-square-end witnesses). The old Fermat `(Y,H)` wheel applies only to the former; a separate N² Pythagorean wheel is supplied and independently verified.

A global theorem proves that some nontrivial integer Pythagorean factor point lies at `Y <= (N^(4/3)-N^(2/3))/2`, with an exact integer-only bound from `ceil(cuberoot(N))`. The threshold where the nearest witness changes kind is `p=N^(1/3)`. An exact monotone staircase `floor(sqrt(N²+Y²))-Y` has closed-form plateau boundaries, allowing huge repeated Y ranges to be skipped. These are honest *geometric* compressions, NOT a claimed faster-than-Strassen or GNFS factor algorithm; the remaining distinct u-values still require arithmetic localization.

Run all three self-contained Python oracles from `vertical_witness_bridge/` to verify 595 semiprimes, 7,925 modular witness tests, 464,800 monotone step checks, and 595 global-bound cases.

### Dyadic averages: connection to original iterated-b and vertical green line

See `iterated_average_vertical_bridge/DYADIC_VERTICAL_THEORY.md` and its `dyadic_scan.cpp` benchmark. The exact rounded-average residual is `r - ((r-N) mod 2^k) - 2^k*j`; under an unknown factor it equals `m*p` when `q=m (mod 2^k)`. A monotone Y-dependent factor-match residual gives rigorous conditional interval certificates for each specified m lane, but not a complete factor locator without covering all possible m. `make -C iterated_average_vertical_bridge test` exercises both GMP and elementary theorems.

## Dyadic lane and Gaussian complex-norm extensions (November 2026)

`dyadic_global_attack/` contains the N-only dyadic multiplier cover for the
three vertical Pythagorean witnesses, a lower bound on residue-only universal
lane covers, exact Newton/Mobius squaring conjugacy, an off-strip Gaussian
norm-one group-order factor search, a separate `gcd(a*a+b*b,N)` norm channel,
and bounded in-process hybrid/rho controls.

The Gaussian method is **conditional known p-1/p+1-type mathematics**, not a
new complexity exponent. The power-of-two Newton path often fails, and close
to the factor strip is not consistently better than fixed Gaussian seeds.
One exact no-go theorem: when the norm H=a*a+b*b is invertible modulo N,
`gcd(N,(Y/H)^2-1)` contains the same prime-divisibility information as the
original CTM residual `E=a*a-b*b-N`. This coordinate change alone does not
create new factor information. The *nonunit* norm H provides an additional
separate factor-witness channel.

Run `cd dyadic_global_attack && make test`. For full proofs and original
benchmark controls read `DYADIC_VERTICAL_SYNTHESIS.md` and `README.md`.
The 120-input comparison found a speedup on deliberately smooth-order input
families, and a slowdown on ordinary controls. Both results are reported.

This collection comprises several **independently testable experiments**,
not a single demonstrated faster-than-QS/GNFS factorization algorithm.
