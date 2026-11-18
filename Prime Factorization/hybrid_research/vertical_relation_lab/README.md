# Squared complex-plane Newton lattice and exact relation experiment

This package studies the *vertical green line* of the supplied `Complex space.md`: squaring `a+bi` with `a²-b²=N` yields `N+iY`, `Y=2ab`, and `H=sqrt(N²+Y²)`. It is a **reproducible research prototype**, not a full GNFS implementation.

## Build

Requires C++17, GCC or Clang, GMP and GMPXX development headers; Python 3 and SymPy for benchmark input generation.

```sh
make
./vertical_lattice sieve 978508015703 420
./vertical_lattice verify 309 3 103
./vertical_lattice count 978508015703 420 561217300000 561217399999 1
./vertical_lattice search 8051 420 100000 6038 0 1
./relation_bench 978508015703 32768 1021 QS --stop-on-factor
./relation_bench 978508015703 32768 1021 CTM --stop-on-factor
./relation_bench 978508015703 32768 1021 CTMWalk --stop-on-factor
./relation_bench 978508015703 32768 1021 GroupingBand --stop-on-factor
```

`count` arguments: `N m Ylow Yhigh [rational_Newton_steps=0]`, using the exact modular Fermat wheel and a Euclidean floor-sum count of the rational tangent/chord lattice envelope. Zero count certifies *no factor* in that `Y` window under the stated conditions. Positive count is only a possible-candidate count. `search N m max_nodes [Yhi [Ylo [steps]]]` recursively refines the same envelope, with a strict region budget. `region_exhausted` applies only to the supplied region, not to all of `N`.

`relation_bench` arguments: `N K B MODE [--stop-on-factor]`. The modes are `QS`, `CTM`, `CTMWalk`, `Grouping`, `GroupingBand`, `SquareY`, `NearPlus`, `NearMinus`. `CTM` is a Fermat-adjacent boundary projection, while `CTMWalk` runs four residue-aligned CTM staircases with the exact addition-only residual update. `SquareY` selects deterministic uniform samples over the balanced `Y` range. The smooth-relation benchmark never assumes that the *exact-factor* Fermat wheel is a valid filter for nonzero residuals.

Run independent mathematical tests:

```sh
python3 tests/test_vertical.py
python3 tests/test_newton.py
python3 tests/test_relations.py
```

Reproduce experimental data:

```sh
python3 tests/benchmark_newton.py
python3 tests/benchmark_relations.py
python3 tests/benchmark_factor_exit.py
python3 tests/benchmark_bounds.py
```

## The proofs and what the experiment means

See **[`math/FORMALIZATION.md`](math/FORMALIZATION.md)** for exact theorem statements and proofs:

1. The vertical curve `H(Y)=sqrt(N²+Y²)` is monotone, but its nearest-integer residual is not.
2. Every genuine factor point satisfies `H-Y=p²`, `H+Y=q²`; mere integer hypotenuse points can be spurious.
3. Corrected Fermat-digit and QR conditions transport to exact `(Y,H) mod m` pairs through inverses modulo `2m`.
4. Exact rational Newton upper endpoints, convex chords, and tangents permit GMP/Euler floor-sum exclusion of whole regions.
5. For CTM's constructed points, `W=H²-N²-Y²=E(2N+E)`, so smooth `E` is not an independent squared-plane smooth relation.
6. CTM and grouping's two-sided modular relations collapse to `N+E ≡ E (mod N)`: duplicate residuals are duplicate relations.
7. A valid square-congruence comparison must smooth both sides, track signs, form separate GF(2) parities, and confirm a nontrivial gcd. Standard QS starts with its left side already a square.

## Results

See `results/` for independent exact-lattice tests, `newton_window_benchmark.json`, the eight-input `relation_benchmark.json`, `factor_exit_benchmark.json`, and `smoothness_bound_sensitivity.json`.

At 40 bits, 80 evenly placed `Y` windows of width 100,000: the integer-ceiling convexity enclosure certified 8 empty, whereas one exact rational Newton step certified 79 empty. A 128-bit factor-bearing validation window retained the actual point while shrinking from 477 enclosure lattice candidates to 1. These are **local** region results, not an efficient exhaustive RSA search.

At 40 bits with `K=32768`, `B=1021`, standard QS recovered a factor after 2,621 candidate evaluations (~5 ms); CTM projection needed 23,534 (~55 ms), CTM staircase 32,431 (~90 ms), and balanced grouping-band 30,998 (~70 ms). At 55 bits the same QS baseline recovered a factor, but none of these new generators produced a factor within this budget. All methods ran deterministic input schedules, and each run is independently validated. **No improvement over conventional QS or full GNFS has been established.**

Because CADO-NFS was not available and the uploaded `gnfs_gmp.cpp` has an algebraic-norm inconsistency (for `alpha²=-c`, the norm is `a²+c b²`), this package does **not** claim a direct GNFS-to-GNFS benchmark. It provides a mathematically valid, common-factor-base first screen for candidate relation-generating advantages, so that smooth residual counts are not confused with useful independent relations.

## Region-level one-sided smoothness mass (new)

`math/PRIME_POWER_MASS.md` proves and `src/prime_power_region.py` computes the exact total valuations of `x²-N` at all primes up to B **over a whole midpoint interval**, via prime-power modular roots and arithmetic-progression counts, without iterating over the midpoints. It provides an upper bound on the number of B-smooth QS relations in that interval and can certify empty regions. It does not identify the winning factor and may be weaker or slower than ordinary QS sieving; see `results/prime_power_region_benchmark.json`.

```
python3 tests/test_prime_power_region.py
python3 tests/benchmark_prime_power_region.py
python3 src/prime_power_region.py 978508015703 989249 989265 7 --bound
```
