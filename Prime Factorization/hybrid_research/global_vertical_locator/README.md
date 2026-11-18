# Deterministic global q-band locator with the direct Y=2ab Fermat wheel

A C++17/GMP research implementation of the squared complex factor-line sieve,
its exact inverse to integer larger-factor bands, and globally batched polynomial
GCD searches. See `math/VERTICAL_Y_WHEEL_AND_GLOBAL_BATCH.md` for precise
proofs, assumptions, and the distinction between Y-line sparsity and true
factor-candidate search complexity.

## Build and tests

After placing this directory alongside the existing `vf-factor-lab` in
`Prime Factorization/hybrid_research/`:

```sh
make
make test
python3 tests/benchmark_global.py    # needs sympy for deterministic test cases
```

The product-tree engine (`gmp_products.hpp`) is taken from the companion
`vf-factor-lab`; do not replace it with an O(K) product loop in the benchmark.

## Commands

```sh
./y_fermat_wheel 978508015703 420 0 733881011777 --classes
./global_vertical sparse 978508015703 561217282496 561217482496
./global_vertical map 978508015703 1299709
./global_vertical wheel 978508015703 2048 100000000
./global_vertical wheel840 978508015703 2048 100000000
./global_vertical batch 978508015703 2048 100000000 404
./global_vertical batch840 978508015703 2048 100000000 306
```

`N` is the only factorization input. The explicit `Y` values are for
illustrating an interval/known example; the `batch` modes search the entire
**assumed-balanced** interval using N only. On other inputs, an exhausted
balanced interval does **not** imply primality.

`m=420` builds exactly six paired `(Y,H)` residue classes for the example. All
candidate `q mod 840` values coprime to 840 remain possible at the modular
level; the direct Y wheel is not presented as a new way to predict q.

## Core distinction

The C++ product engine can rule out an entire q-block with a single batch
GCD. Its approximately fourth-root-type arithmetic scaling is **known**
(Pollard-Strassen) and remains exponential in input bit length; the 40–49-bit
benchmark is not evidence of GNFS-scale complexity. The direct Y wheel
provides reliable vertical-lattice pruning and constant-factor modular
filtering, not a proved new global exponent saving.

The exact vertical density census is reproducible with
`python3 tests/benchmark_y_wheel_density.py`; it requires no third-party Python
packages. The raw six-modulus census is in `results/y_wheel_density.json`.
