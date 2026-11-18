# Unified factorization in C++ with GMP

This is a standalone C++17 port of the existing unified Fermat, VF, and
divisor-product factorizer. Input integers, factors, Fermat state, averages,
and factor bounds use `mpz_class`. The executable has no Python or FLINT runtime
dependency. GMP is linked dynamically; no GMP library source is vendored.

The additional regional worker-pool implementation is documented in
[`PARALLEL.md`](PARALLEL.md). Build it with `make -f Makefile.parallel` and use
`--threads` to select concurrency. The original sequential executable below
remains available as a separately benchmarked implementation.

## Build and run

The following commands are for Linux or WSL, from the project root. Install
GMP development headers/libraries and a C++17 compiler first. The retained
64-bit kernel requires GCC or Clang support for `unsigned __int128`.

```bash
make -f Makefile.gmp
./build/vf_factor_gmp 978508015703 --full --timeout 1
./build/vf_factor_gmp --metadata
```

Equivalently, compile directly:

```bash
mkdir -p build
g++ -O3 -std=c++17 cpp/gmp_factor.cpp vf_factor/fusion_kernel.cpp \
    -lgmpxx -lgmp -o build/vf_factor_gmp
```

`CXX`, `CPPFLAGS`, `CXXFLAGS`, and `LDFLAGS` can be supplied to Make. For example,
set include/library flags when GMP is installed in a nonstandard prefix.
`python3 build_gmp.py` is an optional convenience wrapper; Python is not needed
for the Make or direct compiler build.

The example returns proven prime factors `752867` and `1299709`. Input and
factor values are decimal strings in JSON so consumers cannot round them to
floating-point numbers. The CLI accepts decimal integers only.

## Options and result semantics

Run `./build/vf_factor_gmp --help` for the complete option list.

| Option | Meaning |
|---|---|
| `--timeout SECONDS` | One cooperative time limit; default 10 seconds |
| `--max-steps U64` | One represented-work allowance |
| `--full` | Recursively factor using the same deadline and work allowance |
| `--trace` | Record exact coverage bounds and batch decisions |
| `--batch` | Read one decimal input per line and emit one JSON result per line |
| `--gmp-fermat` | Force GMP Fermat arithmetic below 2^64 as well |
| `--initial-degree INT` | Initial polynomial degree; default 16 |
| `--max-degree INT` | Degree cap; default 512, accepted range up to 65536 |
| `--fermat-ratio FLOAT` | Fermat time allocation / preceding product time; default 2 |
| `--initial-fermat-tests U64` | Initial admissible midpoint tests; default 64 |
| `--no-vf` | Disable VF probes for a controlled comparison |
| `--no-feedback` | Disable reciprocal factor-bound feedback |
| `--no-sharing` | Build separate step-specific polynomials |
| `--no-fermat` | Run the product side alone |
| `--no-probable-prime-check` | Do not stop on an unproved large probable prime |

The work cap counts scalar divisor tests, positions represented by successfully
evaluated polynomial blocks, and exact Fermat square tests. Exceptional recovery
revisits already represented positions and does not charge them twice. Prechecks
do not consume represented search work. Their arithmetic time is included in the
global deadline. An individual GMP operation can finish slightly after a
deadline; timeouts do not interrupt GMP internally.

Single-search results have status `factor`, `prime`, `probable_prime`, or
`budget`. A factor is always a verified proper divisor. Below 2^64, primality
uses a deterministic Miller–Rabin base set. Above 2^64, GMP's return code 1 is
reported as `probable_prime`; only code 2 is reported as proved `prime`.
The probable-prime check uses 25 repetitions as its GMP parameter. Disabling
it leaves an unchecked large prime to the exact exhaustive search, which may
consume the budget.

Full-factorization output separates:

- `prime_factors`: proven prime factors;
- `probable_prime_factors`: unproved probable-prime factors;
- `unresolved`: cofactors whose search exhausted its budget.

Their product is the original input, including multiplicities. Full status is
`complete`, `probable_complete`, or `partial`, respectively. A new recursive
cofactor is not started after the global deadline expires.

Batch mode starts a new independent budget for each input line. Invalid lines
produce an error object and do not stop subsequent lines; the final exit code
is nonzero if any line failed validation. Standalone malformed arguments produce
an error on stderr with a nonzero exit code. Blank input lines are invalid.

```bash
./build/vf_factor_gmp --batch --timeout 1 < inputs.txt > results.jsonl
```

## What was ported

The coordinator preserves the Python algorithm's complete wheel-20/QR64
intersection, period-160 jumps, exact reciprocal factor bounds, VF level
selection, full divisor strips, shared polynomials, and aggregate-GCD recovery.
No unfinished strip is marked excluded. Failed VF probes never certify a
divisor interval.

The native 64-bit Fermat kernel is reused without modifying its source.
Larger inputs use a persistent GMP cursor with a difference recurrence,
`mpz_perfect_square_p`, and exact `mpz_sqrt` on a square hit. The forced-GMP flag
allows a matched comparison of the two cursor implementations below 2^64.

The GMP polynomial engine implements:

1. Balanced product trees for block polynomials.
2. Carry-free coefficient packing into GMP integers (Kronecker substitution)
   for larger polynomial multiplications.
3. Monic remainder computation by reversal and Newton series inversion.
4. Balanced product/remainder trees for multipoint evaluation.
5. Small schoolbook operations where packing overhead would dominate.

For normalized steps, the same polynomial `P_m(X) = product(X+j, j=0..m-1)`
serves both the step-30 divisor wheel and step-five VF requests. After common
factor checks, those steps are invertible modulo N. Remainder-tree divisors are
monic, so the polynomial arithmetic works over composite moduli without assuming
every nonzero coefficient is invertible.

The degree cap is a resource setting. Arbitrary-precision arithmetic does not
imply that balanced integers of arbitrary bit length can be factored quickly.
The port makes no new asymptotic complexity claim.

## C++ API

The coordinator is header-based:

```cpp
#include "cpp/gmp_factor.hpp"
#include <iostream>

int main() {
    vf::Options options;
    options.timeout = 1.0;
    options.max_steps = 1000000;
    const auto result = vf::factor(vf::Z("978508015703"), options);
    std::cout << result.status << " " << result.factor << "\n";
}
```

Compile the caller together with `vf_factor/fusion_kernel.cpp`, and link
`-lgmpxx -lgmp`. `vf::factor_full` returns the three cofactor lists. A caller
managing its own batch can use `vf::find_factor(n, options, budget)` with one
`vf::Budget` instance to share limits explicitly.

## Validation and benchmarks

```bash
make -f Makefile.gmp selftest

# Integration and existing-reference tests (Python test dependencies required):
python3 -m pip install -r requirements.txt
python3 build_fusion.py
python3 build_native.py
python3 build_wheel20.py
python3 -m pytest -q

# Matched benchmark, writing a separate output directory:
python3 gmp_benchmark.py --all-mpz-control --outdir results/gmp-reproduction
```

The C++ arithmetic self-test covers exact roots, deterministic/probable primality,
native/GMP cursor agreement, large-input resume/coverage, and reciprocal bounds.
The polynomial self-test compares packed products and fast monic remainders
against scalar references, including composite rings, duplicate points, signed
progressions, carry bounds, and GCD-N recovery. Integration tests compare exact
fixed-work traces and counts with Python and validate global recursive budgets.

Benchmark inputs are generated before the timed calls. Known factors are used
only for output validation. Internal search time is the primary timing;
external request wall time is retained separately. Persistent batch mode excludes
process startup. Both implementations include per-input polynomial construction.
The forced-GMP control below 2^64 measures cursor arithmetic separately from the
general C++/Python comparison. All methods retain the frozen scheduler settings.

See `../results/gmp_report.html` and `../results/gmp/` for measured outcomes,
settings, compiler/GMP versions, hashes, and every raw timed call. The older
one-second experiments are separate from this matched half-second comparison.

## Primary documentation

- [GMP C++ interface](https://gmplib.org/manual/C_002b_002b-Interface-General)
- [GMP primality and GCD functions](https://gmplib.org/manual/Number-Theoretic-Functions)
- [GMP integer roots](https://gmplib.org/manual/Integer-Roots)
