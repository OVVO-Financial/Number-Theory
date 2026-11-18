# VF Factor Lab

A unified factorizer combining the complete Fermat sieve, iterated-average
GCD searches, and fast polynomial divisor products. The components share
exact search bounds, polynomial constructions, and one global budget.
The earlier standalone implementations and experiments remain as controls.
Built and tested on 17 November 2026.

## Parallel search across regions

The C++ / GMP extension now includes a regional worker pool. Workers search
separate Fermat midpoint windows and divisor-product strips, with VF probes
batched into product jobs. Completed exclusion certificates feed one shared
interval ledger; unfinished gaps remain visible even when jobs finish out of
order. Local workers share one deadline and represented-work allowance.

```bash
make -f Makefile.parallel
./build/vf_factor_gmp_parallel 17315933965806733211 --threads 4 --timeout 2
./build/vf_factor_gmp_parallel 978508015703 --threads 4 --full --timeout 1
```

Explicit `--factor-lo` / `--factor-hi` candidate bands can be assigned manually
to separate machines. This release provides local parallelism and portable
range jobs; network dispatch, cross-host cancellation, and shared cross-host
budgets are not implemented. Product GCDs can recover a divisor below an
assigned band when its multiple occurs inside that band. A fully excluded
partial band reports `region_exhausted`, never a global prime claim.

The frozen 26-input, one-second experiment completed 19, 19, and 20 inputs
with one, two, and four regional workers; the unchanged sequential GMP engine
completed 23. Four workers consumed about 3.84 CPU-seconds per wall-second
on the measured workload,
but their median paired speedup over one worker was 0.958x across 19 joint
successes. This demonstrates functioning parallel search, with mixed
performance rather than a general speedup. Both implementations are retained.

See [`cpp/PARALLEL.md`](cpp/PARALLEL.md) for the mathematics, API, and options,
and `results/parallel_report.html` for the complete measured comparison.
`results/parallel/` contains the frozen protocol, all runs, and an independent
audit. Reproduce with `python3 parallel_benchmark.py --outdir
results/parallel-reproduction`; generate the saved report with
`python3 parallel_report.py`.

## Run the C++ / GMP extension

The standalone C++17 implementation uses GMP arbitrary-precision integers and
its own fast polynomial engine. It requires GMP development headers/libraries
and GCC or Clang with `unsigned __int128` support. Python and FLINT are not
runtime dependencies of this executable.

```bash
make -f Makefile.gmp
./build/vf_factor_gmp 978508015703 --full --timeout 1
./build/vf_factor_gmp --metadata
```

The result is `752867 * 1299709`, with no unresolved cofactors. The port retains
the same complete period-160 sieve, reciprocal bounds, VF progression selection,
normalized shared polynomials, GCD recovery, and global recursive budget.
It keeps the existing fast native cursor below 2^64 and uses GMP for larger
Fermat arguments. `--gmp-fermat` forces the GMP cursor at every size.

Inputs, factors, and large bounds are emitted as decimal JSON strings.
Above 64 bits, a GMP probable-prime result is explicitly `probable_prime`;
`--full` lists these separately in `probable_prime_factors`. They are never
silently placed in the proven `prime_factors` list. Use
`--no-probable-prime-check` to continue an exact search on unchecked large primes.

See [`cpp/README.md`](cpp/README.md) for the C++ API, all options, batch mode,
arithmetic details, and test instructions. `results/gmp_report.html` documents
the matched C++/Python benchmark; `results/gmp/` contains its raw evidence.

## Run the Python reference

```bash
python3 -m pip install -r requirements.txt
python3 build_fusion.py
python3 -m vf_factor.fused 978508015703 --full --timeout 10

# Optional comparison kernels, then the complete test suite:
python3 build_native.py
python3 build_wheel20.py
python3 -m pytest -q
```

Python 3.10 or newer is required. Native acceleration and the optional 64-bit
comparison kernels require a C++17 compiler supporting `unsigned __int128`
(GCC or Clang). The polynomial implementation uses python-flint 0.9.0 / FLINT.
The unified algorithm uses its exact Python Fermat cursor when the native
kernel is unavailable. The supplied archive includes kernel sources and build
scripts; compile them locally for the native benchmark performance.

For `vf_factor.fused`, `--timeout` is one cooperative wall-time budget shared
by every search component. With `--full`, the same deadline also covers all
recursive cofactor splits; `--max-steps` likewise shares one work allowance.
A result with `status: budget` has not proved anything about primality.
Unresolved cofactors are explicitly listed. The older `vf_factor.cli` remains
available for individual comparison methods, with its original per-search
timeout semantics.

## The unified algorithm

The primary API is `vf_factor.fused.factor(n, timeout=...)`. It performs common
prechecks once, then maintains a single interval `lo <= d <= hi` for a possible
smaller divisor. A persistent Fermat cursor intersects the complete wheel-20
condition with QR64 into a period-160 wheel; it resumes where it stopped.

If `B` is the next admissible untested Fermat midpoint, every remaining smaller
factor satisfies

\[
d \le B-\left\lceil\sqrt{B^2-N}\right\rceil.
\]

If every divisor below `lo` has been excluded, Fermat needs no midpoint beyond

\[
R_{\max}=\left\lfloor\frac{\mathrm{lo}^2+N}{2\,\mathrm{lo}}\right\rfloor.
\]

These are reciprocal coverage certificates. Ordinary divisor strips are
truncated by the current upper bound; completed strips lower the Fermat cap.
Every required residue class and final tail must be checked before `lo`
advances. VF witness misses do not exclude any divisor interval.

For degree `m`, one polynomial

\[
P_m(X)=\prod_{j=0}^{m-1}(X+j)
\]

serves both step-30 divisor-wheel blocks and step-five VF blocks. Since both
steps are units after the shared small-factor checks,

\[
\gcd\!\left(\prod_j(\alpha+s j),N\right)
=\gcd\!\left(P_m(\alpha s^{-1}\bmod N),N\right).
\]

The requests share one multipoint batch and normally one aggregate GCD.
If its GCD is N, recovery continues through individual blocks and skips
trivial VF-only hits. The core supports arbitrary-precision integers, with
a Python Fermat cursor when the optional uint64 native cursor cannot be used.

The native 64-bit square-root routine uses a hardware estimate followed by
exact integer corrections. Every square and factor is checked exactly;
larger square arguments use the integer-Newton fallback. An independent
45,389-case root validation is included.

Defaults were chosen on a separate 20-input pilot: maximum degree 512,
Fermat time allocation twice the preceding product-batch cost, with a bounded
slice and a single global deadline. All three components remain enabled.
The final fresh corpus uses seed 926483, distinct from both the pilot and
the initial integer-root development experiment.

```bash
python3 fused_benchmark.py --outdir results/fused-reproduction
python3 fused_report.py
```

`results/fused_report.html` contains the focused synthesis report.
`results/fused/` contains fresh comparisons, original-corpus replay, stress
tests, the frozen protocol, and an exact interaction trace. The report script
regenerates the bundled measurements; the reproduction command writes a
separate directory. `results/fused-development-integer-root/` preserves the
first engineering run and its source snapshot. The older experiments below
are retained as historical controls, not silently overwritten.

## What was built

| Method | Purpose |
|---|---|
| `regional_gmp_parallel` (C++ executable) | Regional worker pool with exact coverage intervals, global work reservations, and explicit candidate bands |
| `fused_gmp` (C++ executable) | Standalone GMP port with shared polynomial arithmetic and arbitrary-precision Fermat state |
| `fused` | Original unified Python/FLINT coordinator with the optional native 64-bit Fermat cursor |
| `published_vf` | Exact-integer reproduction of the later Julia `VF(n, 5)` routine |
| `vf_fast` | Same selected average and neighbouring progressions; 32 scalar rounds followed by fast polynomial block searches |
| `vf_fast_multi` | Shared block polynomials across all retained averaging levels |
| `vf_first` | Ablation: identical fast machinery but always average level one |
| `vf_hybrid` | 64 parity-correct Fermat candidates, followed by `vf_fast` |
| `product_trial` | Fast product search over ordinary odd trial divisors, with no averages |
| `trial_division` | Python trial division with a wheel modulo 30 |
| `fermat` | Python parity-correct Fermat with a modulo-64 square filter |
| `pollard_brent` | Reproducible Python Brent rho with batched GCDs |
| `native_published_vf` | Optimized C++ 64-bit reproduction of the published VF routine |
| `native_trial_division` | C++ 64-bit trial wheel |
| `native_fermat` | C++ 64-bit parity-correct Fermat |
| `native_pollard_brent` | C++ 64-bit Brent rho |
| `flint_native` | FLINT's complete integer-factorization routine as an external performance reference |

All methods share the same checks for factors 2, 3 and 5, perfect squares,
and primality below 2^64. Beyond 64 bits, the core methods retain exact integer
arithmetic and do not turn a probable-prime test into a prime certificate.
Native 64-bit methods explicitly reject values outside their range.

No algorithm receives the known factors stored in the benchmark corpus.
Those factors are used solely to validate outputs outside the timed call.

## Exact mathematical basis

With `a0 = ceil(sqrt(N))`, the k-th average is

\[
A_k = \frac{(2^k-1)N+a_0}{2^k}.
\]

For odd N and integral arguments,

\[
\gcd(A_k-b/2^k,N)=\gcd(a_0-b,N).
\]

The integrality condition is `b = a0 - N (mod 2^k)`. Its intersection with a
compatible decimal-ending class has period `5*2^k`. Thus the corresponding
integer GCD arguments form a step-five arithmetic progression.

The later Julia implementation searches integer neighbours of the average
directly. It does **not** impose the PDF's original finite half-gap interval
or identify every searched offset with a true Fermat half-gap. This distinction
is material: the later routine succeeds on the midpoint-parity regressions
that invalidate an unchanged transfer of the PDF's Fermat digit classes.

### Fast block products

For one stream `y_j = alpha + s*j`, build

\[
F_m(X)=\prod_{i=0}^{m-1}(X+s i).
\]

Then `F_m(alpha + s*m*r)` is the product of block r's m candidates.
The implementation uses a balanced polynomial product tree and fast
multipoint evaluation. It does not individually multiply every represented
candidate. For machine-word moduli it uses FLINT `nmod_poly` with an explicit
monic product/remainder tree; for larger moduli it uses native
`fmpz_mod_poly.multipoint_evaluate`.

Monic divisors make remainder-tree evaluation valid over composite moduli.
No inverse of an arbitrary nonzero residue modulo N is assumed.

One polynomial is reused across streams and levels within a factorization.
There is no cross-input polynomial cache; benchmark times include setup.

### Complete coverage and proper factors

After removing factor 5, a step-five progression has a p-divisible term in
every p consecutive indices for any prime p dividing N. If such a term has
GCD N, the term p positions later differs by 5p and has GCD p. Therefore
`2*floor(sqrt(N))+2` indices suffice for a proper factor when N is composite.

The fast VF implementation retains at least one unrestricted descending
progression through that bound. It respects the source routine's upper bound
on ascending bases. A block GCD equal to N is inspected; if it contains only
trivial hits, searching continues. Partial blocks are handled exactly.

The theoretical acceleration is known Strassen-style product-search
technology. This project claims no new factorization exponent. The practical
question tested is whether the averaging schedules improve that machinery.
The fixed maximum polynomial degree is a memory/work-budget setting, not an
implementation of an unlimited asymptotic degree choice.

## Reproduce the benchmark

```bash
python3 -m vf_factor.benchmark --timeout 1 --repeats 3 --outdir results/reproduction
```

The fixed corpus has 83 cases: 75 composites and 8 prime controls. Generated
semiprimes have exact sizes 24, 32, 40, 48, 56, and 64 bits. The families are
balanced nonclose, close, unbalanced, published examples, midpoint regressions,
squares, cubes, and prime controls. A separate seed provides tuning inputs.

The benchmark pins its process to one available logical CPU and sets FLINT to
one thread. Method order is randomized with a fixed seed. Completed cases have
three timed runs; a budget failure is retained and not rerun. Imports and call
dispatch are warmed, but per-input search/polynomial setup is timed.

Every returned factor is checked by exact division. Prime claims are checked
independently. Timeout rows remain censored failures. The report does not rank
a failed search as faster than a successful one, and paired speedups use only
cases both methods completed. Success counts are shown alongside times.

Python and C++ rho use different deterministic RNG implementations, so their
search trajectories differ. Runtime repetitions repeat the same trajectory;
they measure timing variation, not probability of success over random seeds.
The FLINT reference performs a complete prime factorization, while the other
methods aim to find one proper factor. Its native call cannot be interrupted
mid-operation; this reference is restricted to the <=64-bit main corpus.

Run `python3 report.py` to regenerate the supplied report from `results/main`
and `results/stress`. The reproduction command above writes a separate result
directory so that the original measurements remain available.

### Results files

- `results/main/benchmark-summary.csv`: one median row per case/method.
- `results/main/benchmark-runs.jsonl`: every timed run, factor and counters.
- `results/main/benchmark-corpus.json`: fixed integers and verification metadata.
- `results/main/environment.json`: versions, timing protocol, CPU affinity and overall process peak RSS.
- `results/protocol.json`: parameters frozen before the main experiment.
- `results/test-results.txt`: independent correctness test result.
- `results/native-validation.json`: native-kernel correctness checks.
- `results/benchmark_report.html`: analysis and plots.

`represented_candidates` counts positions represented by computed block
products, including products computed before an early exit. It is not a count
of distinct integers or arithmetic operations. Blocks/streams may overlap.
GCD calls, polynomial evaluations, degree, and construction/evaluation times
are separately recorded. Exception-block recovery revisits represented terms.

For `vf_fast` variants, `max_steps` caps represented progression terms plus
explicit preliminary Fermat tests. A block is evaluated only if its whole
range fits the remaining cap. Other methods document their natural budget
unit (divisor tests, GCD candidates, or rho map evaluations). These counts are
not asserted to be interchangeable measures of CPU work.

## Correctness tests

The automated suite includes exhaustive small inputs, independent direct
products over composite rings, negative starts, non-power-of-two block sizes,
noninvertible differences between evaluation points, moduli larger than 64 bits,
GCD=N recovery, finite final blocks, source-table reproductions, exact source
examples, input validation, and zero time/work budgets. A separate native
validation covers uint64 boundaries and exact source search-order equality.

## Complete Fermat wheel-20 follow-up

The additional experiment implements the exact complete table in the linked
*Complete_Fermat_Sieve_Verified.pdf*. It indexes the pair table by `N mod 20`,
retains the correct parity lane, and generates a cyclic midpoint gap table.
An exact square check enforces the full paired residue condition. Four matched
C++ controls independently switch the wheel and the modulo-64 square filter.

```bash
python3 build_wheel20.py
python3 wheel20_factor.py 10093015271 --qr64 --timeout 10
python3 wheel20_benchmark.py --outdir results/wheel20-reproduction
python3 wheel20_report.py
```

The focused report is `results/wheel20_report.html`; raw measurements and
configuration are under `results/wheel20/`. The benchmark reuses all 83 original
inputs and adds three explicitly separated regressions/structural examples.
The report generator reads the bundled `results/wheel20` measurements, while
the reproduction command writes a separate directory.

The full suite with both native libraries built passes 5,554 tests. The
additional 128 wheel tests include all odd residue classes, the 551 regression,
20,000 small-input searches, exact first-solution and candidate counts,
coupled height residues, and 128-bit square intermediates. The new benchmark
has 1,450 raw calls, with no invalid factors, false prime claims, errors, or
mixed repeat statuses.

Wheel Fermat completed all 24 balanced and 18 close-factor inputs, compared
with 20 and 15 for the product baseline. Product search completed all 18
unbalanced inputs, while wheel Fermat completed nine. Overall composite
completion counts are 65/75 for either wheel variant and 68/75 for product
search. The wheel's extra benefit over a matched native Fermat control is
separate from the advantage of Fermat's search order on close-factor inputs.

## Sources and attribution

- Fred Viole, *Complex Space Factorization*, sections 5.1–5.2 and 7.3–7.4, supplied project PDF.
- Fred Viole, *Iterated Averages.md*, supplied project source; source of the reproduced Julia level-selection and step-five schedule.
- [The Complete Fermat Sieve by Terminal Digit: A corrected and verified reference](https://github.com/OVVO-Financial/Number-Theory/blob/master/Number%20Theory%20Papers/Complete_Fermat_Sieve_Verified.pdf).
- Bostan, Gaudry and Schost, [Linear recurrences and polynomial evaluation](https://mathexp.eu/bostan/publications/BoGaSc04.pdf).
- Hittmeir, [A time-space tradeoff for Lehman's deterministic integer factorization method](https://doi.org/10.1090/mcom/3623), including the Strassen construction.
- [python-flint polynomial documentation](https://python-flint.readthedocs.io/en/latest/fmpz_mod_poly.html).
- [FLINT polynomial evaluation](https://flintlib.org/doc/fmpz_mod_poly.html).
- [SymPy primality-testing documentation](https://docs.sympy.org/latest/modules/ntheory.html): reference generation/validation only, and verification of the shared 64-bit Miller-Rabin precheck's range.

No claim of a new complexity exponent or competitive large-RSA factorization
follows from these finite benchmarks. The data records the actual capability
and limitations of this implementation on the specified hardware and corpus.
