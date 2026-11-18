# Monotone vertical geometry: one-sided square relations and exact region filters

**Status:** exact algebraic theorems, executable GMP implementations and finite benchmarks. This note does not assert improved asymptotic factorization complexity or a GNFS advantage. It extends `FORMALIZATION.md` and the `Complex space.md` construction.

## 1. Exact inverse from squared-green-line coordinate to Fermat midpoint

For an odd input `N>1` define, for a nonnegative **real** `Y`,

    H(Y) = sqrt(N²+Y²)
    A(Y) = sqrt((N + H(Y))/2).

Then `A(Y)` is strictly increasing for `Y>0` because

    A'(Y) = Y/(4 A(Y) H(Y)) > 0.

For an exact odd factor pair `N=pq`, `1<p<=q`, put

    a=(p+q)/2, b=(q-p)/2,
    Y=(q²-p²)/2 = 2ab, H=(q²+p²)/2.

It follows that `A(Y)=a`, and `a²-N=b²`. Thus the *exact* vertical green-line point can be mapped to its Fermat midpoint without rounding error. The deterministic implemented selector for arbitrary integer `Y` uses exact integer square-root ceilings:

    h = ceil(sqrt(N²+Y²)),
    x = ceil(sqrt(ceil((N+h)/2))),
    x = max(x, ceil(sqrt(N))).

This is monotone nondecreasing (not necessarily strictly increasing after rounding). It returns `x=a` at every exact factor point. It is safe as a **candidate map**, not a factor localization certificate: many `Y` values map to the same `x`, and one cannot skip `Y` intervals without an independent proof.

## 2. One-sided smooth relation without two independently smooth factors

Every integer `x>=ceil(sqrt(N))` yields the exact congruence

    x² == x²-N  (mod N).                       (S1)

Therefore, **only** `F(x)=x²-N` needs to factor over a chosen factor base. Given multiple `B`-smooth `F(x_i)`, take a nonzero GF(2) dependency among their signed exponent-parity vectors. For a dependency `D` let

    X = product_{i in D} x_i  (mod N),
    Z² = product_{i in D} F(x_i)  (over integers, sign parity even).

Then `X² == Z² (mod N)` exactly. Recover a nontrivial factor, if any, from `gcd(X-Z,N)` or `gcd(X+Z,N)`; always verify it. The left side is automatically square, unlike `L=N+E` versus `E` in a generic CTM/grouping residual relation.

For a CTM point `(a,b)` with residual `E=a²-b²-N`, we have

    F(a)=a²-N=b²+E.                          (S2)

This is the exact way to turn that CTM point into a **one-sided** quadratic-sieve candidate. The smallness/smoothness of `E` alone does not imply `F(a)` is smooth. The vertical geometry and naive grouping can select candidates `a`; equation (S1) remains identical to standard quadratic sieve. Consequently no new smoothness mechanism has been proved from the coordinate change.

## 3. Three selection schedules tested against the same one-sided QS kernel

`relation_bench` adds these executable modes; each tests smoothness of `x²-N` and uses **the same** parity/dependency extraction as ordinary `QS`.

- **CTMWalkOneSide:** four modulo-20-aligned CTM staircase trajectories, using cached residual increments, offering each distinct `a` only once.
- **GroupingBandOneSide:** under the optional balance premise `p<q<2p`, iterate candidate multiples-of-four `g` in the `N`-derived group band, take `q=g±1`, approximate the partner `floor(N/q)`, and form `a=ceil((q+floor(N/q))/2)`. Each distinct `a` is offered once.
- **VerticalOneSide:** sample deterministic `Y` on `[0,3N/4]`, invert the green line with the GMP ceiling map of Section 1, and offer each distinct `x` once.

`QS` offers `x=ceil(sqrt(N))+j`, every candidate unique. Methods receive the same smoothness base and raw candidate budget, but **do not necessarily offer an equal number of distinct `x`**; the executable reports `unique_one_sided_x`, `valid_relations`, rank, dependencies, factor extraction and wall time. Candidates may duplicate; deduplication is required to avoid mistaking repeated representations for independent information.

These selectors are research baselines, not purportedly optimal or guaranteed faster. In particular, `GroupingBandOneSide` is **conditional on the factor-ratio bound**. Without that premise it does not exhaust the possible factors.

## 4. Exact N-only region exclusion remains separate

`vertical_lattice` implements the necessary square-class transport

    (Y,H) mod m = ((q²-p²)/2,(q²+p²)/2) mod m,

with `p*q=N` and inverses modulo `2m`. It uses the convex chord / tangent and rational Newton upper bounds for `H(Y)=sqrt(N²+Y²)` to count **all integer lattice points inside a certified outer enclosure**, after the modular residues, by Euclidean floor sums. Zero possible points proves that the entire supplied `Y` interval contains no true factor point, without enumerating the `Y` positions. Nonzero points are *not* a factor certificate. The modular factor sieve applies **only to exact factor points**; it must not be imposed on nonzero-residual smooth relations.

An interval-to-midpoint mapping follows from monotonicity of `A(Y)`: if every factor point within `[Y0,Y1]` has a midpoint in `[ceil(A(Y0)),floor(A(Y1))]`, a wheel-20 midpoint test may exclude the entire interval if that *integer* midpoint range contains no admissible residue. This is another necessary `N`-only certificate and is weaker than the full `(Y,H)` envelope in many cases. It does not justify skipping other unexamined intervals.

## 5. Stronger global inference is still the open mathematical step

The exact monotonicity of `H` and `A` does **not** make the nearest-integer residual monotone. Neither does the integer lattice count reveal which among exponentially many surviving regions contains a factor. The currently tested procedures improve local region certification but have not produced a polylogarithmic locator or a subexponential global factorization bound.

An actual breakthrough would be an efficient `N`-only theorem that certifies exclusion of a large *union* of distant `Y` regions without visiting them, or an `N`-specific candidate selector which gives significantly more **independent one-sided smooth relations per total resource unit** than the same QS machinery applied to consecutive `x` values. The benchmark tests this directly and, on the current small corpus, does not demonstrate such an advantage.

## 6. Factor 5 is never a Fermat-wheel lane for the semiprime target

For an odd composite semiprime `N` ending in digit 5, `5 | N`, so `5` is an immediate known nontrivial factor. The general C++ engine performs this precheck **before** `--rsa`'s conditional stage-0 bypass and before the digit wheel. For `gcd(N,10)=1`, the only possible decimal endings of *prime* factors are `1,3,7,9`; the terminal-digit wheel need handle no ending-in-5 semiprime case. The exact `vertical_lattice` program similarly reports precheck factors before using a wheel whose modulus shares a factor with `N`. This policy does not assume that a `5`-ending integer is itself prime or semiprime.

## Reproduction

```
make
python3 tests/test_one_sided.py
python3 tests/benchmark_one_sided.py
./relation_bench 978508015703 32768 1021 QS --stop-on-factor
./relation_bench 978508015703 32768 1021 CTMWalkOneSide --stop-on-factor
./relation_bench 978508015703 32768 1021 GroupingBandOneSide --stop-on-factor
./relation_bench 978508015703 32768 1021 VerticalOneSide --stop-on-factor
```

The output JSON with three-run medians is `results/one_sided_benchmark.json`. It is a finite small-integer experiment and is not a GNFS benchmark.
