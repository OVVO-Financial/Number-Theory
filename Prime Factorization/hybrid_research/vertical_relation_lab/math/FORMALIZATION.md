# Squared complex factor line: exact Newton lattice certificates and relation accounting

**Scope.** This note formalizes mathematical consequences of the user's [`Complex space.md`](https://github.com/OVVO-Financial/Number-Theory/blob/master/Complex%20space.md) and the supplied `complex_space_factorization.cpp`, `gnfs_gmp.cpp`, and earlier monotone-corridor prototype. The original manuscript introduces `z=a+bi`, `N=a²-b²`, its green vertical line `z²=N+i(2ab)`, and Fermat terminal-digit restrictions. Its line stating `a²+b²=N` while explaining the squared construction is inconsistent with its own Fermat formula; the correct equation is `a²-b²=N`, whereas `a²+b²=H` denotes the hypotenuse. The precise lattice counting, rational Newton enclosure, necessary conditions, residual collapse, and relation comparison here are **new derivations**, not theorems attributed to that manuscript. Proofs use exact integer/rational arithmetic throughout.

## Theorem 1. Vertical green line is a monotone *real* curve

Assume `N>0`, `a>=b>=0` are integers, and put

    E = a²-b²-N,
    Y = 2ab,
    H = a²+b².

Then `z² = (N+E)+iY` and, identically,

    H²-Y² = (a²-b²)² = (N+E)².                (1)

For a **true** factor point `E=0`, we obtain

    H = sqrt(N²+Y²),
    (H-Y)(H+Y) = N²,
    H-Y = (a-b)² = p²,
    H+Y = (a+b)² = q²,
    pq = N.                                             (2)

For real `Y>=0`, `h(Y)=sqrt(N²+Y²)` has

    h'(Y) = Y/sqrt(N²+Y²) >=0,
    h''(Y) = N²/(N²+Y²)^(3/2) >0.               (3)

Hence `h` is strictly increasing for `Y>0` and strictly convex. Holding `H=k` *fixed*, the residual `D_k(Y)=k²-N²-Y²` is strictly decreasing for positive `Y`. But choosing `k=ceil(h(Y))` changes `k` discontinuously, so `D_ceil(Y)` is NOT monotone. E.g. for `N=15`, the residual at `Y=0,4,8,12,16` equals `0,15,0,31,3`. A Newton iterate at fixed `Y` determines the real hypotenuse; it does not determine which `Y` yields integer factor coordinates.

The alternative exact conditions are

    a² = (H+N)/2,    b² = (H-N)/2,    Y=2ab.  (4)

Thus checking only whether `h(Y)` is integral is insufficient. The squared conditions `(H±Y)` in (2) are both necessary and sufficient for a positive integer factor pair (with the usual odd-parity assumptions for recovering integer `a,b`).

## Theorem 2. Exact count of integer lattice intersections and genuine factor points

For odd `N>0`, `H,Y` nonnegative integers with `H²-Y²=N²` correspond bijectively to the **unordered** positive divisor pairs `(u,v)` of `N²`, via

    u=H-Y, v=H+Y, uv=N²,
    H=(u+v)/2, Y=(v-u)/2, u<=v.

All divisors of odd `N²` are odd, so `H,Y` are integers. Because `N²` is a perfect square, exactly one divisor pair is central (`u=v=N`). Thus the number of nonnegative Pythagorean points is

    P(N) = (tau(N²)+1)/2.                        (5)

The *factor-bearing* points are exactly those where **both** `u,v` are perfect squares. They have the form `(u,v)=(p²,q²)` with `pq=N, p<=q`. The number is

    F(N) = ceil(tau(N)/2).                        (6)

For distinct-prime `N=pq`, this means 5 Pythagorean points but only 2 actual positive factor pairs, of which `(1,N)` is trivial. Example: `N=309=3*103`. `(H,Y)=(515,412)` obeys `H²-Y²=N²`, but `H-Y=103` and `H+Y=927` are not squares, so this lattice point does **not** recover a factor pair. The actual nontrivial point is `(H,Y)=(5309,5300)`.

## Theorem 3. Exact Fermat/CRT wheel transported to `Y=2ab`

Fix any `m>=1` with `gcd(N,2m)=1`, so possible factors `p,q` are units modulo `2m`. For each unit `r` modulo `2m`, set

    s = N*r^(-1)  (mod 2m),
    y_r = (s²-r²)/2  (mod m),
    h_r = (s²+r²)/2  (mod m).                    (7)

Both numerators are even because `r,s` are odd. Their halves modulo `m` are well-defined independently of the integer representatives modulo `2m`: replacing `r` with `r+2m` changes `r²/2` by a multiple of `m`. Define the finite allowed set

    C_m(N) = {(y_r,h_r) : r in (Z/2mZ)^x}.

For every actual factor point, `(Y mod m,H mod m) in C_m(N)`. This is proved by choosing `r=p mod 2m`, which gives `s=q mod 2m`, and substituting (2). It never rejects an exact factor point, assuming the gcd precondition. If `gcd(N,2m)>1`, factor extraction must precede the wheel; the implementation explicitly rejects such a wheel construction instead of silently excluding that factor.

For `m=420=4*3*5*7` and `N=978508015703` (coprime to `2m`), `C_m(N)` contains 6 `(Y,H)` pairs:

    (36,325), (84,85), (204,325),
    (216,325), (336,85), (384,325).

Only six distinct `Y mod 420` classes remain. These are necessary conditions for **exact factor points**. They are NOT a completeness-preserving sieve for generic `E != 0` smooth relations. A nonzero residual changes the underlying congruence to `pq=N+E`, so its admissible residues must be conditioned on `E` rather than `N`.

## Theorem 4. Certified intervals on the vertical line from rational Newton bounds

Let `0<=A<=Y<=B`, and let `h(Y)=sqrt(N²+Y²)`. Construct rational **upper** bounds `U_A>=h(A)`, `U_B>=h(B)` by starting with `u/v=ceil(sqrt(x))/1` for `x=N²+A²` or `x=N²+B²` and applying an exact Newton step

    (u/v)_new = (u² + x v²)/(2uv).             (8)

Indeed `(u²+xv²)/(2uv)-sqrt(x) = (u-v sqrt(x))²/(2uv)>=0` for positive `u,v`. Successive steps quadratically improve the approximation. All numerators and denominators are GMP integers, so rounding cannot invalidate the upper bound.

Convexity and the tangent inequality give the exact enclosure

    L(Y) = (N²+A Y)/U_A <= h(Y),
    h(Y) <= U(Y) = [U_A(B-Y)+U_B(Y-A)]/(B-A)      (A<B).   (9)

The first inequality uses the true tangent at `A`, equal to `(N²+A Y)/h(A)`, then replaces the denominator by a larger `U_A`. The second inequality uses the chord property for a convex graph, then replaces endpoint values by upper bounds. For `A=B`, the enclosure uses `L=(N²+A²)/U_A` and upper `U=U_A`.

Let `(ry,rh)` be an allowed residue pair modulo `m`. Every possible factor point `(Y,H)` in `[A,B]` lies in one of the intervals

    Y=y0+m j,   H=rh+m k,
    L(Y) <= H <= U(Y),                            (10)

where `y0>=A` is the first `ry mod m` position. Its integer-lattice upper bound count is

    sum_j [ floor((U(y0+mj)-rh)/m)
          - ceil((L(y0+mj)-rh)/m) + 1 ].         (11)

Each term is a floor of a rational *affine* function of `j`; the code evaluates the sum with exact Euclidean floor sums, without iterating over all `Y` in the region. If the total over all allowed pairs is **zero**, there is no exact factor point in that interval. If it is positive, the interval must be further checked; some surviving lattice points may be spurious Pythagorean points, or merely points inside a tangent/chord enclosure.

A bounded recursive search is therefore complete **within its explicitly supplied `Y` domain**: discard intervals with count zero, split nonempty intervals, and at leaves verify the exact perfect-square conditions `(H±Y)` and `p*q=N`. An exhausted *partial* region is never a primality certificate for `N`.

For a balanced semiprime `p<q<2p`, the factor's `Y` satisfies

    0 <= Y = (q²-p²)/2 < 3N/4,                   (12)

since `q²<2N` and `p²>N/2`. Thus `[0,floor(3N/4)]` is an `N`-only search domain **conditionally** covering balanced inputs. Other inputs can have much larger `Y` and cannot be excluded from an incomplete search.

## Theorem 5. Squared-plane residual collapse and GCD equivalence

For *any* integer `(a,b)` and `E=a²-b²-N`, set `H=a²+b²`, `Y=2ab`. Then

    W = H²-N²-Y² = E(2N+E).                    (13)

Hence, modulo `N`, `W = E² (mod N)`, and

    gcd(W,N)=gcd(E²,N).                          (14)

The two gcds `gcd(W,N)` and `gcd(E,N)` detect the same prime divisors of `N`, although the exponents within the gcd can differ. Therefore the squared-plane residual of a **constructed CTM point** does not independently improve the set of prime factors detectable by a gcd. This collapse does NOT hold for arbitrary independently chosen lattice points `(H,Y)` that are not built from integer `a,b`.

Likewise, `E` being smooth does not imply that `W` is smooth: `W` includes the additional potentially huge factor `2N+E`. Conversely a small `|W|` requires a small `E` unless cancellation or special divisibility occurs. Sorting by `W` alone does not yield a monotone integer-factor search.

## Theorem 6. CTM/grouping residual relations collapse to smooth-neighbor congruences

From CTM, write `L=(a-b)(a+b)=N+E`. Therefore

    L congruent E (mod N).                     (15)

A small (even smooth) `E` alone does not give a quadratic-sieve square congruence. One valid but restrictive relation protocol is to require **both** `L` and `|E|` to factor over the same prime base. Record their exponent parity vectors, with a sign bit when `E<0`. A subset whose parity vectors vanish *separately on both sides* gives integer squares `X²=product L_i`, `Y²=product E_i`, whence `X²=Y² (mod N)` and a possible nontrivial gcd of `X±Y` with `N`.

For naive grouping `N=g*floor(N/g)+r` with `0<=r<g`, set `E=-r` and `L=g*floor(N/g)=N+E`. This is **the same** modular relation (15), with a negative residual. Therefore, under this precise relation protocol, the geometry affects which residuals are visited and at what cost, but the relation itself depends **only on `N` and the distinct residual `E`**. Different grouping chains or CTM paths producing the same `E` do NOT create new independent relations. The experiment explicitly deduplicates `E` before smoothness and rank accounting.

For a generic squared-line trial point `(H,Y)`, set `L=H²-Y²`, `E=L-N²`; then `L=E (mod N)`, and the same two-sided protocol applies. When `(H,Y)` arises from integer CTM `(a,b)`, however, `L=(N+E_CTM)²` and the residual is `E_CTM(2N+E_CTM)`: its square-congruence identity is tautological modulo `N`. Such points must not be mislabeled as new independent smoothness information merely because the squared residual is small.

By contrast, ordinary quadratic sieve (QS) uses

    x² = N + f(x),   x² congruent f(x) (mod N).  (16)

Its **left side is already a square**. It therefore requires only `f(x)` to be smooth, and the parity matrix has roughly one factor-base block rather than two. This structural advantage must be charged in any comparison.

## Experiment protocol (implemented, not assumed)

- **Exact wheel/region experiment:** `vertical_lattice`. For `gcd(N,2m)=1`, generate `(Y,H)` classes from (7); use (8)–(11) with 0/1/2 rational Newton steps; report zero-count certified regions; verify all factor candidates with GMP exact squares.
- **Relation experiment:** `relation_bench`. Use the **same** prime base, exact smoothness factorization, explicit modular congruence check, deduplication by residual, two-sided parity as required, exact binary matrix elimination, and nontrivial gcd extraction.
- **Methods:** QS, CTM near-boundary projection (`a=ceil(sqrt N)+j`, `b=floor(sqrt(a²-N))`), actual four-lane modulo-20 CTM staircase (`CTMWalk`), naive odd seed grouping (`Grouping`), balanced group bands (`GroupingBand`), uniform vertical-line `Y` samples (`SquareY`), and direct near-N smooth-neighbor controls (`NearPlus`, `NearMinus`). Only `GroupingBand` requires the supplied `p<q<2p` premise. The exact-factor modular wheel is **not** misapplied to arbitrary residual relations; CTMWalk deliberately prefers one admissible mod-20 region, without claiming to cover all relation candidates.
- **Fixed budget:** eight deterministic semiprimes with 28–63 actual bits, three runs each, 32,768 candidate positions, prime base up to 1,021. Factorization outputs are checked against the original `N` without supplying factors to executable inputs. Additional bound sensitivity uses `B=257,1021,4093` on three cases.
- **Accounting:** record raw smooth `E`, distinct `E`, smoothness of **both** required sides, rank, dependencies, verified factors, candidate evaluations and median CPU-kernel time. Source reports contain raw independent results.
- **Controls/limitations:** QS here is a straightforward C++/GMP *relation-generator* baseline with trial smoothness and a GF(2) dependency engine, not a production QS sieve; `gnfs_gmp.cpp` is not a validated GNFS program and CADO-NFS was unavailable. The benchmark does **not** establish relative performance against full GNFS or cryptographic RSA sizes. In particular, 32,768 samples at a fixed small factor base do not measure subexponential asymptotics.

## Interpretation

The exact Newton vertical lattice enclosure produces **real, certified N-only region exclusions**; one step is dramatically better than an integral-ceiling enclosure on the tested windows. But the full balanced `Y` domain has length `Theta(N)`; to cover it with short, nearly linear windows of width `Delta` would require `Theta(N/Delta)` windows unless a further global argument permits safe skipping. No GNFS-level scaling follows.

For modular relation generation, some CTM and grouping paths do produce valid dependencies and genuine factors. They have NOT shown a consistent advantage over the equally budgeted conventional QS relation generator. CTM's abundant smooth **single residuals** often disappear when the second, necessary smoothness condition is imposed. The strongest possible further development would discover a **new valid one-sided square relation** or an **efficiently discoverable global constraint** on factor-bearing `Y` bands; simply counting smooth `E` values does not suffice.
