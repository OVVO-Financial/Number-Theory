# Prime-power valuation mass: a regionwise one-sided smooth-relation certificate

## Theorem (exact valuation mass without scanning lattice points)

Let `N>1` be odd, let `L<=U` with `L²>N`, and fix a smoothness bound `B>=2` such that `gcd(N,p)=1` for every prime `p<=B` (otherwise a small factor has already been found). Put `F(x)=x²-N>0`.

For each prime `p<=B` and `k>=1`, let `R_(p^k)(N)` be the complete root set of `r² == N (mod p^k)`. Define

    C_(p^k)(L,U) = sum_{r in R_(p^k)} [floor((U-r)/p^k) - floor((L-1-r)/p^k)].

Then exactly

    V_p(L,U) = sum_{x=L}^U v_p(F(x)) = sum_{k>=1} C_(p^k)(L,U).   (P1)

The first equality defines the valuation mass. The second follows from `v_p(z)=sum_{k>=1} 1_(p^k|z)` and counts roots in residue progressions. Stop when `p^k>U²-N`. For odd `p` not dividing `N`, each simple root modulo `p` lifts uniquely to all higher powers via Hensel's congruence. For `p=2`, enumerate the two lifts per root. This computes all `V_p` using only modular roots and floor divisions; **no individual `F(x)` is evaluated**.

## Theorem (certificate and upper bound)

Let `M=L²-N>=2` and define

    P_B(L,U) = product_{p<=B} p^(V_p(L,U)).                  (P2)

This is the product of the `B`-smooth parts of `F(x)` over the entire interval. If `T` values in the interval are completely B-smooth, each contributes its full `F(x)>=M` to this product. Consequently

    M^T <= P_B(L,U).                                       (P3)

Thus the exact upper bound

    T <= min(U-L+1, floor(log(P_B)/log M))                 (P4)

holds whenever `M>=2`. The implementation evaluates the bound with **integer divisions** to avoid logarithmic rounding. In particular,

    P_B(L,U) < M    =>    NO B-smooth F(x) exists in [L,U].  (P5)

The zero certificate is evaluated with multiplication saturated at `M` to avoid constructing huge products. This is a genuine, completely `N`-only **regional relation-exclusion** test, independent of knowing any factor and independent of the two-sided CTM residual.

## Application to the monotone squared green line

The inverse map

    A(Y)=sqrt((N+sqrt(N²+Y²))/2)

is strictly increasing for positive `Y`. Therefore a `Y` window maps to an `A`-axis interval with explicitly computable rational/integer endpoints. The valuation-mass theorem may be applied to that interval to get a certified upper bound on its number of **one-sided quadratic-sieve** smooth relations, without enumerating any individual midpoint `x` in the interval. The exact factor-point Fermat wheel can separately eliminate regions that contain no *exact* factor points, but must **not** be used to discard candidate smooth relations of nonzero residual.

## Complexity and limitations

Computing `V_p` requires enumerating roots modulo `p` (or a faster modular-square-root algorithm), lifting them through `O(log_p(U²-N))` prime powers, and counting each root progression in constant many integer operations. For `p` coprime to `2N`, there are at most two roots modulo each odd prime power and at most four roots for powers of two. The algorithm is sublinear in the window length for **mass counting**.

However, `(P3)` bounds *how many B-smooth relations could occur*, not where they occur or whether they are linearly independent. The bound can be loose for large windows because valuations from different `x` accumulate. Exact computation of `P_B` itself can be expensive for huge windows; use the saturating zero certificate in that regime. No improvement over GNFS or a subexponential factorization method follows from (P1)-(P5).

## Reproducible implementation

`src/prime_power_region.py` implements (P1)-(P5) and the fast zero test. `tests/test_prime_power_region.py` compares valuation sums, bounds, and certified-zero windows with independent brute-force enumeration, covering both `p=2` and odd Hensel lifts.

```
python3 tests/test_prime_power_region.py
python3 src/prime_power_region.py 978508015703 989249 989265 7 --bound
```
