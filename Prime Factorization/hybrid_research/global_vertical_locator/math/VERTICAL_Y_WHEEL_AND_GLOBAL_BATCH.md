# Global factor-band localization from the vertical 2ab coordinate

**Research status (2026-11-17):** exact arithmetic derivations, certified factor-band product GCD, and a matched local benchmark. This is a sound modular wheel and a Pollard–Strassen-type batching specialization, **not** a new subexponential factoring algorithm. The wheel applies to **exact factor points**, not arbitrary CTM or grouping residual relations. The conditional domain `p<q<2p` is not inferred from an arbitrary input.

## 1. Correct squared complex-space identities

For `N=pq` with positive odd `p<=q`, set

```
a=(p+q)/2, b=(q-p)/2, Y=2ab=(q²-p²)/2, H=a²+b²=(q²+p²)/2.
```

Then `a²-b²=N`, `(a+bi)²=N+iY`, `H²=N²+Y²`, `H+Y=q²` and `H-Y=p²`. For the nontrivial factor-bearing point these square identities are essential; **merely having integral H is insufficient**. For odd factor pairs, `Y` is divisible by four because odd squares are 1 modulo 8.

## 2. Spacing and exact inverse

On the continuous factor hyperbola `p=N/q`:

```
Y(q)=(q²-N²/q²)/2,   Y'(q)=q+N²/q³ >0.
```

The derivative is minimal at `q=3^(1/4)*sqrt(N)` with exact minimum `4*sqrt(N)/3^(3/4)`. Therefore any distinct **odd** candidate larger factors, whose q-values differ by at least two, have continuous Y-labels separated by at least `8*sqrt(N)/3^(3/4)`. That is a useful spacing guarantee on the vertical picture, but does **not** reduce the number of candidate q-values.

The inverse is

```
q(Y)=sqrt(Y+sqrt(N²+Y²)).
```

For every nonnegative integer `Y`:

```
floor(q(Y)) = isqrt(Y + isqrt(N²+Y²)).
```

Proof: for integer `k>=0`, `k²<=Y+sqrt(N²+Y²)` is equivalent to `k²-Y<=isqrt(N²+Y²)`. The two integer square roots are exact; no floating point is involved. A ceil is obtained by testing `q^4-2Yq²-N²<0` after taking the floor. Thus `[A,B]` in Y is mapped exactly onto the possible integer q-window `[ceil(q(A)),floor(q(B))]`, and this permits O(1) GMP square-root *window localization*, independent of the number of integer Y positions.

## 3. Full Fermat sieve directly on Y = 2ab

Let `m>=2` and `gcd(N,2m)=1`. Every actual factor point has odd unit residues `r=p (mod 2m)`, `s=q (mod 2m)` satisfying `rs=N (mod 2m)`. For each unit `r` define

```
s=N*r^(-1) (mod 2m)
y_r=((s²-r²)/2) (mod m)
h_r=((s²+r²)/2) (mod m).
```

Because r and s are odd, their squared sum/difference are even. Replacing either chosen odd representative by a multiple of `2m` leaves the half-squares unchanged modulo m. Set

```
C_m(N) = { (y_r,h_r) : r in (Z/2mZ)^× }.
```

Every genuine factor point satisfies `(Y mod m,H mod m) in C_m(N)`. The *projection to Y* is an exact Fermat/CRT sieve on the vertical coordinate, while preserving `(Y,H)` pairs is a strictly stronger necessary condition for lattice-enclosure filtering. For the sample `N=978508015703`, `m=420`, the six paired classes are

```
(36,325), (84,85), (204,325),
(216,325), (336,85), (384,325).
```

Hence only **six of 420 Y classes** survive. With `m=20`, only `(4,5)` and `(16,5)` occur. The class enumeration is deterministic and uses only `N`, extended inverses, and exact modular arithmetic.

**Precheck:** if `gcd(N,2m)>1`, take that divisor first. In particular, a semiprime ending in five is factored by a gcd with 5 and is never passed through the semiprime wheel.

### Critical sharpness theorem: the Y sieve supplies no extra exclusion of unit q-residues

For **each** `q`-residue `s` that is a unit modulo `2m`, choose `r=N*s^(-1) (mod 2m)`. This r is a unit, and the pair `(y_r,h_r)` is in `C_m(N)` by construction. Therefore every unit q-residue can be extended to a *locally valid* `(p,q,Y,H)` tuple modulo `2m`.

So an exact Fermat wheel on `(Y,H)` is valuable for eliminating **Y/H lattice states** and loose tangent/chord enclosures. **Once candidates have already been indexed by q**, transporting that same wheel back cannot eliminate any more q-residues than the ordinary small-prime coprimality wheel. There is no contradiction: most integer Y in the sparse vertical display are not images of ANY integer candidate q; both forms count different objects.

Also, this modular construction is **not admissible for arbitrary nonzero residual E**, since there the congruence is `pq=N+E` rather than `N`.

## 4. Whole-block divisibility certificate and global batching

For any N, consider a finite candidate block `Q={s,s+t,...,s+(d-1)t}` and set

```
P_Q(X)=product_{j=0}^{d-1}(X+t*j) (mod N),
B_Q=P_Q(s) (mod N),
g_Q=gcd(B_Q,N).
```

- If `g_Q=1`, **no candidate q in Q shares any prime factor with N**. The entire block is certified empty of factor-producing q-candidates without checking the individual q's.
- If `1<g_Q<N`, a proper factor is found immediately. A candidate q may be a multiple of a prime factor, not the factor itself.
- If `g_Q=N`, the block may contain distinct prime-factor witnesses or repeated powers; recover by subdividing or checking individual candidates. The C++ implementation does the latter.

All blocks with the same step `t` use **one common polynomial** and evaluate it at many starting points through monic multipoint/product trees over `Z/NZ`. The normalized progression polynomial is valid for `gcd(t,N)=1`; prechecks extract factors 2, 3, 5, and 7 when the `t=840` wheel is used. A fixed small wheel and polynomial degree `d~sqrt(K)` for a total of K possible q-values yield a standard `O-tilde(sqrt(K))` modular-arithmetic batch algorithm, or `O-tilde(N^(1/4))` in the balanced domain K=Theta(sqrt N), subject to the fast-polynomial model and aggregate fallback costs. This is a known Pollard–Strassen class of techniques, **not an improved bound** or a measured GNFS competitor.

The supplied `global_vertical` executable has four matched modes:

```
wheel      : divisibility of unit q residues mod 10
wheel840   : divisibility of unit q residues mod 840
batch      : one shared polynomial in progression step 10
batch840   : one shared polynomial in progression step 840
```

The two `840` modes implement the **exact q-residue image of the Fermat Y-wheel with m=420**. No factor information is supplied to the program, and every output factor is verified against N. For the class `p<q<2p`, `sqrt(N)<=q<sqrt(2N)` is a complete conditional factor range. Exhausting that range is **not a general primality certificate**, and a timeout or incomplete batch is never marked as such.

## 5. Measured result and major limitation

The reproducible three-run C++/GMP benchmark is at `results/global_benchmark.json`. All methods use identical N-only balanced-domain endpoints; a wheel840 run and batch840 run use exactly the same candidate residue classes, as do the wheel10 and batch10 runs. At the largest 49-bit case, batch10 is about 4.76x faster than individual wheel10 GMP checks (44.7ms vs 212.4ms in the measured three-run medians); the 840 wheel is **not consistently faster** than the 10 wheel and does not improve asymptotic complexity. Fixed-deadline/GNFS/RSA claims are not supported.

The relevant research target remains **a stronger N-only region theorem that eliminates a superconstant fraction of q-candidate blocks at lower-than-batch cost** or a demonstrably better smooth-relation generator. Sparsity in the Y coordinate and modular filtering alone are not such a theorem.

## 6. Verification

- `tests/test_y_fermat_wheel.py`: independent Python residue enumeration, preservation of factor points, exact inverse interval counts, every unit q-residue survives in C_m(N), 128-bit test, early factor prechecks.
- `tests/test_global_vertical.py`: C++ executable factor recovery on semiprimes, prime controls, early budget exhaustion, exact (Y,H) reconstruction, and large-integer inverse tests.
- `tests/benchmark_global.py`: deterministic balanced semiprimes, 3 repeats per method, measured runtime, raw JSON. Requires `sympy` only for constructing test-case inputs; the factorizer itself never receives factors.
- `Makefile`: `make`, `make test` (GMP+gmpxx and Python standard library). No GPU measurement or GNFS comparison.

### Whole-domain density demonstration

For `N=978508015703` in the **conditional balanced** domain, `[0,floor(3N/4)]` has 733,881,011,778 integer Y values, while `[ceil(sqrt N),floor(sqrt(2N))]` contains just 409,738 integer q candidates, 163,896 of which are units modulo 10. The direct Fermat wheel has:

| m | allowed Y classes | allowed Y density | allowable unit q candidates (mod 2m) |
|---:|---:|---:|---:|
| 20 | 2 | 10.000% | 163,896 |
| 60 | 2 | 3.333% | 109,264 |
| 420 | 6 | 1.429% | 93,655 |
| 840 | 6 | 0.714% | 93,655 |
| 2520 | 18 | 0.714% | 93,655 |
| 9240 | 30 | 0.325% | 85,142 |

The dramatic Y density reduction is real, but even at m=9240 it leaves roughly 2.38 billion integer Y positions in this domain, compared to tens of thousands of q candidates. The correct order is to use **the exact Y-to-q inverse first**, then filter or batch q. In particular, the unchanged q count at m=840 versus 420 proves that an extra factor of 2 in the Y modulus adds a representational restriction without eliminating any additional q candidates. Source: `results/y_wheel_density.json`.
