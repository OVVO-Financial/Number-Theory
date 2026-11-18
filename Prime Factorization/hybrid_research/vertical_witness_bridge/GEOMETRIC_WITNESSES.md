# Additional factor witnesses and the exact monotone staircase on the vertical green line

**Research status.** Exact elementary theorems and independently tested modular examples. This is not a faster factoring algorithm or an RSA-scale locator. The original factor-specific terminal-digit wheel must be distinguished from the wider Pythagorean witness wheel.

## 1. Two different integral lattices

The original complex-square construction fixes `z=a+bi` with integers `a>=b>=0` and `a²-b²=N`. It gives

- `Y=2ab`, `H=a²+b²`;
- `(H-Y,H+Y)=((a-b)²,(a+b)²)=(p²,q²)`;
- `H²-Y²=N²`.

Its image is **a strict subset** of the positive integer solutions `(Y,H)` to `H²-Y²=N²`. In particular, the additional points below do NOT have `Y=2ab` for integral `a,b`; they are integer points on the *same vertical hyperbola* obtained by allowing the squared-plane coordinates `(H,Y)` directly. This extension is a new search space, not a claim that the source's original square roots were already integral there.

For odd `N`, every `H>=Y>=0` solving `H²-Y²=N²` is uniquely represented by an odd divisor `u|N²`, `1<=u<=N`:

```
v = N²/u; H=(u+v)/2; Y=(v-u)/2; u=H-Y; v=H+Y.
```

As u and v are odd, H and Y are integers. Conversely every solution yields that divisor pair.

## 2. Exactly three nontrivial factor-revealing points for a distinct odd semiprime

Let `N=pq` with odd primes `3<=p<q`. The divisors `u` of `N²=p²q²` with `u<=N` are **exactly**

```
  u in {1, p, q, p², N}.
```

Therefore the nonnegative vertical hyperbola has exactly five integer points. For the three intermediate choices `u=p,q,p²`, `1<gcd(u,N)<N` returns a factor. The other two (`u=1,N`) are trivial.

- `u=p²`: `Y=(q²-p²)/2`, the original integer complex factor point. `gcd(Y,N)=1`.
- `u=q`: `Y=q(p²-1)/2`, a **non-square-end** Pythagorean point, and `gcd(Y,N)=q`.
- `u=p`: `Y=p(q²-1)/2`, a second **non-square-end** point, and `gcd(Y,N)=p`.

Because `Y(u)=(N²/u-u)/2` is strictly decreasing for positive u, the *nearest nontrivial* Y coordinate is obtained from `max(p²,q)`. It is the new `u=q` point iff `q>p²`. In the balanced class `p<q<2p`, we always have `q<p²`, so the original factor point remains nearest.

Examples:

| N | original factor-point Y | nearest nontrivial Pythagorean Y | reason |
|---:|---:|---:|---|
| 309=3*103 | 5300 | **412** | q>p² |
| 187=11*17 | **84** | 84 | q<p² |
| 978508015703=752867*1299709 | **561217382496** | same | q<p² |
| 1009000007063=1009*1000000007 | 500000006999490984 | **509040003563280** | q>p² |

This extra multiplicity does not imply the coordinates can be *found* quickly from N alone. Finding u with nontrivial gcd(N,u) is itself factorization.

## 3. A different Fermat wheel is required to preserve all three witnesses

For a modulus m with gcd(N,2m)=1, the square-factor sieve enumerates odd units r=p (mod 2m), s=N*r^-1 (mod 2m), and stores

```
 ((s²-r²)/2 mod m, (s²+r²)/2 mod m).
```

That preserves ONLY the original integer-complex-factor sublattice; it can reject the extra u=p and u=q Pythagorean witnesses. Example N=187, m=20: the extra witness `u=q=17`, with `(Y,H)=(1020,1037)` and `(Y mod 20,H mod 20)=(0,17)`, is rejected by the square-only wheel.

For the full integer Pythagorean hyperbola, enumerate odd unit residues `r=u (mod 2m)` and `s=N²*r^-1 (mod 2m)` and store the **different** paired classes

```
 ((s-r)/2 mod m, (s+r)/2 mod m).
```

Every genuine Pythagorean point is retained. It is a necessary, not sufficient, modular condition. Example N=187:

| modulus m | square-point paired classes | all Pythagorean paired classes | Pythagorean Y-residue count |
|---:|---:|---:|---:|
| 20 | 2 | 8 | 3 |
| 60 | 2 | 16 | 3 |
| 420 | 6 | 96 | 9 |

All allowed residue classes are *unit-locally* compatible: for every unit r mod 2m, the complementary s exists. Thus finite local congruences alone cannot remove otherwise admissible u or q residue classes after reparameterizing back. More powerful region exclusion must use the simultaneous global size inequalities and congruences, not modular admissibility alone.

## 4. An EXACT monotone integer coordinate and closed-form plateau boundaries

For any nonnegative integer Y, define

```
h(Y)=sqrt(N²+Y²),
U(Y)=floor(h(Y))-Y.
```

For Y>=0, h'(Y)=Y/h(Y) in [0,1), so `floor(h(Y+1))-floor(h(Y))` is either 0 or 1. Therefore `U(Y+1)-U(Y)` is **0 or -1**: U is nonincreasing, with no downward skips. Starting at U(0)=N, it traverses *every* integer level down to 0 as Y tends to infinity. This is the true monotone staircase behind the sparse vertical display.

For each integer k with 1<=k<=N, the threshold is exactly

```
U(Y)>=k  <=>  N² >= 2*k*Y + k²
           <=>  Y <= floor((N²-k²)/(2*k)) = T(k).
```

In particular, all integer Y values attaining the *same* U(Y)=k are

```
T(k+1)+1 <= Y <= T(k).
```

This interval can be jumped over in O(1) large-integer operations without enumerating its Y positions. For `N=978508015703` and `k=1299709`, the corresponding integer u level repeats over **283,404,141,794** consecutive Y positions. That is genuine Y-axis sparsity, but all these positions carry the *same* u-level, so skipping them does not skip independent factor candidates.

For odd k, the rightmost Y=T(k) is an exact Pythagorean point iff k divides N² (then `(N²-k²)/(2k)` is integer). The fractional part of `(N²/k-k)/2` is therefore a **divisor test**, but is a sawtooth in k, not a globally monotone residual. (For even k the numerator is odd and cannot be divisible by 2k when N is odd.)

For small Y satisfying `Y²<2N+1`, `floor(h(Y))=N`, so `U(Y)=N-Y` and `gcd(U(Y),N)=gcd(Y,N)`: a vertical-line GCD search near the origin is *exactly trial division in different coordinates*.

## 5. A sound whole-Y-region batch-GCD certificate

For integer interval `[A,B]`, let

```
lo = max(2, ceil(sqrt(N²+B²))-B),
hi = min(N-1, floor(sqrt(N²+A²))-A).
```

Here `ceil(sqrt(T))` can be computed exactly from `isqrt(T)` and a square check. Because U is nonincreasing, every **nontrivial** Pythagorean witness `(H,Y)` with `A<=Y<=B` has integer `u=H-Y` in `[lo,hi]`. The proof is independent of the factors.

Define (with empty-product convention)

```
G_[A,B] = gcd(N, product_{u=lo}^{hi} u).
```

- `G=1` rigorously excludes all nontrivial integer Pythagorean witnesses in the whole Y interval. No individual Y point is examined.
- `1<G<N` yields a verified nontrivial factor.
- `G=N` is ambiguous; subdivide and retry. Any finite block product can be formed modulo N.

**Critical complexity caveat:** this is exactly a batch-GCD divisor interval after the coordinate inversion. The existing polynomial/product-tree backend can compute it without per-candidate GMP division, but it does **not** provide a better asymptotic bound than standard batch factor search. Some u in the interval may be nonzero multiples of a factor yet not divide N²; this can trigger positive GCD without a Pythagorean point, which is still a valid factor discovery.

For genuine integer-complex factor points, impose u=p² a perfect square, or equivalently additionally check both H-Y and H+Y are squares; doing so returns to the original factor-pair search domain. The broader witness search cannot reuse that square-only filter.

## 6. Why the projected CTM residual does not provide an independent relation

For a constructed integer CTM point `a,b`, let `E=a²-b²-N`, `Y=2ab`, `H=a²+b²`. Then identically

```
H² - N² - Y² = E(2N+E).
```

Thus the squared-plane residual carries the prime divisors of E, not an independent congruence. Near a smooth real curve the residual depends on the rounding; its sawtooth behavior is real. A new algorithm needs an additional **N-only, cross-region constraint** beyond this coordinate identity. In particular, an improvement must exclude large sets of *different* u/q candidates, not only repeated Y states, at lower cost than the existing batch-GCD control.

## 7. Reproducible tests and interpretation

`python3 test_witness_theory.py` checked **595 distinct odd semiprimes**, **7,925** exact Pythagorean modular-witness checks, and counted **2,844 modular extra-witness checks** where the strict square-only paired wheel excludes a non-square witness. `python3 test_monotone_staircase.py` checked **464,800** exact floor-square-root/threshold identities on small and GMP-scale input values. Output in `witness_results.json` and `staircase_results.json`.

No factors were discovered by this test without knowledge of the sample p,q. Enumerating factor-derived witnesses is an **oracle validation of the theorem**, not a demonstrated factor locator. Real global performance must charge search for the unknown u and compare the batch-GCD implementation to an equally informed baseline.

### Suggested next experiment

Test **cross-scale** constraints rather than another fixed-modulus wheel: subdivide the N-only factor-ratio ranges, compute rational tangent/chord or continued-fraction approximations to the reciprocal curves, and combine those with modular restrictions for both `p,q` and `u` over whole boxes; require a rigorous empty-box certificate. Benchmark total candidate-equivalent work and real wall time against the existing product-tree backend and modern factoring baselines. Neither the present theorem nor fixed-modulus sieving proves a superconstant saving beyond the batch baseline.

## 8. Universal N-only Y cutoff: a global restriction (not merely local tiling)

Let `N=pq` with distinct odd primes `p<q`, and define the nearest nontrivial Pythagorean witness

```
Y_min = min( (q²-p²)/2, q*(p²-1)/2 ).
```

The two branches meet exactly when `q=p²`, equivalently `p=N^(1/3)`. The latter is the familiar cube-root scale separating small-factor search from larger-factor methods, here obtained **directly from the relative heights of the two kinds of integer points on the squared complex line**.

**Theorem (global vertical cutoff):**

```
0 < Y_min <= (N^(4/3) - N^(2/3))/2.
```

**Proof.** Put `x=N^(1/3)` and remember `q=N/p`.

- If `p<=x`, then `q>=p²`, and `Y_min=Y_extra=(N/2)(p-1/p)`. Since `p-1/p` is strictly increasing, `Y_min<=(N/2)(x-1/x)`.
- If `p>=x`, then `q<=p²`, and `Y_min=Y_square=(N²/p²-p²)/2`. This is strictly decreasing in `p`, so `Y_min<=(N²/x²-x²)/2`.
- Both upper bounds are the same `(N^(4/3)-N^(2/3))/2`.

This is a factor-independent guarantee: **above this Y bound, we can omit the whole vertical tail if the only goal is to find *some* nontrivial integer Pythagorean factor witness.** Other nontrivial witnesses may still occur there.

One completely integer-only bound is

```
c = ceil_integer_cube_root(N),
B_N = floor(N*(c²-1)/(2*c)),
Y_min <= B_N.
```

This follows because `t-1/t` is increasing for `t>0` and `c>=N^(1/3)`. It uses no floating point. Equivalently the factor-revealing divisor `u=max(q,p²)` always lies in

```
ceil_integer_cube_root(N²) <= u < N,   1<gcd(u,N)<N.
```

The exponent `4/3` is sharp for this three-witness construction in the sense that semiprimes with `q` very near `p²` come arbitrarily close to the continuous extremum in examples. For `p=1009, q=1018091`, the exact nearest Y is `518249042640`, while the integer-only bound `B_N=518762670053` is only about 0.099% larger. See `global_bound_results.json`.

**Complexity limitation:** reducing the *Y-coordinate extent* from `Theta(N²)` to `O(N^(4/3))` is a real geometric restriction but is not, by itself, an asymptotic factoring speedup. Under the exact monotone transform, the surviving distinct u values can still occupy `Theta(N)` integers. Searching the u=q branch is equivalent to searching for small p, and searching u=p² is equivalent to searching candidate p-values. A new complexity improvement must eliminate/aggregate genuinely distinct such candidates beyond ordinary interval batch GCD.

## 9. Practical synthesis: a one-certificate exclusion of the entire high-q tail

The cube-root switch can be made into a **correct two-phase hybrid**, though it reuses known batch-divisibility arithmetic rather than creating a new asymptotic algorithm.

Let `K=floor_integer_cube_root(N)` and construct the product `P_K=prod_{d=2}^K d (mod N)` using the existing block-product GCD engine (this need not be materialized as an enormous integer). Set `g=gcd(N,P_K)`.

- For `N=pq` with distinct odd primes, if `1<g<N`, a small factor is found and the computation ends. The `g=N` ambiguity cannot arise when `K<N^(1/2)` and `p,q` are distinct primes.
- If `g=1`, then `p>K`, hence `p>N^(1/3)` because `p` is an integer. It follows that `q=N/p<N^(2/3)`. **The entire factor-q interval `q>=ceil(N^(2/3))` is now ruled out**, even though all those large q candidates could have passed a small Fermat/CRT wheel individually.
- In this remaining p>N^(1/3) regime, the nearest nontrivial Pythagorean witness is the original **square-end** factor point. Consequently the *stronger* square-factor Fermat `(Y,H)` sieve is sound for a complete search of the remaining conditional semiprime domain; we no longer need to preserve the two non-square-end witnesses to be complete.
- Intersect this with the global Y bound from Section 8 and use Newton/CRT interval filters, followed by the existing polynomial batch GCD over the surviving factor intervals.

This is the desired coordinate synthesis in an exact, but limited, sense: **one global divisibility certificate** permits a large `q` tail to be excluded and changes which vertical-lattice wheel may safely be used afterward. The expensive part is producing the small-factor certificate, and the residual q range still contains many possible factors. The method is a reformulation of the established small-factor/Lehman-style split, not a demonstrated improvement over modern deterministic `N^(1/5+o(1))` algorithms or GNFS.
