# Iterated-average `b` as an exact dyadic residue and a monotone vertical-lattice channel

**Status:** elementary identities proved below, GMP implementation and exhaustive small-semiprime tests. The implementation has not established a new global factoring exponent or a GNFS competitor. This note extends *Iterated Averages.md* and the squared green-line construction in *Complex space.md*. The original `VF()` uses rounded iterated averages with step-5 GCD trials and optional nearby offsets; it is **not** a numerical Newton method for integer factor roots.

## 1. Correcting the contracted-b offset when the Fermat midpoint is unknown

Let `N=pq`, `3<=p<q` odd primes, `a=(p+q)/2`, `b=(q-p)/2`, `r=ceil(sqrt(N))`, and `delta=a-r`. For `w=2^k` define the k-th average

```
A_k = N + (r-N)/w = ((w-1)*N+r)/w.
```

The corrected factor-bearing offset from `r` is `t=r-p=b-delta`, not simply `b` unless `r=a`. Exactly,

```
A_k - t/w = ((w-1)*N + p)/w.
```

This is an **integer** if and only if `q=1 (mod w)` because its numerator equals `p*((w-1)*q+1)` and p is a unit modulo w. If integral, its GCD with N contains p. Thus the apparent shrinking by 2^k is conditional on low bits of the unknown q. In particular, taking k indefinitely does not keep providing useful integer GCD candidates.

The manuscript's N=8051, p=83, q=97 has `r=a=90`, `b=7`, and `q=1 (mod 32)`. Its first three displayed averages therefore admit the correction exactly; at k=3, `(7*N+83)/8=7055` is an integer multiple of 83.

For the manuscript's N=978508015703, p=752867, q=1299709, `r=989196`, `a=1026288`, `b=273421`, `delta=37092`, hence `t=236329`. At k=6, q is not 1 modulo 64, so subtracting `b/64` is not a valid integer-divisibility argument. The actual hit comes from the generalized multiplier below.

## 2. Integer-rounded dyadic residual and its exact GCD equivalence

Let `C_k=(w-1)*N+r`, `s_k=C_k mod w` with `0<=s_k<w`, and `f_k=floor(C_k/w)`. Since `C_k=-N+r (mod w)`, `s_k=(r-N) mod w` and can be computed from N alone. For integer `j`,

```
D_k(j) = r - s_k - w*j,
    w*(f_k-j) = (w-1)*N + D_k(j).
```

Because gcd(w,N)=1,

```
            gcd(f_k-j,N) = gcd(D_k(j),N).      (1)
```

The upward side uses `D_k(-j)=r-s_k+w*j`. No large average needs to be repeatedly constructed. `D_k(j)` changes by exactly `-w` per step, permitting addition-only iteration and polynomial batch GCD over any finite j-block.

Let m be any odd integer with `q=m (mod w)` (including *negative* representatives). The index

```
             j_{k,m} = (r-s_k-m*p)/w          (2)
```

is an integer, and at this index `D_k(j_{k,m})=m*p`. Whenever it lies within the chosen positive/negative scan and `gcd(m,q)=1`, it yields exactly p by a GCD with N. This is the exact connection between iterated averages, q±1 grouping families, and the wider q±3, q±5, ... low-bit families. Since the factor q is unknown, (2) is a **retrospective witness identity**, not a fast oracle for selecting k or m.

**Worked source input:** At N=978508015703, k=6 has `w=64`, `s_6=53`, and `q=61=-3 (mod 64)`. Taking m=-3 gives `j=50746` and `f_6-j=963218792667`; the source's step-5-plus-adjacent check calls the searched offset 50745. At k=8, `w=256`, `s_8=181`, again `q=-3 (mod 256)`, and `j=12686`, with `f_8-j=974685709944`. Both values have GCD 752867 with N.

This illustrates real input-dependent speedups at particular dyadic levels, without proving that an N-only policy selects those levels cheaply.

## 3. Connecting the dyadic residual directly to the vertical green line

For `Y>=0`, define the real functions

```
H(Y) = sqrt(N^2+Y^2),
p(Y) = sqrt(H(Y)-Y),
q(Y) = sqrt(H(Y)+Y) = N/p(Y).
```

At an original **integer-complex factor point** `Y=2ab` with a=(p+q)/2, b=(q-p)/2, these are exactly the integer factors and `H=a^2+b^2`. The derivative is

```
      p'(Y) = -p(Y)/(2*H(Y)) < 0.
```

Fix a multiplier m, a level k and the dyadic residual from (1). Define the *factor-match error*

```
     R_{k,m}(Y,j) = r-s_k-w*j - m*p(Y).
```

For fixed k,m,j, this is strictly monotone on Y>0 when m != 0:

```
     dR/dY = m*p(Y)/(2*H(Y)).
```

It vanishes at the factor point when its q factor satisfies q=m (mod w) and j takes the value in (2). Over any *given* Y-band [Y0,Y1], p(Y) belongs to [p(Y1),p(Y0)]. Consequently all possible j-indices for that fixed m lie between the endpoint values `(r-s_k-m*p(Yi))/w`. If this interval misses a searched j-block, the complete Y-band is excluded **for that specified modular lane m**. Exact rational or integer-square-root enclosures can certify these endpoint bounds.

**Essential limit:** A full factor-free certificate must cover every possible odd q-residue m modulo w, not only promising small m. A short m-set such as {±1,±3} is only a conditional fast path. Without coverage of the other residue classes the monotone error cannot certify that an arbitrary factor is absent. This is the same missing global arithmetic information identified by the earlier direct-Y Fermat wheel.

## 4. Bringing in all three integer Pythagorean witnesses

The **extended integer** green hyperbola `H²-Y²=N²` includes integer points beyond the original square-complex image. Put `u=H-Y`, so `u|N²` at any integer Pythagorean point. For distinct-prime odd N=pq, the three proper-factor witnesses have `u=p`, `u=q`, or `u=p²`. Define the signed dyadic average

```
          Z_{k,eps}(u) = ((w-1)*N + eps*u)/w,
          eps in {-1,+1}.
```

If `eps*u=N (mod w)`, the average is integral, and

```
               gcd(Z_{k,eps}(u),N)=gcd(u,N).  (3)
```

This follows since multiplying by w, a unit modulo N, turns the average into eps*u modulo N. Thus averages preserve the factor-bearing property of any of the three Pythagorean witnesses when their low bits permit integral lifting. In continuous Y-coordinates, `u(Y)=H(Y)-Y` obeys `u'(Y)=-u(Y)/H(Y)<0`. Consequently the real signed average is also monotone, although its integrality/divisibility pattern is discontinuous.

**A universal two-step witness theorem.** At least one `u in {p,q,p²}` satisfies `u=N (mod 4)`:

- If N=1 mod4, choose u=p²=1 mod4.
- If N=3 mod4, p and q have opposite residues mod4, so the one equal to 3 matches N.

Therefore at least one *positive-sign* two-step averaged Pythagorean witness is an integer with nontrivial GCD. This is an existence theorem, **not** an N-only algorithm for locating the witness.

**A universal signed three-step theorem.** At least one u among the same three satisfies `u=±N (mod 8)`:

- If N is ±1 mod8, p²=1 mod8 supplies the appropriate sign.
- If N is ±3 mod8, one of p,q is 3 or 5 mod8, supplying a sign.

This does **not** extend indefinitely. For p=3, q=5, N=15, none of u=3,5,9 is ±15 modulo16. Deeply contracted dyadic witnesses require special q/p low-bit structure and cannot be guaranteed from semiprimality alone.

## 5. What can be excluded without visiting vertical lattice points

- Choose a Y-band. Invert it to a rigorous interval of possible u=H-Y or p(Y)=sqrt(u(Y)). This does not enumerate its Y points.
- Apply the necessary Fermat/witness paired CRT residues from the prior package. Note the square-complex (u=p²) and all-Pythagorean (u=p,q,p²) wheels differ.
- For each selected dyadic m-lane use the **monotone interval for R_{k,m}** to discard the entire band if it cannot intersect the bounded j-search.
- For the remaining lanes batch products of arithmetic-progression residuals `D_k(j)` modulo N and take GCD. A GCD of 1 excludes the block of those specific residuals. GCD N is ambiguous and requires splitting.

This is an exact architecture. But every q unit-residue occurs in one m lane, so *covering every lane* remains necessary for an unconditional completeness claim. It is not a subexponential factor-region locator as presently formulated.

## 6. Independent tests, synthetic census, and honest interpretation

`make test` validates the GMP batch worker against reference GCD calculations. `test_dyadic_vertical.py` checks 5,000 distinct-prime small semiprimes, exact residual identity (1), equation (2), the three-witness GCD property (3), the 2- and 3-step existence claims, and exact vertical monotonicity. In the recorded run it passed **1,098,672 exact checks**. Factor identities in the oracle tests deliberately use known p,q to establish correctness; the GMP worker receives N only.

`benchmark_retrospective.py` constructs 600 deterministic-seed balanced semiprimes with p in approximately [430,000,670,000] and q<2p, then **uses their known factors retrospectively** to calculate the minimum GCD offset at each k. This is not factoring runtime, and the oracle-best k cannot be selected from N by the tested rule.

Representative median offset counts from the fixed-seed census:

| Policy | Median first-hit offset |
|---|---:|
| k=1 only | 56,156 |
| k=6 only | 93,241.5 |
| Source fractional-part selector among k=1..6 | 80,195.5 |
| Best of k=1..6 (factor-informed oracle) | 16,486 |
| Best of k=1..14 (factor-informed oracle) | 9,903 |
| Ordinary Fermat midpoint offset (different operation cost) | 11,517 |

The source fractional-part selector picked the best of six levels in **18.5%** of inputs. Its first-hit offset beat the Fermat offset in **12.5%** of inputs, without comparing GMP GCD time to square-test time. Scanning all K levels incurs roughly K times as many candidate terms at a common depth unless sufficiently shared polynomial evaluation, deduplication, or genuine region exclusion is obtained. The observed oracle improvement is **not** a free algorithmic speedup.

The GMP batch worker reproduced the manuscript N=978508015703, k=8 hit at j=12686. It checked 25,373 signed terms (both directions), used 25,550 modular products and 63 GCD calls with batches of 256, and verified the factor. This is an implementation improvement in GCD count, not a change in candidate asymptotics.

## 7. Key next theorem

The precise missing claim for globally improved complexity is:

> Given only N, can several vertical intervals or many dyadic lanes be ruled out **collectively**, by a certificate cheaper than evaluating the corresponding batch products, in a way that excludes a superconstant fraction of *distinct possible factor candidates* while retaining all three integer Pythagorean witnesses?

Neither the present dyadic contraction nor the monotone geometry alone supplies that theorem. The present work nevertheless proves that iterated averages encode a **testable 2-adic structure** of factor witnesses and supplies an exact, GPU-compatible, batched-linear-residual experiment.
