# A factor-independent dyadic lane cover and two genuinely different complex-plane arithmetic channels

**Status:** exact elementary theorems, finite optimized covers, CPU/GMP implementations and bounded numerical comparisons. **No new general factorization exponent and no result competitive with GNFS on generic RSA semiprimes.** References: the repository's `Complex space.md`, `Iterated Averages.md`, `Naive Grouping.md`, and the v4 research package's `DYADIC_VERTICAL_THEORY.md`.

## 1. Exact lane constraints and a complete cover

Let `N=pq`, where `p<q` are distinct odd primes. Fix `w=2^k`, `k>=3`, and let `n=N mod w`. All inverses here are in `(Z/wZ)^x`. In the prior dyadic notation `D_k(j)=r-s_k-wj`, `s_k=(r-N) mod w`, `r=ceil(sqrt(N))`, so

    D_k(j) == N (mod w).

The integer Pythagorean factors of `N^2` have nontrivial witnesses `u=p, q, p^2`. A multiplier `m` satisfying `m u == N (mod w)` must be, respectively,

    u=p:    m=q=n*p^(-1) (mod w)
    u=q:    m=p
    u=p^2:  m=q*p^(-1)=n*p^(-2) (mod w).

Every odd square is `1 (mod 8)`. Therefore the **N-only complete square-witness lane set**

    M_k(N) = {m (mod 2^k): m==N (mod 8)}

has exactly `2^(k-3)` members, and for every factor pair, at least one genuine `u=p^2` witness has its multiplier in this set. This is a factor-four reduction over all odd multiplier residues. All lane values may be selected from `N` alone; no factor is used by the selector.

For an **arbitrary complete three-witness lane cover** `M`, fix `n` and let `x=p mod w` range over the `phi(w)=2^(k-1)` locally admissible odd residues. Its three possible multipliers are

    x, n*x^(-1), n*x^(-2) mod w.

One fixed `m` hits at most **one** `x` through the first equation, **one** through the second, and **four** through `x^2==n*m^(-1) (mod 2^k)`; there are at most four odd roots of any quadratic congruence at `k>=3`. Consequently

    |M| >= ceil(2^(k-1)/6).

This is a lower bound for purely residue-based, fully universal coverage; it is **not** a lower bound on factoring algorithms using additional arithmetic information from the complete integer `N`, and it is **not** an operation-count lower bound for batched tests. For signed-magnitude lanes, count `+m` and `-m` separately or account for their pair factor.

An exact finite mixed-integer set-cover solver (SciPy/HiGHS) found:

| k | All odd residues | Full square-witness cover | Optimal 3-witness cover for `N==-1 mod 2^k` |
|---|---:|---:|---:|
| 6 | 32 | 8 | 6 |
| 7 | 64 | 16 | 11 |
| 8 | 128 | 32 | 22 |
| 9 | 256 | 64 | 43 |

This is a nontrivial **interaction between witness families** but remains exponential in `k`. Optimization was performed on finite residue tables, not on unknown factors. Exact values for other `N mod 2^k` can differ; the `N==-1` column does not apply to all `N==7 (mod 8)`. The ILP is an analysis tool, not a claimed cheap high-k runtime strategy.

**Important duplication identity.** Independently of how many multipliers are selected, *the actual GCD trial integers at a fixed depth k* are the **single** arithmetic progression

    D_k(j) == N (mod 2^k).

Every `D_(k+1)` value is also a `D_k` value. Multiplier labels explain how a hit could arise, but they do not produce distinct new progressions. To achieve a speedup, the selected set must reduce the **numerical range** searched and that reduction must outweigh setup, multiple-lane, and witness costs. In particular `u=p^2` has scale `N` for balanced inputs, unlike `u=p,q` of scale `sqrt(N)`. A low lane count alone does not yield a fast algorithm.

## 2. A new exact bridge from `2ab` to CRT square roots of one

For *any integer* complex trial point `z=a+bi`, put

    A=a^2-b^2,  E=A-N,  Y=2ab,  H=a^2+b^2.

The Pythagorean identity `H^2=A^2+Y^2` holds regardless of the accuracy of the trial point. Suppose `gcd(H,N)=1` and set `S=Y*H^(-1) mod N`. Then

    S^2-1 == -(E*H^(-1))^2  (mod N).                  (I)

For squarefree `N=pq`, squaring does not change whether each prime divides the residual, so

    gcd(N,S^2-1) = gcd(N,E).                          (II)

**Consequence:** the normalized `2ab` square-root residual itself carries *exactly the same prime-divisibility information* as the original CTM residual `E`; it cannot create a new global factor detector by a monotone rescaling. If `H` is a nonunit modulo `N`, however, `gcd(H,N)` is an **independent detection channel** (Section 5).

At the genuine original factor point `a=(p+q)/2, b=(q-p)/2`, `A=N`, `H-Y=p^2`, `H+Y=q^2`, and `H` is a unit modulo the squarefree `N`. Hence

    S==+1 (mod p),   S==-1 (mod q).

The true vertical factor point determines a nontrivial CRT square root of 1. Recovering it without first locating that point is the difficult task.

**Another no-go identity:** if we normalize the squared-line point `W=N+iY` by its conjugate, then when `Y` is a unit modulo `N`,

    W/conj(W) == -1 (mod N).

That normalization collapses the entire vertical green line to one constant modular element. To extract nontrivial multiplicative-order information we must use an off-line seed and the underlying complex square, not only `W` itself.

## 3. Modular Newton has a conjugacy, not a contraction

Newton on `S^2-1` over `Z/NZ` (where 2 and S are units) reads

    S_new = (S^2+1)/(2*S).

The exact residual update is

    S_new^2-1 = (S^2-1)^2/(4*S^2).                    (III)

Thus as long as divisions remain legal, the set of prime factors of `N` dividing the residual **cannot increase**. If `gcd(S^2-1,N)=1` initially, no number of Newton steps can make the residual zero modulo one of the prime factors. This disproves literal modular Newton convergence as a universal factor locator.

Under the Möbius coordinate `t=(S-1)/(S+1)`, wherever `S+1` is a unit,

    t_new = t^2.

For a complex trial point `S=Y/H`, this simplifies to

    t == - ((a-b)/(a+b))^2  (mod N).                 (IV)

Therefore repeated squaring of `t` tests a power-of-two subgroup condition. If `t^(2^j)==-1` modulo one prime factor but not the other, `gcd(t^(2^j)+1,N)` reveals that factor; equivalently a *pole* of the Möbius/Newton chart appears on one CRT component. It is **conditional multiplicative-order factoring** (a 2-primary form of Pollard `p-1`), not generic real-valued monotone convergence.

`newton2_primary.cpp` implements exact GMP squaring and GCDs. On the structured example `N=65537*65543`, its chosen seed `(a,b)=(2,1)` recovers 65537 after 14 squaring steps. On 120 random/smooth-neighbor semiprimes in the recorded bounded corpus, this pure power-of-two path found no factor. Both statements must be carried together.

## 4. Normalizing a generic complex square gives conditional p-1/p+1 group factoring

For a seed `z=a+bi` whose **norm is a unit modulo N**, define the norm-one Gaussian element

    xi = z / conjugate(z)
       = (a^2-b^2 + i*(2ab))/(a^2+b^2)      (mod N).

This directly reuses the original `2ab` numerator. The seed need not lie on the actual green line. Over an odd prime `r|N`, Gaussian conjugation gives a cyclic norm-one subgroup of order dividing

    r - chi_r(-1) = r-1 if r==1(mod 4), or r+1 if r==3(mod 4).

If one prime's group order divides a chosen smooth exponent `L`, then `xi^L==1` modulo that prime; a GCD of `N` with its **two** coordinate residuals `(Re xi^L-1, Im xi^L)` extracts a factor provided the other prime's coordinate pair is not also 1. This is the established Lucas / Williams `p+1`-type mechanism in Gaussian language, not a new asymptotic theorem.

When `N==3(mod 4)`, one unknown factor is `1(mod 4)` and the other `3(mod 4)`: the two Gaussian group-order types are known to differ **from N alone**. For `N==1(mod 4)`, a Lucas discriminant `D` with Jacobi `(D/N)=-1` can provide an analogous quadratic-character split. Both methods exploit *group structure* rather than enumerating `p` or `q` candidates.

`gaussian_norm_one.cpp` implements this with small deterministic seeds, GMP modular Gaussian products, and prime-power stages. `smooth_order.cpp` independently implements standard stage-one Pollard `p-1`, Lucas `p+1`, and bounded Brent rho. `hybrid_native.cpp` calls these together in one C++ process.

**Crucial distinction:** fixed small Gaussian seeds and near-strip seeds had similar success counts on the tested corpus. There is no measured evidence that choosing an arbitrary seed closer to the green boundary improves smooth-group orders.

## 5. The independent complex norm channel

Checking `gcd(a^2+b^2,N)` alongside `gcd(a^2-b^2-N,N)` targets different modular slopes. For `p==1(mod 4)`, two nonzero roots of `-1 (mod p)` give solutions of `a^2+b^2==0 (mod p)` independent of the slopes `a==+/-b` in the CTM error. Fermat's sum-of-two-squares theorem guarantees `p=a^2+b^2` for some small integer pair when p is such a prime.

If `N==3(mod 4)` with distinct odd prime factors, exactly one prime factor is `1(mod 4)`, so that factor has a norm witness. Under the *balanced assumption* `p<q<2p`, the coordinates of that witness lie within `0<=a,b<= (2N)^(1/4)`, but the two-dimensional box contains `Theta(sqrt(N))` integer candidates. This is an additional target family and useful regional statement, **not** a sub-fourth-root factor search.

The double-channel batch product can use `E*H` modulo N, and a GCD at each batch boundary. It retains the same candidate states and sometimes recovers a norm witness earlier.

An early census script mistakenly omitted the case `E==0` (an exact Fermat factor hit), dramatically overstating the norm benefit. The **corrected** census of 100 cases at each of three sizes found earlier norm hits in **3, 3, and 1 cases**, respectively. `dual_residual.cpp` reproduces an actual independent earlier norm hit for `N=243866858267`: at `a=ceil(sqrt(N))+2392`, the `H` channel extracts 594329, while the `E` exact factor is reached at offset 8497. No universal performance advantage is claimed.

## 6. Measured performance and present limits

The bounded GMP C++ experiment used 120 semiprimes with prime factors of approximately 18–30 bits: 40 ordinary controls, 40 with a planted B-smooth `p-1`, and 40 with a planted B-smooth `p+1`; all factor verification was exact. The workers receive only `N`.

**Matched in-process benchmark (B=256; times are median microseconds across 40 per class):**

| Corpus | Complete research hybrid incl. rho | Brent rho only | Hybrid/rho |
|---|---:|---:|---:|
| Ordinary/random | 1229 | 573 | 2.14 (hybrid slower) |
| `p-1` B-smooth planted | 161 | 542 | 0.30 (hybrid ~3.4x faster) |
| `p+1` B-smooth planted | 159 | 551 | 0.29 (hybrid ~3.5x faster) |

A Gaussian-only B=1024, four-seed run recovered factors on 22/40 ordinary, 37/40 p-1 smooth, and 34/40 p+1 smooth cases in a separate sample from the same corpus. Success is input-family dependent. These are **small to moderate sizes**, not GNFS benchmarks. On bounded 330- and 364-bit RSA-100/RSA-110 checks, the hybrid found no factors.

No module proves universal factorization in a bounded fast path. The `hybrid_native` fallback is **bounded rho**, and an unfinished search must report **unresolved**, never **prime**. For deterministic full balanced-domain coverage, connect this experimental fast path to the existing v4 product-tree engine and charge all of its costs. Its known fourth-root class does not challenge GNFS asymptotically.

## 7. The next mathematically meaningful target

The residue-only dyadic lanes cannot give a fixed-size universal family at growing depth, and modular Newton merely squares an underlying group coordinate. Both restrictions are proved, not implementation failures. The new productive direction is **N-only quadratic-character / norm-one group choice** combined with a global batched fallback. A genuine general exponent gain would require proving that a cheap `N`-derived choice forces small order for one prime and not the other, or produces unusually dense smooth relations compared against equally informed QS/NFS baselines. No such theorem is established by the results here.
