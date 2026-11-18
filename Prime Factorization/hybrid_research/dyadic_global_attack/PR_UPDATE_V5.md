## V5: N-only dyadic covering / modular Newton / Gaussian group-order / independent norm residual

**Base research bundle:** `number_theory_hybrid_pr_ready_v4.zip`. Add the entire `Prime Factorization/hybrid_research/dyadic_global_attack/` directory.

### New exact results

1. Three-witness dyadic lane multipliers are `p, q, q/p (mod 2^k)`. The **N-only square-witness cover** `m==N (mod 8)` is always complete and has `2^(k-3)` lanes, while *every* complete three-witness residue-only cover needs at least `ceil(2^(k-1)/6)` lanes. Offline exact ILP reduced the `N==-1(mod 2^k)` cases to 6, 11, 22, 43 lanes for k=6,7,8,9 respectively.
2. Every dyadic GCD trial at level `k` lies in the one arithmetic progression `D_k==N (mod 2^k)` and higher-depth progressions are nested: lane count does not equal independent new candidate families.
3. With `E=a^2-b^2-N, Y=2ab, H=a^2+b^2` and `H` invertible mod `N`, the normalized vertical residual satisfies `(Y/H)^2-1 == -(E/H)^2 (mod N)`. For squarefree semiprimes, the corresponding GCD equals `gcd(E,N)`. **This is a rigorous limitation on coordinate rescaling.** At the genuine factor point, `Y/H` is a nontrivial CRT square root of 1.
4. Newton iteration on `S^2=1` conjugates exactly to `t->t^2`, with `t=-((a-b)/(a+b))^2 (mod N)`. The root residual's prime support cannot increase while divisions remain legal. The bounded dyadic-GCD fast path thus detects particular 2-primary group orders, not modular convergence.
5. The Gaussian norm-one element `(a^2-b^2+i2ab)/(a^2+b^2)` has order modulo prime `r` dividing `r-chi_r(-1)`: `r-1` if r=1(mod4), `r+1` if r=3(mod4). This is a direct connection to the known Pollard `p-1` / Lucas `p+1` families. When `N=3(mod4)`, the two order types differ **from N alone**. A corrected C++/GMP module uses group exponentiation instead of scanning vertical lattice points.
6. `gcd(N,a^2+b^2)` adds distinct factor witnesses for primes `1(mod4)`. Corrected benchmarks **count exact E=0 Fermat hits**. The extra norm channel won in 3, 3, and 1 of 100 bounded cases at three input scales; early figures that omitted E=0 are superseded.

### Code / test scope

`Makefile`, `hybrid_native.cpp`, `gaussian_norm_one.cpp`, `newton2_primary.cpp`, `smooth_order.cpp`, `dual_residual.cpp`, tests, raw input corpus, ILP optimizer, timings, README and `DYADIC_VERTICAL_SYNTHESIS.md`. `make test` passes; Python reference includes 701,496 lane-cover checks, 48,000 Newton-conjugacy/character checks, 9,840 normalized-vertical identities, 1,445 Gaussian finite-field identities, 86,410 norm-witness checks, and 14 native CLI regression cases. Counts are per-suite; not all represent disjoint tests.

### Measured outcomes

120 semiprimes (40 ordinary, 40 planted p-1-smooth, 40 planted p+1-smooth), common C++17/GMP binary and same rho budget: random median 1229 us hybrid vs 573 us rho; planted p-1 median 161 vs 542 us; planted p+1 median 159 vs 551 us. None of the bounded paths recovered RSA-100/RSA-110 under the recorded budgets.

No GNFS performance claim, general global-search exponent gain, GPU device speedup, or proof that small factor-independent dyadic lanes universally contain an inexpensive early factor. A correct balanced-domain fallback remains in v4's product-tree component; bounded rho and group searches must report *unresolved*, not *prime*.
