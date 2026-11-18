# Exact grouping multiplicity and conditional balanced-band theorem

Let `N=pq` with distinct odd primes `3<=p<q`, and define for each positive integer group count `g`:

    a=floor(N/g), b=N mod g, c=g-b,
    N=ga+b=ca+b(a+1)=g(a+1)-c.

A **direct** success has `g=p` or `g=q`. When `1<b<g`, `b|a` implies `N=b((a/b)g+1)`, hence `b=p`, `g|(q-1)` and `g>p`. Conversely if `g|(q-1)`, `g>p` and `q-1=kg`, then `N=p(kg+1)`, `a=pk`, `b=p`, so the exact regrouping succeeds. Likewise `c|a+1`, `c>1` is equivalent to `g|(q+1)` and `g>p` (the group count is admissible and `c=p`). Thus the full mathematical success set is

    {p,q} union {g|q-1:g>p} union {g|q+1:g>p}.

It has cardinality `2+D_{>p}(q-1)+D_{>p}(q+1)`, since the latter two families do not intersect for `g>p>=3` (common divisor divides 2). This counts successful **group sizes**, not independent doubling chains. Each chain is indexed by the *odd part* of `g`.

Writing `q±1=2^h u` with odd `u`, the successful odd roots arising from this neighbor are exactly divisors `s|u` satisfying `2^h s>p`. Direct roots `p` and `q` must be counted separately, removing any overlap with `p` when `p|q²-1`.

For **balanced semiprimes** `p<q<2p`, every proper divisor of `q-1` or `q+1` is `<=p`; therefore only `q-1` and `q+1` appear as non-direct successful group sizes. Exactly one of those is divisible by four. Its odd part is the unique successful odd root below `sqrt(N)/2`, provided seed 1 and quotient-one regrouping are admitted.

The N-only conditional band is the union, over `h>=2`, of odd `s` with

    ceil((ceil(sqrt N)-1)/2^h) <= s <= floor((floor(sqrt(2N))+1)/2^h),
    4*s*s < N.

It contains that unique root. Its asymptotic size is `((sqrt(2)-1)/4)*sqrt(N)+O(1)`, retaining approximately `41.421%` of the naive small odd-root range. This is a **constant-factor search-space reduction**, not a proof of GNFS-level speed or subexponential factorization.

A broader `gcd(a,b)` or `gcd(a+1,c)` test detects additional factors in some unbalanced cases (e.g. `798607`), but should be analyzed separately: it is not the same predicate as exact regrouping and costs additional GCD work. The executable `ng_gmp.cpp` implements both with `--gcd` as an explicit option. Search budgets and unsuccessful results are never advertised as a proof of primality.

**Important scope:** these theorems assume `N` is a product of distinct odd primes; prime powers and general composites require separate treatment. Computing the complete successful set from `q±1` presupposes knowing `q`, so the formula is an explanation of multiplicity, not a direct factorization algorithm.
