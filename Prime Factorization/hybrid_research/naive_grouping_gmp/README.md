# Corrected GMP Naive Grouping worker

Standalone bounded, deterministic exact-regrouping worker; the optional `--gcd` flag broadens it to GCD-detectable factor chains. All hits are verified; failures to split are explicitly `budget_or_region_exhausted`, **not primality claims**.

```sh
make && make check
./ng_gmp 187 seed 1 1000
./ng_gmp 798607 seed 59 1000
./ng_gmp 798607 seed 13 1000 --gcd
./ng_gmp 978508015703 groups 500000 500000
```

For N=pq with distinct odd primes p<q, the original exact regrouping succeeds at total group size `g` precisely for direct factors `g=p,q`, or for `g>p` dividing `q-1` or `q+1`. Distinct chains are indexed by the odd part of `g`. `groups` is a *conditional balanced* search, using the known fact that if `p<q<2p` then one of `q-1,q+1` is divisible by 4. It checks `g == 0 (mod 4)` in `[floor(sqrt(N))-2, floor(sqrt(2N))+2]`. It does not assert a factorization theorem for arbitrary unbalanced semiprimes.

The `seed` and `seeds` modes initialize `g` using an exact quotient/remainder; every doubling then obeys `a'=floor(a/2)`, `b'=b+(a mod 2)g`, requiring no further full division by `N`. `--gcd` compares `gcd(a,b)` and `gcd(a+1,g-b)`; it has extra costs and need not help balanced cases. Every mode accepts finite work budgets. A prior C++ helper rejected valid quotient-1 cross-group successes; this implementation allows them.

**Factor-5 precheck:** Every nontrivial semiprime with last decimal digit 5 has 5 as a factor. Such inputs are handled immediately rather than routed through a dedicated Fermat residue lane.
