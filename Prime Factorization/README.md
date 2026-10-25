# Prime Factorization Routines

This section presents the methods of factorization derived from the complex space insights detailed [here](../Complex%20space.md).

* The method of [Complex Trial Multiplication](Complex%20Trial%20Multiplication.md)

* The method of [Iterated Averages](Iterated%20Averages.md)

* The method of [Naive Grouping](Naive%20Grouping.md)

* The method of [Simultaneous Complex Factorization](Simultaneous%20Complex%20Factorization.md)

---

## Unified engine

[**Unified Complex Space Factorization**](Unified%20Complex%20Space%20Factorization.md) combines all
of the above with the corrected Fermat terminal-digit sieve into a single deterministic
`O(N^(1/3))` routine — [`julia/complex_space_factorization.cpp`](../julia/complex_space_factorization.cpp).

It shows that the corrected terminal-digit sieve *is* the quadratic-residue sieve at
modulus 20, and generalises it to a CRT wheel worth ~420,000x instead of 5x.
