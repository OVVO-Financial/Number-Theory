# Vertical green-line: three Pythagorean witnesses, exact monotone staircase, global Y bound

This standalone research addition addresses the original complex-space `Y=2ab` green line **and its larger integer Pythagorean lattice**. Not every Pythagorean point is the square of an integer Gaussian complex factor point. The module differentiates these two lattices and proves that the larger one has three nontrivial factor-revealing points for `N=pq` (p and q distinct odd primes). It supplies an N-only `O(N^(4/3))` **global Y cutoff** and an exact monotone integer step function in Y.

Run independent Python-only exact tests (no third-party dependencies):

```sh
python3 test_witness_theory.py
python3 test_monotone_staircase.py
python3 test_global_witness_bound.py
```

Proof: `GEOMETRIC_WITNESSES.md`. Raw reproducible outputs: `witness_results.json`, `staircase_results.json`, `global_bound_results.json`.

**Important limitations:** These are arithmetical certificates and oracle verification, not a new efficient general factor locator. Knowing that a factor-bearing u lies inside `[N^(2/3),N)` does not reveal which u. The monotone vertical-coordinate staircase is a compression of the same divisor search, and the batch-GCD certificate is exactly a known product-tree factor test after coordinate inversion. The previous Fermat-wheel implementation is correct for *square-factor points* but must NOT be applied unmodified to the additional non-square Pythagorean points: an independent wheel is derived and tested.
