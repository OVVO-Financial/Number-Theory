# Real-pair complex and quaternion arithmetic

Companion implementations and executable tests for Fred Viole's paper
*Real-Pair Representations of Complex Numbers and Quaternions: Algebra,
Arithmetic Coordinates, and Implementations in R and Python*.

All arithmetic is implemented in real pairs or quadruples. Built-in
complex values are used for conversion/reference checks only.

## Files

- `Real_Pair_Representations.R`: base R implementation
- `Real_Pair_Representations_tests.R`: base R assertions and 1,000 seeded randomized cases
- `Real_Pair_Representations.py`: standard-library Python implementation
- `Real_Pair_Representations_tests.py`: Python unittest suite with 10,000 seeded randomized quaternion cases
- `Real_Pair_Representations_examples.py`: worked Python examples
- `Real_Pair_Representations_independent_check.py`: separate reference formula cross-check

## Run from this directory

```sh
Rscript Real_Pair_Representations_tests.R
python3 Real_Pair_Representations_tests.py
python3 Real_Pair_Representations_examples.py
python3 Real_Pair_Representations_independent_check.py
```

Python requires 3.9+ without third-party packages. The R code uses base R.
The principal conventions are `phi(a+bi) = (a-b, a+b)` and
`psi(a+bi+cj+dk) = (a-b, a+b, c-d, c+d)`. Quaternion division is
right division: `p/q = p * inverse(q)`.

The scripts include exact basis cases (`ij=k`, `ji=-k`, `j^-1=-j`),
quaternion norm multiplicativity, both-sided inverse identities, and
numerical comparisons against the standard coefficient-by-coefficient
quaternion multiplication formula. The R suite is included for
reproducibility; run it in an R installation before reporting its results.
