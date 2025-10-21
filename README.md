# Number Theory

This project investigates prime factorisation and analytic number theory by expressing complex numbers as **simultaneous pairs of real values**.  Treating the imaginary unit in this way lets us reason about complex arithmetic, visualise factor searches, and implement sieves while staying entirely within the real plane.

The repository now lives on a single branch and the documents referenced below are kept alongside the code.  All of the links have been updated to use relative paths so they work both locally and on GitHub.

## Overview

* [Factors in the complex space](Complex%20space.md) introduces the simultaneous-real mapping and illustrates how trial division, Fermat's method, and complex trial multiplication interact on the complex lattice.
* The `Prime Factorization` directory gathers algorithms that combine these insights into practical routines.
* The `R` and `julia` directories contain the code used to generate the visualisations and to prototype the factorisation methods.

## Mapping simultaneous reals

The simultaneous representation of a complex value `a + bi` is defined as `[(a - b), (a + b)]`.  A few examples:

| Complex number | Simultaneous real pair |
| -------------- | ---------------------- |
| `i`            | `[-1, 1]`              |
| `2i`           | `[-2, 2]`              |
| `-i`           | `[1, -1]`              |
| `8 + 5i`       | `[3, 13]`              |
| `-2 - 3i`      | `[1, -5]`              |

These pairs are the foundation for the geometric interpretations used throughout the project.

## Repository structure

| Path | Description |
| ---- | ----------- |
| [`Complex space.md`](Complex%20space.md) | Narrative overview of the complex-space factor strip and how different algorithms traverse it. |
| [`Prime Factorization/`](Prime%20Factorization/README.md) | Detailed write-ups for each factorisation method together with Julia source files. |
| [`Prime Factorization/README.md`](Prime%20Factorization/README.md) | Entry point for the factorisation notes shown below. |
| [`R/`](R/) | R scripts for generating the diagrams plus a gallery of their output. |
| [`julia/`](julia/) | Julia implementations of the algorithms, including supporting utilities such as `Newton_Square_Root.jl`. |
| [`Number Theory Papers/`](Number%20Theory%20Papers/) | Research papers that motivate the approach and provide mathematical background. |
| [`Images/`](Images/) | Static images embedded throughout the documentation. |

## Factorisation methods

The factorisation routines derived from the complex-space analysis are documented here:

* [Complex Trial Multiplication](Prime%20Factorization/Complex%20Trial%20Multiplication.md)
* [Iterated Averages](Prime%20Factorization/Iterated%20Averages.md)
* [Naive Grouping](Prime%20Factorization/Naive%20Grouping.md)
* [Simultaneous Complex Factorization](Prime%20Factorization/Simultaneous%20Complex%20Factorization.md)

Each note explains the intuition, provides worked examples, and links to the corresponding Julia implementation.

## Working with the code

* **Julia prototypes:** The `julia/` directory mirrors the documentation above.  For example, `julia/Simultaneous_Complex_Factorization.jl` implements the algorithm described in the accompanying markdown file.
* **R visualisations:** R scripts such as `R/Complex Space Generator.R` reproduce the plots shown in the documentation.  The `R/README.md` file contains sample output and screenshots of the generated graphics.

## References

Additional background, proofs, and derivations can be found in the [Number Theory Papers](Number%20Theory%20Papers/) collection.  These papers formalise the simultaneous-real mapping, explore the resulting factor properties, and describe supporting sieves.
